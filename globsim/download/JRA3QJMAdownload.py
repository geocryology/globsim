import logging
import netCDF4 as nc
import numpy as np
import pandas as pd
import pygrib
import xarray as xr

from abc import ABC, abstractmethod
from datetime import datetime, timedelta
from os import rename
from pathlib import Path
from typing import Optional, Union
from urllib.error import URLError

from globsim.download.GenericDownload import GenericDownload
from globsim.download.JRA3Q_helpers import GribSubsetter, download_daily_gribs, download_constant
from globsim.download.JRAdownload import getDate
from globsim.download.JRA3Q_dl import GetAccessor
from globsim.download.jra_dict_formatters import getPressureLevels


logger = logging.getLogger('globsim.download')

CHUNK_LAT = 4
CHUNK_LON = 4


class J3QD(GenericDownload):

    _tunits = "hours since 1947-01-01"
    _tcal = "standard"
    REANALYSIS = 'jra3q'

    def __init__(self, pfile):
        super().__init__(pfile)
        self.retry_delay_min = 1
        par = self.par
        self._set_input_directory("jra3q")
        levels = getPressureLevels(GribSubsetter.DEFAULT_LEV_HPA, self.elevation['min'], self.elevation['max'])
        self._subsetter = GribSubsetter(round((float(par['bbW']) - 1.25) / 1.25) * 1.25,
                                        round(float(par['bbE']) / 1.25) * 1.25,
                                        round((float(par['bbS']) - 1.25) / 1.25) * 1.25,
                                        round(float(par['bbN']) / 1.25) * 1.25,
                                        levels=levels)

        # time bounds
        self.date  = getDate(par)

        # Connect
        self.credential = Path(par['credentials_directory'], ".netrc")
        self.connect()

        # chunk size for downloading and storing data [days]
        logger.info("Chunk size ignored for JRA3Q download. Using monthly")
        # self.chunk_size = par['chunk_size'] * 2000

        # variables
        logger.info(f"Variable list ignored for {self.REANALYSIS}")

        if par.get("download_only") is not None:
            self.download_only = par.get("download_only")
        else:
            self.download_only = False

    def file_in_progress(self, filename):
        return str(Path(filename).with_suffix(".inp"))
    
    def complete(self, filename):
        f = self.file_in_progress(filename)
        rename(f, filename)

    def connect(self):
        try:
            self.access = GetAccessor(netrc_file=self.credential)
        except ValueError:
            logger.error("Could not find credentials. Aborting")
            raise ValueError(f"Could not find credentials file {self.credential}")
        except URLError:
            logger.error("Could not connect to JRA server. This aint gonna work.")

    @staticmethod
    def date2num(dt:datetime) -> float:
        t = nc.date2num(dt, J3QD._tunits, J3QD._tcal)
        t = np.round(t, 5)  # catch any floating-point errors
        return t

    def nc_filename(self, year, month, kind="sa") -> str:
        return f"{self.REANALYSIS}_{kind}_{year}{month}01_to_{year}{month}31.nc"

    def retrieve(self):
        days = pd.date_range(self.date['beg'], self.date['end'] + timedelta(days=1))  # add extra day but end when we get to it

        # Constant surface data
        to_grib = [download_constant(self.access, self.directory)]
        t = days[0]
        to_ncfile = Path(self.output_directory, f'{self.REANALYSIS}_to_{t.strftime(r"%Y%m%d")}_to_{t.strftime(r"%Y%m%d")}.nc')
        gribs_to_netcdf(*to_grib,
                        netcdf_file=to_ncfile,
                        tstart=self.date2num(t),
                        tstop=self.date2num(days[1]),
                        kind="to",
                        subsetter=self._subsetter)

        # Temporal data
        running_sa, running_sf, running_pl, running_dates = [], [], [], []

        for date in days:
            end_of_days = (date == days[-1])
            
            # Aggregation subroutine
            if (not self.download_only) and ((date.day == 1 and len(running_dates) > 1) or end_of_days):
                t1 = running_dates[0]
                t2 = running_dates[-1]
                tstart = J3QD.date2num(t1)
                tstop = J3QD.date2num(date)  # t2 + timedelta(days=1)

                sa_ncfile = Path(self.output_directory, f'{self.REANALYSIS}_sa_{t1.strftime(r"%Y%m%d")}_to_{t2.strftime(r"%Y%m%d")}.nc')
                sf_ncfile = Path(self.output_directory, f'{self.REANALYSIS}_sf_{t1.strftime(r"%Y%m%d")}_to_{t2.strftime(r"%Y%m%d")}.nc')
                pl_ncfile = Path(self.output_directory, f'{self.REANALYSIS}_pl_{t1.strftime(r"%Y%m%d")}_to_{t2.strftime(r"%Y%m%d")}.nc')

                # aggregate gribs then delete them
                logger.info(f"Aggregating gribs to netcdf: {t1} to {t2}")

                if not sa_ncfile.is_file():
                    gribs_to_netcdf(*running_sa, netcdf_file=self.file_in_progress(sa_ncfile),
                                    tstart=tstart, tstop=tstop, kind="sa", subsetter=self._subsetter, engine="pygrib")
                    self.complete(sa_ncfile)
                
                if not sf_ncfile.is_file():
                    gribs_to_netcdf(*running_sf, netcdf_file=self.file_in_progress(sf_ncfile), 
                                    tstart=tstart, tstop=tstop, kind="sf", subsetter=self._subsetter, engine="pygrib")
                    self.complete(sf_ncfile)
                
                if not pl_ncfile.is_file():
                    gribs_to_netcdf(*running_pl, netcdf_file=self.file_in_progress(pl_ncfile),
                                    tstart=tstart, tstop=tstop, kind="pl", subsetter=self._subsetter, engine="xarray")
                    self.complete(pl_ncfile)
                
                #  delete grib files
                logger.info("Not deleting grib files")
                #for grib in running_sa + running_sf + running_pl:
                    #Path(grib).unlink()

                # clear running
                running_sa, running_sf, running_pl, running_dates = [], [], [], []

            if end_of_days:
                break

            # Download
            logger.info(f"Downloading gribs for {date}")
            sa, sf, pl = download_daily_gribs(self.access, self.directory, date.year, date.month, date.day)
            running_sa += sa
            running_sf += sf
            running_pl += pl
            running_dates.append(date)


class VariableError(ValueError):
    pass


class ConversionHandler(ABC):
    """
    Handles the translation of grib files into a netcdf file
    """
    _include = []
    _dims = ()

    def __init__(self,
                 tstart: float,
                 tstop: float,
                 subsetter:Optional[GribSubsetter] = None):

        self._tstart = tstart
        self._tstop = tstop
        self.subsetter = subsetter if subsetter else GribSubsetter(0, 360, -90, 90)

    @property
    def tstep(self):
        raise NotImplementedError("implement in child class")

    @property
    def times(self):
        t = np.arange(self._tstart, self._tstop, self.tstep)
        return t

    @abstractmethod
    def empty_file(self, *args, **kwargs):
        raise NotImplementedError("implement in child class")

    def create_variable(self, shortName, ncd, ncvar_dict):
        var = ncd.createVariable(shortName, "f8", self._dims, chunksizes=self.chunk_sizes(ncd))
        for attr in ncvar_dict:
            var.setncattr(attr, ncvar_dict[attr])

    def chunk_sizes(self, ncd):
        return None
    
    def valid_record(self, record):
        is_valid = record.shortName in self._include
        return is_valid
    
    def get_valid_vars(self, dataset) -> "list[str]":
        return [v for v in dataset.variables if v in self._include.keys()]

    def xread_grib(self, gribfile) -> "xr.Dataset":
        d = xr.open_dataset(gribfile, engine='cfgrib')
        return d

    def xgrib_to_nc(self, ncfile, gribfile):
        logger.debug(f"de-gribbing {Path(gribfile).name} with xarray")
        dataset = self.xread_grib(gribfile)
        for v in self.get_valid_vars(dataset):
            var = dataset[v]
            self.xvar_to_nc(var, ncfile)
    
    def xvar_to_nc(self, var, ncfile):
        if self.subsetter is not None:
            vals, lats, lons = self.subsetter.xsubset(var)
        else:
            vals, lats, lons = var.values, var.latitude.values, var.longitude.values

        # check if ncfile exists
        if not Path(ncfile).exists():
            logger.warning(f"missing file {ncfile}. Creating it")
            self.empty_file(ncfile, lats, lons, self.times)

        with nc.Dataset(ncfile, 'a') as ncd:
            # check if variable is in ncfile
            # if not, add it
            ncvar_dict = {"units" : var.units,
                          "long_name" : var.name}

            if var.name not in ncd.variables:
                self.create_variable(var.name, ncd, ncvar_dict)

            # check if timestep is in ncfile
            timestep = J3QD.date2num(pd.Timestamp(var.valid_time.values).to_pydatetime())
            try:
                t_i = np.where(self.times == timestep)[0][0]
            except IndexError:
                t_i = self.handle_missing_timestep(var.valid_time.values, ncfile)

            # write data
            if t_i is not None:
                self.xwrite_record(ncd=ncd,
                                   varname=var.name,
                                   t_i=t_i,
                                   vals=vals,
                                   var=var)

    def grib_to_nc(self, ncfile, record):
        varname = record.shortName

        if not self.valid_record(record):
            raise VariableError(f"{varname}")

        if self.subsetter is not None:
            vals, lats, lons = self.subsetter.subset(record)
        else:
            vals, lats, lons = record.data()

        # check if ncfile exists
        if not Path(ncfile).exists():
            logger.warning(f"missing file {ncfile}. Creating it")
            lats = lats[:,0]
            lons = lons[0,:]

            self.empty_file(ncfile, lats, lons, self.times)  # if not, create it

        with nc.Dataset(ncfile, 'a') as ncd:
            # check if variable is in ncfile
            # if not, add it
            ncvar_dict = {"units" : record.units,
                          "long_name" : record.name,
                          "comment": record.__repr__()}

            if varname not in ncd.variables:
                self.create_variable(record.shortName, ncd, ncvar_dict)

            # check if timestep is in ncfile
            timestep = J3QD.date2num(record.validDate)
            try:
                t_i = np.where(self.times == timestep)[0][0]
            except IndexError:
                t_i = self.handle_missing_timestep(record, ncfile)
            
            # write temporary array
            # self.write_mem(ncd, t_i, vals, record)

            # write data
            if t_i is not None:
                self.write_record(ncd=ncd,
                                  varname=varname,
                                  t_i=t_i,
                                  vals=vals,
                                  record=record)

    def handle_missing_timestep(self, date, ncfile):
        print(f"err: timestamp {date} not in ncfile {ncfile}")
        return None

    def write_record(self, ncd, varname, t_i, vals, record):
        ncd[varname][t_i, :, :] = vals

    def xwrite_record(self, ncd, varname, t_i, vals, var):
        ncd[varname][t_i, :, :] = vals

    def valid_messages(self):
        m = [d["messagenumber"] for d in list(self._include.values()) if "messagenumber" in d]
        return m


class ToConverter(ConversionHandler):
    _include = {"z": None}
    _dims = ("time", "latitude", "longitude",)

    @property
    def tstep(self):
        return 24

    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)

    def write_record(self, ncd, varname, t_i, vals, record):
        ncd[varname][0, :, :] = vals

    def xwrite_record(self, ncd, varname, t_i, vals, var):
        ncd[varname][0, :, :] = vals

    def handle_missing_timestep(self, date, ncfile):
        return 0


class SaConverter(ConversionHandler):
    _include = {"sp": None,
                "2t": None,
                "2sh": None,
                "2r": None,
                "10u": None,
                "10v": None}
    _dims = ("time", "latitude", "longitude",)

    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = J3QD.date2num(datetime(1980, 3, 1, 0))
        t2 = J3QD.date2num(datetime(1980, 3, 1, 6))
        return t2 - t1

    def chunk_sizes(self, ncd):
        return (ncd['time'].shape[0], CHUNK_LAT, CHUNK_LON,)
    
    def xread_grib(self, gribfile) -> xr.Dataset:
        d = xr.open_dataset(gribfile, engine='cfgrib', backend_kwargs={"filter_by_keys":{'typeOfLevel': 'surface'}})
        return d

    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)


class FcstConverter(ConversionHandler):
    _include = {"tprate": {"messagenumber": 6},
                "sp": {"messagenumber": 12},
                'dswrf': {"messagenumber": 13},
                "dlwrf": {"messagenumber": 17}
                }
    _dims = ("time", "latitude", "longitude",)
   
    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = J3QD.date2num(datetime(1980, 3, 1, 0))
        t2 = J3QD.date2num(datetime(1980, 3, 1, 1))
        return t2 - t1

    def valid_record(self, record):
        return record.messagenumber in self.valid_messages()
        
    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)

    def xread_grib(self, gribfile) -> xr.Dataset:
        d = xr.open_dataset(gribfile, engine='cfgrib', backend_kwargs={"filter_by_keys":{'typeOfLevel': 'surface'}})
        return d

    def chunk_sizes(self, ncd):
        return (ncd['time'].shape[0], CHUNK_LAT, CHUNK_LON,)


class PlConverter(ConversionHandler):
    _dims = ("time", "level", "latitude", "longitude",)
    _include = {"t":None,  # temperature
                "u":None,  # u-wind
                "v":None,  # v-wind
                "q":None,  # specific humidity
                "r":None,   # relative humidity
                "gh":None    # geopotential
                }

    def __init__(self, tstart, tstop, subsetter):
        super().__init__(tstart, tstop, subsetter)
        # self.mem_array = np.empty(len(self.times), len(self.subsetter.levels), len(self.subsetter.lats), len(self.subsetter.lons))

    def valid_record(self, record):
        good_variable = record.shortName in self._include
        lv_hPa = record.level if record.pressureUnits == "hPa" else record.level / 10
        good_level = lv_hPa in self.subsetter.levels

        return (good_variable and good_level)

    def xread_grib(self, gribfile) -> xr.Dataset:
        d = xr.open_dataset(gribfile, engine='cfgrib', backend_kwargs={"filter_by_keys":{'typeOfLevel': 'isobaricInhPa'}})
        return d
    
    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = J3QD.date2num(datetime(1980, 3, 1, 0))
        t2 = J3QD.date2num(datetime(1980, 3, 1, 6))
        return t2 - t1

    def empty_file(self, filename:str, lats, lons, times):
        levels = self.subsetter.levels
        empty_pressure_levels_file(filename, lats, lons, times, levels)

    def write_record(self, ncd, varname, t_i, vals, record):
        lev = record.level if record.pressureUnits == "hPa" else record.level / 10

        lev_i = np.where(ncd['level'][:] == lev)[0][0]

        ncd[varname][t_i, lev_i, :, :] = vals

    def xwrite_record(self, ncd, varname, t_i, vals, var):
        ncd[varname][t_i, :, :, :] = vals[::-1, :, :]  # put levels in backwards to match JRA55
    
    def chunk_sizes(self, ncd):
        return (ncd['time'].shape[0], ncd['level'].shape[0], CHUNK_LAT, CHUNK_LON,)


def empty_constant_file(filename: str, lats, lons):
    with nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC') as ncd:
        _ = ncd.createDimension('latitude', len(lats))
        _ = ncd.createDimension('longitude', len(lons))

        var_lat = ncd.createVariable('latitude', 'f8', ('latitude',))
        var_lat.units = "degrees_north"
        var_lon = ncd.createVariable('longitude', 'f8', ('longitude',))
        var_lon.units = "degrees_east"

        var_lat[:] = lats
        var_lon[:] = lons
        

def empty_surface_file(filename:str, lats, lons, times):
    with nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC') as ncd:
        _ = ncd.createDimension('time', len(times))
        _ = ncd.createDimension('latitude', len(lats))
        _ = ncd.createDimension('longitude', len(lons))

        var_time = ncd.createVariable('time', 'f8', ('time',))
        var_time.units = J3QD._tunits
        var_time.calendar = J3QD._tcal
        var_lat = ncd.createVariable('latitude', 'f8', ('latitude',))
        var_lat.units = "degrees_north"
        var_lon = ncd.createVariable('longitude', 'f8', ('longitude',))
        var_lon.units = "degrees_east"

        var_time[:] = times
        var_lat[:] = lats
        var_lon[:] = lons


def empty_pressure_levels_file(filename:str, lats, lons, times, levs):
    with nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC') as ncd:
        _ = ncd.createDimension('level', len(levs))
        _ = ncd.createDimension('time', len(times))
        _ = ncd.createDimension('latitude', len(lats))
        _ = ncd.createDimension('longitude', len(lons))

        var_time = ncd.createVariable('time', 'f8', ('time',))
        var_time.units = J3QD._tunits
        var_time.calendar = J3QD._tcal
        var_lat = ncd.createVariable('latitude', 'f8', ('latitude',))
        var_lat.units = "degrees_north"
        var_lon = ncd.createVariable('longitude', 'f8', ('longitude',))
        var_lon.units = "degrees_east"
        var_lev = ncd.createVariable('level', 'f8', ('level',))
        var_lev.units = "hPa"

        var_time[:] = times
        var_lat[:] = lats
        var_lon[:] = lons
        var_lev[:] = levs


def gribs_to_netcdf(*grib_files:"Union[str,Path]", netcdf_file:"Union[str,Path]", tstart:float, tstop:float, kind:str, subsetter=None, engine="xarray"):
    logger.info(f"creating {netcdf_file}")
    notified = []
    Converter = get_converter(kind)
    converter = Converter(tstart, tstop, subsetter)

    for gribfile in grib_files:
        if engine == "xarray":
            converter.xgrib_to_nc(netcdf_file, gribfile)

        elif engine == "pygrib":
            logger.debug(f"de-gribbing {Path(gribfile).name} with pygrib")

            with pygrib.open(str(gribfile)) as f:  # pygrib doesn't like Path objects
                for record in f:
                    if not converter.valid_record(record):
                        continue
                    try:
                        converter.grib_to_nc(netcdf_file, record)
                    except VariableError as e:
                        if e.args[0] not in notified:
                            logger.warning(f"warning: skip variable {e}")
                            notified.append(e.args[0])

        else:
            raise ValueError("engine must be 'pygrib' or 'xarray'")


def get_converter(kind:str) -> "type[ConversionHandler]":
    if kind == "sf":
        return FcstConverter
    elif kind == "sa":
        return SaConverter
    elif kind == "pl":
        return PlConverter
    elif kind == "to":
        return ToConverter
    else:
        raise ValueError(f"unknown kind {kind}")
