import netCDF4 as nc
import numpy as np
import pandas as pd
import pygrib
import logging
from urllib.error import URLError
from pathlib import Path
from typing import Optional, Union

from abc import ABC, abstractmethod
from globsim.download.GenericDownload import GenericDownload
from globsim.download.JRA3Q_helpers import GribSubsetter, download_daily_gribs, download_constant
from globsim.download.JRAdownload import getDate
from globsim.download.JRA3Q_dl import GetAccessor
from datetime import datetime

logger = logging.getLogger('globsim.download')


class J3QD(GenericDownload):

    _tunits = "hours since 1947-01-01"
    _tcal = "standard"

    def __init__(self, pfile):
        super().__init__(pfile)
        self.retry_delay_min = 1
        par = self.par
        self._set_input_directory("jra3q")
        self._subsetter = GribSubsetter(round((float(par['bbW']) - 1.25) / 1.25) * 1.25,
                                        round(float(par['bbE']) / 1.25) * 1.25,
                                        round((float(par['bbS']) - 1.25) / 1.25) * 1.25,
                                        round(float(par['bbN']) / 1.25) * 1.25)

        # time bounds
        self.date  = getDate(par)

        # Connect
        self.credential = Path(par['credentials_directory'], ".netrc")
        self.connect()

        # chunk size for downloading and storing data [days]
        logging.info("Chunk size ignored for JRA3Q download. Using monthly")
        # self.chunk_size = par['chunk_size'] * 2000

        # variables
        logging.info("Variable list ignored for JRA3Q")

    def connect(self):
        try:
            self.access = GetAccessor(netrc_file=self.credential)
        except ValueError:
            logging.error("Could not find credentials. Aborting")
            raise ValueError(f"Could not find credentials file {self.credential}")
        except URLError:
            logging.error("Could not connect to JRA server. This aint gonna work.")

    @staticmethod
    def date2num(dt:datetime) -> float:
        return nc.date2num(dt, J3QD._tunits, J3QD._tcal)

    def nc_filename(self, year, month, kind="sa") -> str:
        return f"jra3q_{kind}_{year}{month}01_to_{year}{month}31.nc"

    def retrieve(self):
        days = pd.date_range(self.date['beg'], self.date['end'])

        # Constant surface data
        to_grib = [download_constant(self.access, self.directory)]
        t = days[0]
        to_ncfile = Path(self.directory, f'jra3q_to_{t.strftime(r"%Y%m%d")}_to_{t.strftime(r"%Y%m%d")}.nc')
        gribs_to_netcdf(*to_grib,
                        netcdf_file=to_ncfile,
                        tstart=self.date2num(t),
                        tstop=self.date2num(days[1]),
                        kind="to",
                        subsetter=self._subsetter)

        # Temporal data
        running_sa, running_sf, running_pl, running_dates = [], [], [], []

        for date in days:
            # Aggregation subroutine
            if date.day == 1 and len(running_dates) > 1:
                t1 = running_dates[0]
                t2 = running_dates[-1]
                tstart = J3QD.date2num(t1)
                tstop = J3QD.date2num(date)  # t2 + timedelta(days=1)

                sa_ncfile = Path(self.output_directory, f'jra3q_sa_{t1.strftime(r"%Y%m%d")}_to_{t2.strftime(r"%Y%m%d")}.nc')
                sf_ncfile = Path(self.output_directory, f'jra3q_sf_{t1.strftime(r"%Y%m%d")}_to_{t2.strftime(r"%Y%m%d")}.nc')
                pl_ncfile = Path(self.output_directory, f'jra3q_pl_{t1.strftime(r"%Y%m%d")}_to_{t2.strftime(r"%Y%m%d")}.nc')

                # aggregate gribs then delete them
                logger.info("Aggregating gribs to netcdf")
                gribs_to_netcdf(*running_sa, netcdf_file=sa_ncfile, tstart=tstart, tstop=tstop, kind="sa", subsetter=self._subsetter)
                gribs_to_netcdf(*running_sf, netcdf_file=sf_ncfile, tstart=tstart, tstop=tstop, kind="sf", subsetter=self._subsetter)
                gribs_to_netcdf(*running_pl, netcdf_file=pl_ncfile, tstart=tstart, tstop=tstop, kind="pl", subsetter=self._subsetter)

                #  delete grib files
                logger.info("Not deleting grib files")
                #for grib in running_sa + running_sf + running_pl:
                    #Path(grib).unlink()

                # clear running
                running_sa, running_sf, running_pl, running_dates = [], [], [], []

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

    def create_variable(self, record, ncd, ncvar_dict):
        varname = record.shortName
        var = ncd.createVariable(varname, "f8", self._dims)
        for attr in ncvar_dict:
            var.setncattr(attr, ncvar_dict[attr])

    def valid_record(self, record):
        is_valid = record.shortName in self._include
        return is_valid

    def grib_to_nc(self, ncfile, record, subsetter=None):
        varname = record.shortName

        subsetter = self.subsetter if subsetter is None else subsetter

        if not self.valid_record(record):
            raise VariableError(f"{varname}")

        if subsetter is not None:
            vals, lats, lons = subsetter.subset(record)
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
                self.create_variable(record, ncd, ncvar_dict)

            # check if timestep is in ncfile
            timestep = J3QD.date2num(record.validDate)
            try:
                t_i = np.where(self.times == timestep)[0][0]
            except IndexError:
                t_i = self.handle_missing_timestep(record, ncfile)
            
            # write temporary array 
            self.write_mem(ncd, t_i, vals, record)

            # write data
            if t_i is not None:
                self.write_record(ncd=ncd,
                                  varname=varname,
                                  t_i=t_i,
                                  vals=vals,
                                  record=record)

    def handle_missing_timestep(self, record, ncfile):
        print(f"err: timestamp {record.validDate} not in ncfile {ncfile}")
        return None

    def write_record(self, ncd, varname, t_i, vals, record):
        ncd[varname][t_i, :, :] = vals


class ToConverter(ConversionHandler):
    _include = ["z"]
    _dims = ("time", "latitude", "longitude",)

    @property
    def tstep(self):
        return 24

    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)

    def write_record(self, ncd, varname, t_i, vals, record):
        ncd[varname][0, :, :] = vals

    def handle_missing_timestep(self, record, ncfile):
        return 0


class SaConverter(ConversionHandler):
    _include = ["sp", "2t", "2sh", "2r", "10u", "10v"]
    _dims = ("time", "latitude", "longitude",)

    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = J3QD.date2num(datetime(1980, 3, 1, 0))
        t2 = J3QD.date2num(datetime(1980, 3, 1, 6))
        return t2 - t1

    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)


class FcstConverter(ConversionHandler):
    _include = {6: "tprate",
                12: "sp",
                13: 'dswrf',
                17: "dlwrf"}
    _dims = ("time", "latitude", "longitude",)

    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = J3QD.date2num(datetime(1980, 3, 1, 0))
        t2 = J3QD.date2num(datetime(1980, 3, 1, 1))
        return t2 - t1

    def valid_record(self, record):
        return record.messagenumber in self._include.keys()

    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)


class PlConverter(ConversionHandler):
    _dims = ("time","level", "latitude", "longitude",)
    _include = ["t", "u", "v"]

    def __init__(self):
        super().__init__()
        # self.mem_array = np.empty((len(self.times), len(self.subsetter.levels), len(self.subsetter.lats), len(self.subsetter.lons))
                                  
    def valid_record(self, record):
        good_variable = record.shortName in self._include
        lv_hPa = record.level if record.pressureUnits == "hPa" else record.level / 10
        good_level = lv_hPa in self.subsetter.levels

        return (good_variable and good_level)

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

        print(f"writing record at t[{t_i}] and lev[{lev_i}]")
        ncd[varname][t_i, lev_i, :, :] = vals


def empty_constant_file(filename: str, lats, lons):
    with nc.Dataset(filename, 'w', format='NETCDF4') as ncd:
        _ = ncd.createDimension('latitude', len(lats))
        _ = ncd.createDimension('longitude', len(lons))

        var_lat = ncd.createVariable('latitude', 'f8', ('latitude',))
        var_lat.units = "degrees_north"
        var_lon = ncd.createVariable('longitude', 'f8', ('longitude',))
        var_lon.units = "degrees_east"

        var_lat[:] = lats
        var_lon[:] = lons
        

def empty_surface_file(filename:str, lats, lons, times):
    with nc.Dataset(filename, 'w', format='NETCDF4') as ncd:
        _ = ncd.createDimension('time', len(times))
        _ = ncd.createDimension('latitude', len(lats))
        _ = ncd.createDimension('longitude', len(lons))

        var_time = ncd.createVariable('time', 'i8', ('time',))
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
    with nc.Dataset(filename, 'w', format='NETCDF4') as ncd:
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


def gribs_to_netcdf(*grib_files:"Union[str,Path]", netcdf_file:"Union[str,Path]", tstart:float, tstop:float, kind:str, subsetter=None):
    notified = []
    Converter = get_converter(kind)
    converter = Converter(tstart, tstop, subsetter)

    for grib_file in grib_files:
        with pygrib.open(str(grib_file)) as f:  # pygrib doesn't like Path objects
            for record in f:
                try:
                    converter.grib_to_nc(netcdf_file, record, subsetter=subsetter)
                except VariableError as e:
                    if e.args[0] not in notified:
                        logger.warning(f"warning: skip variable {e}")
                        notified.append(e.args[0])


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
