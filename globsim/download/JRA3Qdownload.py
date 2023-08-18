import netCDF4 as nc
import numpy as np
import pygrib

from pathlib import Path
from abc import ABC, abstractmethod
from globsim.download.GenericDownload import GenericDownload
from datetime import datetime, timedelta


class J3QD(GenericDownload):

    _tunits = "hours since 1947-01-01"
    _tcal = "standard"

    def __init__(self):
        pass


class VariableError(ValueError):
    pass


class ConversionHandler(ABC):
    """
    Handles the translation of grib files into a netcdf file
    """
    _include = []
    _dims = ()

    def __init__(self, tstart: int, tstop: int, subsetter=None):
        self._tstart = tstart
        self._tstop = tstop
        self.subsetter = subsetter

    @property
    def tstep(self):
        raise NotImplementedError("implement in child class")

    @property
    def times(self):
        t = np.arange(self._tstart, self._tstop + self.tstep, self.tstep)
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
            print(f"missing file {ncfile}. Creating it")
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
            timestep = nc.date2num(record.validDate, J3QD._tunits, J3QD._tcal)
            try:
                t_i = np.where(self.times == timestep)[0][0]
            except IndexError:
                print(f"err: timestamp {record.validDate} not in ncfile")
                return
            
            # write data
            self.write_record(ncd=ncd,
                              varname=varname,
                              t_i=t_i,
                              vals=vals,
                              record=record)

    def write_record(self, ncd, varname, t_i, vals, record):
        ncd[varname][t_i, :, :] = vals


class SaConverter(ConversionHandler):
    _include = ["sp", "2t", "2sh", "2r", "10u", "10v"]
    _dims = ("time", "latitude", "longitude",)

    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = nc.date2num(datetime(1980, 3, 1, 0), J3QD._tunits, J3QD._tcal)
        t2 = nc.date2num(datetime(1980, 3, 1, 6), J3QD._tunits, J3QD._tcal)
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
        t1 = nc.date2num(datetime(1980, 3, 1, 0), J3QD._tunits, J3QD._tcal)
        t2 = nc.date2num(datetime(1980, 3, 1, 1), J3QD._tunits, J3QD._tcal)
        return t2 - t1
    
    def valid_record(self, record):
        return record.messagenumber in self._include.keys()
    
    def empty_file(self, filename:str, lats, lons, times):
        empty_surface_file(filename, lats, lons, times)


class PlConverter(ConversionHandler):
    _dims = ("time","level", "latitude", "longitude",)
    _include = ["t"]#, "u", "v"]

    def valid_record(self, record):
        good_variable = record.shortName in self._include
        lv_hPa = record.level if record.pressureUnits == "hPa" else record.level / 10
        good_level = lv_hPa in self.subsetter.levels
        
        return (good_variable and good_level)

    @property
    def tstep(self):
        # 6-hourly, calculated on the fly in case we change time units in the future
        t1 = nc.date2num(datetime(1980, 3, 1, 0), J3QD._tunits, J3QD._tcal)
        t2 = nc.date2num(datetime(1980, 3, 1, 6), J3QD._tunits, J3QD._tcal)
        return t2 - t1
    
    def empty_file(self, filename:str, lats, lons, times):
        levels = self.subsetter.levels
        empty_pressure_levels_file(filename, lats, lons, times, levels)

    def write_record(self, ncd, varname, t_i, vals, record):
        lev = record.level if record.pressureUnits == "hPa" else record.level / 10

        lev_i = np.where(ncd['level'][:] == lev)[0][0]
        
        print(f"writing record at t[{t_i}] and lev[{lev_i}]")
        ncd[varname][t_i, lev_i, :, :] = vals


def empty_surface_file(filename:str, lats, lons, times):
    with nc.Dataset(filename, 'w', format='NETCDF4') as ncd:
        dim_time = ncd.createDimension('time', len(times))
        dim_lat = ncd.createDimension('latitude', len(lats))
        dim_lon = ncd.createDimension('longitude', len(lons))

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
        dim_lev = ncd.createDimension('level', len(levs))
        dim_time = ncd.createDimension('time', len(times))
        dim_lat = ncd.createDimension('latitude', len(lats))
        dim_lon = ncd.createDimension('longitude', len(lons))

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


def gribs_to_netcdf(*grib_files, netcdf_file, tstart, tstop, kind, subsetter=None):
    Converter = get_converter(kind)
    converter = Converter(tstart, tstop, subsetter)

    for grib_file in grib_files:
        with pygrib.open(grib_file) as f:
            for record in f:
                try:
                    converter.grib_to_nc(netcdf_file, record, subsetter=subsetter)
                except VariableError as e:
                    print(f"debug: skip variable {e}")


def get_converter(kind:str) -> "type[ConversionHandler]":
    if kind == "sa":
        return SaConverter
    elif kind == "sf":
        return FcstConverter
    elif kind == "pl":
        return PlConverter
    else:
        raise ValueError(f"unknown kind {kind}")

