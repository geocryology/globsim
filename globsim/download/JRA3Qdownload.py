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
from globsim.download.JRAdownload import JRAdownload
from globsim.download.jra_dict_formatters import getPressureLevels
from globsim.download.RDA import Rdams 
from globsim.download.jra_dict_formatters import J3QDictFormatter
from globsim.common_utils import variables_skip

logger = logging.getLogger('globsim.download')


class J3QD(JRAdownload):

    JRA_VERSION = 'jra3q'
    API = Rdams
    DICT_FORMATTER = J3QDictFormatter
    dsID = 'd640000'
    timeName = 'time'

    def getVars(self, dsi):
        '''get all variables in tar file'''

        flist = list(Path(self.directory).glob(f'*{dsi}*'))
        varlist = []
        for f in flist:
            fname = f.name
            for part in fname.split("."):
                if len(part.split("-")) == 4:
                    var = part.split("-")[0]
                    varlist.append(var)

        variables = np.unique(varlist)
        return variables
    
    def get_nc_files(self, variable, dsi):
        return np.sort(list(Path(self.directory).glob(f"*{dsi}*{variable}*.nc")))
    
    def ncfVar(self, param):
        varname = param.split('-')[0]
        d = {'rh2m': 'Relative humidity',
            'tmp2m': 'Temperature',
            'ugrd10m': 'u-component of wind',
            'vgrd10m': 'v-component of wind',
            'lat': 'latitude',
            'lon': 'longitude',
            'time': 'time',
            'hgt': 'Geopotential height',
            'rh': 'Relative humidity',
            'tmp': 'Temperature',
            'ugrd': 'u-component of wind',
            'vgrd': 'v-component of wind',}
        return d.get(varname, 'skip')

    
    def makeNCF(self, dsi):        
        self.extract_downloaded_tar_files(dsi)
        variables = self.getVars(dsi)
        dataLev = self.getDataLev(dsi)
        self.getDimName(dataLev)
        
        nc_template_files = self.get_nc_files(variables[0], dsi)
        
        ncf = xr.open_mfdataset(nc_template_files.tolist(), decode_times=False)

        if dataLev == 'pl':
            Levs = ncf['pressure_level']
        else: 
            Levs = None

        Times = ncf['time']
        Lats  = ncf['lat']
        Lons  = ncf['lon']

        ncf.close()

        file_new = self.getOutFile(ncf, dataLev)
        ncn = new_jra_download_file(file_new, Times, Lats, Lons, Levs)
        
        # assign dimensions
        ncn['time'][:] = Times
        ncn['longitude'][:] = Lons
        ncn['latitude'][:] = Lats
        '''
        self.formatter = self.DICT_FORMATTER(date=4,
                                    area=self.area,
                                    elevation=self.elevation,
                                    variables=self.variables)
        '''
        import pdb;pdb.set_trace()
        for vari in variables:
            flist = self.get_nc_files(vari, dsi)
            ncf = xr.open_mfdataset(flist.tolist(), decode_times=False)
            for n, var in enumerate(ncf.variables.keys()):
                if variables_skip(self.ncfVar(var)):
                    print("Skipping variable: ", var)
                    continue
                logger.info(f"Creating variable: {var}")
                if dataLev == 'pl':
                    vari = ncn.createVariable(self.ncfVar(var), 'f4',
                                              ('time', 'level',
                                               'latitude', 'longitude',))
                    vari[:,:,:,:] = ncf[var][:,:,:,:]
                else:
                    vari = ncn.createVariable(self.ncfVar(var),'f4',
                                              ('time',
                                               'latitude', 'longitude'))
                    vari[:,:,:] = ncf[var][:,:,:]
                vari.long_name = ncf[var].long_name
                vari.jra_name = self.formatter.lookup_param(var)
                vari.units     = ncf[var].units

            ncf.close()
            for f in flist:
                print("TODO: UNLINK FILE ")

        ncn.close()


# initialize new data file and create group
def new_jra_download_file(filename, Times, Lats, Lons, Levs=None) -> "nc.Dataset":
    ncn = nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC')

    # make dimensions
    if Levs is not None:
        ncn.createDimension('level', len(Levs))
        levels     = ncn.createVariable('level', 'i4',('level',))
        levels.long_name  = 'pressure level'
        levels.units      = 'mbar'
        levels[:] = Levs

    ncn.createDimension('time', len(Times))
    ncn.createDimension('latitude', len(Lats))
    ncn.createDimension('longitude', len(Lons))

    # make dimension variables
    times      = ncn.createVariable('time', 'd',('time',))
    latitudes  = ncn.createVariable('latitude', 'f8', ('latitude',))
    longitudes = ncn.createVariable('longitude', 'f8', ('longitude',))

    times.standard_name = 'time'
    times.units     = Times.units
    times.calendar  = 'standard'
    latitudes.standard_name  = Lats.long_name
    latitudes.units      = Lats.units
    longitudes.standard_name = Lons.long_name
    longitudes.units     = Lons.units

    return ncn
   
class JraDownloadHandler:

    def __init__(self):
        