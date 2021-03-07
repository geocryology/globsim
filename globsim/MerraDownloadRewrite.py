#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright Xiaojing Quan & Stephan Gruber
# =============================================================================
# REVISION HISTORY
# 20170510 -- Initial Version Created
# 20171208 -- First Draft Completed
#
# ==============================================================================
# A scripts for downloading MERRA-2 reanalysis data:
# -- Geoppotential Height at pressure levels [time*level*lat*lon] (Unit:m) (6-hourly/day)
# -- Air Temperature at pressure levels [time*level*lat*lon] (Unit:K) (6-hourly/day)
# -- Relative Humidity at Pressure Levels[time*level*lat*lon] (Unit:1) (3-hourly/day)
# -- Easteward Wind at Pressure Levels [time*level*lat*lon] (Unit:m/s) (6-hourly/day)
# -- Northward Wind at Pressure Levels [time*level*lat*lon] (Unit:m/s) (6-hourly/day)
# -- Air Temperature at 2 Meter [time*lat*lon] (Unit:K) (1-hourly/day)
# -- Eastward Wind at 2 Meter [time*lat*lon] (Unit:K) (1-hourly/day)
# -- Northward Wind at 2 Meter [time*lat*lon] (Unit: m/s) (1-hourly/day)
# -- Eastward Wind at 10 Meter  [time*lat*lon] (Unit: m/s) (1-hourly/day)
# -- Northward Wind at 10 Meter [time*lat*lon] (Unit: m/s) (1-hourly/day)
# -- Precipitation Flux [time*lat*lon] (Unit: kg/m2/s) (1-hourly/day)
# -- Surface Incoming Shoertwave Flux [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Incoming Shortwave Flux Assuming Clear Sky [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Net Downward Longwave Flux [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Net Downward Longwave Flux Assuming Clear Sky [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Longwave Flux Emitted from Surface [time*lat*lon] (Unit:W/m2) (1-hourly/day)
#
# -- Surface Absorbed Longwave Flux [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Absorbed Longwave Flux Assuming Clear Sky [time*lat*lon] (Unit:W/m2) (1-hourly/day)
#
# Saved as netCDF 4
# ====================HOW TO RUN THIS ==========================================
#
# (1) Register a New User in Earthdata Login:
#  https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login
#
# (2) Authorize NASA GESDISC DATA ARCHIVE in Earthdata Login
# https://disc.gsfc.nasa.gov/registration/authorizing-gesdisc-data-access-in-earthdata_login
#
# (3) Adapt the script below with:
#    - Authrized Username and Password (setup in .merrarc file),
#    - Input parameters: Date, Area, Elevation, Chunk_size, Variables, etc.
#      (setup in Globsim download parameter file )
#
# (4) Obtaining the URL addresses of the objected datasets at:
#     https://disc.sci.gsfc.nasa.gov/daac-bin/FTPSubset2.pl
#
# (5) Obtianing the mutiple datasets with spefici spacial and temporal)
#
# (6) Get all varialbes which are needed, and saved in NetCDF files
#
# ==============================================================================
# IMPORTANT Notes:

# 1. Samples of Selected URLs list:

# 3d,6-hourly,Instantaneous,Pressure-Level, Analyzed Meteorological Fields
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4')

# 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4/2016/01/MERRA2_400.inst3_3d_asm_Np.20160101.nc4')

# 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')

# 2d,1-hourly, single-level, full horizontal resolution, Surface Flux Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXFLX.5.12.4/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4')

# 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXRAD.5.12.4/1981/01/MERRA2_100.tavg1_2d_rad_Nx.19810101.nc4')

# 2d 1-Hourly,Time-Averaged,Single-Level,Assimilation,Single-Level Diagnostics V5.12.4
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXSLV.5.12.4/2016/01/MERRA2_400.tavg1_2d_slv_Nx.20160101.nc4')

# 2. Radiation Variables Processing:
# Considering the variables below as equal utilization:
#
# downwelling_shortwave_flux_in_air_assuming_clear_sky = Surface_Incoming_Shoertwave_Flux
#
# downwelling_shortwave_flux_in_air_assuming_clear_sky = Surface_Incoming_Shortwave_Flux_Assuming_Clear_Sky
#
# downwelling_longwave_flux_in_air = Longwave Flux Emitted from Surface + Surface Net Downward Longwave Flux
#
# downwelling_longwave_flux_in_air_assuming_clear_sky = Longwave Flux Emitted from Surface + Surface Net Downward Longwave Flux Assuming Clear Sky

# ==============================================================================

from __future__        import print_function

import numpy as np
import netCDF4 as nc
import itertools
import pandas as pd
import time as tc
import sys

from pydap.client      import open_url, open_dods
from pydap.cas.urs     import setup_session
from datetime          import datetime, timedelta
from os                import path, makedirs
from netCDF4           import Dataset
from math              import floor, atan2, pi

from scipy.interpolate import interp1d, griddata
from time              import sleep

from globsim.common_utils import str_encode, series_interpolate, variables_skip, GenericDownload, GenericScale, GenericInterpolate
from globsim.meteorology import LW_downward, pressure_from_elevation
from globsim.nc_elements import netcdf_base, new_scaled_netcdf

import warnings

warnings.filterwarnings("ignore", category=UserWarning, module='urllib2')


class MERRAgeneric():
    """
    Parent class for other merra classes.
    """

    def netCDF_empty(self, ncfile_out, stations, nc_in):
        # TODO: change date type from f4 to f8 for lat and lon
        '''
        Creates an empty station file to hold interpolated reults. The number of
        stations is defined by the variable stations, variables are determined by
        the variable list passed from the gridded original netCDF.

        ncfile_out: full name of the file to be created
        stations:   station list read with common_utils.StationListRead()
        variables:  variables read from netCDF handle
        lev:        list of pressure levels, empty is [] (default)
        '''
        rootgrp = netcdf_base(nc_in, ncfile_out, len(stations), None,
                              'hours since 1980-01-01 00:00:00')

        station = rootgrp["station"]
        latitude = rootgrp["latitude"]
        longitude = rootgrp["longitude"]
        height = rootgrp["height"]

        # assign station characteristics
        station[:]   = list(stations['station_number'])
        latitude[:]  = list(stations['latitude_dd'])
        longitude[:] = list(stations['longitude_dd'])
        height[:]    = list(stations['elevation_m'])

        # extra treatment for pressure level files
        try:
            lev = nc_in.variables['level'][:]
            print("== 3D: file has pressure levels")
            level           = rootgrp.createDimension('level', len(lev))
            level           = rootgrp.createVariable('level','i4',('level'))
            level.long_name = 'pressure_level'
            level.units     = 'hPa'
            level[:] = lev
        except Exception:
            print("== 2D: file without pressure levels")
            lev = []

        # remove extra variables
        varlist_merra = [str_encode(x) for x in nc_in.variables.keys()]
        self.MERRA_skip(varlist_merra)

        # create and assign variables based on input file
        for n, var in enumerate(varlist_merra):
            if variables_skip(var):
                continue
            print("VAR: ", str_encode(var))
            # extra treatment for pressure level files
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4', ('time', 'level', 'station'))
            else:
                tmp = rootgrp.createVariable(var,'f4', ('time', 'station'))
            tmp.long_name = str_encode(nc_in.variables[var].long_name)  # for merra2
            tmp.units     = str_encode(nc_in.variables[var].units)

        # close the file
        rootgrp.close()


class SaveNCDF_pl_3dm():
    """ write output netCDF file for abstracted variables from original
    meteorological data at pressure levels
    """

    def safety_checks(self, data_3dmasm, data_3dmana):
        assert(np.array_equal(data_3dmasm[0].lat.data, data_3dmana[0].lat.data))
        assert(np.array_equal(data_3dmasm[0].lev.data, data_3dmana[0].lev.data))
        assert(np.array_equal(data_3dmasm[0].lon.data, data_3dmana[0].lon.data))

    def saveData(self, data_3dmasm, data_3dmana, chunk_size, dir_data, elevation):
        """
        build NetCDF file for saving output variables (T, U, V, H, RH, lat, lon, levels, time)
        """
        self.safety_checks(data_3dmasm, data_3dmana)

        # set up File
        file_ncdf  = path.join(dir_data, f"merra_pl_{data_3dmana[0].time.begin_date}_to_{data_3dmana[-1].time.begin_date}.nc")
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        rootgrp.source      = 'MERRA-2, meteorological variables from metadata at pressure levels'

        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('level', len(data_3dmana[0].lev))
        _  = rootgrp.createDimension('lats', len(data_3dmana[0].lat))
        _  = rootgrp.createDimension('lons', len(data_3dmana[0].lon))

        # Create coordinate variables
        Latitudes               = rootgrp.createVariable('latitude', 'f4',('lats'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = data_3dmana[0].lat.data    # pass the values of latitude

        Longitudes               = rootgrp.createVariable('longitude', 'f4',('lons'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = data_3dmana[0].lon.data

        Level = rootgrp.createVariable('level','i4', ('level'))
        Level.standard_name = "air_pressure"
        Level.long_name = "vertical level"
        Level.units = "hPa"
        Level.positive = "down"
        Level.axis = "Z"
        Level[:] = data_3dmana[0].lev.data

        Time  = rootgrp.createVariable('time', 'i4', ('time'))
        Time.standard_name = "time"
        Time.units  = "hours since 1980-1-1 00:00:00"
        Time.calendar      = "gregorian"

        netCDFTime = []
        for x in data_3dmana:
            time_datetime = nc.num2date(x.time.data, units=x.time.units, calendar=Time.calendar)  # Convert units
            time_nc = nc.date2num(time_datetime, units=Time.units, calendar=Time.calendar)
            netCDFTime.extend(time_nc)

        Time[:] = netCDFTime

        data_variables_3dmasm = [v for v in list(data_3dmasm[0]) if v not in ['time', 'lat', 'lon', 'lev']]
        data_variables_3dmana = [v for v in list(data_3dmana[0]) if v not in ['time', 'lat', 'lon', 'lev']]

        # Treat 6-hourly data
        for x in data_variables_3dmana:
            out_var = rootgrp.createVariable(x, 'f4', ('time', 'level','lats','lons'), fill_value=9.9999999E14)
            out_var.long_name = data_3dmana[0][x].long_name
            out_var.units         = data_3dmana[0][x].units
            out_var.missing_value = data_3dmana[0][x].missing_value

            # stack data arrays
            all_data = np.concatenate([dataset[x].data[0] for dataset in data_3dmana], axis=0)

            # TODO: [NB] What's going on with the extrapolation here? Is this necessary?
            if x in ["U", "V"]:
                all_data = MERRAdownload.wind_rh_Extrapolate(all_data)

            elif x in ["T"]:
                h_data = np.concatenate([dataset["H"].data[0] for dataset in data_3dmana], axis=0)  # could also get from rootgrp["H"]
                all_data = MERRAdownload.tempExtrapolate(t_total=all_data, h_total=h_data, elevation=elevation)

            out_var[:] = all_data

        # Treat 3-hourly data
        for x in data_variables_3dmasm:
            out_var = rootgrp.createVariable(x, 'f4', ('time', 'level','lats','lons'), fill_value=9.9999999E14)
            out_var.long_name = data_3dmasm[0][x].long_name
            out_var.units         = data_3dmasm[0][x].units
            out_var.missing_value = data_3dmasm[0][x].missing_value

            # stack data arrays
            all_data = np.concatenate([dataset[x].data[0] for dataset in data_3dmasm], axis=0)

            # TODO: [NB] What's going on with the extrapolation here? Is this necessary?
            if x in ["RH"]:
                all_data = MERRAdownload.wind_rh_Extrapolate(all_data)

            out_var[:] = all_data[::2,:,:,:]  # downsample to 6-hourly

        rootgrp.close()


class SaveNCDF_sa():
    """ write output netCDF file for abstracted variables from original 2D
        meteorological Diagnostics dataset and suface flux Diagnostics datasets
    """
    def saveData(self, data_2dm, dir_data):
        """build a NetCDF file for saving output variables
        (T2M,U2M, V2M, U10M,V10M, PRECTOT, PRECTOTCORR,T2MDEW,QV2M,
        lat, lon, time)

        Parameters
        ----------
        data_2m : list
            list of DatasetTypes retrieved from
        dir_data : str
            Path to directory in which files will be saved
        """

        # get the varable list and time indices

        # set up File
        file_ncdf  = path.join(dir_data, f"merra_sa_{data_2dm[0].time.begin_date}_to_{data_2dm[-1].time.begin_date}.nc")
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        rootgrp.source      = 'MERRA-2, meteorological variables from metadata at surface level'

        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('lats', len(data_2dm[0].lat))
        _  = rootgrp.createDimension('lons', len(data_2dm[0].lon))

        # Create coordinate variables
        Latitudes               = rootgrp.createVariable('latitude', 'f4',('lats'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = data_2dm[0].lat.data    # pass the values of latitude

        Longitudes               = rootgrp.createVariable('longitude', 'f4',('lons'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = data_2dm[0].lon.data

        data_variables = [v for v in list(data_2dm[0]) if v not in ['time', 'lat', 'lon']]

        # Fill in data variables
        for x in data_variables:
            out_var = rootgrp.createVariable(x, 'f4', ('time','lats','lons'), fill_value=9.9999999E14)
            out_var.long_name = data_2dm[0][x].long_name
            out_var.units         = data_2dm[0][x].units
            out_var.missing_value = data_2dm[0][x].missing_value

            # stack data arrays
            all_data = np.concatenate([dataset[x].data[0] for dataset in data_2dm], axis=0)
            out_var[:] = all_data

        # Fill in the time
        Time  = rootgrp.createVariable('time', 'i4', ('time'))
        Time.standard_name = "time"
        Time.units  = "hours since 1980-1-1 00:00:00"
        Time.calendar      = "gregorian"

        netCDFTime = []
        for x in data_2dm:
            time_datetime = nc.num2date(x.time.data, units=x.time.units, calendar=Time.calendar)  # Convert units
            time_nc = nc.date2num(time_datetime, units=Time.units, calendar=Time.calendar)
            netCDFTime.extend(time_nc)

        Time[:] = netCDFTime

        rootgrp.close()


class SaveNCDF_sf():
    """ write output netCDF file for abstracted variables from original 2D
        radiation Diagnostics datasets  datasets
    """
    def saveData(self, data_2dr, data_2ds, data_2dv, dir_data):
        """Build a NetCDF file for saving output variables (SWGDN, SWGDNCLR,
           LWRNT, LWGEM, LWGNTCLR, LWGAB, LWGABCLR, lat, lon, time)
           """
        self.safety_checks(data_2dr, data_2ds, data_2dv)

        # set up File
        file_ncdf  = path.join(dir_data, f"merra_sf_{data_2dr[0].time.begin_date}_to_{data_2dr[-1].time.begin_date}.nc")
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        rootgrp.source      = 'MERRA-2, radiation variables from metadata at surface level'

        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('lats', len(data_2dr[0].lat))
        _  = rootgrp.createDimension('lons', len(data_2dr[0].lon))

        # Create coordinate variables
        Latitudes               = rootgrp.createVariable('latitude', 'f4',('lats'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = data_2dr[0].lat.data    # pass the values of latitude

        Longitudes               = rootgrp.createVariable('longitude', 'f4',('lons'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = data_2dr[0].lon.data

        data_variables_2dr = [v for v in list(data_2dr[0]) if v not in ['time', 'lat', 'lon']]
        data_variables_2ds = [v for v in list(data_2ds[0]) if v not in ['time', 'lat', 'lon']]
        data_variables_2dv = [v for v in list(data_2dv[0]) if v not in ['time', 'lat', 'lon']]

        # Fill in data variables from the three sources
        for source, namelist in zip([data_2dr, data_2ds, data_2dv],
                                    [data_variables_2dr, data_variables_2ds, data_variables_2dv]):

            for x in namelist:
                out_var = rootgrp.createVariable(x, 'f4', ('time','lats','lons'), fill_value=9.9999999E14)
                out_var.long_name = source[0][x].long_name
                out_var.units         = source[0][x].units
                out_var.missing_value = source[0][x].missing_value

                # stack data arrays
                all_data = np.concatenate([dataset[x].data[0] for dataset in source], axis=0)
                out_var[:] = all_data

        # Fill in the time
        Time  = rootgrp.createVariable('time', 'i4', ('time'))
        Time.standard_name = "time"
        Time.units  = "hours since 1980-1-1 00:00:00"
        Time.calendar      = "gregorian"

        netCDFTime = []
        for x in data_2dr:
            time_datetime = nc.num2date(x.time.data, units=x.time.units, calendar=Time.calendar)  # Convert units
            time_nc = nc.date2num(time_datetime, units=Time.units, calendar=Time.calendar)
            netCDFTime.extend(time_nc)

        Time[:] = netCDFTime

        # Calculate variables
        all_variables = list(rootgrp.variables)

        if 'LWGEM' in all_variables and 'LWGNT' in all_variables:
            self.derive_lwgdn(rootgrp)

        if 'LWGEM' in all_variables and 'LWGNTCLR' in all_variables:
            self.derive_lwgdnclr(rootgrp)

        rootgrp.close()

    def derive_lwgdn(self, rootgrp):
        out_var = rootgrp.createVariable("LWGDN", 'f4', ('time','lats','lons'), fill_value=9.9999999E14)
        out_var.long_name = "downwelling_longwave_flux_in_air"
        out_var.units = 'W/m2'
        out_var[:] = rootgrp['LWGEM'][:].data + rootgrp['LWGNT'][:].data

    def derive_lwgdnclr(self, rootgrp):
        out_var = rootgrp.createVariable("LWGDNCLR", 'f4', ('time','lats','lons'), fill_value=9.9999999E14)
        out_var.long_name = "downwelling_longwave_flux_in_air_assuming_clear_sky"
        out_var.units = 'W/m2'
        out_var[:] = rootgrp['LWGEM'][:].data + rootgrp['LWGNTCLR'][:].data

    def safety_checks(self, data_2dr, data_2ds, data_2dv):
        """ ensure assumptions are true of datasets """
        assert(np.array_equal(data_2dr[0].time.data, data_2ds[0].time.data))
        assert(np.array_equal(data_2dr[0].time.data, data_2dv[0].time.data))
        assert(data_2dr[0].time.units == data_2dv[0].time.units)


class SaveNCDF_sc():
    """ write output netCDF file for abstracted variables from original 2D
        Constant Model Parameters
    """

    def saveData(self, just_the_data, dir_data):
        """ Create a NetCDF file for saving output variables """

        #set up file path and names
        file_ncdf  = path.join(dir_data,("merra_sc" + ".nc"))
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        print("Saved File Type:", rootgrp.file_format)
        rootgrp.source      = 'MERRA2 constant model parameters'

        # Create dimensions
        time  = rootgrp.createDimension('time', None)
        lats  = rootgrp.createDimension('lats', len(just_the_data.lat))
        lons  = rootgrp.createDimension('lons', len(just_the_data.lon))

        # Output the results of extracted variables
        data_variables = [v for v in list(just_the_data) if v not in ['time', 'lat', 'lon']]

        for x in list(data_variables):
            out_var = rootgrp.createVariable(x, 'f4',
                                                ('time','lats','lons'),
                                                fill_value=9.9999999E14)

            out_var.long_name = just_the_data[x].long_name
            out_var.units = just_the_data[x].units
            out_var.missing_value = just_the_data[x].missing_value

            out_var[:] = just_the_data[x].data[0]

        # Create dimension variables
        Time  = rootgrp.createVariable('time', 'i4', ('time'))

        Time.standard_name = "time"
        Time.units  = "hours since 1992-01-02 03:00:00"
        Time.calendar      = "gregorian"
        Time[:] = just_the_data.time.data

        Latitudes  = rootgrp.createVariable('latitude', 'f4', ('lats'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = just_the_data.lat.data

        Longitudes  = rootgrp.createVariable('longitude', 'f4', ('lons'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = just_the_data.lon.data

        rootgrp.close()

# TODO: update functions of checking available files
class MERRAdownload(GenericDownload):
    """
    Class for MERRA-2 data that has methods for querying NASA GES DISC server,
    and returning all variables usually wanted.
    """
    
    def __init__(self, pfile):
        super().__init__(pfile)
        par = self.par

        self._set_data_directory("merra2")

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size']

        # time bounds
        init_date  = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                      'end': datetime.strptime(par['end'], '%Y/%m/%d')}

        self.date = self.update_time_bounds(init_date)

        # start connection session
        self.credential = path.join(par['credentials_directory'], ".merrarc")
        self.account = open(self.credential, "r")
        self.inf = self.account.readlines()
        self.username = ''.join(self.inf[0].split())
        self.password = ''.join(self.inf[1].split())
        self.start_session()

        # Create Subsetter objects
        self.build_subsetters()
    
    def update_time_bounds(self, date):
        """ Update time bounds for downloading so that files aren't re-downloaded
        date : dict
            Dictionary with keys 'beg' and 'end' and values as datetime objects
        Returns
        -------
        date : dict
            Updated dictionary with keys 'beg' and 'end' and values as datetime objects
        """
        import glob
        import re
        files = glob.glob(f"{self.directory}/*.nc") 
        
        pl = re.compile(r"merra_pl_(\d{8})_to_(\d{8})")
        sa = re.compile(r"merra_sa_(\d{8})_to_(\d{8})")
        sf = re.compile(r"merra_sf_(\d{8})_to_(\d{8})")

        last_dates = [[F.search(x).group(2) for x in files if F.search(x)] for F in [pl, sa, sf]]
        
        if last_dates == [[],[],[]]:  # No matching files
            return date
        
        else:
            last_completed = [datetime.strptime(max(D), "%Y%m%d") for D in last_dates]

            if min(last_completed) == date['end']:
                print("All data files have already been downloaded. Exiting.")
                sys.exit(0)

            if min(last_completed) != max(last_completed):  # Partially finished chunk
                date['beg'] = min(last_completed) - timedelta(days=self.chunk_size - 1)
            
            else:  # Completely finished chunk
                date['beg'] = max(last_completed) + timedelta(days=1)

            print(f"Previously downloaded files detected. Beginning download at {date['beg']}")
            
        return date

    def getURLs(self, date):
        """ Set up urls by given range of date and type of data to get objected
        url address

        Args:
        date : dict

        Returns
        -------
        list
        """
        # Setup the based url strings
        baseurl_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/')
        baseurl_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/')

        baseurl_3dn = (baseurl_3d + 'M2I6NPANA.5.12.4/{YM}/MERRA2_{FN}.inst6_3d_ana_Np.{YMD}.nc4')
        baseurl_3da = (baseurl_3d + 'M2I3NPASM.5.12.4/{YM}/MERRA2_{FN}.inst3_3d_asm_Np.{YMD}.nc4')
        baseurl_2dm = (baseurl_2d + 'M2I1NXASM.5.12.4/{YM}/MERRA2_{FN}.inst1_2d_asm_Nx.{YMD}.nc4')
        baseurl_2dr = (baseurl_2d + 'M2T1NXRAD.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_rad_Nx.{YMD}.nc4')
        baseurl_2ds = (baseurl_2d + 'M2T1NXFLX.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_flx_Nx.{YMD}.nc4')
        baseurl_2dv = (baseurl_2d + 'M2T1NXSLV.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_slv_Nx.{YMD}.nc4')

        # build the urls list
        urls_3dmana = []
        urls_3dmasm = []
        urls_2dm = []
        urls_2dr = []
        urls_2ds = []
        urls_2dv = []

        for d in pd.date_range(date['beg'], date['end']):
            ym = d.strftime("%Y/%m")
            ymd = d.strftime("%Y%m%d")
            fn = self.get_file_number(d.year)

            urls_3dmana.append(baseurl_3dn.format(YM=ym, FN=fn, YMD=ymd))
            urls_3dmasm.append(baseurl_3da.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dm.append(baseurl_2dm.format(YM=ym, FN=fn, YMD=ymd))
            urls_2ds.append(baseurl_2ds.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dr.append(baseurl_2dr.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dv.append(baseurl_2dv.format(YM=ym, FN=fn, YMD=ymd))

        # Setup URL for getting constant model parameters (2D, single-level, full horizontal resolution)
        url_2dc = ['https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4']

        return urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv

    @staticmethod
    def get_file_number(year):
        """
        The file names consist of a number and a meta data string.
        The number changes over the years. 1980 until 1991 it is 100,
        1992 until 2000 it is 200, 2001 until 2010 it is  300
        and from 2011 until now it is 400.
        """
        file_number = ''

        if year >= 1980 and year < 1992:
            file_number = '100'
        elif year >= 1992 and year < 2001:
            file_number = '200'
        elif year >= 2001 and year < 2011:
            file_number = '300'
        elif year >= 2011:
            file_number = '400'
        else:
            raise Exception('The specified year is out of range.')
        return file_number

    @staticmethod
    def tempExtrapolate(t_total, h_total, elevation):
        """ Processing 1D vertical extrapolation for Air Temperature, at where
            the values are lacking  (marked by 9.9999999E14) from merra-2 3D
            Analyzed Meteorological Fields datasets.
            IMPORTANT TIP:
            To set up 'ele_max = 2500' (meter) or higher
            Reason: to make sure. get enough levels of geopotential height for
            conducting 1dinterp (linear) (2 points of values needed at least)
        """

        # restructure t_total [time*lev*lat*lon] to [lat*lon*time*lev]
        t_total = t_total[:,:,:,:].transpose((2,3,0,1))
        h_total = h_total[:,:,:,:].transpose((2,3,0,1))

        # find the value gap and conduct 1d extrapolation
        for i in range(0, len(t_total)):
            for j in range(0, len(t_total[0])):
                t_time = t_total[i][j][:]
                h_time = h_total[i][j][:]
                for k in range(0, len(t_time)):
                    t_lev = t_time[k][:]
                    h_lev = h_time[k][:]
                    id_interp = []
                    for z in range(0, len(t_lev)):
                        # find the indices of levels with missing values
                        if t_lev[z] > 99999:
                            id_interp.append(z)

                            if id_interp != []:
                                # get the levels of geopotential heights with missing values
                                lev_interp = h_lev[id_interp]
                                # pass the index of first found level with existing value to z_top
                                z_top = id_interp[-1] + 1
                                # get values at the lowest 3 levels of geopotential heights with existed values
                                lev_3p = h_lev[z_top:z_top + 3]
                                # get values at the lowest 3 levels of air temperature with existed values
                                t_3p = t_lev[z_top:z_top + 3]
                                # Using spicy.interpolate.interp1d function-------------------------
                                # Require >= 2 points of levs and t in minimum
                                if len(lev_3p) >= 2:
                                    # build linear function based on given values at lowest 3 levels of air temperature and geopotential heights
                                    f = interp1d(lev_3p, t_3p, kind='linear', fill_value='extrapolate')
                                    # use built function to calculate the values of air temperature at the found missing-values levels
                                    t_interp = f(lev_interp)
                                    # fill the calculated values into missing-values levels
                                    t_lev[id_interp] = t_interp
                                else:
                                    print('Numbers of points for extrapolation are too low (less then 2):', len(lev_3p))
                                    print('Failed to conduct extrapolation at some points in the output')
                                    print('Current ele_max =', elevation['max'])
                                    print('Higher Value of "ele_max" is needed to reset: > 2500')
                                    sys.exit(0)

                        else:
                            t_lev[z] = t_lev[z]
                        h_lev[z] = h_lev[z]

                    # assign back
                    t_time[k][:] = t_lev
                    h_time[k][:] = h_lev

                # replace the extrapolated value [time * level] to each individual cell
                t_total[i][j][:] = t_time
                h_total[i][j][:] = h_time

        # restructure back
        t_total = t_total[:,:,:,:].transpose((2,3,0,1))
        h_total = h_total[:,:,:,:].transpose((2,3,0,1))

        return t_total

    @staticmethod
    def wind_rh_Extrapolate(data_total):
        """Processing 1D vertical extrapolation for UV and RH, at where the
           values are lacking (marked by 9.9999999E14) from merra-2 3d Analyzed
           Meteorological Fields datasets.
           Wind (U,V) and Relative Humidity (RH) are utilized the value of at
           lowest pressure levels to the ones with value gaps
        """
        # restructure u_total,v_total [time*lev*lat*lon] to [lat*lon*time*lev]
        data_total = data_total[:,:,:,:].transpose((2,3,0,1))

        # find and fill the value gap
        for i in range(0, len(data_total)):
            for j in range(0,len(data_total[0])):
                data_time = data_total[i][j][:]
                for k in range(0,len(data_time)):
                    data_lev = data_time[k][:]
                    id_interp = []
                    for z in range(0, len(data_lev)):
                        if data_lev[z] > 99999:
                            id_interp.append(z)

                        if id_interp != []:
                            z_top = id_interp[-1] + 1
                            data_lev[id_interp] = data_lev[z_top]
                        else:
                            data_lev[z] = data_lev[z]

                    data_time[k][:] = data_lev

                # replace the interpolation value to each single pixel
                data_total[i][j][:] = data_time

        # restructure back
        data_total = data_total[:,:,:,:].transpose((2,3,0,1))

        return data_total

    # build full dictionary between variable names from input parameter
    # file and original merra2 data products
    full_variables_dic = {
            'air_temperature': ['air_temperature',
                                '2-meter_air_temperature'],

            'relative_humidity' : ['relative_humidity',
                                   '2-metre_dewpoint_temperature',
                                   '2-metre_specific_humidity'],

            'precipitation_amount': ['total_precipitation',
                                     'total_precipitation_corrected'],

            'wind_from_direction': ['eastward_wind',
                                    'northward_wind',
                                    '2-meter_eastward_wind',
                                    '2-meter_northward_wind',
                                    '10-meter_eastward_wind',
                                    '10-meter_northward_wind'],

            'wind_speed':['eastward_wind',
                          'northward_wind',
                          '2-meter_eastward_wind',
                          '2-meter_northward_wind',
                          '10-meter_eastward_wind',
                          '10-meter_northward_wind'],

            'downwelling_shortwave_flux_in_air': ['surface_incoming_shortwave_flux'],

            'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['surface_incoming_shortwave_flux_assuming_clear_sky'],

            'downwelling_longwave_flux_in_air': ['surface_net_downward_longwave_flux',
                                                 'longwave_flux_emitted_from_surface'],

            'downwelling_longwave_flux_in_air_assuming_clear_sky':['surface_net_downward_longwave_flux_assuming_clear_sky',
                                                                   'longwave_flux_emitted_from_surface']}

    # build variables Standards Names and referenced Names for downloading
    # from orginal MERRA-2 datasets

    # 3D Analyzed Meteorological fields data
    full_variables_pl_ana = {'geopotential_height':'H',
                             'air_temperature':'T',
                             'eastward_wind':'U',
                             'northward_wind': 'V'}
    # 3D Assimilated Meteorological fields data
    full_variables_pl_asm = {'relative_humidity': 'RH'}

    # 2D Single-Level Diagnostics data (instantaneous)
    full_variables_sm = {'2-meter_air_temperature': 'T2M',
                         '2-meter_eastward_wind': 'U2M',
                         '2-meter_northward_wind':'V2M',
                         '10-meter_eastward_wind':'U10M',
                         '10-meter_northward_wind':'V10M',
                         '2-metre_specific_humidity':'QV2M'}

    # 2D surface flux diagnostics data
    full_variables_sf = {'total_precipitation': 'PRECTOT',
                         'total_precipitation_corrected': 'PRECTOTCORR'}

    # 2D single-level diagnostics (time-averageed)
    full_variables_sv = {'2-metre_dewpoint_temperature': 'T2MDEW'}

    # 2D radiation diagnostics data
    full_variables_sr = {'surface_incoming_shortwave_flux': 'SWGDN',
                         'surface_incoming_shortwave_flux_assuming_clear_sky': 'SWGDNCLR',
                         'surface_net_downward_longwave_flux':'LWGNT',
                         'longwave_flux_emitted_from_surface': 'LWGEM',
                         'surface_net_downward_longwave_flux_assuming_clear_sky': 'LWGNTCLR'}

    def start_session(self):
        self.session = setup_session(self.username, self.password, 
        check_url="https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4")

    def build_subsetters(self):
        self.subsetters = {
            "3dmana": MERRASubsetter('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4', self.session),
            "3dmasm": MERRASubsetter('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4/2016/01/MERRA2_400.inst3_3d_asm_Np.20160101.nc4', self.session),
            "2dm": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I1NXASM.5.12.4/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4', self.session),
            "2ds": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXFLX.5.12.4/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4', self.session),
            "2dr": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXRAD.5.12.4/1981/01/MERRA2_100.tavg1_2d_rad_Nx.19810101.nc4', self.session),
            "2dv": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXSLV.5.12.4/2016/01/MERRA2_400.tavg1_2d_slv_Nx.20160101.nc4', self.session)
            }

        for s in self.subsetters.values():
            s.set_lon_range(self.area['west'], self.area['east'])
            s.set_lat_range(self.area['south'], self.area['north'])
            s.set_elev_range(self.elevation['min'], self.elevation['max'])

        # set target variables for subsetter objects
        self.subsetters['3dmasm'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_pl_asm))  # Add H manually to match old code
        self.subsetters['3dmana'].set_variables(["H"] + self.getVariables(self.full_variables_dic, self.full_variables_pl_ana))

        self.subsetters['2dm'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sm))
        self.subsetters['2ds'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sf))
        self.subsetters['2dr'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sr))
        self.subsetters['2dv'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sv))

    def getVariables(self, full_variables_dic, full_variables_type):
        """
        build the major variable list for retrieving between lists from
        download parameter file and the product
        """
        get_variables = []

        for i in range(0, len(self.variables)):
            for var in full_variables_dic.keys():
                if  var == self.variables[i]:
                    var_names = full_variables_dic[var]
                    #  Set up the variables list for accessing original type
                    # of MERRA-2 datasets (3D and 2D)
                    for j in range(0, len(var_names)):
                        for var1 in full_variables_type.keys():
                            if var1 == var_names[j]:
                                get_variables.append(full_variables_type[var1])
        # set the list
        get_variables = list(set(get_variables))
        """
        #add the basic variables into the list
        if 'relative_humidity' in list(full_variables_type.keys()) :
            get_variables.extend(['lat','lon','lev','time'])
        elif 'air_temperature' in list(full_variables_type.keys()):
            # !ADD Geopotential Height in the first element of downloading
            # list. Must be the first one
            get_variables.insert(0,'H')
            # add the variables names of latitude, longitude, levels and time
            get_variables.extend(['lat','lon','lev','time'])
        else:
            # add the variables names of latitude, longitude and time
            get_variables.extend(['lat','lon','time'])
        """
        return get_variables

    def retrieve(self):
        """
        Retrive all required MERRA-2 data from NASA Goddard Earth Sciences Data
        and Information Services Center
        """

        # Chunk size for spliting files and download [days], Format:Integer
        chunk_size = int(self.chunk_size)

        # Get merra-2 2d Constant Model Parameters (outside of time & date looping!)
        print('''========== Get Variables From MERRA2 2d Time-Invariant ==========''')
        print('''========== Single-level, Constant Model Parameters ==========''')

        self.download_merra_to()

        print("----------- Result set NO.0: completed -----------")

        # build chunks
        all_dates = pd.date_range(self.date['beg'], self.date['end'])
        chunks = [{'beg': x, 'end': y} for (x,y) in zip(all_dates[0::chunk_size], all_dates[chunk_size - 1::chunk_size])]

        if not chunks[-1]['end'] == self.date['end']:  # if chunks not an even multiple
            chunks.append({'beg': chunks[-1]['beg'] + timedelta(days=chunk_size),
                           'end': self.date['end']})

        for date_range in chunks:
            print(f"=== Downloading chunk {date_range['beg']} to {date_range['end']}===")

            # Build the urls for the chunk range
            urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = self.getURLs(date_range)

            # Access, subset (server-side) and download the data (in memory)
            print("M2I6NPANA", end="\r")
            data_3dmana = [self.subsetters['3dmana'].subset_dataset(url, metadata=True) for url in urls_3dmana]

            print("M2I3NPASM", end="\r")
            data_3dmasm = [self.subsetters['3dmasm'].subset_dataset(url, metadata=True) for url in urls_3dmasm]

            print("M2I1NXASM", end="\r")
            data_2dm = [self.subsetters['2dm'].subset_dataset(url, metadata=True) for url in urls_2dm]

            print("M2T1NXFLX", end="\r")
            data_2ds = [self.subsetters['2ds'].subset_dataset(url, metadata=True) for url in urls_2ds]

            print("M2T1NXRAD", end="\r")
            data_2dr = [self.subsetters['2dr'].subset_dataset(url, metadata=True) for url in urls_2dr]

            print("M2T1NXSLV", end="\r")
            data_2dv = [self.subsetters['2dv'].subset_dataset(url, metadata=True) for url in urls_2dv]

            print("Writing to NC file", end="\r")

            SaveNCDF_sa().saveData(data_2dm, self.directory)
            SaveNCDF_sf().saveData(data_2dr, data_2ds, data_2dv, self.directory)
            SaveNCDF_pl_3dm().saveData(data_3dmasm, data_3dmana, chunk_size, self.directory, self.elevation)

    def download_merra_to(self):
        """ Download time-invariant data for the specified lat-lon boundaries """
        # check if it already exists
        if path.isfile(path.join(self.directory, ("merra_sc" + ".nc"))):
            print("WARNING:  file 'merra_to.nc' already exists and is being skipped")
            return()

        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = self.getURLs(self.date)
        S = MERRASubsetter(url_2dc[0], self.session)

        get_variables_2dc = ['PHIS','FRLAKE','FRLAND','FRLANDICE',
                             'FROCEAN','SGH','lat','lon','time']

        S.set_elev_range(self.elevation['min'], self.elevation['max'])
        S.set_lat_range(self.area['south'], self.area['north'])
        S.set_lon_range(self.area['west'], self.area['east'])
        S.set_variables(get_variables_2dc[:-3])

        data = S.subset_dataset(url_2dc[0], metadata=True)

        SaveNCDF_sc().saveData(data, self.directory)


class MERRASubsetter:

    LEVS = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750,
                     725, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250,
                     200, 150, 100, 70, 50, 40, 30, 20, 10, 7.0, 5.0, 4.0,
                     3.0, 2.0, 1.0, 0.7, 0.5, 0.4, 0.3, 0.1])

    lat_values = np.arange(-90, 90.5, 0.5)
    lon_values = np.arange(-180, 180, 0.625)

    def __init__(self, url, session):
        """ For subsetting datasets on the MERRA OpenDAP server

        url : str
            path to pydap-accessible URL endpoint for example dataset (same variables and dimensions)

        session : pydap access session object created from setup_session

        2d VAR[lat:lat][lon:lon][time:time] """
        self.session = session
        self.dataset = open_url(url, session=session)
        self.set_time_step(self.dataset)

    def set_time_step(self, dataset_type):
        try:
            self.n_timesteps = dataset_type['time'].shape[0]
        except KeyError:
            raise KeyError("The selected dataset does not have a time dimension, or the time dimension is not named 'time'.")

    def subset_time(self):
        """ Get the full day """
        return "[0:{}]".format(self.n_timesteps - 1)

    def subset_lat(self, lat_min, lat_max):
        indices = np.where((self.lat_values >= lat_min) & (self.lat_values <= lat_max))
        return "[{}:{}]".format(np.min(indices), np.max(indices))

    def subset_lon(self, lon_min, lon_max):
        indices = np.where((self.lon_values >= lon_min) & (self.lon_values <= lon_max))
        return "[{}:{}]".format(np.min(indices), np.max(indices))

    def subset_lev(self, elev_min, elev_max):
        Pmax = pressure_from_elevation(elev_min) + 55
        Pmin = pressure_from_elevation(elev_max) - 55

        indices = np.where((self.LEVS >= Pmin) & (self.LEVS <= Pmax))
        return "[{}:{}]".format(np.min(indices), np.max(indices))

    def subset_2d_variable(self, variable_name, lat_min, lat_max, lon_min, lon_max):
        var_string = (variable_name
                      + self.subset_time()
                      + self.subset_lat(lat_min, lat_max)
                      + self.subset_lon(lon_min, lon_max))

        return var_string

    def subset_3d_variable(self, variable_name, elev_min, elev_max, lat_min, lat_max, lon_min, lon_max):
        var_string = (variable_name
                      + self.subset_time()
                      + self.subset_lev(elev_min, elev_max)
                      + self.subset_lat(lat_min, lat_max)
                      + self.subset_lon(lon_min, lon_max))

        return var_string

    def set_lat_range(self, lat_min, lat_max):
        self.lat_min = lat_min
        self.lat_max = lat_max

    def set_lon_range(self, lon_min, lon_max):
        self.lon_min = lon_min
        self.lon_max = lon_max

    def set_elev_range(self, elev_min, elev_max):
        self.elev_min = elev_min
        self.elev_max = elev_max

    def set_variables(self, variables):
        self.variables = variables

    def subset_variable(self, variable):
        """ Generate subset slice string for a variable

        Returns
        -------
        str : variable name with slices, e.g. "H[0:3][4:59][23:44]"
        """
        is_3d = ('lev' in self.dataset[variable].dimensions)

        if is_3d:
            slice_string = self.subset_3d_variable(variable, self.elev_min, self.elev_max,
                                                   self.lat_min, self.lat_max,
                                                   self.lon_min, self.lon_max)
        else:
            slice_string = self.subset_2d_variable(variable,
                                                   self.lat_min, self.lat_max,
                                                   self.lon_min, self.lon_max)

        return slice_string

    def create_dods_url(self, dataset_url):
        """ Create a DODS uri link to access a subset of the data by adding variables and slices """
        uri_parameters = []

        for var in self.variables:
            uri_parameters.append(self.subset_variable(var))

        uri_parameters.append('lon' + self.subset_lon(self.lon_min, self.lon_max))
        uri_parameters.append('lat' + self.subset_lat(self.lat_min, self.lat_max))
        uri_parameters.append('time' + self.subset_time())

        if 'lev' in self.dataset:
            uri_parameters.append('lev' + self.subset_lev(self.elev_min, self.elev_max))

        dods_url = dataset_url + ".dods?" + ",".join(uri_parameters)

        return dods_url

    def subset_dataset(self, url, metadata=False):
        """ Return a subset dataset as

        Parameters
        ----------
        url : [type]
            base URL for a dataset that can be opened with pydap.open_url
        variables : list
            List of MERRA variable names to be downloaded

        Returns
        -------
        DatasetType : pydap DatasetType
        """

        dods_url = self.create_dods_url(url)
        ds = open_dods(dods_url, session=self.session, metadata=metadata)

        return ds
