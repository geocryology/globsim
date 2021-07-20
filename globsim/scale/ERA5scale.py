# -*- coding: utf-8 -*-
#
# Methods for downloading ERA5 data from the ECMWF server for limited
# areas and limited times.
#
#  OVERALL WORKFLOW
#
#  See: https://github.com/geocryology/globsim/wiki/Globsim
#       and the examples directory on github
#
# ECMWF and netCDF information:
# https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
#
# For variable codes and units, see:
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
# Check ECMWF job status: http://apps.ecmwf.int/webmars/joblist/
#
#  CF-Convention: this is how you check netcdf file conformity:
#  http://pumatest.nerc.ac.uk/cgi-bin/cf-checker-dev.pl
#
# ===============================================================================

from __future__ import print_function

import netCDF4 as nc
import numpy as np
import urllib3

from datetime import timedelta
from math import floor, atan2, pi
from os import path
from scipy.interpolate import interp1d

from globsim.common_utils import series_interpolate
from globsim.meteorology import spec_hum_kgkg
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale

urllib3.disable_warnings()


class ERA5scale(GenericScale):
    """
    Class for ERA5 data that has methods for scaling station data to
    better resemble near-surface fluxes.

    Processing kernels have names in UPPER CASE.

    Args:
        sfile: Full path to a Globsim Scaling Parameter file.

    Example:
        ERAd = ERA5scale(sfile)
        ERAd.process()
    """

    def __init__(self, sfile, era5type='reanalysis'):
        super().__init__(sfile)
        par = self.par

        if era5type == 'reanalysis':
            self.ens = False
            self.src = 'era5_'
        elif era5type == 'ensemble_members':
            self.ens = True
            self.src = 'era5_ens'

        # input file handles
        self.nc_pl = nc.Dataset(path.join(self.intpdir,
                                          '{}_pl_'.format(self.src) +
                                          self.list_name + '_surface.nc'), 'r')
        self.nc_sa = nc.Dataset(path.join(self.intpdir,
                                          '{}_sa_'.format(self.src) +
                                self.list_name + '.nc'), 'r')
        self.nc_sf = nc.Dataset(path.join(self.intpdir,
                                          '{}_sf_'.format(self.src) +
                                self.list_name + '.nc'), 'r')
        self.nc_to = nc.Dataset(path.join(self.intpdir,
                                          '{}_to_'.format(self.src) +
                                self.list_name + '.nc'), 'r')
        self.nstation = len(self.nc_to.variables['station'][:])

        # time vector for output data
        # get time and convert to datetime object
        nctime = self.nc_pl.variables['time'][:]
        self.t_unit = self.nc_pl.variables['time'].units  # "hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_pl.variables['time'].calendar
        time = nc.num2date(nctime, units=self.t_unit, calendar=self.t_cal)

        # interpolation scale factor
        self.time_step = par['time_step'] * 3600    # [s] scaled file
        self.interval_in = (time[1] - time[0]).seconds  # interval in seconds
        self.interpN = floor(self.interval_in / self.time_step)

        # number of time steps for output, include last value
        self.nt = int(floor((max(time) - min(time)).total_seconds()
                            / 3600 / par['time_step'])) + 1

        # vector of output time steps as datetime object
        # 'seconds since 1900-01-01 00:00:0.0'
        mt = min(time)
        self.times_out = [mt + timedelta(seconds=(x * self.time_step))
                          for x in range(0, self.nt)]

        # vector of output time steps as written in ncdf file [s]
        self.scaled_t_units = 'seconds since 1900-01-01 00:00:0.0'
        self.times_out_nc = nc.date2num(self.times_out,
                                        units=self.scaled_t_units,
                                        calendar=self.t_cal)

    def getValues(self, ncf, varStr, ni=10):
        if self.ens:
            values = ncf.variables[varStr][:, ni, :]
        else:
            values = ncf.variables[varStr][:]

        return values

    def indProcess(self, ni=10):
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                print(kernel_name)
                getattr(self, kernel_name)(ni)

    def process(self, ni=10):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """

        stations = self.stations['station_name']
        # iterate thorugh kernels and start process
        if self.ens:
            for ni in self.nc_sa['number']:
                src = '{}_{}'.format(self.src, ni)
                self.output_file = self.getOutNCF(self.par, src)
                self.rg = new_scaled_netcdf(self.output_file,
                                            self.nc_pl, self.times_out_nc,
                                            t_unit=self.scaled_t_units,
                                            station_names=stations)
                self.indProcess(ni)

        else:
            self.output_file = self.getOutNCF(self.par, self.src)
            self.rg = new_scaled_netcdf(self.output_file,
                                        self.nc_pl, self.times_out_nc,
                                        t_unit=self.scaled_t_units,
                                        station_names=stations)
            self.indProcess(ni=10)

        # close netCDF files
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()

    def PRESS_Pa_pl(self, ni=10):
        """
        Surface air pressure from pressure levels.
        """
        # add variable to ncdf file
        vn = 'PRESS_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_pressure ERA-5 pressure levels only'
        var.units     = 'Pa'
        var.standard_name = 'surface_air_pressure'

        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_pl, 'air_pressure', ni)
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            # scale from hPa to Pa
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                          time_in * 3600, values[:, n]) * 100

    def AIRT_C_pl(self, ni=10):
        """
        Air temperature derived from pressure levels, exclusively.
        """
        # add variable to ncdf file
        vn = 'AIRT_pl'  # variable name
        var           = self.rg.createVariable(vn, 'f4', ('time', 'station'))
        var.long_name = 'air_temperature ERA-5 pressure levels only'
        var.units     = 'degrees_C'
        var.standard_name = 'air_temperature'

        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_pl, 't', ni)  # self.nc_pl.variables['t'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600, values[:, n] - 273.15)

    def AIRT_C_sur(self, ni=10):
        """
        Air temperature derived from surface data, exclusively.
        """
        # add variable to ncdf file
        vn = 'AIRT_sur'   # variable name
        var           = self.rg.createVariable(vn, 'f4', ('time', 'station'))
        var.long_name = '2_metre_temperature ERA-5 surface only'
        var.units     = self.nc_sa.variables['t2m'].units.encode('UTF8')

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sa, 't2m', ni)  # self.nc_sa.variables['t2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
                                                             values[:, n] - 273.15)

    def AIRT_redcapp(self, ni=10):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP.
        """
        print("AIRT_redcapp")

    def PREC_mm_sur(self, ni=10):
        """
        Precipitation sum in mm for the time step given.
        """
        # add variable to ncdf file
        vn  = 'PREC_sur'  # variable name
        var = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Total precipitation ERA-5 surface only'
        var.units     = 'kg m-2 s-1'
        var.standard_name = 'precipitation_amount'

        # interpolation scale factor
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)

        # total prec [mm] in 1 second
        values  = self.getValues(self.nc_sf, 'tp', ni)  # self.nc_sf.variables['tp'][:]*1000/self.interval_in
        values  = values * 1000 / self.interval_in

        # interpolate station by station
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.time_step

    def RH_per_sur(self, ni=10):
        """
        Relative humdity derived from surface data, exclusively. Clipped to
        range [0.1,99.9]. Kernel AIRT_C_sur must be run before.
        """
        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sa, 'd2m', ni)  # self.nc_sa.variables['d2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            dewp[:, n] = series_interpolate(self.times_out_nc,
                                            time_in * 3600, values[:, n] - 273.15)

        # add variable to ncdf file
        vn = 'RH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity ERA-5 surface only'
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        # simple: https://doi.org/10.1175/BAMS-86-2-225
        RH = 100 - 5 * (self.rg.variables['AIRT_sur'][:, :] - dewp[:, :])
        self.rg.variables[vn][:, :] = RH.clip(min=0.1, max=99.9)

    def WIND_sur(self, ni=10):
        """
        Wind speed and direction temperature derived from surface data,
        exclusively.
        """
        # temporary variable, interpolate station by station
        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sa, 'u10', ni)  # self.nc_sa.variables['u10'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            U[:, n] = series_interpolate(self.times_out_nc,
                                         time_in * 3600, values[:, n])

        # temporary variable, interpolate station by station
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sa, 'v10', ni)  # self.nc_sa.variables['v10'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            V[:, n] = series_interpolate(self.times_out_nc,
                                         time_in * 3600, values[:, n])

        # wind speed, add variable to ncdf file, convert
        vn_spd = 'WSPD_sur'  # variable name
        var           = self.rg.createVariable(vn_spd,'f4',('time', 'station'))
        var.long_name = '10 wind speed ERA-5 surface only'
        var.units     = 'm s-1'
        var.standard_name = 'wind_speed'

        # wind direction, add variable to ncdf file, convert, relative to North
        vn_dir = 'WDIR_sur'  # variable name
        var           = self.rg.createVariable(vn_dir,'f4',('time', 'station'))
        var.long_name = '10 wind direction ERA-5 surface only'
        var.units     = 'degree'
        var.standard_name = 'wind_from_direction'

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            WS = np.sqrt(np.power(V,2) + np.power(U,2))
            WD = [atan2(V[i, n], U[i, n]) * (180 / pi) +
                  180 for i in np.arange(V.shape[0])]
            self.rg.variables[vn_spd][:, n] = WS
            self.rg.variables[vn_dir][:,n] = WD

    def SW_Wm2_sur(self, ni=10):
        """
        Short-wave downwelling radiation derived from surface data, exclusively.
        This kernel only interpolates in time.
        """

        # add variable to ncdf file
        vn = 'SW_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Surface solar radiation downwards ERA-5 surface only'
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_shortwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sf, 'ssrd', ni)  # self.nc_sf.variables['ssrd'][:]/3600/self.interval_in#[w/m2/s]
        values  = values / 3600 / self.interval_in  # [w/m2/s]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.time_step

    def LW_Wm2_sur(self, ni=10):
        """
        Long-wave downwelling radiation derived from surface data, exclusively.
        This kernel only interpolates in time.
        """

        # add variable to ncdf file
        vn = 'LW_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Surface thermal radiation downwards ERA-5 surface only'
        var.standard_name = 'surface_downwelling_longwave_flux'
        var.units     = 'W m-2'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sf, 'strd', ni)  # self.nc_sf.variables['strd'][:]/3600/self.interval_in #[w m-2 s-1]
        values  = values / 3600 / self.interval_in  # [w m-2 s-1]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.time_step

    def SH_kgkg_sur(self, ni=10):
        '''
        Specific humidity [kg/kg]
        https://crudata.uea.ac.uk/cru/pubs/thesis/2007-willett/2INTRO.pdf
        '''

        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_sa, 'd2m', ni)  # self.nc_sa.variables['d2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            dewp[:, n] = series_interpolate(self.times_out_nc,
                                            time_in * 3600, values[:, n] - 273.15)

        # compute
        SH = spec_hum_kgkg(dewp[:, :],
                           self.rg.variables['PRESS_pl'][:, :])

        # add variable to ncdf file
        vn = 'SH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Specific humidity ERA-5 surface only'
        var.units     = '1'
        var.standard_name = 'specific_humidity'
        self.rg.variables[vn][:, :] = SH

