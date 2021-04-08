#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Methods for downloading ERA-Interim data from the ECMWF server for limited
# areas and limited times.
#
#
# (C) Copyright Stephan Gruber (2013–2017)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

import logging
import netCDF4 as nc
import numpy as np

from datetime import timedelta
from math import floor, atan2, pi
from os import path
from scipy.interpolate import interp1d

from globsim.common_utils import series_interpolate, cummulative2total
from globsim.meteorology import spec_hum_kgkg, LW_downward
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale

logger = logging.getLogger()


class ERAIscale(GenericScale):
    NAME = "ERA-I"

    """
    Class for ERA-Interim data that has methods for scaling station data to
    better resemble near-surface fluxes.

    Processing kernels have names in UPPER CASE.

    Args:
        sfile: Full path to a Globsim Scaling Parameter file.

    Example:
        ERAd = ERAIscale(sfile)
        ERAd.process()
    """

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file handles
        self.nc_pl = nc.Dataset(path.join(self.intpdir, 'erai_pl_' +
                                          self.list_name + '_surface.nc'), 'r')
        self.nc_sa = nc.Dataset(path.join(self.intpdir,'erai_sa_' +
                                          self.list_name + '.nc'), 'r')
        self.nc_sf = nc.Dataset(path.join(self.intpdir, 'erai_sf_' +
                                          self.list_name + '.nc'), 'r')
        self.nc_to = nc.Dataset(path.join(self.intpdir, 'erai_to_' +
                                          self.list_name + '.nc'), 'r')
        self.nstation = len(self.nc_to.variables['station'][:])

        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(par, 'erai')

        # time vector for output data
        # get time and convert to datetime object
        nctime = self.nc_sf.variables['time'][:]
        self.t_unit = self.nc_sf.variables['time'].units #"hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_sf.variables['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)

        # interpolation scale factor
        self.time_step = par['time_step'] * 3600    # [s] scaled file
        interval_in = (time[1]-time[0]).seconds #interval in seconds
        self.interpN = floor(interval_in/self.time_step)

        #number of time steps for output, include last value
        self.nt = int(floor((max(time) - min(time)).total_seconds()
                      / 3600 / par['time_step']))+1

        # vector of output time steps as datetime object
        # 'seconds since 1900-01-01 00:00:0.0'
        mt = min(time)
        self.times_out = [mt + timedelta(seconds = (x*self.time_step))
                          for x in range(0, self.nt)]

        # vector of output time steps as written in ncdf file [s]
        self.scaled_t_units = 'seconds since 1900-01-01 00:00:00'
        self.times_out_nc = nc.date2num(self.times_out,
                                        units = self.scaled_t_units,
                                        calendar = self.t_cal)


    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """
        if path.isfile(self.output_file):
            print("Warning, output file already exists. This may cause problems")
        self.rg = new_scaled_netcdf(self.output_file, self.nc_pl,
                                 self.times_out_nc,
                                 t_unit = self.scaled_t_units,
                                 station_names = self.stations['station_name'])

        # iterate through kernels and start process
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                print(kernel_name)
                getattr(self, kernel_name)()

        # close netCDF files
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()

    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """
        # add variable to ncdf file
        vn = 'PRESS_pl' # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_pressure {} pressure levels only'.format(self.NAME)
        var.units     = 'Pa'.encode('UTF8')
        var.standard_name = 'surface_air_pressure'

        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = self.nc_pl.variables['air_pressure'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            #scale from hPa to Pa
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                             time_in * 3600, values[:, n]) * 100

    def AIRT_C_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """
        # add variable to ncdf file
        vn = 'AIRT_pl' # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_temperature {} pressure levels only'.format(self.NAME)
        var.units     = self.nc_pl.variables['t'].units.encode('UTF8')

        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = self.nc_pl.variables['t'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                        time_in*3600, values[:, n]-273.15)

    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """
        # add variable to ncdf file
        vn = 'AIRT_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '2_metre_temperature {} surface only'.format(self.NAME)
        var.units     = self.nc_sa.variables['t2m'].units.encode('UTF8')

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['t2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                    time_in*3600,
                                                    values[:, n]-273.15)

    def AIRT_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP.
        """
        print("AIRT_redcapp")

    def PREC_mm_sur(self):
        """
        Precipitation sum in mm for the time step given.
        """
        # add variable to ncdf file
        vn = 'PREC_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Total precipitation {} surface only'.format(self.NAME)
        var.units     = 'kg m-2 s-1'
        var.standard_name = 'precipitation_amount'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['tp'][:]*1000 #[mm]

        nctime = self.nc_sf['time'][:]
        self.t_unit = self.nc_sf['time'].units #"hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_sf['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)

        interval_in = (time[1]-time[0]).seconds

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in*3600,
                         cummulative2total(values[:, n], time)/interval_in,
                         kind = 'linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.time_step

    def RH_per_sur(self):
        """
        Relative humdity derived from surface data, exclusively. Clipped to
        range [0.001, 0.999]. Kernel AIRT_ERAI_C_sur must be run before.
        """
        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['d2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            dewp[:, n] = series_interpolate(self.times_out_nc,
                                            time_in*3600, values[:, n]-273.15)

        # add variable to ncdf file
        vn = 'RH_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity {} surface only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        # simple: https://doi.org/10.1175/BAMS-86-2-225
        RH = 100 - 5 * (self.rg.variables['AIRT_sur'][:, :]-dewp[:, :])
        self.rg.variables[vn][:, :] = RH.clip(min=0.1, max=99.9)/100


    def WIND_sur(self):
        """
        Wind speed and direction temperature derived from surface data,
        exclusively.
        """
        # temporary variable, interpolate station by station
        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['u10'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            U[:, n] = series_interpolate(self.times_out_nc,
                                         time_in*3600, values[:, n])

        # temporary variable, interpolate station by station
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['v10'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            V[:, n] = series_interpolate(self.times_out_nc,
                                         time_in*3600, values[:, n])

        # wind speed, add variable to ncdf file, convert
        vna = 'WSPD_sur' # variable name
        var           = self.rg.createVariable(vna,'f4',('time', 'station'))
        var.long_name = '10 wind speed {} surface only'.format(self.NAME)
        var.units     = 'm s-1'
        var.standard_name = 'wind_speed'

        # wind direction, add variable to ncdf file, convert, relative to North
        vnb = 'WDIR_sur' # variable name
        var           = self.rg.createVariable(vnb,'f4',('time', 'station'))
        var.long_name = '10 wind direction {} surface only'.format(self.NAME)
        var.units     = 'degree'
        var.standard_name = 'wind_from_direction'

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            WS = np.sqrt(np.power(V,2) + np.power(U,2))
            WD = [atan2(V[i, n], U[i, n])*(180/pi) +
                  180 for i in np.arange(V.shape[0])]
            self.rg.variables[vna][:, n] = WS
            self.rg.variables[vnb][:,n] = WD

    def SW_Wm2_sur(self):
        """
        Short-wave downwelling radiation derived from surface data, exclusively.
        This kernel only interpolates in time.
        """

        # add variable to ncdf file
        vn = 'SW_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Surface solar radiation downwards {} surface only'.format(self.NAME)
        var.units     = self.nc_sf.variables['ssrd'].units.encode('UTF8')
        var.standard_name = 'surface_downwelling_shortwave_flux'

        # interpolation scale factor
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['ssrd'][:]/3600 #w/m2

        nctime = self.nc_sf['time'][:]
        self.t_unit = self.nc_sf['time'].units #"hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_sf['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)

        interval_in = (time[1]-time[0]).seconds

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in*3600,
                         cummulative2total(values[:, n], time)/interval_in,
                         kind = 'linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.time_step


    def LW_Wm2_sur(self):
        """
        Long-wave downwelling radiation derived from surface data, exclusively.
        This kernel only interpolates in time.
        """

        # add variable to ncdf file
        vn = 'LW_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Surface thermal radiation downwards {} surface only'.format(self.NAME)
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_longwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['strd'][:]/3600 #[w m-2]

        nctime = self.nc_sf['time'][:]
        self.t_unit = self.nc_sf['time'].units #"hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_sf['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)

        # interpolation scale factor
        interval_in = (time[1]-time[0]).seconds

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in*3600,
                         cummulative2total(values[:, n], time)/interval_in,
                         kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)*self.time_step

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg]
        https://crudata.uea.ac.uk/cru/pubs/thesis/2007-willett/2INTRO.pdf
        '''
        # add variable to ncdf file
        vn = 'SH_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Specific humidity {} surface only'.format(self.NAME)
        var.units     = '1'
        var.standard_name = 'specific_humidity'

        # temporary variable, interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['d2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            dewp[:, n] = series_interpolate(self.times_out_nc,
                                            time_in*3600, values[:, n]-273.15)

        # compute
        SH = spec_hum_kgkg(dewp[:, :],
                           self.rg.variables['PRESS_pl'][:, :])


        self.rg.variables[vn][:, :] = SH

    def LW_Wm2_topo(self):
        """
        Long-wave radiation downwards [W/m2]
        https://www.geosci-model-dev.net/7/387/2014/gmd-7-387-2014.pdf
        """
        # get sky view, and interpolate in time
        N = np.asarray(list(self.stations['sky_view'][:]))

        # add variable to ncdf file
        vn = 'LW_topo' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Incoming long-wave radiation {} surface only'.format(self.NAME)
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_longwave_flux'

        # compute
        for i in range(0, len(self.rg.variables['RH_sur'][:])):
            for n, s in enumerate(self.rg.variables['station'][:].tolist()):
                LW = LW_downward(self.rg.variables['RH_sur'][i, n],
                     self.rg.variables['AIRT_sur'][i, n]+273.15, N[n])
                self.rg.variables[vn][i, n] = LW

