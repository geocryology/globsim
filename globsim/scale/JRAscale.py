#!/usr/bin/env python
# -*- coding: utf-8 -*-
import netCDF4 as nc
import numpy as np
import logging

from pathlib import Path
from math import atan2, pi
from scipy.interpolate import interp1d

from globsim.common_utils import str_encode, series_interpolate
from globsim.meteorology import LW_downward
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale
from globsim import __version__ as globsim_version

logger = logging.getLogger('globsim.scale')


class JRAscale(GenericScale):
    """
    Class for JRA-55 data that has methods for scaling station data to
    better resemble near-surface fluxes.

    Processing kernels have names in UPPER CASE.

    Args:
        sfile: Full path to a Globsim Scaling Parameter file.

    Example:
        JRAd = JRAscale(sfile)
        JRAd.process()
    """
    NAME = "JRA-55"

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file names
        self.nc_pl = nc.Dataset(Path(self.intpdir, f'jra_pl_{self.list_name}_surface.nc'),
                                'r')
        self.nc_sa = nc.Dataset(Path(self.intpdir, f'jra_sa_{self.list_name}.nc'),
                                'r')
        self.nc_sf = nc.Dataset(Path(self.intpdir, f'jra_sf_{self.list_name}.nc'),
                                'r')

        for dataset in [self.nc_pl, self.nc_sa, self.nc_sf]:
            dataset.globsim_version = globsim_version

        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(par, 'jra55')

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl.variables['time'], par['time_step'])

        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.time_step,
                                                      num_times=self.nt,
                                                      output_units=self.t_unit,
                                                      output_calendar=self.t_cal)

    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """
        self.rg = new_scaled_netcdf(self.output_file, self.nc_pl,
                                    self.times_out_nc, self.nc_pl['time'].units)

        # add station names to netcdf
        # first convert to character array
        names_out = nc.stringtochar(np.array(self.stations['station_name'], 'S32'))

        # create space in the netcdf
        _            = self.rg.createDimension('name_strlen', 32)
        st           = self.rg.createVariable('station_name', "S1",
                                              ('station', 'name_strlen'))
        st.standard_name = 'platform_name'
        st.units     = ''

        # add data
        st[:] = names_out

        # iterate through kernels and start process
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                logger.info(f"running scaling kernel: '{kernel_name}'")
                getattr(self, kernel_name)()

        # self.conv_geotop()

        # close netCDF files
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()

    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """
        # add variable to ncdf file
        vn = 'PRESS_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_pressure {} pressure levels only'.format(self.NAME)
        var.units     = 'Pa'
        var.standard_name = 'surface_air_pressure'

        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = self.nc_pl.variables['air_pressure'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            # scale from hPa to Pa
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values[:, n]) * 100

    def AIRT_C_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """
        vn = 'AIRT_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_temperature {} pressure levels only'.format(self.NAME)
        var.units     = 'degrees_C'
        var.standard_name = 'air_temperature'

        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:]
        values  = self.nc_pl.variables['Temperature'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n]) - 273.15

    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """

        # add variable to ncdf file
        vn = 'AIRT_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '2_metre_temperature {} surface only'.format(self.NAME)
        var.units     = 'degrees_C'
        var.standard_name = 'air_temperature'

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['Temperature'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n]) - 273.15

    def RH_per_sur(self):
        """
        Relative Humidity derived from surface data, exclusively.
        """
        # add variable to ncdf file
        vn = 'RH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'relative humidity {} surface only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['Relative humidity'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

    def WIND_sur(self):
        """
        Wind at 10 metre derived from surface data, exclusively.
        """

        # add variable to ncdf file
        vn = '10 metre U wind component'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '10 metre U wind component'
        var.units     = str_encode(self.nc_sa.variables['u-component of wind'].units)

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['u-component of wind'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

        # add variable to ncdf file
        vn = '10 metre V wind component'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '10 metre V wind component'
        var.units     = str_encode(self.nc_sa.variables['v-component of wind'].units)

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['v-component of wind'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

        # add variable to ncdf file
        vn = 'WSPD_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '10 metre wind speed {} surface only'.format(self.NAME)
        var.units     = 'm s-1'
        var.standard_name = 'wind_speed'

        # add variable to ncdf file
        vn = 'WDIR_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '10 metre wind direction {} surface only'.format(self.NAME)
        var.units     = 'degree'
        var.standard_name = 'wind_from_direction'

        # convert
        # u is the ZONAL VELOCITY, i.e. horizontal wind TOWARDS EAST.
        # v is the MERIDIONAL VELOCITY, i.e. horizontal wind TOWARDS NORTH.
        V = self.rg.variables['10 metre V wind component'][:]
        U = self.rg.variables['10 metre U wind component'][:]

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            WS = np.sqrt(np.power(V, 2) + np.power(U, 2))
            WD = [atan2(V[i, n], U[i, n]) * (180 / pi) +
                  180 for i in np.arange(V.shape[0])]

            self.rg.variables['WSPD_sur'][:, n] = WS
            self.rg.variables['WDIR_sur'][:,n] = WD

    def SW_Wm2_sur(self):
        """
        solar radiation downwards derived from surface data, exclusively.
        """

        # add variable to ncdf file
        vn = 'SW_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Surface solar radiation downwards {} surface only'.format(self.NAME)
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_shortwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:]
        values  = self.nc_sf.variables['Downward solar radiation flux'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

    def LW_Wm2_sur(self):
        """
        Long-wave radiation downwards derived from surface data, exclusively.
        """

        # add variable to ncdf file
        vn = 'LW_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Surface thermal radiation downwards {} surface only'.format(self.NAME)
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_longwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:]
        values  = self.nc_sf.variables['Downward longwave radiation flux'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

    def PREC_mm_sur(self):
        """
        Precipitation derived from surface data, exclusively.
        Convert unit: mm/day to mm/s (kg m-2 s-1)
        """

        # add variable to ncdf file
        vn = 'PREC_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Total precipitation {} surface only'.format(self.NAME)
        var.units     = 'kg m-2 s-1'
        var.comment = "units [kg m-2 s-1] corresponds to [mm/s] for water (density 1000 [kg m-3])"
        var.standard_name = 'precipitation_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:]
        values  = self.nc_sf.variables['Total precipitation'][:] / (24 * 3600)  # [mm/s]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.scf

    def LW_Wm2_topo(self):
        """
        Long-wave radiation downwards [W/m2]
        https://www.geosci-model-dev.net/7/387/2014/gmd-7-387-2014.pdf
        """
        # get sky view, and interpolate in time
        N = np.asarray(list(self.stations['sky_view'][:]))

        # add variable to ncdf file
        vn = 'LW_topo'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Incoming long-wave radiation {} surface only'.format(self.NAME)
        var.units     = str_encode('W m-2')
        var.standard_name = 'surface_downwelling_longwave_flux'

        # compute
        for i in range(0, len(self.rg.variables['RH_sur'][:])):
            for n, s in enumerate(self.rg.variables['station'][:].tolist()):
                LW = LW_downward(self.rg.variables['RH_sur'][i, n],
                                 self.rg.variables['AIRT_sur'][i, n] + 273.15, N[n])
                self.rg.variables[vn][i, n] = LW

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg] derived from surface data, exclusively
        '''
        # add variable to ncdf file
        vn = 'SH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'specific humidity {} surface only'.format(self.NAME)
        var.units     = 'kg/kg'
        var.standard_name = 'specific_humidity'

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['Specific humidity'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])
