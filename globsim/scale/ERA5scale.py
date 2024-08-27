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
import netCDF4 as nc
import numpy as np
import urllib3
import logging
import pytz

from math import atan2, pi
from scipy.interpolate import interp1d
from pathlib import Path
from pysolar.solar import get_azimuth_fast

from globsim.common_utils import series_interpolate
from globsim.meteorology import spec_hum_kgkg
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.toposcale import lw_down_toposcale, elevation_corrected_sw, sw_partition, solar_zenith, sw_toa, shading_corrected_sw_direct, illumination_angle
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
import globsim.constants as const
import globsim.redcapp as redcapp

urllib3.disable_warnings()
logger = logging.getLogger('globsim.scale')


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
    src = 'era5'

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file handles
        self.nc_pl_sur = nc.Dataset(Path(self.intpdir, f"{self.src}_pl_{self.list_name}_surface.nc"), 'r')
        self.nc_pl = nc.Dataset(Path(self.intpdir, f"{self.src}_pl_{self.list_name}.nc"), 'r')
        self.nc_sa = nc.Dataset(Path(self.intpdir, f'{self.src}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.intpdir, f'{self.src}_sf_{self.list_name}.nc'), 'r')
        self.nc_to = nc.Dataset(Path(self.intpdir, f'{self.src}_to_{self.list_name}.nc'), 'r')

        # Check data integrity
        _check_timestep_length(self.nc_pl_sur.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        if not (self.nc_sf['time'].shape == self.nc_pl_sur['time'].shape == self.nc_sa['time'].shape):
            logger.critical(f"Different number of time steps in input data: sf ({self.nc_sf['time'].shape[0]}) pl ({self.nc_pl_sur['time'].shape[0]}) sa ({self.nc_sa['time'].shape[0]})")

        try:
            self.nc_to.set_auto_scale(True)
        except Exception:
            logger.error("Could not set autoscaling! Values using geopotential (era_to) may be incorrect")

        self.nstation = len(self.nc_to.variables['station'][:])

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl_sur.variables['time'], par['time_step'])
        self.scaled_t_units = 'seconds since 1900-01-01 00:00:0.0'

        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.time_step,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.t_cal)

    def getValues(self, ncf, varStr):
        return ncf.variables[varStr][:]

    def indProcess(self):
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                logger.info(f"running scaling kernel: '{kernel_name}'")
                kernel = getattr(self, kernel_name)
                _ = kernel()
            else:
                logger.error(f"Missing kernel {kernel_name}")

    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """

        stations = self.stations['station_name']
        # iterate thorugh kernels and start process

        self.output_file = self.getOutNCF(self.par, self.src)
        self.rg = new_scaled_netcdf(self.output_file,
                                    self.nc_pl_sur, self.times_out_nc,
                                    t_unit=self.scaled_t_units,
                                    station_names=stations)
        # add surface height
        self.add_grid_elevation(self.rg, self.getValues(self.nc_to, 'z')[0, :] / const.G)

        self.indProcess()

        logger.info(f"Created scaled output file {self.output_file}")
        # close netCDF files
        self.rg.close()
        self.nc_pl_sur.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()

    def PRESS_Pa_pl(self):
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
        time_in = self.nc_pl_sur.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_pl_sur, 'air_pressure')
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            # scale from hPa to Pa
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values[:, n]) * 100

    def AIRT_C_pl(self):
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
        time_in = self.input_times_in_output_units(self.nc_pl_sur)
        values  = self.getValues(self.nc_pl_sur, 't')  # self.nc_pl_sur.variables['t'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in, values[:, n] - 273.15)

    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """
        # add variable to ncdf file
        vn = 'AIRT_sur'   # variable name
        var           = self.rg.createVariable(vn, 'f4', ('time', 'station'))
        var.long_name = '2_metre_temperature ERA-5 surface only'
        var.units     = self.nc_sa.variables['t2m'].units.encode('UTF8')
        var.standard_name = 'air_temperature' 

        # interpolate station by station
        time_in = self.input_times_in_output_units(self.nc_sa)
        values  = self.getValues(self.nc_sa, 't2m')  # self.nc_sa.variables['t2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values[:, n] - 273.15)

    def AIRT_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP Cao et al. (2017) 10.5194/gmd-10-2905-2017
        """
        logger.warning(f"Globsim implementation of REDCAPP only provides Delta_T_c")

        # add variable to ncdf file
        var = redcapp.add_var_delta_T(self.rg)

        # get T from surface level
        T_sa  = self.getValues(self.nc_sa, 't2m')  
        # get grid surface elevation from geopotential  (Cao: elev. @ coarse-scale topography)
        h_sur = self.getValues(self.nc_to, 'z')[0, :] / const.G  # [m]
        # get pressure-level temperatures
        airT_pl = self.getValues(self.nc_pl, 't')
        # get pressure-level elevations from geopotential
        elevation = self.getValues(self.nc_pl, 'z') / const.G  # [m]

        Delta_T_c = redcapp.delta_T_c(T_sa=T_sa, 
                                      airT_pl=airT_pl, 
                                      elevation=elevation,
                                      h_sur=h_sur)  

        time_in = self.input_times_in_output_units(self.nc_sa)
        values  = Delta_T_c
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            var[:, n] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values[:, n])

    def PREC_mm_sur(self):
        """
        Precipitation sum in mm for the time step given.
        """
        # add variable to ncdf file
        vn  = 'PREC_sur'  # variable name
        var = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Total precipitation ERA-5 surface only'
        var.units     = 'mm s-1'
        var.standard_name = 'precipitation_flux'

        # interpolation scale factor
        time_in = self.input_times_in_output_units(self.nc_sf)

        # total prec [mm] in 1 second
        values  = self.getValues(self.nc_sf, 'tp')  # self.nc_sf.variables['tp'][:]*1000/self.interval_in
        # https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
        # We assume interval_in is in hours (1 for regular, 3 for ensemble)
        values  = values * 1000 / (self.interval_in)  # [m] -> [mm s-1]
        # interpolate station by station
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.scf

    def RH_per_sur(self):
        """
        Relative humdity derived from surface data, exclusively. Clipped to
        range [0.1,99.9]. Kernel AIRT_C_sur must be run before.
        """
        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.input_times_in_output_units(self.nc_sa)
        values  = self.getValues(self.nc_sa, 'd2m')  # self.nc_sa.variables['d2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            dewp[:, n] = series_interpolate(self.times_out_nc,
                                            time_in, values[:, n] - 273.15)

        # add variable to ncdf file
        vn = 'RH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity ERA-5 surface only'
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        RH = self._rh()(self.rg.variables['AIRT_sur'][:, :], dewp[:, :])
        self.rg.variables[vn][:, :] = RH.clip(min=0.1, max=99.9)

    def RH_per_pl(self):
        """
        Relative humdity derived from pressure-level.
        """
        # temporary variable,  interpolate station by station
        rh = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_pl_sur.variables['time'][:].astype(np.int64)

        values  = self.getValues(self.nc_pl_sur, 'r')

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            rh[:, n] = series_interpolate(self.times_out_nc,
                                            time_in, values[:, n])

        # add variable to ncdf file
        vn = 'RH_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity ERA-5 pressure-levels'
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        self.rg.variables[vn][:, :] = rh

    def WIND_sur(self):
        """
        Wind speed and direction temperature derived from surface data,
        exclusively.
        """
        # temporary variable, interpolate station by station
        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.input_times_in_output_units(self.nc_sa)
        values  = self.getValues(self.nc_sa, 'u10')  # self.nc_sa.variables['u10'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            U[:, n] = series_interpolate(self.times_out_nc,
                                         time_in, values[:, n])

        # temporary variable, interpolate station by station
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.input_times_in_output_units(self.nc_sa)
        values  = self.getValues(self.nc_sa, 'v10')  # self.nc_sa.variables['v10'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            V[:, n] = series_interpolate(self.times_out_nc,
                                         time_in, values[:, n])


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

        WS = np.sqrt(np.power(V,2) + np.power(U,2))
        self.rg.variables[vn_spd][:, :] = WS

        # Wind direction: 
        # U - easterly component  V - northerly component
        # convert to "clockwise from north" from "anti-clockwise from x-axis" : 90 - angle
        # convert "from direction" : + 180
        WD = 90 - (np.arctan2(V, U) * (180 / np.pi)) + 180
        WD = np.mod(WD, 360)
        self.rg.variables[vn_dir][:, :] = WD

    def SW_Wm2_sur(self):
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
        time_in = self.input_times_in_output_units(self.nc_sf)
        values  = self.getValues(self.nc_sf, 'ssrd')  # self.nc_sf.variables['ssrd'][:]/3600/self.interval_in#[w/m2/s]
        values  = values / (self.interval_in)  # [J m-2] -> [w m-2]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)

    def SW_Wm2_topo(self):
        """
        Short-wave downwelling radiation corrected using a modified version of TOPOscale.
        Partitions into direct and diffuse
        """

        # add variable to ncdf file
        vn_diff = 'SW_topo_diffuse'  # variable name
        var           = self.rg.createVariable(vn_diff,'f4',('time', 'station'))
        var.long_name = 'TOPOscale-corrected diffuse solar radiation'
        var.units     = 'W m-2'
        var.standard_name = 'surface_diffuse_downwelling_shortwave_flux_in_air'

        vn_dir = 'SW_topo_direct'  # variable name
        var           = self.rg.createVariable(vn_dir,'f4',('time', 'station'))
        var.long_name = 'TOPOscale-corrected direct solar radiation'
        var.units     = 'W m-2'
        var.standard_name = 'surface_direct_downwelling_shortwave_flux_in_air'

        # interpolate station by station
        nc_time = self.input_times_in_output_units(self.nc_sf)
        py_time = nc.num2date(nc_time[:], self.scaled_t_units, self.t_cal, only_use_cftime_datetimes=False)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        lat = self.getValues(self.nc_pl_sur, 'latitude')
        lon = self.getValues(self.nc_pl_sur, 'longitude')
        sw = self.getValues(self.nc_sf, 'ssrd') / self.interval_in  # [J m-2] --> [W m-2]
        grid_elev = self.getValues(self.nc_to, 'z')[0,:] / const.G  # z has 2 dimensions from the scaling step
        station_elev = self.getValues(self.nc_pl_sur, 'height')

        svf = self.get_sky_view()
        slope = self.get_slope()
        aspect = self.get_aspect()

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            zenith = solar_zenith(lat=lat[n], lon=lon[n], time=py_time)
    
            diffuse, corrected_direct = elevation_corrected_sw(zenith=zenith,
                                                               grid_sw=sw[:,n],
                                                               lat=np.ones_like(sw[:,n]) * lat[n],
                                                               lon=np.ones_like(sw[:,n]) * lon[n],
                                                               time=py_time,
                                                               grid_elevation=np.ones_like(sw[:,n]) * grid_elev[n],
                                                               sub_elevation=np.ones_like(sw[:,n]) * station_elev[n])

            diffuse = diffuse * svf[n]  # apply sky-view factor

            if not np.all(slope == 0):
                azimuth = get_azimuth_fast(lat[n], lon[n], py_time)
                cos_i_sub = illumination_angle(zenith, azimuth, slope[n], aspect[n])
                cos_i_grid = np.cos(np.radians(zenith))
                corrected_direct = shading_corrected_sw_direct(corrected_direct, cos_i_sub, cos_i_grid)
                
                sensible_values_mask = np.where(cos_i_grid < 0.001, 0, 1) * np.where(corrected_direct > 1366, 0, 1)
                corrected_direct *= sensible_values_mask

            f = interp1d(nc_time, corrected_direct, kind='linear')
            self.rg.variables[vn_dir][:, n] = f(self.times_out_nc)

            f = interp1d(nc_time, diffuse, kind='linear')
            self.rg.variables[vn_diff][:, n] = f(self.times_out_nc)

    def LW_Wm2_sur(self):
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
        time_in = self.input_times_in_output_units(self.nc_sf)
        values  = self.getValues(self.nc_sf, 'strd')  # self.nc_sf.variables['strd'][:]/3600/self.interval_in #[w m-2 s-1]
        values  = values / (self.interval_in)  # [J m-2] -> [w m-2]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)

    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        # add variable to ncdf file
        vn = 'LW_topo'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'TOPOscale-corrected thermal radiation downwards ERA-5'
        var.standard_name = 'surface_downwelling_longwave_flux'
        var.units     = 'W m-2'

        # interpolate station by station
        time_in = self.input_times_in_output_units(self.nc_sf)
        t_sub = self.getValues(self.nc_pl_sur, 't')  # [K]
        rh_sub = self.getValues(self.nc_pl_sur, 'r')  # [%]
        t_grid = self.getValues(self.nc_sa, 't2m')  # [K]
        dewp_grid = self.getValues(self.nc_sa, 'd2m')  # [K]
        lw_grid  = self.getValues(self.nc_sf, 'strd') / self.interval_in  # [w m-2 s-1]

        rh_grid = self._rh()(t_grid, dewp_grid)

        lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)
        
        svf = self.get_sky_view()

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            data = lw_sub[:, n] * svf[n]
            f = interp1d(time_in, data, kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) 

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg]
        https://crudata.uea.ac.uk/cru/pubs/thesis/2007-willett/2INTRO.pdf
        '''

        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.input_times_in_output_units(self.nc_sa)
        values  = self.getValues(self.nc_sa, 'd2m')  # self.nc_sa.variables['d2m'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            dewp[:, n] = series_interpolate(self.times_out_nc,
                                            time_in, values[:, n] - 273.15)

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


    def input_times_in_output_units(self, ncf):
        """
        Convert time in input data to time in Globsim.
        """
        raw = ncf['time'][:].astype(np.int64)
        time = nc.num2date(raw, units=ncf['time'].units, calendar=ncf['time'].calendar)
        converted = nc.date2num(time, units=self.scaled_t_units, calendar=self.t_cal)
        return converted.astype(np.int64)