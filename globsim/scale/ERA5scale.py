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
from globsim.meteorology import spec_hum_kgkg, relhu_approx_lawrence
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.toposcale import lw_down_toposcale, elevation_corrected_sw, sw_partition, solar_zenith, sw_toa, shading_corrected_sw_direct, illumination_angle
from globsim.scale.GenericScale import GenericScale, _check_timestep_length

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

    def __init__(self, sfile, era5type='reanalysis'):
        super().__init__(sfile)
        par = self.par

        if era5type == 'reanalysis':
            self.ens = False
            self.src = 'era5'
        elif era5type == 'ensemble_members':
            self.ens = True
            self.src = 'era5_ens'

        # input file handles
        self.nc_pl = nc.Dataset(Path(self.intpdir, f"{self.src}_pl_{self.list_name}_surface.nc"), 'r')
        self.nc_sa = nc.Dataset(Path(self.intpdir, f'{self.src}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.intpdir, f'{self.src}_sf_{self.list_name}.nc'), 'r')
        self.nc_to = nc.Dataset(Path(self.intpdir, f'{self.src}_to_{self.list_name}.nc'), 'r')

        # Check data integrity
        _check_timestep_length(self.nc_pl.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        if not (self.nc_sf['time'].shape == self.nc_pl['time'].shape == self.nc_sa['time'].shape):
            logger.critical(f"Different number of time steps in input data: sf ({self.nc_sf['time'].shape[0]}) pl ({self.nc_pl['time'].shape[0]}) sa ({self.nc_sa['time'].shape[0]})")

        try:
            self.nc_to.set_auto_scale(True)
        except Exception:
            logger.error("Could not set autoscaling! Values using geopotential (era_to) may be incorrect")

        self.nstation = len(self.nc_to.variables['station'][:])

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl.variables['time'], par['time_step'])
        self.scaled_t_units = 'seconds since 1900-01-01 00:00:0.0'

        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.time_step,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.t_cal)

    def getValues(self, ncf, varStr, ni=10):
        if self.ens:
            values = ncf.variables[varStr][:, ni, :]
        else:
            values = ncf.variables[varStr][:]

        return values

    def indProcess(self, ni=10):
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                logger.info(f"running scaling kernel: '{kernel_name}'")
                getattr(self, kernel_name)(ni)
            else:
                logger.error(f"Missing kernel {kernel_name}")

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

        logger.info(f"Created scaled output file {self.output_file}")
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
                                                             time_in * 3600,
                                                             values[:, n]) * 100

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
        var.units     = 'mm s-1'
        var.standard_name = 'precipitation_flux'

        # interpolation scale factor
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)

        # total prec [mm] in 1 second
        values  = self.getValues(self.nc_sf, 'tp', ni)  # self.nc_sf.variables['tp'][:]*1000/self.interval_in
        # https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
        # We assume interval_in is in hours (1 for regular, 3 for ensemble)
        values  = values * 1000 / (self.interval_in)  # [m] -> [mm s-1]

        # interpolate station by station
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.scf

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

        RH = relhu_approx_lawrence(self.rg.variables['AIRT_sur'][:, :], dewp[:, :])
        self.rg.variables[vn][:, :] = RH.clip(min=0.1, max=99.9)

    def RH_per_pl(self, ni=10):
        """
        Relative humdity derived from pressure-level. 
        """
        import pdb;pdb.set_trace()
        # temporary variable,  interpolate station by station
        rh = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = self.getValues(self.nc_pl, 'r', ni)  
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            rh[:, n] = series_interpolate(self.times_out_nc,
                                            time_in * 3600, values[:, n])

        # add variable to ncdf file
        vn = 'RH_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity ERA-5 pressure-levels'
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        self.rg.variables[vn][:, :] = rh

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

        WS = np.sqrt(np.power(V,2) + np.power(U,2))
        self.rg.variables[vn_spd][:, :] = WS

        # Wind direction: 
        # U - easterly component  V - northerly component
        # convert to "clockwise from north" from "anti-clockwise from x-axis" : 90 - angle
        # convert "from direction" : + 180
        WD = 90 - (np.arctan2(V, U) * (180 / np.pi)) + 180
        self.rg.variables[vn_dir][:, :] = WD

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
        values  = values / (self.interval_in)  # [J m-2] -> [w m-2]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)

    def SW_Wm2_topo(self, ni=10):
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
        nc_time = self.nc_sf.variables['time']
        py_time = nc.num2date(nc_time[:], nc_time.units, nc_time.calendar, only_use_cftime_datetimes=False)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        lat = self.getValues(self.nc_pl, 'latitude', ni)
        lon = self.getValues(self.nc_pl, 'longitude', ni)
        sw = self.getValues(self.nc_sf, 'ssrd', ni) / self.interval_in  # [J m-2] --> [W m-2]
        grid_elev = self.getValues(self.nc_to, 'z', ni)[0,:] / 9.80665  # z has 2 dimensions from the scaling step
        station_elev = self.getValues(self.nc_pl, 'height', ni)

        svf = self.get_sky_view()
        slope = self.get_slope()
        aspect = self.get_aspect()


        interpolation_time = nc_time[:].astype(np.int64)
        
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

            f = interp1d(interpolation_time * 3600, corrected_direct, kind='linear')
            self.rg.variables[vn_dir][:, n] = f(self.times_out_nc)

            f = interp1d(interpolation_time * 3600, diffuse, kind='linear')
            self.rg.variables[vn_diff][:, n] = f(self.times_out_nc)

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
        values  = values / (self.interval_in)  # [J m-2] -> [w m-2]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)

    def LW_Wm2_topo(self, ni=10):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        # add variable to ncdf file
        vn = 'LW_topo'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'TOPOscale-corrected thermal radiation downwards ERA-5'
        var.standard_name = 'surface_downwelling_longwave_flux'
        var.units     = 'W m-2'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        t_sub = self.getValues(self.nc_pl, 't', ni)  # [K]
        rh_sub = self.getValues(self.nc_pl, 'r', ni)  # [%]
        t_grid = self.getValues(self.nc_sa, 't2m', ni)  # [K]
        dewp_grid = self.getValues(self.nc_sa, 'd2m', ni)  # [K]
        lw_grid  = self.getValues(self.nc_sf, 'strd', ni) / self.interval_in  # [w m-2 s-1]
        rh_grid = relhu_approx_lawrence(t_grid, dewp_grid)

        lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)
        
        svf = self.get_sky_view()

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            data = lw_sub[:, n] * svf[n]
            f = interp1d(time_in * 3600, data, kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) 

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
