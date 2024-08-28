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
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4'
#        '/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4')

# 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4'
#        '/2016/01/MERRA2_400.inst3_3d_asm_Np.20160201.nc4')

# 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4'
#         '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')

# 2d,1-hourly, single-level, full horizontal resolution, Surface Flux Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXFLX.5.12.4'
#        '/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4')

# 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T3NPRAD.5.12.4'
#       '/2016/01/MERRA2_400.tavg3_3d_rad_Np.20160102.nc4')

# 2d 1-Hourly,Time-Averaged,Single-Level,Assimilation,Single-Level Diagnostics V5.12.4
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXSLV.5.12.4'
# '       '/2016/01/MERRA2_400.tavg1_2d_slv_Nx.20160101.nc4')

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

import logging
import netCDF4 as nc
import numpy as np
import pytz
import warnings

from os import path, makedirs
from math import atan2, pi
from pysolar.solar import get_azimuth_fast
from scipy.interpolate import interp1d

from globsim.common_utils import series_interpolate
from globsim.scale.toposcale import lw_down_toposcale, solar_zenith, elevation_corrected_sw, illumination_angle, shading_corrected_sw_direct
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
import globsim.constants as const
import globsim.redcapp as redcapp

warnings.filterwarnings("ignore", category=UserWarning, module='netCDF4')

logger = logging.getLogger('globsim.scale')


class MERRAscale(GenericScale):
    """
    Class for MERRA data that has methods for scaling station data to
    better resemble near-surface fluxes.

    Processing kernels have names in UPPER CASE.

    Args:
        sfile: Full path to a Globsim Scaling Parameter file.

    Example:
        MERRAd = MERRAscale(sfile)
        MERRAd.process()
    """
    NAME = "MERRA-2"

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file names
        self.nc_pl_sur = nc.Dataset(path.join(self.intpdir, f'merra2_pl_{self.list_name}_surface.nc'), 'r')
        self.nc_pl = nc.Dataset(path.join(self.intpdir, f'merra2_pl_{self.list_name}.nc'), 'r')
        self.nc_sa = nc.Dataset(path.join(self.intpdir, f'merra2_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(path.join(self.intpdir, f'merra2_sf_{self.list_name}.nc'), 'r')
        self.nc_sc = nc.Dataset(path.join(self.intpdir, f'merra2_sc_{self.list_name}.nc'), 'r')

        # Check data integrity
        _check_timestep_length(self.nc_pl_sur.variables['time'], "pl_sur")
        _check_timestep_length(self.nc_pl.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        # self.nc_sc = nc.Dataset(path.join(self.intpdir,
        #                                  'merra2_to_' +
        #                        self.list_name + '.nc'), 'r')
        self.nstation = len(self.nc_pl_sur.variables['station'][:])

        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(par, 'merra2')

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl_sur.variables['time'], par['time_step'])
        self.scaled_t_units = 'seconds since 1980-01-01 00:00:00'

        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.time_step,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.t_cal)

    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """

        if not path.isdir(path.dirname(self.output_file)):
            makedirs(path.dirname(self.outfile))

        self.rg = new_scaled_netcdf(self.output_file,
                                    self.nc_pl_sur,
                                    self.times_out_nc,
                                    t_unit=self.scaled_t_units)
        
        # add grid elevation to netCDF
        self.add_grid_elevation(self.rg, self.nc_sc['PHIS'][0, :] / const.G)
        
        # add station names to netcdf
        # first convert to character array
        names_out = nc.stringtochar(np.array(self.stations['station_name'], 'S32'))

        # create space in the netcdf
        _        = self.rg.createDimension('name_strlen', 32)
        st           = self.rg.createVariable('station_name', "S1",
                                              ('station', 'name_strlen'))
        st.standard_name = 'platform_name'
        st.units     = ''

        # add data
        st[:] = names_out

        # iterate through kernels and start process
        self.run_kernels()

        # close netCDF files
        self.rg.close()
        self.nc_pl_sur.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_pl.close()

    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """
        # add variable to ncdf file
        vn = 'PRESS_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_pressure MERRA-2 pressure levels only'
        var.units     = 'Pa'
        var.standard_name = 'surface_air_pressure'

        # interpolate station by station
        time_in = self.nc_pl_sur.variables['time'][:].astype(np.int64)
        values  = self.nc_pl_sur.variables['air_pressure'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            # scale from hPa to Pa
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
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
        time_in = self.nc_pl_sur.variables['time'][:].astype(np.int64)
        values  = self.nc_pl_sur.variables['T'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
                                                             values[:, n] - 273.15)

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
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['T2M'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
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
        T_sa  = self.nc_sa.variables['T2M'][::6, :]  # every 6th hour to match pl data 
        # TODO: interpolate pl data to 24h instead of downsampling sa data.

        # get grid surface elevation from geopotential  (Cao: elev. @ coarse-scale topography)
        h_sur = self.nc_sc['PHIS'][0, :] / const.G  # remove time dimension. Time-invariant.
        
        # get pressure-level temperatures
        airT_pl = self.nc_pl.variables['T'][:]
        # get pressure-level elevations from geopotential
        elevation = self.nc_pl['H'][:]  # 

        Delta_T_c = redcapp.delta_T_c(T_sa=T_sa, 
                                      airT_pl=airT_pl, 
                                      elevation=elevation,
                                      h_sur=h_sur)  

        time_in = time_in = self.nc_pl.variables['time'][:].astype(np.int64)
        values  = Delta_T_c
        
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            var[:, n] = series_interpolate(self.times_out_nc, 
                                           time_in * 3600,
                                           values[:, n])

    def RH_per_pl(self):
        """
        Relative Humdity derived from pressure level data, exclusively.Clipped to
        range [0.1,99.9].
        """

        # temporary variable,  interpolate station by station
        time_in = self.nc_pl_sur.variables['time'][:].astype(np.int64)
        values = self.nc_pl_sur.variables['RH'][:]

        # add variable to ncdf file
        vn = 'RH_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity {} surface only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            rh = series_interpolate(self.times_out_nc, time_in * 3600,
                                    values[:, n])
            rh *= 100  # Convert to % 
            self.rg.variables[vn][:, n] = rh

    def RH_per_sur(self):
        """
        Relative Humdity derived from surface data, exclusively.Clipped to
        range [0.1,99.9].
        """

        # temporary variable,  interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        dewp = self.nc_sf.variables['T2MDEW'][:]
        t = self.nc_sa.variables['T2M'][:]
        relhu = self._rh()(t, dewp)

        # add variable to ncdf file
        vn = 'RH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Relative humidity {} surface only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            rh = series_interpolate(self.times_out_nc, time_in * 3600,
                                    relhu[:, n])
            self.rg.variables[vn][:, n] = rh

    def WIND_sur(self):
        """
        Wind speed and direction at 10 metre derived from surface data,
        exclusively.
        """

        # temporary variable, interpolate station by station
        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['U10M'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            U[:, n] = series_interpolate(self.times_out_nc,
                                         time_in * 3600,
                                         values[:, n])

        # temporary variable, interpolate station by station
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['V10M'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            V[:, n] = series_interpolate(self.times_out_nc,
                                         time_in * 3600,
                                         values[:, n])

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

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            WS = np.sqrt(np.power(V[:, n], 2) + np.power(U[:, n], 2))
            WD = [atan2(V[i, n], U[i, n]) * (180 / pi) + 180 for i in np.arange(V.shape[0])]
            WD = np.mod(WD, 360)
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
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['SWGDN'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
                                                             values[:, n])

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
        nc_time = self.nc_sf.variables['time']
        py_time = nc.num2date(nc_time[:], nc_time.units, nc_time.calendar, only_use_cftime_datetimes=False)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        lat = self.nc_pl_sur['latitude'][:]
        lon = self.nc_pl_sur['longitude'][:]
        sw = self.nc_sf['SWGDN'][:]  # [W m-2]
        grid_elev = self.nc_sc["PHIS"][0, :] / const.G  # [m]
        station_elev = self.nc_pl_sur["height"][:]  # [m]

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
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['LWGDN'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
                                                             values[:, n])

    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        # add variable to ncdf file
        vn = 'LW_topo'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'TOPOscale-corrected thermal radiation downwards ERA-5'
        var.standard_name = 'surface_downwelling_longwave_flux'
        var.units     = 'W m-2'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        """ I'm cutting corners here by repeating variables with longer timesteps. This should be improved [NB]"""
        t_sub = np.repeat(self.nc_pl_sur['T'][:], 6, axis=0)  # [K]
        rh_sub = np.repeat(self.nc_pl_sur['RH'][:], 6, axis=0) * 100  # [%]
        t_grid = self.nc_sa['T2M'][:]  # [K]
        dewp_grid = self.nc_sf['T2MDEW'][:]  # [K]
        lw_grid  = self.nc_sf["LWGDN"]  # [w m-2 s-1]
        rh_grid = self._rh()(t_grid, dewp_grid)

        lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)
        
        svf = self.get_sky_view()

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            data = lw_sub[:, n] * svf[n]
            f = interp1d(time_in * 3600, data, kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) 

    def PREC_mm_sur(self):
        """
        Precipitation derived from surface data, exclusively.
        Convert units: kg/m2/s to kg/m2/s
        1 kg/m2 = 1mm
        """

        # add variable to ncdf file
        vn = 'PREC_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Total precipitation {} surface only'.format(self.NAME)
        var.units     = 'kg m-2 s-1'
        var.comment = "units [kg m-2 s-1] corresponds to [mm/s] for water (density 1000 [kg m-3])"
        var.standard_name = 'precipitation_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['PRECTOT'][:]  # mm s-1 aka kg/m2/s
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in * 3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.scf

    def PRECCORR_mm_sur(self):
        """
        Corrected Precipitation derived from surface data, exclusively.
        Convert units: kg/m2/s to mm/time_step (hours)
        1 kg/m2 = 1mm
        """

        # add variable to ncdf file
        vn = 'PRECCORR_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Corrected Total precipitation {} surface only'.format(self.NAME)
        var.units     = 'kg m-2 s-1'
        var.comment = "units [kg m-2 s-1] corresponds to [mm/s] for water (density 1000 [kg m-3])"
        var.standard_name = 'precipitation_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        values  = self.nc_sf.variables['PRECTOTCORR'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
                                                             values[:, n]) * self.scf

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg] derived from surface data, exclusively.
        '''

        # add variable to ncdf file
        vn = 'SH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Specific humidity {} surface only'.format(self.NAME)
        var.units     = '1'
        var.standard_name = 'specific_humidity'

        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)
        values  = self.nc_sa.variables['QV2M'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in * 3600,
                                                             values[:, n])
'''
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
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_longwave_flux'

        # compute
        for i in range(0, len(self.rg.variables['RH_sur'][:])):
            for n, s in enumerate(self.rg.variables['station'][:].tolist()):
                LW = LW_downward(self.rg.variables['RH_sur'][i, n],
                                 self.rg.variables['AIRT_sur'][i, n] + 273.15,
                                 N[n])

                self.rg.variables[vn][i, n] = LW

'''