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
# (6) Get all variables which are needed, and saved in NetCDF files
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
import logging
import netCDF4 as nc
import numpy as np
import pytz
import warnings

from os import path, makedirs
from math import atan2, pi
from pysolar.solar import get_azimuth_fast

from globsim.common_utils import series_interpolate
from globsim.scale.toposcale import lw_down_toposcale, solar_zenith, elevation_corrected_sw, illumination_angle, shading_corrected_sw_direct
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
import globsim.scale.kernel_templates as kt
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
        self.nc_to = self.nc_sc  # alias 

        # Check data integrity
        _check_timestep_length(self.nc_pl_sur.variables['time'], "pl_sur")
        _check_timestep_length(self.nc_pl.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        # self.nc_sc = nc.Dataset(path.join(self.intpdir,
        #                                  'merra2_to_' +
        #                        self.list_name + '.nc'), 'r')

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
        
    def get_grid_elevation_m(self, station_ix=None, preserve_dims=False):
        return self.get_station_values("to", "PHIS", station_ix, preserve_dims=preserve_dims) / const.G

    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """

        if not path.isdir(path.dirname(self.output_file)):
            makedirs(path.dirname(self.outfile))

        self.set_valid_stations(self.nc_pl_sur)
        valid_indices = self.valid_stations['nc_index']
        self.rg = new_scaled_netcdf(self.output_file,
                                    self.nc_pl_sur,
                                    self.times_out_nc,
                                    t_unit=self.scaled_t_units,
                                    valid_indices=valid_indices,
                                    station_names=self.valid_stations['station_name_scale'])
        
        # add grid elevation to netCDF
        self.add_grid_elevation(self.rg, self.nc_sc['PHIS'][0, valid_indices.values] / const.G)
        
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
        vn = kt.PRESS_Pa_pl(self.rg, self.NAME)
    
        time_in = self.input_times_in_output_units(self.nc_pl_sur)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("pl_sur", "air_pressure", interp_ix, units="Pa") 
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values) 

    def AIRT_C_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """
        vn = kt.AIRT_C_pl(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_pl_sur)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("pl_sur", "T", interp_ix, units="degree_C")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values)

    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """
        vn = kt.AIRT_C_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sa)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sa", "T2M", interp_ix, units="degree_C")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)
    
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

        time_in = self.input_times_in_output_units(self.nc_pl)
        values  = Delta_T_c
        
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            var[:, n] = series_interpolate(self.times_out_nc, 
                                           time_in,
                                           values[:, n])

    def RH_per_pl(self):
        """
        Relative Humdity derived from pressure level data, exclusively.Clipped to
        range [0.1,99.9].
        """
        vn = kt.RH_per_pl(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_pl_sur)

        for siteslist_ix, interp_ix in self.iterate_stations():
            rh  = self.get_station_values("pl_sur", "RH", interp_ix, units="percent")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, rh)

    def RH_per_sur(self):
        """
        Relative Humdity derived from surface data, exclusively.Clipped to
        range [0.1,99.9].
        """
        vn = kt.RH_per_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            d2m_values  = self.get_station_values("sf", "T2MDEW", interp_ix, units="degree_C")
            t2m_values = self.get_station_values("sa", "T2M", interp_ix, units="degree_C")
            rh_values = self._rh()(t2m_values, d2m_values).clip(min=0.1, max=99.9)
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, rh_values )

    def WIND_sur(self):
        """
        Wind speed and direction at 10 metre derived from surface data,
        exclusively.
        """
        vn_u, vn_v, vn_spd, vn_dir = kt.WIND_sur(self.rg, self.NAME)

        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)
        
        time_in = self.input_times_in_output_units(self.nc_sa)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values_u  = self.get_station_values('sa', 'U10M', interp_ix, units="m s-1")
            values_v  = self.get_station_values('sa', 'V10M', interp_ix, units="m s-1")
            U[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values_u)
            V[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values_v)

        self.rg.variables[vn_u][:, :] = U
        self.rg.variables[vn_v][:, :] = V

        WS = np.sqrt(np.power(V,2) + np.power(U,2))
        self.rg.variables[vn_spd][:, :] = WS

        WD = 90 - (np.arctan2(V, U) * (180 / np.pi)) + 180
        WD = np.mod(WD, 360)
        self.rg.variables[vn_dir][:, :] = WD

    def SW_Wm2_sur(self):
        """
        solar radiation downwards derived from surface data, exclusively.
        """
        vn = kt.SW_Wm2_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "SWGDN", interp_ix)
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def SW_Wm2_topo(self):
        """
        Short-wave downwelling radiation corrected using a modified version of TOPOscale.
        Partitions into direct and diffuse
        """
        vn_dir, vn_diff, vn_glob = kt.SW_Wm2_topo(self.rg, self.NAME)
        
        # interpolate station by station
        nc_time = self.input_times_in_output_units(self.nc_sf)
        py_time = nc.num2date(nc_time[:], self.scaled_t_units, self.t_cal, only_use_cftime_datetimes=False)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        
        lat = self.get_values('pl_sur', 'latitude')
        lon = self.get_values('pl_sur', 'longitude')
        
        svf = self.get_sky_view()
        slope = self.get_slope()
        aspect = self.get_aspect()
        
        grid_elev = self.get_values('to', 'PHIS', (0, slice(None,None,1))) / const.G  # z has 2 dimensions from the scaling step
        station_elev = self.get_values("pl_sur","height")
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            zenith = solar_zenith(lat=lat[interp_ix], lon=lon[interp_ix], time=py_time)
            sw = self.get_station_values('sf', 'SWGDN', interp_ix, preserve_dims=False)

            diffuse, corrected_direct = elevation_corrected_sw(zenith=zenith,
                                                               grid_sw=sw,
                                                               lat=np.ones_like(sw) * lat[interp_ix],
                                                               lon=np.ones_like(sw) * lon[interp_ix],
                                                               time=py_time,
                                                               grid_elevation=np.ones_like(sw) * grid_elev[interp_ix],
                                                               sub_elevation=np.ones_like(sw) * station_elev[interp_ix])

            diffuse = diffuse * svf[siteslist_ix]  # apply sky-view factor

            if slope[siteslist_ix] != 0:
                azimuth = get_azimuth_fast(lat[interp_ix], lon[interp_ix], py_time)
                cos_i_sub = illumination_angle(zenith, azimuth, slope[siteslist_ix], aspect[siteslist_ix])
                cos_i_grid = np.cos(np.radians(zenith))
                corrected_direct = shading_corrected_sw_direct(corrected_direct, cos_i_sub, cos_i_grid)
                
                sensible_values_mask = np.where(cos_i_grid < 0.001, 0, 1) * np.where(corrected_direct > 1366, 0, 1)
                corrected_direct *= sensible_values_mask

            global_sw = diffuse + corrected_direct

            self.rg.variables[vn_dir][:, siteslist_ix] = series_interpolate(self.times_out_nc, nc_time, corrected_direct)
            self.rg.variables[vn_diff][:, siteslist_ix] = series_interpolate(self.times_out_nc, nc_time, diffuse)
            self.rg.variables[vn_glob][:, siteslist_ix] = series_interpolate(self.times_out_nc, nc_time, global_sw)

    def LW_Wm2_sur(self):
        """
        Long-wave radiation downwards derived from surface data, exclusively.
        """
        vn = kt.LW_Wm2_sur(self.rg, self.NAME) 

        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "LWGDN", interp_ix)
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        vn = kt.LW_Wm2_topo(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)
        svf = self.get_sky_view()

        for siteslist_ix, interp_ix in self.iterate_stations():
            """ I'm cutting corners here by repeating variables with longer timesteps. This should be improved [NB]"""
            t_sub = self.get_station_values("pl_sur", 'T', interp_ix)
            t_sub = np.repeat(t_sub, 6, axis=0)
            rh_sub = self.get_station_values("pl_sur", 'RH', interp_ix)  # [%]
            rh_sub = np.repeat(rh_sub, 6, axis=0) * 100  # [%]
            t_grid = self.get_station_values("sa", 'T2M', interp_ix)  # [K]
            dewp_grid = self.get_station_values("sf", 'T2MDEW', interp_ix)  # [K]
            rh_grid = self._rh()(t_grid - 273.15, dewp_grid - 273.15) 
            lw_grid  = self.get_station_values("sf", 'LWGDN', interp_ix) # [w m-2 s-1]
            lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)
            values = lw_sub * svf[siteslist_ix]
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def PREC_mm_sur(self):
        """
        Precipitation derived from surface data, exclusively.
        Convert units: kg/m2/s to kg/m2/s
        1 kg/m2 = 1mm
        """
        vn  = kt.PREC_mm_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "PRECTOT", interp_ix) * self.scf
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def PRECCORR_mm_sur(self):
        """
        Corrected Precipitation derived from surface data, exclusively.
        Convert units: kg/m2/s to mm/time_step (hours)
        1 kg/m2 = 1mm
        """
        vn  = kt.PRECCORR_mm_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "PRECTOTCORR", interp_ix) * self.scf
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg] derived from surface data, exclusively.
        '''
        vn = kt.SH_kgkg_sur(self.rg, self.NAME) 
        
        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sa", "QV2M", interp_ix, units="kg kg-1")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)
    
    