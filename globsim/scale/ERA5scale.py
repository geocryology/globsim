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

from cfunits import Units
from pathlib import Path
from pysolar.solar import get_azimuth_fast

from globsim.common_utils import series_interpolate
from globsim.meteorology import spec_hum_kgkg
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.toposcale import lw_down_toposcale, elevation_corrected_sw, solar_zenith, shading_corrected_sw_direct, illumination_angle
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
import globsim.scale.kernel_templates as kt
import globsim.constants as const
import globsim.redcapp as redcapp
import globsim.dreamit as dreamit

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
    REANALYSIS = 'era5'
    NAME = "ERA-5"
    SCALING = {"sf": {"ssrd": (1/3600, 0),  # [J m-2] -> [W m-2]
                      "strd": (1/3600, 0),  # [J m-2] -> [W m-2]
                      "tp": (1000/3600, 0)},  # [m] -> [mm s-1]
               "sa": {},
               "pl": {},
               "to": {},
               "pl_sur": {}}
    

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file handles
        self.nc_pl_sur = nc.Dataset(Path(self.intpdir, f"{self.REANALYSIS}_pl_{self.list_name}_surface.nc"), 'r')
        self.nc_pl = nc.Dataset(Path(self.intpdir, f"{self.REANALYSIS}_pl_{self.list_name}.nc"), 'r')
        self.nc_sa = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_sf_{self.list_name}.nc'), 'r')
        self.nc_to = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_to_{self.list_name}.nc'), 'r')

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

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl_sur.variables['time'], par['time_step'])
        self.scaled_t_units = 'seconds since 1900-01-01 00:00:0.0'

        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.time_step,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.t_cal)
        
    def get_grid_elevation_m(self, station_ix=None, preserve_dims=False):
        return self.get_station_values("to", "z", station_ix, preserve_dims=preserve_dims) / const.G

    def get_pressure_level_height_m(self, station_ix, preserve_dims=False):
        return self.get_station_values("pl", "z", station_ix, 
                                       preserve_dims=preserve_dims) / const.G

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

        self.set_valid_stations(self.nc_pl_sur)
        # iterate thorugh kernels and start process

        self.output_file = self.getOutNCF(self.par, self.REANALYSIS)
        valid_indices = self.valid_stations['nc_index']
        self.rg = new_scaled_netcdf(self.output_file,
                                    self.nc_pl_sur, self.times_out_nc,
                                    t_unit=self.scaled_t_units,
                                    valid_indices=valid_indices,
                                    station_names=self.valid_stations['station_name_scale'])
        # add surface height
        self.add_grid_elevation(self.rg, self.getValues(self.nc_to, 'z')[0, valid_indices.values] / const.G)

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
            values  = self.get_station_values("pl_sur", "t", interp_ix, units="degree_C")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)
        
    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """
        vn = kt.AIRT_C_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sa)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sa", "t2m", interp_ix, units="degree_C")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values) - 273.15

    def AIRT_DReaMIT(self):
        """
        Air temperature derived from surface data, pressure level data, and
        dynamically-computed inversion metrics as shown by the method DReaMIT
        """
        kt.AIRT_DReaMIT(self.rg)

        nc_time = self.nc_sa.variables['time']
        time_in = self.input_times_in_output_units(self.nc_sa)

        hypsometry = self.get_hypsometry()
        list_params = dreamit.get_model_params('era5')
        time_frac_year = dreamit.time_frac_year(nc_time)
        pl_height = self.get_values("pl_sur", "height")

        for siteslist_ix, interp_ix in self.iterate_stations():
            # get T from pressure level
            T_pl_in = self.get_station_values("pl", "t", interp_ix, preserve_dims=True)
            # get height from pressure level
            h_pl_in = self.get_station_values("pl", "z", interp_ix, preserve_dims=True)
            # get station temperature from pressure level surface
            T_pl_surface_in = self.get_station_values("pl_sur", "t", interp_ix, preserve_dims=True)
            # get station height from pressure level surface
            h_pl_surface_in = pl_height[interp_ix: interp_ix+1]
            # get grid height from to
            grid_elev_in = self.get_station_values("to", "z", interp_ix, preserve_dims=True)
            # get T from surface
            T_sur_in = self.get_station_values("sa", "t2m", interp_ix, preserve_dims=True)            

            mtrcs = dreamit.dreamit_metrics(reanalysis=self.NAME,
                                            T_pl_in=T_pl_in,
                                            h_pl_in=h_pl_in,
                                            T_pl_surface_in=T_pl_surface_in,
                                            h_pl_surface_in=h_pl_surface_in,
                                            grid_elev_in=grid_elev_in)
            
            z_top_inversion_m, T_lapse_grid_C, T_lapse_station_C, lapse_Cperm = mtrcs



            AIRT_DReaMIT_C, beta_t_C = dreamit.dreamit_air_T(T_lapse_grid=T_lapse_grid_C,
                                                            T_lapse_station=T_lapse_station_C,
                                                            T_sur=T_sur_in,
                                                            time_frac_year=time_frac_year,
                                                            hyps=np.atleast_1d(hypsometry),
                                                            params=list_params)

            self.rg.variables['z_top_inversion_m'][:, siteslist_ix] = np.interp(self.times_out_nc,
                                                                     time_in, z_top_inversion_m)
            self.rg.variables['T_lapse_grid_C'][:, siteslist_ix] = np.interp(self.times_out_nc,
                                                                     time_in, T_lapse_grid_C) - 273.15
            self.rg.variables['T_lapse_station_C'][:, siteslist_ix] = np.interp(self.times_out_nc,
                                                                     time_in, T_lapse_station_C) - 273.15
            self.rg.variables['lapse_Cperm'][:, siteslist_ix] = np.interp(self.times_out_nc,
                                                                     time_in, lapse_Cperm)
            self.rg.variables['AIRT_DReaMIT_C'][:, siteslist_ix] = np.interp(self.times_out_nc,
                                                                     time_in, AIRT_DReaMIT_C) - 273.15
        self.rg.variables['beta_t_C'][:] = np.interp(self.times_out_nc,
                                                     time_in, beta_t_C[:])

    def AIRT_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP Cao et al. (2017) 10.5194/gmd-10-2905-2017
        """
        logger.warning(f"Globsim implementation of REDCAPP only provides Delta_T_c")

        # add variable to ncdf file
        var = redcapp.add_var_delta_T(self.rg)

        time_in = self.input_times_in_output_units(self.nc_sa)

        for siteslist_ix, interp_ix in self.iterate_stations():
            # get T from surface level
            T_sa  = self.get_station_values("sa", 't2m', interp_ix, preserve_dims=True)  
            # get grid surface elevation from geopotential  (Cao: elev. @ coarse-scale topography)
            h_sur = self.get_station_values("to", 'z', interp_ix, preserve_dims=False) / const.G  # [m]
            # get pressure-level temperatures
            airT_pl = self.get_station_values("pl", 't', interp_ix, preserve_dims=True)
            # get pressure-level elevations from geopotential
            elevation = self.get_station_values("pl", 'z', interp_ix, preserve_dims=True) / const.G  # [m]

            Delta_T_c = redcapp.delta_T_c(T_sa=T_sa, 
                                        airT_pl=airT_pl, 
                                        elevation=elevation,
                                        h_sur=h_sur)  

            values  = Delta_T_c[:, 0]
            var[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)                              

    def PREC_mm_sur(self):
        """
        Precipitation sum in mm for the time step given.
        """
        vn  = kt.PREC_mm_sur(self.rg, self.NAME)

        # interpolation scale factor
        time_in = self.input_times_in_output_units(self.nc_sf)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "tp", interp_ix)  
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values) * self.scf

    def RH_per_sur(self):
        """
        Relative humdity derived from surface data, exclusively. Clipped to
        range [0.1,99.9]. Kernel AIRT_C_sur must be run before.
        """
        vn = kt.RH_per_sur(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sa)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            d2m_values  = self.get_station_values("sa", "d2m", interp_ix, units="degree_C")
            t2m_values = self.get_station_values("sa", "t2m", interp_ix, units="degree_C")
            rh_values = self._rh()(t2m_values, d2m_values).clip(min=0.1, max=99.9)
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, rh_values )

    def RH_per_pl(self):
        """
        Relative humdity derived from pressure-level.
        """
        vn = kt.RH_per_pl(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_pl_sur)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("pl_sur", "r", interp_ix, units="percent")
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def WIND_sur(self):
        """
        Wind speed and direction temperature derived from surface data,
        exclusively.
        """
        vn_u, vn_v, vn_spd, vn_dir = kt.WIND_sur(self.rg, self.NAME)

        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)

        time_in = self.input_times_in_output_units(self.nc_sa)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values_u  = self.get_station_values("sa", "u10", interp_ix, units="m s-1")
            values_v  = self.get_station_values("sa", "v10", interp_ix, units="m s-1")
            U[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values_u)
            V[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values_v)            

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
        vn = kt.SW_Wm2_sur(self.rg, self.NAME)
        
        time_in = self.input_times_in_output_units(self.nc_sf)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "ssrd", interp_ix)
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)
    
    def SW_Wm2_topo(self):
        """
        Short-wave downwelling radiation corrected using a modified version of TOPOscale.
        Partitions into direct and diffuse
        """
        vn_dir, vn_diff, vn_glob = kt.SW_Wm2_topo(self.rg, self.NAME)

        nc_time = self.input_times_in_output_units(self.nc_sf)
        py_time = nc.num2date(nc_time[:], self.scaled_t_units, self.t_cal, only_use_cftime_datetimes=False)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        
        lat = self.get_values('pl_sur', 'latitude')
        lon = self.get_values('pl_sur', 'longitude')

        svf = self.get_sky_view()
        slope = self.get_slope()
        aspect = self.get_aspect()

        grid_elev = self.get_values('to', 'z', (0, slice(None,None,1))) / const.G  # z has 2 dimensions from the scaling step
        station_elev = self.get_values("pl_sur", "height")

        #for n, s in enumerate(self.rg.variables['station'][:].tolist()):
        for siteslist_ix, interp_ix in self.iterate_stations():
            zenith = solar_zenith(lat=lat[interp_ix], lon=lon[interp_ix], time=py_time)
            sw = self.get_station_values('sf', 'ssrd', interp_ix, preserve_dims=False)
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
        Long-wave downwelling radiation derived from surface data, exclusively.
        This kernel only interpolates in time.
        """
        vn = kt.LW_Wm2_sur(self.rg, self.NAME) 

        time_in = self.input_times_in_output_units(self.nc_sf)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", "strd", interp_ix)
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        vn = kt.LW_Wm2_topo(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)
        svf = self.get_sky_view()

        for siteslist_ix, interp_ix in self.iterate_stations():
            t_sub = self.get_station_values("pl_sur", 't', interp_ix, units="degree_K")  # [K]
            rh_sub = self.get_station_values("pl_sur", 'r', interp_ix, units="percent")  # [%]
            t_grid = self.get_station_values("sa", 't2m', interp_ix, units="degree_K")  # [K]
            dewp_grid = self.get_station_values("sa", 'd2m', interp_ix, units="degree_K")  # [K]
            rh_grid = self._rh()(t_grid - 273.15, dewp_grid - 273.15)
            lw_grid  = self.get_station_values("sf", 'strd', interp_ix)

            lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)

            values = lw_sub * svf[siteslist_ix]
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg]
        https://crudata.uea.ac.uk/cru/pubs/thesis/2007-willett/2INTRO.pdf
        '''
        vn = kt.SH_kgkg_sur(self.rg, self.NAME) 

        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.input_times_in_output_units(self.nc_sa)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values = self.get_station_values("sa", "d2m", interp_ix, units="degree_C")
            dewp[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

        SH = spec_hum_kgkg(dewp[:, :], self.rg.variables['PRESS_pl'][:, :])
        self.rg.variables[vn][:, :] = SH


