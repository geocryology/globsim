#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import netCDF4 as nc

from pathlib import Path

from globsim.common_utils import series_interpolate
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
from globsim.scale.scalenames import ScaleNames as SN
import globsim.scale.kernel_templates as kt
import globsim.dreamit as dreamit

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
    NAME = "JRA-generic"
    REANALYSIS = "jra-generic"

    VARNAMES = {}
    CONVERTERS = {}
    

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file names
        self.nc_pl_sur = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_pl_{self.list_name}_surface.nc'), 'r')
        self.nc_pl = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_pl_{self.list_name}.nc'), 'r')
        self.nc_sa = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_sf_{self.list_name}.nc'), 'r')
        try:
            self.nc_to = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_to_{self.list_name}.nc'), 'r')
        except AttributeError:
            logger.error("Missing invariant ('*_to') file. Some scaling kernels may fail. ")

        # Check data integrity
        _check_timestep_length(self.nc_pl_sur.variables['time'], "pl_sur")
        _check_timestep_length(self.nc_pl.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl_sur.variables['time'], par['time_step'])
        
        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.scaled_tstep_h,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.scaled_t_cal)
    
    def AIRT_DReaMIT(self):
        """
        Air temperature derived from surface data, pressure level data, and
        dynamically-computed inversion metrics as shown by the method DReaMIT
        """
        kt.AIRT_DReaMIT(self.rg)
        
        # Time netCDF file
        nc_time = self.nc_sa.variables['time']
        time_in = self.get_values("sa", SN.time)
       
        hypsometry = self.get_hypsometry()
        list_params = dreamit.get_model_params('jra3qg')
        time_frac_year = dreamit.time_frac_year(nc_time)
        pl_height = self.get_values("pl_sur", SN.elevation)

        for siteslist_ix, interp_ix in self.iterate_stations():
            # get T from pressure level
            T_pl_in = self.get_station_values("pl", "Temperature", interp_ix, preserve_dims=True)
            # get height from pressure level
            h_pl_in = self.get_station_values("pl", "Geopotential height", interp_ix, preserve_dims=True)
            # get station temperature from pressure level surface
            T_pl_surface_in = self.get_station_values("pl_sur", "Temperature", interp_ix, preserve_dims=True)
            # get station height from pressure level surface
            h_pl_surface_in = pl_height[interp_ix: interp_ix+1]  
            # get grid height from to
            grid_elev_in = self.get_station_values("to", "Geopotential", interp_ix, preserve_dims=True)
            # get T from surface
            T_sur_in = self.get_station_values("sa", "Temperature", interp_ix, preserve_dims=True, units="degree_C")
            
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
                                                            hyps=hypsometry[siteslist_ix],
                                                            params=list_params)

            self.rg.variables['z_top_inversion_m'][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                                     time_in, z_top_inversion_m)
            self.rg.variables['T_lapse_grid_C'][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                                     time_in, T_lapse_grid_C) - 273.15
            self.rg.variables['T_lapse_station_C'][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                                     time_in, T_lapse_station_C) - 273.15
            self.rg.variables['lapse_Cperm'][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                                     time_in, lapse_Cperm)
            self.rg.variables['AIRT_DReaMIT_C'][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                                     time_in, AIRT_DReaMIT_C) - 273.15
        self.rg.variables['beta_t_C'][:] = series_interpolate(self.times_out_nc,
                                                              time_in, beta_t_C[:])