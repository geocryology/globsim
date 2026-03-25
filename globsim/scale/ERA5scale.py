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

from cfunits import Units
from pathlib import Path

from globsim.meteorology import spec_hum_kgkg
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
from globsim.scale.scalenames import ScaleNames as SN
import globsim.scale.kernel_templates as kt
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
    SCALING = {"sf": {}, 
               "sa": {},
               "pl": {},
               "to": {},
               "pl_sur": {}}
    
    VARNAMES = {
    "sa":     { SN.time:        "time",
                SN.longitude:     "longitude",
                SN.latitude:      "latitude",
                SN.temperature: "t2m",
                SN.dewpoint:    "d2m",
                SN.rh:           "d2m", # calculate from dewpoint and temperature
                SN.specific_humidity: "d2m", # calculate from dewpoint and temperature
                SN.u_wind:      "u10",
                SN.v_wind:      "v10",
                SN.ozone:       "tco3",
                SN.water_vapour:"tcwv"},
    "sf":     {SN.time:        "time",
               SN.longitude:     "longitude",
               SN.latitude:      "latitude",
               SN.sw_down_flux:       "ssrd",
               SN.sw_down_accumulated:"ssrd",
               SN.lw_down_flux:       "strd",
               SN.lw_down_accumulated: "strd",
               SN.precipitation_total: "tp",
               SN.precipitation_rate: "tp"},
    "pl":     {SN.time:        "time",
               SN.longitude:     "longitude",
               SN.latitude:      "latitude",
               SN.temperature:   "t",
               SN.rh:            "r",
               SN.geopotential:  "z",
               SN.elevation:     "z"},
    "pl_sur": {SN.time:        "time",   
               SN.longitude:     "longitude",
               SN.latitude:      "latitude",
               SN.temperature:   "t",
               SN.rh:            "r",
               SN.elevation:     "height",
               SN.pressure:      "air_pressure"},
    "to":     {SN.longitude:     "longitude",
               SN.latitude:      "latitude",
               SN.geopotential:  "z",
               SN.elevation:     "z"},
}
    
    CONVERTERS = {("sf", SN.sw_down_flux): "_radiation_to_flux",
                  ("sf", SN.lw_down_flux): "_radiation_to_flux",
                  ("sf", SN.precipitation_rate): "_precip_tot_to_flux",
                  ("sa", SN.rh): "_dewpoint_to_rh",
                  ("sa", SN.specific_humidity): "_dewpoint_to_sh",
                  ("to", SN.elevation): "_geopotential_to_m",
                  ("pl", SN.elevation): "_geopotential_to_m"
                  }

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file handles
        self.nc_pl_sur = nc.Dataset(Path(self.interp_dir, f"{self.REANALYSIS}_pl_{self.list_name}_surface.nc"), 'r')
        self.nc_pl = nc.Dataset(Path(self.interp_dir, f"{self.REANALYSIS}_pl_{self.list_name}.nc"), 'r')
        self.nc_sa = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_sf_{self.list_name}.nc'), 'r')
        self.nc_to = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_to_{self.list_name}.nc'), 'r')

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
        
        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.scaled_tstep_h,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.scaled_t_cal)

    def getValues(self, ncf, varStr):
        return ncf.variables[varStr][:]
    
    def _precip_tot_to_flux(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        """Convert total precipitation to precipitation rate by dividing by time step."""
        input_units = Units(nc_var.units)
        water_density = Units("kg m-3")
        converted_data = data / self.get_time_step("sf")
        converted_units = input_units * water_density / Units("s") 

        return converted_data, converted_units.units
    
    def _dewpoint_to_sh(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        dewp = Units.conform(data, Units(nc_var.units), Units("degree_C"))
        pressure = self.get_values("pl_sur", SN.pressure, _slice, units="Pa")
        sh = spec_hum_kgkg(dewp, pressure)
        
        return sh, "kg kg-1"

    def _dewpoint_to_rh(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        d2m_values = Units.conform(data, Units(nc_var.units), Units("degree_C"))
        t2m_values = self.get_values("sa", SN.temperature, _slice, units="degree_C")
        rh_values = self._rh()(t2m_values, d2m_values).clip(min=0.1, max=99.9)
        return rh_values, "percent"
    
    def _radiation_to_flux(self, data, nc_var, _slice):
        """Convert accumulated radiation to flux by dividing by time step."""
        input_units = Units(nc_var.units)
        converted_data = data / self.get_time_step("sf")
        converted_units = input_units / Units("s")

        return converted_data, converted_units.units
'''    
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
        pl_height = self.get_values("pl_sur", SN.elevation)

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
    
'''