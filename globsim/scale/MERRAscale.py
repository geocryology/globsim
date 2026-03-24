#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import netCDF4 as nc
import numpy as np
import pytz
import warnings

from os import path, makedirs
from math import atan2, pi
from pysolar.solar import get_azimuth_fast
from cfunits import Units

from globsim.common_utils import series_interpolate
from globsim.scale.toposcale import lw_down_toposcale, solar_zenith, elevation_corrected_sw, illumination_angle, shading_corrected_sw_direct
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
import globsim.scale.kernel_templates as kt
import globsim.constants as const
import globsim.redcapp as redcapp
from globsim.scale.scalenames import ScaleNames as SN

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
    REANALYSIS = "merra2"
    VARNAMES = {
        "sa":     {SN.time:        "time",
                   SN.temperature: "T2M",
                   SN.specific_humidity: "QV2M",
                   SN.rh:          "T2M",
                   SN.u_wind:      "U10M",
                   SN.v_wind:      "V10M"},
        "sf":     {SN.time:        "time",
                   SN.dewpoint:      "T2MDEW",
                   SN.sw_down_flux:       "SWGDN",
                   SN.lw_down_flux:       "LWGDN",
                   SN.precipitation_rate: "PRECTOT"},
        "pl":     {SN.time:        "time",
                   SN.temperature:   "T",
                   SN.rh:            "RH",
                   SN.elevation:     "H"},
        "pl_sur": {SN.temperature:   "T",
                   SN.elevation:     "height",
                   SN.pressure:      "air_pressure",
                   SN.rh:            "RH"},
        "to":     {SN.time:          "time",
                   SN.geopotential:  "PHIS",
                   SN.elevation:     "PHIS"},
    }

    CONVERTERS = {("to", SN.elevation): "_geopotential_to_m",
                  ("sa", SN.rh): "_temp_to_rh",
                  }

    def _temp_to_rh(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        """Convert temperature to relative humidity using dewpoint."""
        t2m = Units.conform(data, Units(nc_var.units), Units("degree_C"))
        d2m = self.get_values("sf", SN.dewpoint, _slice, units="degree_C")
        rh = self._rh()(t2m, d2m).clip(min=0.1, max=99.9)
        
        return rh, "percent"
        
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

        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.scaled_tstep_h,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.scaled_t_cal)

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
        self.nc_sc.close()

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

