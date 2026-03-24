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
                   SN.pressure:      "air_pressure",
                   SN.rh:            "RH"},
        "to":     {SN.time:        "time",
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

