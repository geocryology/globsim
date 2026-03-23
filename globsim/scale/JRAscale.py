#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import netCDF4 as nc
import numpy as np
import pytz

from cfunits import Units
from math import atan2, pi
from pathlib import Path
from pysolar.solar import get_azimuth_fast

from globsim.common_utils import series_interpolate
from globsim.scale.toposcale import (lw_down_toposcale, illumination_angle,
                                     shading_corrected_sw_direct, elevation_corrected_sw, 
                                     solar_zenith)
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale, _check_timestep_length
from globsim.scale.scalenames import ScaleNames as SN
import globsim.scale.kernel_templates as kt
import globsim.constants as const
import globsim.redcapp as redcapp
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

    
    VARNAMES = {
        "sa":     {SN.time:        "time",
                   SN.temperature: "Temperature",
                   SN.rh:          "Relative humidity",
                   SN.u_wind:      "u-component of wind",
                   SN.v_wind:      "v-component of wind",
                   SN.pressure:     "air_pressure",
                   SN.specific_humidity: "Specific humidity"},
        "sf":     {SN.time:          "time",
                   SN.sw_down:       "Downward solar radiation flux",
                   SN.lw_down:       "Downward longwave radiation flux",
                   SN.precipitation_total: "Total precipitation",
                   SN.precipitation_rate: "Total precipitation"},
        "pl":     {SN.time:          "time",
                   SN.temperature:   "Temperature",
                   SN.geopotential:  "Geopotential height"},
        "pl_sur": {SN.time:          "time",
                   SN.temperature:   "Temperature",
                   SN.rh:            "Relative humidity",
                   SN.pressure:      "air_pressure"},
        "to":     {SN.time:          "time",
                   SN.geopotential:  "Geopotential"},
    }

    CONVERTERS = {
        ("sf", SN.precipitation_rate): "_daily_precip_to_rate",
    }
    

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

        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(par, f'{self.REANALYSIS}')

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl_sur.variables['time'], par['time_step'])

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
        self.set_valid_stations(self.nc_pl_sur)
        valid_indices = self.valid_stations['nc_index']
        self.rg = new_scaled_netcdf(ncfile_out=self.output_file, 
                                    nc_interpol=self.nc_pl_sur,
                                    times_out=self.times_out_nc, 
                                    t_unit=self.nc_pl_sur['time'].units,
                                    valid_indices=valid_indices,
                                    station_names=self.valid_stations['station_name_scale'],)

        # add station names to netcdf
        elev = self.get_values("to", "Geopotential", (0, slice(None,None,1)))
        self.add_grid_elevation(self.rg, elev[valid_indices.values] / const.G)  # [m]

        # iterate through kernels and start process
        self.run_kernels()

        # close netCDF files
        self.rg.close()
        self.nc_pl_sur.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
    
    def AIRT_DReaMIT(self):
        """
        Air temperature derived from surface data, pressure level data, and
        dynamically-computed inversion metrics as shown by the method DReaMIT
        """
        kt.AIRT_DReaMIT(self.rg)
        
        # Time netCDF file
        nc_time = self.nc_sa.variables['time']
        time_in = self.get_values("sa", "time")
       
        hypsometry = self.get_hypsometry()
        list_params = dreamit.get_model_params('jra3qg')
        time_frac_year = dreamit.time_frac_year(nc_time)
        pl_height = self.get_values("pl_sur", "height")

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

    def AIRT_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP Cao et al. (2017) 10.5194/gmd-10-2905-2017
        """
        logger.warning(f"Globsim implementation of REDCAPP only provides Delta_T_c")
        
        # add variable to ncdf file
        var = redcapp.add_var_delta_T(self.rg)
        time_in = self.get_values("sa", "time")

        for siteslist_ix, interp_ix in self.iterate_stations():
            # get T from surface level
            T_sa  = self.get_station_values("sa", "Temperature", interp_ix, preserve_dims=True)
            # get grid surface elevation from geopotential  (Cao: elev. @ coarse-scale topography)
            h_sur = self.get_station_values("to", "Geopotential", interp_ix, preserve_dims=False) / const.G  # [m]
            # get pressure-level temperatures
            airT_pl = self.get_station_values("pl", "Temperature", interp_ix, preserve_dims=True)
            # get pressure-level elevations from geopotential
            elevation = self.get_station_values("pl", "Geopotential height", interp_ix, preserve_dims=True)  # JRA units say [gpm] but range [-10 to 4600] suggests its already [m]
            # [m] but we need to convert to [m] from [gpm] (geopotential meters)
   
            Delta_T_c = redcapp.delta_T_c(T_sa=T_sa, 
                                          airT_pl=airT_pl, 
                                          elevation=elevation,
                                          h_sur=h_sur) 
            
            values = Delta_T_c[:, 0]  # remove station dimension for interpolation
            var[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)                                                

    def SW_Wm2_topo(self):
        """
        Short-wave downwelling radiation corrected using a modified version of TOPOscale.
        Partitions into direct and diffuse
        """     
        vn_dir, vn_diff, vn_glob = kt.SW_Wm2_topo(self.rg, self.NAME)

        nc_time = self.nc_sf.variables['time']
        py_time = nc.num2date(nc_time[:], nc_time.units, nc_time.calendar,
                               only_use_cftime_datetimes=False)# , only_use_python_datetime=True)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        interpolation_time = nc_time[:].astype(np.int64)

        lat = self.get_values("pl_sur","latitude")
        lon = self.get_values("pl_sur","longitude")

        svf = self.get_sky_view()
        slope = self.get_slope()
        aspect = self.get_aspect()

        grid_elev = self.get_values("to", "Geopotential", (0, slice(None,None,1))) / const.G  # [m]
        station_elev = self.get_values("pl_sur","height")  # [m]

        for siteslist_ix, interp_ix in self.iterate_stations():
            zenith = solar_zenith(lat=lat[interp_ix], lon=lon[interp_ix], time=py_time)
            sw = self.get_station_values("sf", "Downward solar radiation flux", interp_ix, preserve_dims=False)  # [W m-2]

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

            self.rg.variables[vn_dir][:, siteslist_ix] = series_interpolate(self.times_out_nc, interpolation_time, corrected_direct)
            self.rg.variables[vn_diff][:, siteslist_ix] = series_interpolate(self.times_out_nc, interpolation_time, diffuse)
            self.rg.variables[vn_glob][:, siteslist_ix] = series_interpolate(self.times_out_nc, interpolation_time, global_sw)


    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        vn = kt.LW_Wm2_topo(self.rg, self.NAME)

        time_sf = self.get_values("sf", "time").astype(np.int64)
        time_pl = self.get_values("pl", "time").astype(np.int64)
        svf = self.get_sky_view()

        for siteslist_ix, interp_ix in self.iterate_stations():
            t_sub = self.upscale(time_pl, self.get_station_values("pl_sur", "Temperature", interp_ix), time_sf)  # [K]
            rh_sub = self.upscale(time_pl, self.get_station_values("pl_sur", "Relative humidity", interp_ix), time_sf)  # [%]
            t_grid = self.upscale(time_pl, self.get_station_values("sa", "Temperature", interp_ix), time_sf)  # [K]
            rh_grid = self.upscale(time_pl, self.get_station_values("sa", "Relative humidity", interp_ix), time_sf)  # [%]
            lw_grid  = self.get_station_values("sf", "Downward longwave radiation flux", interp_ix)  # [w m-2 s-1]

            lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)
            
            values = lw_sub * svf[siteslist_ix]
            self.rg.variables[vn][:, siteslist_ix] = self.upscale(time_sf, values, self.times_out_nc)

