#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import netCDF4 as nc
import numpy as np
import pytz


from math import atan2, pi
from pathlib import Path
from pysolar.solar import get_azimuth_fast
from scipy.interpolate import interp1d

from globsim.common_utils import str_encode, series_interpolate
from globsim.scale.toposcale import (lw_down_toposcale, illumination_angle,
                                     shading_corrected_sw_direct, elevation_corrected_sw, 
                                     solar_zenith)
from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.GenericScale import GenericScale, _check_timestep_length

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
    NAME = "JRA-55"
    REANALYSIS = "jra55"
    SCALING = {"sf": {"Total precipitation": (24 * 3600, 0)},
               "sa": {},
               "pl": {}}

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file names
        self.nc_pl = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_pl_{self.list_name}_surface.nc'), 'r')
        self.nc_sa = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_sf_{self.list_name}.nc'), 'r')
        try:
            self.nc_to = nc.Dataset(Path(self.intpdir, f'{self.REANALYSIS}_to_{self.list_name}.nc'), 'r')
        except AttributeError:
            logger.error("Missing invariant ('*_to') file. Some scaling kernels may fail. ")

        # Check data integrity
        _check_timestep_length(self.nc_pl.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(par, f'{self.REANALYSIS}')

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl.variables['time'], par['time_step'])

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
        self.rg = new_scaled_netcdf(self.output_file, self.nc_pl,
                                    self.times_out_nc, self.nc_pl['time'].units)

        # add station names to netcdf
        # first convert to character array
        names_out = nc.stringtochar(np.array(self.stations['station_name'], 'S32'))

        # create space in the netcdf
        _            = self.rg.createDimension('name_strlen', 32)
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
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()

    def get_file(self, file:str) -> "nc.Dataset":
        if file == "sa":
            f = self.nc_sa
        elif file == "sf":
            f = self.nc_sf
        elif file == "pl":
            f = self.nc_pl
        elif file == "to":
            f = self.nc_to
        else:
            raise ValueError("sa, sf, to, or pl")
        
        return f
    
    def get_name(self, file:str, jra55name:str) -> str:
        return jra55name
    
    def get_values(self, file:str, jra55name:str, _slice=None, attr=None):
        f = self.get_file(file)
        n = self.get_name(file, jra55name)
        
        if attr is not None:
            v = f.variables[n].getncattr(attr)
        
        else:
            if _slice is None:
                v = f.variables[n][:]
            else:
                v = f.variables[n][_slice]
            
            if jra55name in self.SCALING.get(file).keys():
                scale, offset = self.SCALING.get(file).get(jra55name)
                v *= scale
                v += offset

        return v
    
    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """
        # add variable to ncdf file
        vn = 'PRESS_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))
        var.long_name = 'air_pressure {} pressure levels only'.format(self.NAME)
        var.units     = 'Pa'
        var.standard_name = 'surface_air_pressure'

        # interpolate station by station
        time_in = self.get_values("pl","time").astype(np.int64)
        values  = self.get_values("pl", "air_pressure")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            # scale from hPa to Pa
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc,
                                                             time_in,
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
        time_in = self.get_values("pl","time")
        values  = self.get_values("pl","Temperature")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n]) - 273.15

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
        time_in = self.get_values("sa","time")
        values  = self.get_values("sa", "Temperature")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n]) - 273.15

    def RH_per_sur(self):
        """
        Relative Humidity derived from surface data, exclusively.
        """
        # add variable to ncdf file
        vn = 'RH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'relative humidity {} surface only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        # interpolate station by station
        time_in = self.get_values("sa","time")
        values  = self.get_values("sa", "Relative humidity")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

    def RH_per_pl(self):
        """
        Relative Humidity derived from pressure-level data, exclusively.
        """
        # add variable to ncdf file
        vn = 'RH_pl'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'relative humidity {} pressure-level only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'

        # interpolate station by station
        time_in = self.get_values("pl","time")
        values  = self.get_values("pl", "Relative humidity")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

    def WIND_sur(self):
        """
        Wind at 10 metre derived from surface data, exclusively.
        """

        # add variable to ncdf file
        vn = '10 metre U wind component'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '10 metre U wind component'
        var.units     = self.get_values("sa","u-component of wind", attr="units")

        # interpolate station by station
        time_in = self.get_values("sa","time")
        values  = self.get_values("sa", "u-component of wind")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

        # add variable to ncdf file
        vn = '10 metre V wind component'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = '10 metre V wind component'
        var.units     = self.get_values("sa","v-component of wind", attr="units")

        # interpolate station by station
        time_in = self.get_values("sa","time")
        values  = self.get_values("sa","v-component of wind")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

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

        # convert
        # u is the ZONAL VELOCITY, i.e. horizontal wind TOWARDS EAST.
        # v is the MERIDIONAL VELOCITY, i.e. horizontal wind TOWARDS NORTH.
        V = self.rg.variables['10 metre V wind component'][:]
        U = self.rg.variables['10 metre U wind component'][:]

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            WS = np.sqrt(np.power(V, 2) + np.power(U, 2))
            WD = [atan2(V[i, n], U[i, n]) * (180 / pi) +
                  180 for i in np.arange(V.shape[0])]

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
        time_in = self.get_values("sf","time")
        values  = self.get_values("sf","Downward solar radiation flux")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

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
        lat = self.get_values("pl","latitude")
        lon = self.get_values("pl","longitude")
        sw = self.get_values("sf","Downward solar radiation flux")  # [W m-2]
        grid_elev = self.get_values("to", "Geopotential", (0, slice(None,None,1))) / 9.80665  # [m]
        station_elev = self.get_values("pl","height")  # [m]

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

            f = interp1d(interpolation_time, corrected_direct, kind='linear')
            self.rg.variables[vn_dir][:, n] = f(self.times_out_nc)

            f = interp1d(interpolation_time, diffuse, kind='linear')
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
        time_in = self.get_values("sf","time")
        values  = self.get_values("sf","Downward longwave radiation flux")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])

    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        # add variable to ncdf file
        vn = 'LW_topo'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'TOPOscale-corrected thermal radiation downwards'
        var.standard_name = 'surface_downwelling_longwave_flux'
        var.units     = 'W m-2'

        # interpolate station by station
        time_in = self.get_values("sf","time").astype(np.int64)
        
        t_sub = self.get_values("pl","Temperature")  # [K]
        rh_sub = self.get_values("pl","Relative humidity")  # [%]
        t_grid = self.get_values("sa","Temperature")  # [K]
        rh_grid = self.get_values("sa","Relative humidity")  # [%]
        lw_grid  = self.nc_sf["Downward longwave radiation flux"]  # [w m-2 s-1]

        lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)
        
        svf = self.get_sky_view()

        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            values = lw_sub[:, n] * svf[n]
            f = interp1d(time_in, values, kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)

    def PREC_mm_sur(self):
        """
        Precipitation derived from surface data, exclusively.
        Convert unit: mm/day to mm/s (kg m-2 s-1)
        """

        # add variable to ncdf file
        vn = 'PREC_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'Total precipitation {} surface only'.format(self.NAME)
        var.units     = 'kg m-2 s-1'
        var.comment = "units [kg m-2 s-1] corresponds to [mm/s] for water (density 1000 [kg m-3])"
        var.standard_name = 'precipitation_flux'

        # interpolate station by station
        time_in = self.get_values("sf","time")
        values  = self.get_values("sf","Total precipitation") / (24 * 3600)  # [mm/s]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            f = interp1d(time_in, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * self.scf

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg] derived from surface data, exclusively
        '''
        # add variable to ncdf file
        vn = 'SH_sur'  # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))
        var.long_name = 'specific humidity {} surface only'.format(self.NAME)
        var.units     = 'kg/kg'
        var.standard_name = 'specific_humidity'

        # interpolate station by station
        time_in = self.get_values("sa","time")
        values  = self.get_values("sa","Specific humidity")
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc,
                                                    time_in, values[:, n])
