#!/usr/bin/env python
# -*- coding: utf-8 -*-
import netCDF4 as nc
import numpy as np
import logging
from datetime import datetime
from os import path

from globsim.common_utils import str_encode, variables_skip
from globsim.interpolate.GenericInterpolate import GenericInterpolate
from globsim.nc_elements import netcdf_base, new_interpolated_netcdf
from globsim.interp import calculate_weights, ele_interpolate

logger = logging.getLogger('globsim.interpolate')


class JRAinterpolate(GenericInterpolate):
    """
    Algorithms to interpolate JRA55 netCDF files to station coordinates.
    All variables retains their original units and time-steps.

    Referenced from era_interim.py (Dr.Stephan Gruber): Class ERAinterpolate()

    Args:
        ifile: Full path to a Globsim Interpolate Paramter file
        JRAinterpolate(ifile)


    Example:
        ifile = '/home/xquan/src/globsim/examples/par/examples.globsim_interpolate'
        JRAinterpolate(ifile)

    """
    REANALYSIS = "jra55"
    SA_INTERVAL = 6
    PL_INTERVAL = 6
    SF_INTERVAL = 6
    T_UNITS = 'hours since 1800-01-01 00:00:0.0'

    dpar_sa = {'air_temperature'   : ['surface_temperature'],  # [K] 2m values
                'relative_humidity' : ['relative_humidity'],  # [%]
                'wind_speed' : ['eastward_wind', 'northward_wind']}  # [m s-1] 2m & 10m values
    
    dpar_sf = dpar = {'downwelling_shortwave_flux_in_air': ['downwelling_shortwave_flux_in_air'],  # [W/m2] short-wave downward
                'downwelling_longwave_flux_in_air' : ['downwelling_longwave_flux_in_air'],  # [W/m2] long-wave downward
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['downwelling_shortwave_flux_in_air_assuming_clear_sky'],  # [W/m2] short-wave downward assuming clear sky
                'downwelling_longwave_flux_in_air_assuming_clear_sky': ['downwelling_longwave_flux_in_air_assuming_clear_sky'],
                'precipitation_amount' : ['total_precipitation']}  # [W/m2] long-wave downward assuming clear sky

    dpar_pl = {'air_temperature'   : ['air_temperature'],  # [K]
               'relative_humidity' : ['relative_humidity'],  # [%]
               'wind_speed'        : ['eastward_wind', 'northward_wind']}  # [m s-1]
    
    GEOPOTENTIAL = "Geopotential height"

    def __init__(self, ifile, **kwargs):
        super().__init__(ifile, **kwargs)
        par = self.par

        self.input_dir = path.join(par['project_directory'], self.REANALYSIS)

        # Override inherited chunk size
        self.cs *= 200

        # Load MF Datasets
        self.mf_to = nc.MFDataset(path.join(self.input_dir, f'{self.REANALYSIS}_to_*.nc'), 'r', aggdim="time")
        self.mf_sa = nc.MFDataset(path.join(self.input_dir, f'{self.REANALYSIS}_sa_*.nc'), 'r', aggdim="time")
        self.mf_sf = nc.MFDataset(path.join(self.input_dir, f'{self.REANALYSIS}_sf_*.nc'), 'r', aggdim="time")
        self.mf_pl = nc.MFDataset(path.join(self.input_dir, f'{self.REANALYSIS}_pl_*.nc'), 'r', aggdim="time")
        
        # Check dataset integrity
        logger.info("Check data integrity (sa)")
        self.ensure_datset_integrity(self.mf_sa['time'], self.SA_INTERVAL)
        logger.info("Check data integrity (sf)")
        self.ensure_datset_integrity(self.mf_sf['time'], self.SF_INTERVAL)
        logger.info("Check data integrity (pl)")
        self.ensure_datset_integrity(self.mf_pl['time'], self.PL_INTERVAL)
        logger.info("Data integrity ok")

    def JRA2station(self, ncf_in: "nc.MFDataset", ncfile_out, points,
                    variables=None, date=None):

        """
        Biliner interpolation from fields on regular grid (latitude, longitude)
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all JRA55 files). The type
        of variable and file structure are determined from the input.

        This function creates an empty of netCDF file to hold the interpolated
        results, by calling self.netCDF_empty(). Then, data is
        interpolated in temporal chunks and appended. The temporal chunking can
        be set in the interpolation parameter file.

        Args:
        ncf_in: a MFDataset of globsim-downloaded files

        ncfile_out: Full path to the output netCDF file to write.

        points: A dataframe of locations. See method StationListRead in
                common_utils.py for more details.

        variables:  List of variable(s) to interpolate such as
                    [air_temperature, easteard_wind, northward_wind,
                    relative_humidy, surface_temperature,
                    downwelling_shortwave_flux_in_air,
                    downwelling_longwave_flux_in_air,
                    downwelling_shortwave_flux_in_air_assuming_clear_sky,
                    downwelling_longwave_flux_in_air_assuming_clear_sky].
                    Defaults to using all variables available.

        date: Directory to specify begin and end time for the derived time
                series. Defaluts to using all times available in ncf_in.

        cs: chunk size, i.e. how many time steps to interpolate at once. This
            helps to manage overall memory usage (small cs is slower but less
            memory intense).
        """
        # Check station bounds
        self.validate_stations_extent(ncf_in)

        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # build the output of empty netCDF file
        rootgrp = new_interpolated_netcdf(ncfile_out, self.stations, ncf_in,
                                          time_units=self.T_UNITS)
        rootgrp.source = f'{self.REANALYSIS}, interpolated bilinearly to stations'
        rootgrp.close()

        # open the output netCDF file, set it to be appendable ('a')
        ncf_out = nc.Dataset(ncfile_out, 'a')

        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]
        # "hours since 1900-01-01 00:00:0.0"
        t_unit = ncf_in.variables['time'].units
        try:
            t_cal = ncf_in.variables['time'].calendar
        except AttributeError:  # attribute doesn't exist
            t_cal = u"gregorian"  # standard
        time = nc.num2date(nctime, units=t_unit, calendar=t_cal)

        # detect invariant files (topography etc.)
        invariant = True if len(np.unique(time)) == 1 else False

        # restrict to date/time range if given
        if date is None:
            tmask = time < datetime(3000, 1, 1)
        else:
            tmask = (time < date['end']) * (time >= date['beg'])

        # get time indices
        time_in = nctime[tmask]

        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter = len(time_in) // self.cs
        niter += ((len(time_in) % self.cs) > 0)

        # Create source grid
        sgrid = self.create_source_grid(ncf_in)
        subset_grid, lon_slice, lat_slice = self.create_subset_source_grid(sgrid, self.stations_bbox)
        
        # loop over chunks
        for n in range(niter):
            # indices
            beg = n * self.cs
            # restrict last chunk to lenght of tmask plus one (to get last time)
            if invariant:
                end = beg
            else:
                end = min(n * self.cs + self.cs, len(time_in)) - 1

            # time to make tmask for chunk
            beg_time = nc.num2date(time_in[beg], units=t_unit, calendar=t_cal)
            if invariant:
                # allow topography to work in same code, len(nctime) = 1
                end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                # end = 1
            else:
                end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)

            # '<= end_time', would damage appending
            tmask_chunk = (time <= end_time) * (time >= beg_time)
            if invariant:
                # allow topography to work in same code
                tmask_chunk = np.array([True])

            # get the interpolated variables
            dfield, variables = self.interp2D(ncf_in,
                                              self.stations, tmask_chunk,
                                              subset_grid, lon_slice, lat_slice,
                                              variables=None, date=None)

            # append time
            ncf_out.variables['time'][:] = np.append(ncf_out.variables['time'][:],
                                                     time_in[beg:end + 1])

            # append variables
            for i, var in enumerate(variables):
                if variables_skip(var):
                    continue

                if pl:
                    lev = ncf_in.variables['level'][:]

                    ncf_out.variables[var][beg:end + 1,:,:] = dfield.data[:,i,:,:].transpose((1,2,0))
                else:
                    ncf_out.variables[var][beg:end + 1,:] = dfield.data[:,i,:].transpose((1,0))

        # close the file
        ncf_in.close()
        ncf_out.close()

    def levels2elevation(self, ncfile_in, ncfile_out):
        """
        Linear 1D interpolation of pressure level data available for individual
        stations to station elevation. Where and when stations are below the
        lowest pressure level, they are assigned the value of the lowest
        pressure level.

        """
        # open file

        ncf = nc.MFDataset(ncfile_in, 'r', aggdim='time')
        height = ncf.variables['height'][:]
        nt = len(ncf.variables['time'][:])
        nl = len(ncf.variables['level'][:])

        # list variables
        varlist = [str_encode(x) for x in ncf.variables.keys()]
        for V in ['time', 'station', 'latitude', 'longitude', 'level', 'height']:
            varlist.remove(V)

        # === open and prepare output netCDF file =============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = netcdf_base(ncfile_out, len(height), nt,
                              self.T_UNITS)
        rootgrp.source = f'{self.REANALYSIS}, interpolated (bi)linearly to stations'

        # access variables
        time = rootgrp['time']
        station = rootgrp['station']
        latitude = rootgrp['latitude']
        longitude = rootgrp['longitude']
        height = rootgrp['height']

        # assign base variables
        time[:]      = ncf.variables['time'][:]
        station[:]   = ncf.variables['station'][:]
        latitude[:]  = ncf.variables['latitude'][:]
        longitude[:] = ncf.variables['longitude'][:]
        height[:]    = ncf.variables['height'][:]

        # create and assign variables from input file
        for var in varlist:
            vname = str_encode(ncf.variables[var].long_name)
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))
            tmp.long_name = str_encode(ncf.variables[var].long_name)
            tmp.units     = str_encode(ncf.variables[var].units)

        # add air pressure as new variable
        var = 'air_pressure'
        varlist.append(var)
        tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))
        tmp.long_name = str_encode(var)
        tmp.units     = str_encode('hPa')
        # end file preparation ================================================

        # loop over stations
        for n, h in enumerate(height):
            # convert geopotential [millibar] to height [m]
            # shape: (time, level)
            elevation = ncf.variables[self.GEOPOTENTIAL][:,:,n]

            elev_diff, va, vb = ele_interpolate(elevation, h, nl)
            wa, wb = calculate_weights(elev_diff, va, vb)

            # loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                if var == 'air_pressure':
                    # pressure [Pa] variable from levels, shape: (time, level)
                    data = np.repeat([ncf.variables['level'][:]],
                                     len(time),axis=0).ravel()
                else:
                    # read data from netCDF
                    data = ncf.variables[var][:,:,n].ravel()

                multvawa = np.multiply(data[va], wa)
                multvbwb = np.multiply(data[vb], wb)
                ipol = multvawa + multvbwb
                rootgrp.variables[var][:,n] = ipol  # assign to file

        rootgrp.close()
        ncf.close()
        # closed file =========================================================

    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to
        the more generically JRA-like interpolation functions.
        """
        try:
            self.JRA2station(self.mf_to,
                             path.join(self.output_dir,f'{self.REANALYSIS}_to_{self.list_name}.nc'),
                             self.stations, ['Geopotential'], date=None)
        except OSError:
            logger.error("Could not find invariant ('*_to') geopotential files for JRA. These were not downloaded in earlier versions of globsim. You may need to download them."
                         "  . Some scaling kernels may not work. Future versions of globsim may be less accepting of missing files.")

        # === 2D Interpolation for Surface  Data ===
        # dictionary to translate CF Standard Names into JRA55
        # pressure level variable keys.
        varlist = self.TranslateCF2short(self.dpar_sa)
        self.JRA2station(self.mf_sa,
                         path.join(self.output_dir,f'{self.REANALYSIS}_sa_' + self.list_name + '.nc'),
                         self.stations,
                         varlist, date=self.date)

        # 2D Interpolation for Radiation Data
        # dictionary to translate CF Standard Names into JRA55
        # pressure level variable keys.
        varlist = self.TranslateCF2short(self.dpar_sf)
        self.JRA2station(self.mf_sf,
                         path.join(self.output_dir,f'{self.REANALYSIS}_sf_' + self.list_name + '.nc'),
                         self.stations,
                         varlist,
                         date=self.date)
        
        varlist = self.TranslateCF2short(self.dpar_pl).append('geopotential_height')
        self.JRA2station(self.mf_pl,
                         path.join(self.output_dir,f'{self.REANALYSIS}_pl_' + self.list_name + '.nc'),
                         self.stations,
                         varlist,
                         date=self.date)

        # 1D Interpolation for Pressure Level Data
        self.levels2elevation(path.join(self.output_dir,f'{self.REANALYSIS}_pl_' + self.list_name + '.nc'),
                              path.join(self.output_dir,f'{self.REANALYSIS}_pl_' + self.list_name + '_surface.nc'))
