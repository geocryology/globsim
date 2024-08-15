"""
functions for creating netcdf files
"""
import logging
import netCDF4 as nc
import numpy as np
from os import path
import xarray as xr

from globsim.common_utils import variables_skip, str_encode
from globsim import __version__ as globsim_version


logger = logging.getLogger('globsim.nc_elements')


def nc_new_file(ncfile_out, featureType="timeSeries", fmt='NETCDF4_CLASSIC'):
    rootgrp = nc.Dataset(ncfile_out, 'w', format=fmt)
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.featureType = featureType
    rootgrp.globsim_version = globsim_version

    return rootgrp


def ncvar_add_latitude(rootgrp, dimensions=('station')):
    latitude  = rootgrp.createVariable('latitude', 'f8', dimensions)
    latitude.long_name  = 'latitude'
    latitude.units  = 'degrees_north'
    latitude.standard_name = 'latitude'
    latitude.axis  = 'Y'

    return latitude


def ncvar_add_longitude(rootgrp, dimensions=('station')):
    longitude           = rootgrp.createVariable('longitude', 'f8', dimensions)
    longitude.long_name = 'longitude'
    longitude.units     = 'degrees_east'
    longitude.standard_name = 'longitude'
    longitude.axis  = 'X'

    return longitude


def ncvar_add_time(rootgrp, units, calendar, dimensions=('time'), dtype='i4'):
    time           = rootgrp.createVariable('time', dtype, dimensions)
    time.long_name = 'time'
    time.units     = units
    time.calendar  = calendar
    time.standard_name = 'time'
    time.axis = 'T'

    return time


def ncvar_add_station(rootgrp, dimensions=('station')):
    station            = rootgrp.createVariable('station', 'i4', dimensions)
    station.long_name  = 'station for time series data'
    station.units      = '1'

    return station


def ncvar_add_number(rootgrp, dimensions=('number')):
    number             = rootgrp.createVariable('number', 'i4', dimensions)
    number.long_name   = 'ensemble_member'
    number.units       = '1'

    return number


def ncvar_add_ellipsoid_height(rootgrp, dimensions=('station')):
    height           = rootgrp.createVariable('height', 'f4', dimensions)
    height.long_name = 'Elevation of station'
    height.units     = 'm'
    height.axis      = 'Z'
    height.standard_name = 'height_above_reference_ellipsoid'
    height.positive  = 'up'

    return height


def new_scaled_netcdf(ncfile_out, nc_interpol, times_out,
                      t_unit, station_names=None):
    """
    Create netCDF file for scaled results (same for all reanalyses)
    Returns the file object so that kernel functions can
    successively write variables to it.

    """

    if path.isfile(ncfile_out):  # raise exception if file exists
        raise FileExistsError("File already exists: {}".format(ncfile_out))

    # make netCDF outfile, variables are written in kernels
    rootgrp = nc_new_file(ncfile_out, fmt='NETCDF4')
    rootgrp.source      = 'Reanalysis data interpolated and scaled to stations'

    # dimensions
    n_station = len(nc_interpol.variables['station'][:])
    station = rootgrp.createDimension('station', n_station)
    time    = rootgrp.createDimension('time', len(times_out))

    # base variables
    time           = ncvar_add_time(rootgrp, units=t_unit,
                                    calendar='gregorian', dimensions=('time'),
                                    dtype='i8')
    station        = ncvar_add_station(rootgrp)
    latitude       = ncvar_add_latitude(rootgrp)
    longitude      = ncvar_add_longitude(rootgrp)
    height         = ncvar_add_ellipsoid_height(rootgrp)

    crs           = rootgrp.createVariable('crs','i4')
    crs.long_name = 'coordinate system'
    crs.grid_mapping_name = 'latitude_longitude'
    crs.longitude_of_prime_meridian = 0.0
    crs.semi_major_axis = 6378137
    crs.inverse_flattening = 298.2572236
    crs.wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

    # assign base variables
    time[:]      = times_out
    station[:]   = nc_interpol.variables['station'][:]
    latitude[:]  = nc_interpol.variables['latitude'][:]
    longitude[:] = nc_interpol.variables['longitude'][:]
    height[:]    = nc_interpol.variables['height'][:]

    # add station names to netcdf
    if station_names is not None:
        # first convert to character array
        names_out = nc.stringtochar(np.array(station_names, 'S32'))

        # create space in the netcdf
        _            = rootgrp.createDimension('name_strlen', 32)
        st           = rootgrp.createVariable('station_name', "S1", ('station', 'name_strlen'))
        st.standard_name = 'platform_name'
        st.units     = ''

        # add data
        st[:] = names_out

    return rootgrp


def new_interpolated_netcdf(ncfile_out, stations, nc_in, time_units, calendar=None):
    """
    Creates an empty station file to hold interpolated reults. The number of
    stations is defined by the variable stations, variables are determined by
    the variable list passed from the gridded original netCDF.

    ncfile_out: full name of the file to be created
    stations:   station list read with common_utils.StationListRead()
    variables:  variables read from netCDF handle
    lev:        list of pressure levels, empty is [] (default)
    """
    logger.info(f"Creating new file {ncfile_out} from ")

    rootgrp = netcdf_base(ncfile_out, len(stations), None, time_units, nc_in, calendar)

    station = rootgrp['station']
    latitude = rootgrp['latitude']
    longitude = rootgrp['longitude']
    height = rootgrp['height']

    # assign station characteristics
    station[:]   = list(stations['station_number'])
    latitude[:]  = list(stations['latitude_dd'])
    longitude[:] = list(stations['longitude_dd'])
    height[:]    = list(stations['elevation_m'])

    # extra treatment for pressure level files
    try:
        lev = nc_in.variables['level'][:]
        logger.info(f"Source dataset is 3D (has pressure levels)")
        level           = rootgrp.createDimension('level', len(lev))
        level           = rootgrp.createVariable('level', 'i4', ('level'))
        level.long_name = 'pressure_level'
        level.units     = 'hPa'
        level[:] = lev
    except Exception:
        logger.info(f"Source dataset is 2D (without pressure levels)")
        lev = []

    try:
        num = rootgrp['number'][:]
    except Exception:
        num = []

    # create and assign variables based on input file
    for n, var in enumerate(nc_in.variables):
        if variables_skip(var):
            continue
        
        # extra treatment for pressure level files
        if len(num):
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4',('time', 'number',
                                                       'level', 'station'))
            else:
                tmp = rootgrp.createVariable(var,'f4',('time','number',
                                                       'station'))
        else:
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4', ('time','level',
                                                        'station'))
            else:
                tmp = rootgrp.createVariable(var,'f4', ('time','station'))

        # copy attributes
        input_var = nc_in.variables[var]

        if isinstance(input_var, xr.Variable):  # XARRAY
            for key, val in input_var.attrs.items():
                if key in ['_FillValue']:
                    continue
                tmp.setncattr(key, val)
        else:
            for key in input_var.ncattrs():  # NETCDF4
                if key in ['_FillValue']:
                    continue
                tmp.setncattr(key, getattr(input_var, key))
        units = tmp.units if hasattr(tmp, 'units') else '??'
        logger.info(f"Created new empty variable: {str_encode(var)} [{units}]")
    
    return rootgrp


def netcdf_base(ncfile_out, n_stations, n_time, time_units, nc_in=None, calendar=None):
    # Build the netCDF file
    rootgrp = nc_new_file(ncfile_out)

    # dimensions
    _ = rootgrp.createDimension('station', n_stations)
    _ = rootgrp.createDimension('time', None)

    # base variables
    if calendar is None:
        calendar = 'gregorian'
    _ = ncvar_add_time(rootgrp, units=time_units,
                       calendar=calendar, dimensions=('time'))
    _ = ncvar_add_station(rootgrp)
    _ = ncvar_add_latitude(rootgrp)
    _ = ncvar_add_longitude(rootgrp)
    _ = ncvar_add_ellipsoid_height(rootgrp)

    # for ERA5 ensemble member only
    try:
        num    = nc_in.variables['number'][:]
        number = rootgrp.createDimension('number', len(num))
        number = ncvar_add_number(rootgrp)
        number = rootgrp['number']
        number[:] = num
    except Exception:
        pass

    return rootgrp
