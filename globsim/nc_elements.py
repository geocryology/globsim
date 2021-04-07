"""
functions for creating netcdf files
"""
from datetime import datetime
import netCDF4 as nc
import numpy as np

from globsim.common_utils import variables_skip, str_encode
from globsim._version import __version__
from os import path

GLOBAL_METADATA = {
    'standard_name_vocabulary': 'CF Standard Name Table v27',
}

def nc_new_file(ncfile_out, featureType="timeSeries", fmt='NETCDF4_CLASSIC'):
    rootgrp = nc.Dataset(ncfile_out, 'w', format=fmt)
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.featureType = featureType
    rootgrp.date_created = datetime.now().isoformat()[:19]  # ACDD attribute
    rootgrp.globsim_version = __version__

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
    height.long_name = 'Elevation relative to ellipsoid'
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
    add_history(rootgrp, 'scale', nc_interpol)

    return rootgrp


def new_interpolated_netcdf(ncfile_out, stations, nc_in, time_units):
    """
    Creates an empty station file to hold interpolated reults. The number of
    stations is defined by the variable stations, variables are determined by
    the variable list passed from the gridded original netCDF.

    ncfile_out: full name of the file to be created
    stations:   station list read with common_utils.StationListRead()
    variables:  variables read from netCDF handle
    lev:        list of pressure levels, empty is [] (default)
    """
    rootgrp = netcdf_base(ncfile_out, len(stations), None, time_units, nc_in)

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
        print("== 3D: file has pressure levels")
        level           = rootgrp.createDimension('level', len(lev))
        level           = rootgrp.createVariable('level', 'i4', ('level'))
        level.long_name = 'pressure_level'
        level.units     = 'hPa'
        level[:] = lev
    except Exception:
        print("== 2D: file without pressure levels")
        lev = []

    try:
        num = rootgrp['number'][:]
    except Exception:
        num = []

    # create and assign variables based on input file
    for n, var in enumerate(nc_in.variables):
        if variables_skip(var):
            continue
        print("VAR: ", str_encode(var))
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

        tmp.long_name = nc_in.variables[var].long_name.encode('UTF8')
        tmp.units     = nc_in.variables[var].units.encode('UTF8')

    # add history
    add_history(rootgrp, "interpolate", nc_in)

    return rootgrp


def netcdf_base(ncfile_out, n_stations, n_time, time_units, nc_in=None):
    # Build the netCDF file
    rootgrp = nc_new_file(ncfile_out)

    # dimensions
    _ = rootgrp.createDimension('station', n_stations)
    _ = rootgrp.createDimension('time', None)

    # base variables
    _ = ncvar_add_time(rootgrp, units=time_units,
                       calendar='gregorian', dimensions=('time'))
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

def add_history(rootgrp, globsim_command, nc_in=None):
    if nc_in and hasattr(nc_in, 'history'):
        print("adding history from nc_in")
        rootgrp.history = nc_in.history
    
    if not hasattr(rootgrp, 'history'):
        print("no history found")
        rootgrp.history = ""

    add_newline = (rootgrp.history != "") and (rootgrp.history[-1] != "\n")
    newline = "\n" if add_newline else ""

    if globsim_command.lower() in ['download', 'interpolate', 'scale']:
        new_history = f"{rootgrp.history}{newline}{datetime.now().isoformat()[:19]} globsim-{__version__} {globsim_command}"
        rootgrp.history = new_history
