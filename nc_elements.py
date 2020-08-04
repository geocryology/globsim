"""
functions for creating netcdf files
"""
import netCDF4 as nc

def nc_create_basic_timeseries(ncfile_out):
    rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4_CLASSIC')
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.featureType = "timeSeries"

    return rootgrp


def ncvar_add_latitude(rootgrp, dimensions=('station')):
    latitude  = rootgrp.createVariable('latitude', 'f4', dimensions)
    latitude.long_name  = 'latitude'
    latitude.units  = 'degrees_north'
    latitude.standard_name = 'latitude'
    latitude.axis  = 'Y'

    return latitude


def ncvar_add_longitude(rootgrp, dimensions=('station')):
    longitude           = rootgrp.createVariable('longitude', 'f4', dimensions)
    longitude.long_name = 'longitude'
    longitude.units     = 'degrees_east'
    longitude.standard_name = 'longitude'
    longitude.axis  = 'X'

    return longitude


def ncvar_add_time(rootgrp, units, calendar):
    time           = rootgrp.createVariable('time', 'i4', ('time'))
    time.long_name = 'time'
    time.units     = units
    time.calendar  = calendar
    time.standard_name = 'time'

    return time


def ncvar_add_station(rootgrp, dimensions=('station')):
    station             = rootgrp.createVariable('station', 'i4', dimensions)
    station.long_name   = 'station for time series data'
    station.units       = '1'

    return station


def ncvar_add_geoid_height(rootgrp, dimensions=('station')):
    height           = rootgrp.createVariable('height', 'f4', dimensions)
    height.long_name = 'height_above_reference_ellipsoid'
    height.units     = 'm'
    height.standard_name = 'geoid_height_above_reference_ellipsoid'

    return height