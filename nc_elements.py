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

def ScaledFileOpen(ncfile_out, nc_interpol, times_out, t_unit, station_names=None):
    """
    Create netCDF file for scaled results (same for all reanalyses)
    Returns the file object so that kernel functions can
    successively write variables to it.

    """

    if path.isfile(ncfile_out):  # raise exception if file exists
        raise FileExistsError("File already exists: {}".format(ncfile_out))


    # make netCDF outfile, variables are written in kernels
    rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.source      = 'Reanalysis data interpolated and scaled to stations'
    rootgrp.featureType = "timeSeries"

    name = ncfile_out[-9:-3]

    # dimensions
    station = rootgrp.createDimension('station',
                    len(nc_interpol.variables['station'][:]))
    time    = rootgrp.createDimension('time', len(times_out))

    # base variables
    time           = rootgrp.createVariable('time', 'i8',('time'))
    time.long_name = 'time'
    time.axis      = 'T'

    time.units = t_unit
    time.calendar  = 'gregorian'

    station             = rootgrp.createVariable('station', 'i4',('station'))
    station.long_name   = 'station index for time series data'
    station.units       = ''
    station.standard_name = 'platform_id'

    latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
    latitude.standard_name  = 'latitude'
    latitude.units      = 'degrees_N'
    latitude.long_name  = 'WGS84 latitude'
    latitude.axis       = 'Y'

    longitude           = rootgrp.createVariable('longitude','f4',('station'))
    longitude.standard_name = 'longitude'
    longitude.long_name = 'WGS84 longitude'
    longitude.units     = 'degrees_E'
    longitude.axis      = 'X'

    height           = rootgrp.createVariable('height','f4',('station'))
    height.standard_name = 'height_above_reference_ellipsoid'
    height.long_name = 'height_above_reference_ellipsoid'
    height.units     = 'm'
    height.axis      = 'Z'
    height.positive  = 'up'

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
        nchar        = rootgrp.createDimension('name_strlen', 32)
        st           = rootgrp.createVariable('station_name', "S1", ('station', 'name_strlen'))
        st.standard_name = 'platform_name'
        st.units     = ''

        # add data
        st[:] = names_out

    return rootgrp


def create_empty_netcdf(ncfile_out, stations, nc_in, time_units):
    """
    Creates an empty station file to hold interpolated reults. The number of
    stations is defined by the variable stations, variables are determined by
    the variable list passed from the gridded original netCDF.

    ncfile_out: full name of the file to be created
    stations:   station list read with generic.StationListRead()
    variables:  variables read from netCDF handle
    lev:        list of pressure levels, empty is [] (default)
    """
    rootgrp = netcdf_base(ncfile_out, stations, time_units)

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

    # create and assign variables based on input file
    for n, var in enumerate(nc_in.variables):
        if variables_skip(var):
            continue
        print("VAR: ", str_encode(var))
        # extra treatment for pressure level files
        if len(lev):
            tmp = rootgrp.createVariable(var, 'f4', ('time', 'level', 'station'))
        else:
            tmp = rootgrp.createVariable(var, 'f4', ('time', 'station'))
        tmp.long_name = nc_in.variables[var].long_name.encode('UTF8')
        tmp.units     = nc_in.variables[var].units.encode('UTF8')

    return rootgrp

def netcdf_base(ncfile_out, stations, time_units):
    # TODO change date type from f4 to f8 for lat and lon
    # Build the netCDF file
    rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4_CLASSIC')
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.featureType = "timeSeries"

    # dimensions
    station = rootgrp.createDimension('station', len(stations))
    time    = rootgrp.createDimension('time', None)

    # base variables
    time           = rootgrp.createVariable('time', 'i4', ('time'))
    time.long_name = 'time'
    time.units     = time_units
    time.calendar  = 'gregorian'
    station             = rootgrp.createVariable('station', 'i4', ('station'))
    station.long_name   = 'station for time series data'
    station.units       = '1'
    latitude            = rootgrp.createVariable('latitude', 'f4', ('station'))
    latitude.long_name  = 'latitude'
    latitude.units      = 'degrees_north'
    longitude           = rootgrp.createVariable('longitude', 'f4', ('station'))
    longitude.long_name = 'longitude'
    longitude.units     = 'degrees_east'
    height           = rootgrp.createVariable('height', 'f4', ('station'))
    height.long_name = 'height_above_reference_ellipsoid'
    height.units     = 'm'

    return rootgrp