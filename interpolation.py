#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Methods for interpolating 2D (longitude, latitude) and 3D (level, longitude, latitude) of 
# climate and other data from their original regular grid to specific points (longitude, 
# latitude, elevation) 
#
# These methods accept data arrays as input and return arrays as output. For each data source
# such as ERA-Interim, MERRA, etc. a separate script is needed to read teh data and provide it in
# the right form for interpolation.  
#
#
# (C) Copyright Stephan Gruber (2017)
#         
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# The highly efficient ESMPy package is used for accomplishing interpolation
# https://www.earthsystemcog.org/projects/esmpy/.
#
#  CF-Convention: this is how you check netcdf file conformity: 
#  http://pumatest.nerc.ac.uk/cgi-bin/cf-checker-dev.pl 
#
#===============================================================================
import csv
import ESMF
import numpy as np
import netCDF4 as nc
from datetime import datetime
    
    
def csv_stations(filename):
    """
    Read point data from csv following template below. Returns list of dictionaries.
    ---------------    
    station_number,longitude_deg,latitude_deg,elevation_m
    1,249.32,64.12,278.1
    2,249.95,64.92,318.8
    ---------------     
    """
    # as list of dictionaries
    with open(filename, 'rb') as f:
	reader = csv.DictReader(f)
	points = []
	for row in reader:
	   points.append(row)    
	return points
    

def int2Dpoint(ncfile_in, ncfile_out, points, variables=None, date=None):    
    """
    Biliner 2D interpolation from fields on a regular grid (latitude, longitude) to individual 
    point stations (latitude, longitude).
    
          
    Args:
        ncfile_in: Full path to am ERA-Interim derived 2D netCDF file.
              
        ncfile_out: Full path to the output netCDF file to write.  
              
        points: A dictionary of locations ().
        
        variables:  List of variable(s) to interpolate such as ['airt', 'rh', 'geop', 'wind'].
                    Defaults to using all variables available.
        
        date: Directory to specify begin and end time for the derived time series. Defaluts to 
              using the entire time available in ncfile_in.
              
    Example:
        from datetime import datetime
        date  = {'beg' : datetime(2008, 1, 1),
                 'end' : datetime(2008,12,31)}
        variables  = ['ssrd','tp']       
        stations = csv_stations("points.csv")      
        int2Dpoint('big/era_sa_19790101_to_20161231.nc', 'era_sa_inter.nc', stations, variables=variables, date=date)        
    """   
        
    # open netcdf file handle
    ncf = nc.Dataset(ncfile_in, 'r')

    # get spatial dimensions
    lat  = ncf.variables['latitude'][:]
    lon  = ncf.variables['longitude'][:]
    
    # get time and convert to datetime object
    nctime = ncf.variables['time'][:]
    t_unit = ncf.variables['time'].units # get unit  "hours since 1900-01-01 00:00:0.0"
    try :
        t_cal = ncf.variables['time'].calendar
    except AttributeError : # Attribute doesn't exist
        t_cal = u"gregorian" # or standard
    time = nc.num2date(nctime, units = t_unit, calendar = t_cal)

    # restrict to date/time range if given
    if date is None:
        tmask = time < datetime(3000, 1, 1)
    else:
        tmask = (time <= date['end']) * (time >= date['beg'])

    # test if time steps to interpolate remain
    nt = sum(tmask)
    if nt == 0:
        raise ValueError('No time steps from netCDF file selected.')
    
    # get variables
    varlist = [x.encode('UTF8') for x in ncf.variables.keys()]
    varlist.remove('time')
    varlist.remove('latitude')
    varlist.remove('longitude')
    
    #list variables that should be interpolated
    if variables is None:
        variables = varlist
    #test is variables given are available in file
    if (set(variables) < set(varlist) == 0):
        raise ValueError('One or more variables not available in netCDF file.')
        
    # create source grid from a SCRIP formatted file
    sgrid = ESMF.Grid(filename=ncfile_in, filetype=ESMF.FileFormat.GRIDSPEC)

    # create source field on source grid
    sfield = ESMF.Field(sgrid, name='sgrid',
                       staggerloc=ESMF.StaggerLoc.CENTER,
                       ndbounds=[len(variables), nt])

    # assign data from ncdf: (variale, time, latitude, longitude) 
    for n, var in enumerate(variables):
        sfield.data[n,:,:,:] = ncf.variables[var][tmask,:,:] 

    # create locstream structure, CANNOT have third dimension for 2D interpolation
    locstream = ESMF.LocStream(len(points), coord_sys=ESMF.CoordSys.SPH_DEG)
    locstream["ESMF:Lon"] = [float(e['longitude_deg']) for e in points]
    locstream["ESMF:Lat"] = [float(e['latitude_deg']) for  e in points]
    
    # test if points are inside netCDF file grid
    geom = (min(locstream["ESMF:Lon"]) < min(lon) +
            max(locstream["ESMF:Lon"]) > max(lon) +
            min(locstream["ESMF:Lat"]) < min(lat) +
            max(locstream["ESMF:Lat"]) > max(lat))
    if geom > 0:
        raise ValueError('Points outside netCDF file.')

    # create destination field
    dfield = ESMF.Field(locstream, name='dfield', ndbounds=[len(variables), nt])

    # regridding function
    regrid2D = ESMF.Regrid(sfield, dfield,
                           regrid_method=ESMF.RegridMethod.BILINEAR,
                           unmapped_action=ESMF.UnmappedAction.ERROR,
                           dst_mask_values=None)
                          
    # regrid operation, create destination field (variables, times, points)
    dfield = regrid2D(sfield, dfield)                          
    
    # === write output netCDF file ==============================================================
    # dimensions: station, time
    # variables: latitude(station), longitude(station), elevation(station), ...(time, station) 
    # stations are integer numbers
    # create a file (Dataset object, also the root group).
    rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.source = 'ERA-Interim, interpolated bilinearly to stations'
    rootgrp.featureType = "timeSeries"

    # dimensions
    station = rootgrp.createDimension('station', len(points))
    time    = rootgrp.createDimension('time', nt)

    # base variables
    time           = rootgrp.createVariable('time',     'i4',('time'))
    time.long_name = 'time'
    time.units     = 'hours since 1900-01-01 00:00:0.0'
    time.calendar  = 'gregorian'
    station             = rootgrp.createVariable('station',  'i4',('station'))
    station.long_name   = 'station for time series data'
    station.units       = '1'
    latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
    latitude.long_name  = 'latitude'
    latitude.units      = 'degrees_north'    
    longitude           = rootgrp.createVariable('longitude','f4',('station'))
    longitude.long_name = 'longitude'
    longitude.units     = 'degrees_east'  
    height           = rootgrp.createVariable('height','f4',('station'))
    height.long_name = 'height_above_reference_ellipsoid'
    height.units     = 'm'  
    
    # assign base variables
    time[:] = nctime[tmask]
    station[:]   = [float(e['station_number']) for e in points]
    latitude[:]  = [float(e['latitude_deg'])   for e in points]
    longitude[:] = [float(e['longitude_deg'])  for e in points]
    height[:]    = [float(e['elevation_m'])    for e in points]
    
    # create and assign variables from input file
    for n, var in enumerate(variables):
        vname = ncf.variables[var].name.encode('UTF8')
        tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))
        tmp.long_name = ncf.variables[var].long_name.encode('UTF8')
        tmp.units     = ncf.variables[var].units.encode('UTF8')  
        # assign values
        print dfield.data.shape
        print n
        tmp[:] = dfield.data[n,:,:]
    
    # close file
    rootgrp.close()
                 
 
def ERA2point(ncfile_in, ncfile_out, points, variables=None, date=None):    
    """
    Biliner 2D interpolation from fields on a regular grid (latitude, longitude) to individual 
    point stations (latitude, longitude).
    
          
    Args:
        ncfile_in: Full path to am ERA-Interim derived 2D netCDF file.
              
        ncfile_out: Full path to the output netCDF file to write.  
              
        points: A dictionary of locations ().
        
        variables:  List of variable(s) to interpolate such as ['airt', 'rh', 'geop', 'wind'].
                    Defaults to using all variables available.
        
        date: Directory to specify begin and end time for the derived time series. Defaluts to 
              using the entire time available in ncfile_in.
              
    Example:
        from datetime import datetime
        date  = {'beg' : datetime(2008, 1, 1),
                 'end' : datetime(2008,12,31)}
        variables  = ['ssrd','tp']       
        stations = csv_stations("points.csv")      
        int2Dpoint('big/era_sa_19790101_to_20161231.nc', 'era_sa_inter.nc', stations, variables=variables, date=date)        
    """   
        
    # open netcdf file handle
    ncf = nc.Dataset(ncfile_in, 'r')
    # is it a file with pressure levels?
    pl = 'level' in ncf.dimensions.keys()

    # get spatial dimensions
    lat  = ncf.variables['latitude'][:]
    lon  = ncf.variables['longitude'][:]
    if pl: #only for pressure level files
        lev  = ncf.variables['level'][:]
        nlev = len(lev)
    
    # get time and convert to datetime object
    nctime = ncf.variables['time'][:]
    t_unit = ncf.variables['time'].units # get unit  "hours since 1900-01-01 00:00:0.0"
    try :
        t_cal = ncf.variables['time'].calendar
    except AttributeError : # Attribute doesn't exist
        t_cal = u"gregorian" # or standard
    time = nc.num2date(nctime, units = t_unit, calendar = t_cal)

    # restrict to date/time range if given
    if date is None:
        tmask = time < datetime(3000, 1, 1)
    else:
        tmask = (time <= date['end']) * (time >= date['beg'])

    # test if time steps to interpolate remain
    nt = sum(tmask)
    if nt == 0:
        raise ValueError('No time steps from netCDF file selected.')
    
    # get variables
    varlist = [x.encode('UTF8') for x in ncf.variables.keys()]
    varlist.remove('time')
    varlist.remove('latitude')
    varlist.remove('longitude')
    if pl: #only for pressure level files
        varlist.remove('level')
    
    #list variables that should be interpolated
    if variables is None:
        variables = varlist
    #test is variables given are available in file
    if (set(variables) < set(varlist) == 0):
        raise ValueError('One or more variables not available in netCDF file.')
        
    # create source grid from a SCRIP formatted file
    sgrid = ESMF.Grid(filename=ncfile_in, filetype=ESMF.FileFormat.GRIDSPEC)

    # create source field on source grid
    if pl: #only for pressure level files
        sfield = ESMF.Field(sgrid, name='sgrid',
                            staggerloc=ESMF.StaggerLoc.CENTER,
                            ndbounds=[len(variables), nt, nlev])
    else: # 2D files
        sfield = ESMF.Field(sgrid, name='sgrid',
                            staggerloc=ESMF.StaggerLoc.CENTER,
                            ndbounds=[len(variables), nt])
                            
    # assign data from ncdf: (variale, time, latitude, longitude) 
    for n, var in enumerate(variables):
        if pl: # only for pressure level files
            sfield.data[n,:,:,:,:] = ncf.variables[var][tmask,:,:,:] 
        else:
            sfield.data[n,:,:,:] = ncf.variables[var][tmask,:,:] 

    # create locstream structure, CANNOT have third dimension for 2D interpolation
    locstream = ESMF.LocStream(len(points), coord_sys=ESMF.CoordSys.SPH_DEG)
    locstream["ESMF:Lon"] = [float(e['longitude_deg']) for e in points]
    locstream["ESMF:Lat"] = [float(e['latitude_deg']) for  e in points]
    
    # test if points are inside netCDF file grid
    geom = (min(locstream["ESMF:Lon"]) < min(lon) +
            max(locstream["ESMF:Lon"]) > max(lon) +
            min(locstream["ESMF:Lat"]) < min(lat) +
            max(locstream["ESMF:Lat"]) > max(lat))
    if geom > 0:
        raise ValueError('Points outside netCDF file.')

    # create destination field
    if pl: # only for pressure level files
        dfield = ESMF.Field(locstream, name='dfield', ndbounds=[len(variables), nt, nlev])
    else:
        dfield = ESMF.Field(locstream, name='dfield', ndbounds=[len(variables), nt])    

    # regridding function
    regrid2D = ESMF.Regrid(sfield, dfield,
                           regrid_method=ESMF.RegridMethod.BILINEAR,
                           unmapped_action=ESMF.UnmappedAction.ERROR,
                           dst_mask_values=None)
                          
    # regrid operation, create destination field (variables, times, points)
    dfield = regrid2D(sfield, dfield)                          
    
    # === write output netCDF file ==============================================================
    # dimensions: station, time OR station, time, level
    # variables: latitude(station), longitude(station), elevation(station), ...(time, level, station) 
    # stations are integer numbers
    # create a file (Dataset object, also the root group).
    rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.source      = 'ERA-Interim, interpolated bilinearly to stations'
    rootgrp.featureType = "timeSeries"

    # dimensions
    station = rootgrp.createDimension('station', len(points))
    time    = rootgrp.createDimension('time', nt)
    if pl: # only for pressure level files
        level = rootgrp.createDimension('level', nlev)

    # base variables
    time           = rootgrp.createVariable('time',     'i4',('time'))
    time.long_name = 'time'
    time.units     = 'hours since 1900-01-01 00:00:0.0'
    time.calendar  = 'gregorian'
    station             = rootgrp.createVariable('station',  'i4',('station'))
    station.long_name   = 'station for time series data'
    station.units       = '1'
    latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
    latitude.long_name  = 'latitude'
    latitude.units      = 'degrees_north'    
    longitude           = rootgrp.createVariable('longitude','f4',('station'))
    longitude.long_name = 'longitude'
    longitude.units     = 'degrees_east'  
    height           = rootgrp.createVariable('height','f4',('station'))
    height.long_name = 'height_above_reference_ellipsoid'
    height.units     = 'm'  
    if pl: # only for pressure level files
        level           = rootgrp.createVariable('level','i4',('level'))
        level.long_name = 'pressure_level'
        level.units     = 'millibars'  
       
    # assign base variables
    time[:] = nctime[tmask]
    if pl: # only for pressure level files
        level[:] = lev
    station[:]   = [float(e['station_number']) for e in points]
    latitude[:]  = [float(e['latitude_deg'])   for e in points]
    longitude[:] = [float(e['longitude_deg'])  for e in points]
    height[:]    = [float(e['elevation_m'])    for e in points]
    
    # create and assign variables from input file
    for n, var in enumerate(variables):
        vname = ncf.variables[var].name.encode('UTF8')
        if pl: # only for pressure level files
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'level', 'station'))
        else:
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))    
        tmp.long_name = ncf.variables[var].long_name.encode('UTF8')
        tmp.units     = ncf.variables[var].units.encode('UTF8')  
        # assign values
        if pl: # only for pressure level files
            tmp[:] = dfield.data[n,:,:,:]
        else:
            tmp[:] = dfield.data[n,:,:]    
    
    # close file
    rootgrp.close()
                 
       
                                                                                                                 
                                                                                        