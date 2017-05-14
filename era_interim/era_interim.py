#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Methods for downloading ERA-Interim data from the ECMWF server for limited
# areas and limited times.
#
#
# (C) Copyright Stephan Gruber (2013–2017)
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
#  OVERALL WORKFLOW
# 
#  from era_interim import *
#  from datetime import datetime
#
#  date  = {'beg' : datetime(1979, 1, 1),
#           'end' : datetime(2017, 1, 1)}
#  area  = {'north' :  40.0,
#           'south' :  41.0,
#           'west'  :  60.0,
#           'east'  :  61.0}
#  elevation = {'min' :   50, 
#               'max' : 2000}           
#  directory = '/home/sgruber/storage/Workplace/Stephan/Ekati'   
#          
#  ERAbat = ERAbatch(date, area, elevation, directory, 5) 
#  ERAbat.retrieve()
#
#  station_file = '/home/sgruber/storage/Workplace/Stephan/ESMF/points.csv'
#  ERAint = ERAinterp()
#  ERAint.all_files(directory, station_file)
#
# ECMWF and netCDF information:
# https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
#
# For variable codes and units, see: 
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
# Check ECMWF job status: http://apps.ecmwf.int/webmars/joblist/
#
#  CF-Convention: this is how you check netcdf file conformity: 
#  http://pumatest.nerc.ac.uk/cgi-bin/cf-checker-dev.pl 
#
#===============================================================================

from datetime import datetime, timedelta
from ecmwfapi.api import ECMWFDataServer
from math     import exp, floor
from os       import path
import numpy   as np
import csv
import ESMF
import netCDF4 as nc



class ERAgeneric(object):
    """
    Parent class for other ERA-Interim classes.
    """
        
    def areaString(self, area):
        """Converts numerical coordinates into string: North/West/South/East"""
        res  = str(round(area['north'],2)) + "/"
        res += str(round(area['west'], 2)) + "/"
        res += str(round(area['south'],2)) + "/"
        res += str(round(area['east'], 2))        
        return(res)
        
    def dateString(self, date):
        """Converts datetime objects into string"""
        res  = (date['beg'].strftime("%Y-%m-%d") + "/to/" +
                date['end'].strftime("%Y-%m-%d"))       
        return(res)    
        
    def getPressure(self, elevation):
        """Convert elevation into air pressure using barometric formula"""
        g  = 9.80665   #Gravitational acceleration [m/s2]
        R  = 8.31432   #Universal gas constant for air [N·m /(mol·K)]    
        M  = 0.0289644 #Molar mass of Earth's air [kg/mol]
        P0 = 101325    #Pressure at sea level [Pa]
        T0 = 288.15    #Temperature at sea level [K]
        #http://en.wikipedia.org/wiki/Barometric_formula
        return P0 * exp((-g * M * elevation) / (R * T0)) / 100 #[hPa] or [bar]
    
    def getPressureLevels(self, elevation):
        """Restrict list of ERA-interim pressure levels to be downloaded"""
        Pmax = self.getPressure(elevation['min']) + 55
        Pmin = self.getPressure(elevation['max']) - 55 
        levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 
                         800, 825, 850, 875, 900, 925, 950, 975, 1000])
        mask = (levs >= Pmin) * (levs <= Pmax) #select
        levs = '/'.join(map(str, levs[mask]))
        return levs

    def getDictionaryGen(self, area, date):
        """
        Makes dictionary of generic variables for a server call
        """
        dictionary_gen = {
          'area'    : self.areaString(area),
          'date'    : self.dateString(date),
          'dataset' : "interim",
          'stream'  : "oper",
          'class'   : "ei",
          'format'  : "netcdf",
          'grid'    : "0.75/0.75"}
        return dictionary_gen
 
    def getDstring(self):
        return ('_' + self.date['beg'].strftime("%Y%m%d") + "_to_" +
                      self.date['end'].strftime("%Y%m%d"))
    
    def download(self):
        #TODO test for file existence
        server = ECMWFDataServer()
        print server.trace('=== ERA Interim: START ====')
        server.retrieve(self.getDictionary())
        print server.trace('=== ERA Interim: STOP =====')

    def getNCDF(self):
        return self.file_ncdf                       
        
    def __str__(self):
        string = ("List of generic variables to query ECMWF server for "
                  "ERA-Interim data: {0}")
        return string.format(self.getDictionary) 
                 
        
class ERApl(ERAgeneric):
    """Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server.
       
    Args:
        date: A dictionary specifying the time period desired with a begin 
              and an end date given as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg].  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].
        
        directory: Directory to hold output files
              
    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(2013, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  15.0,
                 'west'  :  60.0,
                 'east'  : 105.0}
        elevation = {'min' :    0, 
                     'max' : 8850}
        variable  = ['airt','rh','wind','geop']             
        directory = '/Users/stgruber/Desktop'             
        ERApl = ERApl(date, area, elevation, variable, directory) 
        ERApl.download()
    """
    def __init__(self, date, area, elevation, variable, directory):
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'era_pl.nc')
        dpar = {'airt' : '130.128',           # [K]
                'rh'   : '157.128',           # [%]
                'geop' : '129.128',           # [m2 s-2]
                'wind' : '131.128/132.128'}   # [m s-1]
        self.param = ''        
        for var in variable:
            self.param += dpar.get(var)+'/'
        self.param = self.param.rstrip('/') #fix last  
    
    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "pl",
           'levellist': self.getPressureLevels(self.elevation),
           'time'     : "00/06/12/18",
           'step'     : "0",
           'type'     : "an",
           'param'    : self.param,
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary    
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary) 
        
        
        
class ERAsa(ERAgeneric):
    """
    Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server for surface analysis variables (airt2,dewp2, ozone, vapor,
    wind10).
       
    Args:
        
        target:    File name of the netcdf file to be created.
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind']
              
    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(2013, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  15.0,
                 'west'  :  60.0,
                 'east'  : 105.0}
        variable  = ['airt2', 'dewp2', 'wind10', 'ozone', 'vapor']             
        directory = '/Users/stgruber/Desktop'             
        ERAsa = ERAsa(date, area, variable, directory) 
        ERAsa.download()      
    """
    def __init__(self, date, area, variable, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'era_sa.nc')
        dpar = {'airt2'  : '167.128',           # [K] 2m values
                'dewp2'  : '168.128',           # [K] 2m values
                'ozone'  : '206.128',           # [kg m-2] Total column ozone
                'vapor'  : '137.128',           # [kg m-2] Total column W vapor                                                                
                'wind10' : '165.128/166.128'}   # [m s-1] 10m values
        self.param = ''        
        for var in variable:
            self.param += dpar.get(var)+'/'
        self.param = self.param.rstrip('/') #fix last


    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "00/06/12/18",
           'step'     : "0",
           'type'     : "an",
           'param'    : self.param,
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary  
        
    def __str__(self):
        string = ("Class for ERA-Interim surface analysis data: {0}")
        return string.format(self.getDictionary)         




class ERAsf(ERAgeneric):
    """
    Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server for surface forecast variables (prec, swin, lwin).
       
    Args:        
        target:    File name of the netcdf file to be created.
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].
        
        ERAgen:    ERAgeneric() object with generic information on the area and
                   time span to be
              
    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(1994, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  40.0,
                 'west'  :  60.0,
                 'east'  :  60.0}
        variable  = ['prec','swin','lwin']             
        directory = '/Users/stgruber/Desktop'             
        ERAsf = ERAsf(date, area, variable, directory) 
        ERAsf.download()   
    """
    def __init__(self, date, area, variable, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'era_sf.nc')
        dpar = {'prec'   : '228.128',   # [m] total precipitation
                'swin'   : '169.128',   # [J m-2] short-wave downward
                'lwin'   : '175.128'}   # [J m-2] long-wave downward
        self.param = ''        
        for var in variable:
            self.param += dpar.get(var)+'/'
        self.param = self.param.rstrip('/') #fix last


    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "00/12",
           'step'     : "3/6/9/12",
           'type'     : "fc",
           'param'    : self.param,
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary  
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary)      


                
class ERAto(ERAgeneric):
    """
    Returns an object for downloading and handling ERA-Interim 
    topography (invariant).
       
    Args:      
              
    Example:
        area  = {'north' :  40.0,
                 'south' :  45.0,
                 'west'  :  60.0,
                 'east'  :  65.0}            
        directory = '/Users/stgruber/Desktop'             
        ERAto = ERAto(area, directory) 
        ERAto.download()       
    """
    def __init__(self, area, directory):
        self.area       = area
        self.date       = {'beg' : datetime(1979, 1, 1),
                           'end' : datetime(1979, 1, 1)}
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'era_to.nc')

    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "12",
           'step'     : "0",
           'type'     : "an",
           'param'    : "129.128",
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary   
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary) 
                 

class ERAinterp(object):
    """
    Collection of methods to interpolate ERA-Interim netCDF files to station
    coordinates. All variables retain theit original units and time stepping.
    """
    
    def csv_stations(self, filename):
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
    
    
    def ERA2station(self, ncfile_in, ncfile_out, points, 
                    variables=None, date=None):    
        """
        Biliner interpolation from fields on regular grid (latitude, longitude) 
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all ERA-Interim files).
          
        Args:
            ncfile_in: Full path to am ERA-Interim derived netCDF file.
              
            ncfile_out: Full path to the output netCDF file to write.  
              
            points: A dictionary of locations. See method csv_stations for more
                    details.
        
            variables:  List of variable(s) to interpolate such as 
                        ['airt', 'rh', 'geop', 'wind'].
                        Defaults to using all variables available.
        
            date: Directory to specify begin and end time for the derived time 
                  series. Defaluts to using all times available in ncfile_in.
              
        Example:
            from datetime import datetime
            date  = {'beg' : datetime(2008, 1, 1),
                     'end' : datetime(2008,12,31)}
            variables  = ['ssrd','tp']       
            stations = csv_stations("points.csv")      
            int2Dpoint('era_sa.nc', 'era_sa_inter.nc', stations, 
                       variables=variables, date=date)        
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
        t_unit = ncf.variables['time'].units #"hours since 1900-01-01 00:00:0.0"
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
            raise ValueError('One or more variables not in netCDF file.')
        
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
                sfield.data[n,:,:,:,:] = ncf.variables[var][tmask,:,:,:].transpose((0,1,3,2)) 
            else:
                sfield.data[n,:,:,:] = ncf.variables[var][tmask,:,:].transpose((0,2,1)) 

        # create locstream, CANNOT have third dimension!!!
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
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt, nlev])
        else:
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt])    

        # regridding function
        regrid2D = ESMF.Regrid(sfield, dfield,
                               regrid_method=ESMF.RegridMethod.BILINEAR,
                               unmapped_action=ESMF.UnmappedAction.ERROR,
                               dst_mask_values=None)
                          
        # regrid operation, create destination field (variables, times, points)
        dfield = regrid2D(sfield, dfield)        
        sfield.destroy() #f ree memory                  
    
        # === write output netCDF file =========================================
        # dimensions: station, time OR station, time, level
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, level, station) or (time, station)
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
                tmp   = rootgrp.createVariable(vname,
                                               'f4',('time', 'level', 'station'))
            else:
                tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))   
                 
            tmp.long_name = ncf.variables[var].long_name.encode('UTF8')
            tmp.units     = ncf.variables[var].units.encode('UTF8')  
            # assign values
            if pl: # only for pressure level files
                tmp[:] = dfield.data[n,:,:,:]
            else:
                tmp[:] = dfield.data[n,:,:]    
    
        rootgrp.close()
        ncf.close()
        
        # closed file ==========================================================
        
    def levels2elevation(self, ncfile_in, ncfile_out):    
        """
        Linear 1D interpolation of pressure level data available for individual
        stations to station elevation. Where and when stations are below the 
        lowest pressure level, they are assigned the value of the lowest 
        pressure level.
        
        """
        # open file 
        ncf = nc.Dataset(ncfile_in, 'r')
        height = ncf.variables['height'][:]
        nt = len(ncf.variables['time'][:])
        nl = len(ncf.variables['level'][:])
        
        # list variables
        varlist = [x.encode('UTF8') for x in ncf.variables.keys()]
        varlist.remove('time')
        varlist.remove('station')
        varlist.remove('latitude')
        varlist.remove('longitude')
        varlist.remove('level')
        varlist.remove('height')
        varlist.remove('z')

        # === open and prepare output netCDF file ==============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'ERA-Interim, interpolated (bi)linearly to stations'
        rootgrp.featureType = "timeSeries"

        # dimensions
        station = rootgrp.createDimension('station', len(height))
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
        time[:] = ncf.variables['time'][:]
        station[:]   = ncf.variables['station'][:]
        latitude[:]  = ncf.variables['latitude'][:]
        longitude[:] = ncf.variables['longitude'][:]
        height[:]    = ncf.variables['height'][:]
        
        # create and assign variables from input file
        for var in varlist:
            vname = ncf.variables[var].name.encode('UTF8')
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))    
            tmp.long_name = ncf.variables[var].long_name.encode('UTF8')
            tmp.units     = ncf.variables[var].units.encode('UTF8')  
        # end file prepation ===================================================
    
                                                                                                
        # loop over stations
        for n, h in enumerate(height): 
            # convert geopotential [millibar] to height [m]
            # shape: (time, level)
            ele = ncf.variables['z'][:,:,n] / 9.80665
            # TODO: check if height of stations in data range (+50m at top, lapse r.)
            
            # difference in elevation. 
            # level directly above will be >= 0
            dele = ele - h
            # vector of level indices that fall directly above station. 
            # Apply after ravel() of data.
            va = np.argmin(dele + (dele < 0) * 100000, axis=1) 
            # mask for situations where station is below lowest level
            mask = va < (nl-1)
            va += np.arange(ele.shape[0]) * ele.shape[1]
            
            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask # +1 when OK, +0 when below lowest level
            
            # weights
            wa = np.absolute(dele.ravel()[vb]) 
            wb = np.absolute(dele.ravel()[va])
            
            wt = wa + wb
            
            wa /= wt # Apply after ravel() of data.
            wb /= wt # Apply after ravel() of data.
            
            #loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                #read data from netCDF
                data = ncf.variables[var][:,:,n].ravel()
                ipol = data[va]*wa + data[vb]*wb   # interpolated value                    
                rootgrp.variables[var][:,n] = ipol # assign to file   
    
        rootgrp.close()
        # closed file ==========================================================    

    def all_files(self, directory, stations_csv):
        #TODO
        """
        Interpolate point time series from downloaded data.
        """
        
        #read station points      
        stations = self.csv_stations(stations_csv)      
        
        #interpolation
        self.ERA2station(path.join(directory,'era_pl.nc'), 
                         path.join(directory,'era_pl_stations.nc'), stations)
        self.levels2elevation(path.join(directory,'era_pl_stations.nc'), 
                              path.join(directory,'era_pl_stations_surface.nc'))
        self.ERA2station(path.join(directory,'era_sa.nc'), 
                         path.join(directory,'era_sa_stations.nc'), stations)
        self.ERA2station(path.join(directory,'era_sf.nc'), 
                        path.join(directory,'era_sf_stations.nc'), stations)
        self.ERA2station(path.join(directory,'era_to.nc'), 
                         path.join(directory,'era_to_stations.nc'), stations)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        
        
class ERAbatch(object):
    """
    Returns an object for ERA-Interim data that has methods for querying 
    the ECMWF server, returning all variables usually needed.
       
    Args:
        date: A dictionary specifying the time period desired with a begin 
              and an end date given as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east.  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
                
        directory: Directory to hold output files
        
        increment_days: How many days should be requested from ECMWF 
                        server per increment? This helps to keep the 
                        processing faster.
        
        n_outfile: Over how many nc files should the final output be
                   distributed? The dafault is 1, but for very large
                   areas it may be better to use several files. 
              
    Example:
        import era_download as ed
        from datetime import datetime
        date  = {'beg' : datetime(1980, 1, 1),
                 'end' : datetime(2015, 1, 2)}
        area  = {'north' :  40.0,
                 'south' :  41.0,
                 'west'  :  60.0,
                 'east'  :  61.0}
        elevation = {'min' :   50, 
                     'max' : 2000}           
        directory = '/Users/stgruber/Desktop/aaa'             
        batch = ed.ERAbatch(date, area, elevation, directory, 5) 
        batch.retrieve()
    """
        
    def __init__(self, date, area, elevation, directory):#, 
#                 increment_days, n_outfile=1):
        self.date      = date
        self.area      = area
        self.elevation = elevation
        self.directory = directory
#        self.increment = increment_days
        self.nc_files  = ''
 #       self.n_outfile = n_outfile
        #TODO: check directory
        #TODO ensure increments is smaller or equal than chosen time window

    def getFileNames(self):  
        return self.nc_files
        #TODO: add methods to get file names with * in them
              
    def retrieve(self):
        #define variables

        var_pl = ['airt', 'rh', 'wind', 'geop']  
        var_sa = ['airt2', 'dewp2', 'wind10']
        var_sf = ['prec', 'swin', 'lwin']

#        #enter time loop, assign dummy date_i to start with
#        date_i = {'beg' : datetime(1994, 1, 1), 'end' : datetime(1999, 1, 2)}
#        slices = floor(float((self.date['end'] - self.date['beg']).days)/
#                       self.increment)+1
#
#        for ind in range (0, int(slices)): 
#            #prepare time slices   
#            date_i['beg'] = self.date['beg'] + timedelta(days = 
#                            self.increment * ind)
#            date_i['end'] = self.date['beg'] + timedelta(days = 
#                            self.increment * (ind+1) - 1)
#            if ind == (slices -1):
#                date_i['end'] = self.date['end']
#            
#            #actual functions                                                                           
#            pl = ERApl(date_i, self.area, self.elevation, 
#                       var_pl, self.directory) 
#            sa = ERAsa(date_i, self.area, var_sa, self.directory) 
#            sf = ERAsf(date_i, self.area, var_sf, self.directory) 
#            self.ERAli   = [pl, sa, sf] #combine in list
#        
#            #download from ECMWF server convert to netCDF  
#            for era in self.ERAli:
#                era.download()

            
        #actual functions                                                                           
        pl = ERApl(self.date, self.area, self.elevation, 
                   var_pl, self.directory) 
        sa = ERAsa(self.date, self.area, var_sa, self.directory) 
        sf = ERAsf(self.date, self.area, var_sf, self.directory) 
        self.ERAli   = [pl, sa, sf] #combine in list
        
        #download from ECMWF server convert to netCDF  
        for era in self.ERAli:
            era.download()            
                                         
        #topography
        top = ERAto(self.area, self.directory)
        top.download()
                                       
    def append(self):
        #TODO
        '''
        Append data to files in a directory.
        '''                                               
                    
    def __str__(self):
        return "Object for ERA-Interim data download and conversion"                

