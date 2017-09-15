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
#  See: https://github.com/geocryology/globsim/wiki/Globsim
#       and the examples directory on github
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
from os       import path, listdir
from generic import ParameterIO, StationListRead, ScaledFileOpen
from fnmatch import filter

import numpy   as np
import netCDF4 as nc

try:
    import ESMF
except ImportError:
    print("*** ESMF not imported, interpolation not possible. ***")
    pass   
            

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

    def TranslateCF2ERA(self, variables, dpar):
        """
        Translate CF Standard Names into ERA-Interim code numbers.
        """
        self.param = ''
        for var in variables:
            try:
                self.param += dpar.get(var)+'/'
            except TypeError:
                pass    
        self.param = self.param.rstrip('/') #fix last                                        
                
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
        
        variables:  List of variable(s) to download based on CF Standard Names.
        
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
        variables  = ['air_temperature', 'relative_humidity']             
        directory = '/Users/stgruber/Desktop'             
        ERApl = ERApl(date, area, elevation, variables, directory) 
        ERApl.download()
    """
    def __init__(self, date, area, elevation, variables, directory):
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        outfile = 'era_pl' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)
 
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable names. 
        dpar = {'air_temperature'   : '130.128',           # [K]
                'relative_humidity' : '157.128',           # [%]
                'wind_speed'        : '131.128/132.128'}   # [m s-1]
                
        # translate variables into those present in ERA pl data        
        self.TranslateCF2ERA(variables, dpar)
        self.param += '/129.128' # geopotential always needed [m2 s-2] 
    
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
                  "ERA-Interim data: {0}")
        return string.format(self.getDictionary) 
        
        
        
class ERAsa(ERAgeneric):
    """
    Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server for surface analysis variables (airt2,dewp2, ozone, vapor,
    wind10).
       
    Args:
        
        target:    File name of the netcdf file to be created.
        
        variable:  List of variable(s) to download based on CF Standard Names 
              
    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(2013, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  15.0,
                 'west'  :  60.0,
                 'east'  : 105.0}
        variables  = ['air_temperature', 'relative_humidity']             
        directory = '/Users/stgruber/Desktop'             
        ERAsa = ERAsa(date, area, variables, directory) 
        ERAsa.download()      
    """
    def __init__(self, date, area, variables, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile = 'era_sa' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)
        
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable names. 
        dpar = {'air_temperature'   : '167.128',  # [K] 2m values
                'relative_humidity' : '168.128',  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky' : 
                    '206.128/137.128',  # [kg m-2] Total column ozone 
                                        # [kg m-2] Total column W vapor                                                             
                'wind_speed' : '165.128/166.128'}   # [m s-1] 10m values
                
        # translate variables into those present in ERA pl data        
        self.TranslateCF2ERA(variables, dpar)
        print(self.param)


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
        
        variables:  List of variable(s) to download that can include one, several
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
        variables  = ['prec','swin','lwin']             
        directory = '/Users/stgruber/Desktop'             
        ERAsf = ERAsf(date, area, variables, directory) 
        ERAsf.download()   
    """
    def __init__(self, date, area, variables, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile = 'era_sf' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable names. 
        dpar = {'precipitation_amount'              : '228.128',   # [m] total precipitation
                'downwelling_shortwave_flux_in_air' : '169.128',   # [J m-2] short-wave downward
                'downwelling_longwave_flux_in_air'  : '175.128'}   # [J m-2] long-wave downward
                
        # translate variables into those present in ERA pl data        
        self.TranslateCF2ERA(variables, dpar)
        print(self.param)


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
           'param'    : "129.128/172.128", # geopotential and land-sea mask
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary   
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary) 
                 

class ERAinterpolate(object):
    """
    Collection of methods to interpolate ERA-Interim netCDF files to station
    coordinates. All variables retain theit original units and time stepping.
    """
        
    def __init__(self, ifile):
        #read parameter file
        self.ifile = ifile
        par = ParameterIO(self.ifile)
        self.dir_inp = path.join(par.project_directory,'eraint') 
        self.dir_out = path.join(par.project_directory,'station')
        self.variables = par.variables
        self.list_name = par.list_name
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        
        #read station points 
        self.stations = StationListRead(self.stations_csv)  
        #convert longitude to ERA notation if using negative numbers  
        self.stations['longitude_dd'] = self.stations['longitude_dd'] % 360             
        
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}
    
    def ERA2station(self, ncfile_in, ncfile_out, points,
                    variables=None, date=None):    
        """
        Biliner interpolation from fields on regular grid (latitude, longitude) 
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all ERA-Interim files).
          
        Args:
            ncfile_in: Full path to am ERA-Interim derived netCDF file. This can
                       contain wildcards to point to multiple files if temporal
                       chunking was used.
              
            ncfile_out: Full path to the output netCDF file to write.  
              
            points: A dictionary of locations. See method StationListRead in
                    generic.py for more details.
        
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
            stations = StationListRead("points.csv")      
            ERA2station('era_sa.nc', 'era_sa_inter.nc', stations, 
                       variables=variables, date=date)        
        """   
        # open netcdf file handle, can be one file of several with wildcards
        ncf = nc.MFDataset(ncfile_in, 'r')
        
        # is it a file with pressure levels?
        pl = 'level' in ncf.dimensions.keys()

        # get spatial dimensions
        lat  = ncf.variables['latitude'][:]
        lon  = ncf.variables['longitude'][:]
        if pl: # only for pressure level files
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
        
        # Create source grid from a SCRIP formatted file. As ESMF needs one
        # file rather than an MFDataset, give first file in directory.
        ncsingle = filter(listdir(self.dir_inp), path.basename(ncfile_in))[0]
        ncsingle = path.join(self.dir_inp, ncsingle)
        sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)

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
        locstream = ESMF.LocStream(len(self.stations), coord_sys=ESMF.CoordSys.SPH_DEG)
        locstream["ESMF:Lon"] = list(self.stations['longitude_dd'])
        locstream["ESMF:Lat"] = list(self.stations['latitude_dd'])

        # create destination field
        if pl: # only for pressure level files
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt, nlev])
        else:
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt])    

        # regridding function, consider ESMF.UnmappedAction.ERROR
        regrid2D = ESMF.Regrid(sfield, dfield,
                               regrid_method=ESMF.RegridMethod.BILINEAR,
                               unmapped_action=ESMF.UnmappedAction.IGNORE,
                               dst_mask_values=None)
                          
        # regrid operation, create destination field (variables, times, points)
        dfield = regrid2D(sfield, dfield)        
        sfield.destroy() #free memory                  
		
        # === write output netCDF file =========================================
        # dimensions: station, time OR station, time, level
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, level, station) or (time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4_CLASSIC')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'ERA-Interim, interpolated bilinearly to stations'
        rootgrp.featureType = "timeSeries"

        # dimensions
        station = rootgrp.createDimension('station', len(self.stations))
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
        station[:]   = list(self.stations['station_number'])
        latitude[:]  = list(self.stations['latitude_dd'])
        longitude[:] = list(self.stations['longitude_dd'])
        height[:]    = list(self.stations['elevation_m'])
    
        # create and assign variables from input file
        for n, var in enumerate(variables):
            vname = ncf.variables[var].long_name.encode('UTF8')
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
        ncf = nc.MFDataset(ncfile_in, 'r', aggdim='time')
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
        varlist.remove('Geopotential')

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
            vname = ncf.variables[var].long_name.encode('UTF8')
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))    
            tmp.long_name = ncf.variables[var].long_name.encode('UTF8')
            tmp.units     = ncf.variables[var].units.encode('UTF8')  
        # end file prepation ===================================================
    
                                                                                                
        # loop over stations
        for n, h in enumerate(height): 
            # convert geopotential [millibar] to height [m]
            # shape: (time, level)
            ele = ncf.variables['Geopotential'][:,:,n] / 9.80665
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


    def TranslateCF2short(self, dpar):
        """
        Map CF Standard Names into short codes used in ERA-Interim netCDF files.
        """
        varlist = [] 
        for var in self.variables:
            varlist.append(dpar.get(var))
        # drop none
        varlist = [item for item in varlist if item is not None]      
        # flatten
        varlist = [item for sublist in varlist for item in sublist]         
        return(varlist) 
    
    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to 
        the more generically ERA-like interpolation functions.
        """                       

        # === 2D Interpolation for Surface Analysis Data ===    
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['t2m'],  # [K] 2m values
                'relative_humidity' : ['d2m'],  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky' : 
                    ['tco3', 'tcwv'],   # [kg m-2] Total column ozone 
                                        # [kg m-2] Total column W vapor                                                             
                'wind_speed' : ['u10', 'v10']}   # [m s-1] 10m values   
        varlist = self.TranslateCF2short(dpar)                      
        self.ERA2station(path.join(self.dir_inp,'era_sa_*.nc'), 
                         path.join(self.dir_out,'era_sa_' + 
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)          
        
        # 2D Interpolation for Surface Forecast Data    'tp', 'strd', 'ssrd' 
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys.       
        dpar = {'precipitation_amount'              : ['tp'],   # [m] total precipitation
                'downwelling_shortwave_flux_in_air' : ['ssrd'], # [J m-2] short-wave downward
                'downwelling_longwave_flux_in_air'  : ['strd']} # [J m-2] long-wave downward
        varlist = self.TranslateCF2short(dpar)                           
        self.ERA2station(path.join(self.dir_inp,'era_sf_*.nc'), 
                        path.join(self.dir_out,'era_sf_' + 
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)          
                        
        # 2D Interpolation for Invariant Data      
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys.            
        dummy_date  = {'beg' : datetime(1979, 1, 1, 12, 0),
                       'end' : datetime(1979, 1, 1, 12, 0)}        
        self.ERA2station(path.join(self.dir_inp,'era_to.nc'), 
                         path.join(self.dir_out,'era_to_' + 
                                   self.list_name + '.nc'), self.stations,
                                   ['z', 'lsm'], date = dummy_date)    
        
        
        # === 2D Interpolation for Pressure Level Data ===
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['t'],           # [K]
                'relative_humidity' : ['r'],           # [%]
                'wind_speed'        : ['u', 'v']}    # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('z')
        self.ERA2station(path.join(self.dir_inp,'era_pl_*.nc'), 
                         path.join(self.dir_out,'era_pl_' + 
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        # 1D Interpolation for Pressure Level Data
        self.levels2elevation(path.join(self.dir_out,'era_pl_' + 
                                        self.list_name + '.nc'), 
                              path.join(self.dir_out,'era_pl_' + 
                                        self.list_name + '_surface.nc'))
        
class ERAdownload(object):
    """
    Class for ERA-Interim data that has methods for querying 
    the ECMWF server, returning all variables usually needed.
       
    Args:
        pfile: Full path to a Globsim Download Parameter file. 
              
    Example:          
        ERAd = ERAdownload(pfile) 
        ERAd.retrieve()
    """
        
    def __init__(self, pfile):
        # read parameter file
        self.pfile = pfile
        par = ParameterIO(self.pfile)
        
        # assign bounding box
        self.area  = {'north':  par.bbN,
                      'south':  par.bbS,
                      'west' :  par.bbW,
                      'east' :  par.bbE}
                 
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}

        # elevation
        self.elevation = {'min' : par.ele_min, 
                          'max' : par.ele_max}
        
        # data directory for ERA-Interim  
        self.directory = path.join(par.project_directory, "eraint")  
        if path.isdir(self.directory) == False:
            raise ValueError("Directory does not exist: " + self.directory)   
     
        # variables
        self.variables = par.variables
            
        # chunk size for downloading and storing data [days]        
        self.chunk_size = par.chunk_size            
              
    def retrieve(self):
        """
        Retrieve all required ERA-Interim data from MARS server.
        """
        #TODO: append to dataset while keeping chunk size
        
        # prepare time loop
        date_i = {}
        slices = floor(float((self.date['end'] - self.date['beg']).days)/
                       self.chunk_size)+1

        for ind in range (0, int(slices)): 
            #prepare time slices   
            date_i['beg'] = self.date['beg'] + timedelta(days = 
                            self.chunk_size * ind)
            date_i['end'] = self.date['beg'] + timedelta(days = 
                            self.chunk_size * (ind+1) - 1)
            if ind == (slices-1):
                date_i['end'] = self.date['end']
            
            #actual functions                                                                           
            pl = ERApl(date_i, self.area, self.elevation, 
                       self.variables, self.directory) 
            sa = ERAsa(date_i, self.area, self.variables, self.directory) 
            sf = ERAsf(date_i, self.area, self.variables, self.directory) 
        
            #download from ECMWF server convert to netCDF  
            ERAli = [pl, sa, sf]
            for era in ERAli:
                era.download()          
                                         
        # topography
        top = ERAto(self.area, self.directory)
        top.download()
        
        # report inventory
        self.inventory()     
                                                                                                                                                                                                     
    def inventory(self):
        """
        Report on data avaialbe in directory: time slice, variables, area 
        """
        print("\n\n\n")
        print("=== INVENTORY FOR GLOBSIM ERA-INTERIM DATA === \n")
        print("Download parameter file: \n" + self.pfile + "\n")
        # loop over filetypes, read, report
        file_type = ['era_pl_*.nc', 'era_sa_*.nc', 'era_sf_*.nc', 'era_t*.nc']
        for ft in file_type:
            infile = path.join(self.directory, ft)
            nf = len(filter(listdir(self.directory), ft))
            print(str(nf) + " FILE(S): " + infile)
            
            if nf > 0:
                # open dataset
                ncf = nc.MFDataset(infile, 'r')
                
                # list variables
                keylist = [x.encode('UTF8') for x in ncf.variables.keys()]
                print("    VARIABLES:")
                print("        " + str(len(keylist)) + 
                      " variables, inclusing dimensions")
                for key in keylist:
                    print("        " + ncf.variables[key].long_name)
                
                # time slice
                time = ncf.variables['time']
                tmin = '{:%Y/%m/%d}'.format(nc.num2date(min(time[:]), 
                                     time.units, calendar=time.calendar))
                tmax = '{:%Y/%m/%d}'.format(nc.num2date(max(time[:]), 
                                     time.units, calendar=time.calendar))
                print("    TIME SLICE")                     
                print("        " + str(len(time[:])) + " time steps")
                print("        " + tmin + " to " + tmax)
                      
                # area
                lon = ncf.variables['longitude']
                lat = ncf.variables['latitude']
                nlat = str(len(lat))
                nlon = str(len(lon))
                ncel = str(len(lat) * len(lon))
                print("    BOUNDING BOX / AREA")
                print("        " + ncel + " cells, " + nlon + 
                      " W-E and " + nlat + " S-N")
                print("        N: " + str(max(lat)))
                print("        S: " + str(min(lat)))
                print("        W: " + str(min(lon)))
                print("        E: " + str(max(lon)))
                            
                ncf.close()
        
    def __str__(self):
        return "Object for ERA-Interim data download and conversion"                


class ERAscale(object):
    """
    Class for ERA-Interim data that has methods for scaling station data to
    better resemble near-surface fluxes.
    
    Processing kernels have names in UPPER CASE.
       
    Args:
        sfile: Full path to a Globsim Scaling Parameter file. 
              
    Example:          
        ERAd = ERAscale(sfile) 
        ERAd.process()
    """
        
    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        par = ParameterIO(self.sfile)
        
        # read kernels
        self.kernels = par.kernels
        if not isinstance(self.kernels, list):
            self.kernels = [self.kernels]
            
        # input file names
        self.nc_pl = nc.Dataset(path.join(par.project_directory,'eraint/era_pl_' + 
                                par.list_name + '_surface.nc'), 'r')
        self.nc_sa = nc.Dataset(path.join(par.project_directory,'eraint/era_sa_' + 
                                par.list_name + '.nc'), 'r')
        self.nc_sf = nc.Dataset(path.join(par.project_directory,'eraint/era_sf_' + 
                                par.list_name + '.nc'), 'r')
        self.nc_to = nc.Dataset(path.join(par.project_directory,'eraint/era_to_' + 
                                par.list_name + '.nc'), 'r')
                               
        # output file 
        self.outfile = par.output_file  
        
        # time vector for output data
        # get time and convert to datetime object
        nctime = self.nc_pl.variables['time'][:]
        self.t_unit = self.nc_pl.variables['time'].units #"hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_pl.variables['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)
        
        #number of time steps
        nt = int(floor((max(time) - min(time)).total_seconds() 
                       / 3600 / par.time_step))
        
        # vector of output time steps as datetime object
        self.times_out    = [min(time) + timedelta(hours=x) for x in range(0, nt)]
        # vector of output time steps as written in ncdf file
        self.times_out_nc = nc.date2num(self.times_out, units = self.t_unit, 
                                        calendar = self.t_cal)

        
    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one 
        variable and adds it to the netCDF file.
        """    
        self.rg = ScaledFileOpen(self.outfile, self.nc_pl, self.times_out_nc)
        
        # iterate thorugh kernels and start process
        for kernel_name in self.kernels:
            getattr(self, kernel_name)()
            
        # close netCDF files   
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()
        
    def AIRT_ERA_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """        
        print("AIRT_ERA_pl")
        
    def AIRT_ERA_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """   
        
        # add variable to ncdf file
        vn = 'AIRT_ERA_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '2_metre_temperature ERA-I surface only'
        var.units     = self.nc_sa.variables['2 metre temperature'].units.encode('UTF8')  
        
        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['2 metre temperature'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n])            
        
        
    def AIRT_ERA_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP.
        """       
        print("AIRT_ERA_redcapp")            