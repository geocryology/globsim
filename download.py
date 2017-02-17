#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# (C) Copyright Stephan Gruber & Bin Cao
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# For variable codes and units, see:
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#cc
#===============================================================================

from datetime import datetime, timedelta
from ecmwfapi import ECMWFDataServer
from os       import path, remove
from math     import exp, floor

import numpy   as np
import pygrib  as pg
import netCDF4 as nc
import glob    as gl

# Check ECMWF job status: http://apps.ecmwf.int/webmars/joblist/
#TODO multiprocessing
#TODO -interpolate per ERA-gridcell nearest to center of points given
#     -stratify list of points to do that



class ERAgeneric(object):
    """Parent class for other ERA-Interim classes.
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
        levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
                         775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])
        mask = (levs >= Pmin) * (levs <= Pmax) #select
        levs = '/'.join(map(str, levs[mask]))
        return levs

    def getDictionaryGen(self, area, date):
        """Makes dictionary of generic variables for a server call"""
        dictionary_gen = {
          'area'    : self.areaString(area),
          'date'    : self.dateString(date),
          'dataset' : "interim",
          'stream'  : "oper",
          'class'   : "ei",
          'grid'    : "0.75/0.75"}
        return dictionary_gen

    def getDstring(self):
        return ('_' + self.date['beg'].strftime("%y%m%d") + "_to_" +
                      self.date['end'].strftime("%y%m%d"))

    def download(self):
        #TODO test for file existence
        server = ECMWFDataServer()
        print server.trace('=== ERA Interim: START ====')
        server.retrieve(self.getDictionary())
        print server.trace('=== ERA Interim: STOP =====')

    def toNCDF(self):
        gribFile(self.file_grib).toNCDF(self.file_ncdf)
        self.deleteGrib()

    def getNCDF(self):
        return self.file_ncdf

    def deleteGrib(self):
        remove(self.file_grib)

    def __str__(self):
        string = ("List of generic variables to query ECMWF server for "
                  "ERA-Interim data: {0}")
        return string.format(self.getDictionary)


class ERApl(ERAgeneric):
    """Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server, for converting grib to ncdf.

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
        ERApl.toNCDF()
    """
    
    def __init__(self, date, area, elevation, variable, directory):
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        self.file_grib  = path.join(self.directory,'ecmwf_erai_pl'+
                                    self.getDstring()+'.grib')
        self.file_ncdf  = path.join(self.directory,'ecmwf_erai_pl'+
                                    self.getDstring()+'.nc')
        dpar = {'airt' : '130.128',           # [K]
                'geop' : '129.128'}   
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
           'target'   : self.file_grib
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary)



class ERAsa(ERAgeneric):
    """Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server, for converting grib to ncdf.

    Args:

        target:    File name of the grib file to be created.

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
        ERAsa.toNCDF()
    """
    def __init__(self, date, area, variable, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        self.file_grib  = path.join(self.directory,'ecmwf_erai_sa'+
                                    self.getDstring()+'.grib')
        self.file_ncdf  = path.join(self.directory,'ecmwf_erai_sa'+
                                    self.getDstring()+'.nc')
        dpar = {'airt2'  : '167.128'}          # [K] 2m values
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
           'target'   : self.file_grib
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("Class for ERA-Interim surface analysis data: {0}")
        return string.format(self.getDictionary)







class ERAto(ERAgeneric):
    """Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server, for converting grib to ncdf.

    Args:
        area: download era-interim area
        directory: directory to save era-interim d

    Example:
        area  = {'north' :  40.0,
                       'south' :  45.0,
                       'west'  :  60.0,
                       'east'  :   65.0}
        directory = '/Users/stgruber/Desktop'
        ERAto = ERAto(area, directory)
        ERAto.download()
        ERAto.toNCDF()
    """
    def __init__(self, area, directory):
        self.area       = area
        self.date       = {'beg' : datetime(1979, 1, 1),
                                    'end' : datetime(1979, 1, 1)}
        self.directory  = directory
        self.file_grib  = path.join(self.directory,'ecmwf_erai_to.grib')
        self.file_ncdf  = path.join(self.directory,'ecmwf_erai_to.nc')

    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "12",
           'step'     : "0",
           'type'     : "an",
           'param'    : "129.128",
           'target'   : self.file_grib
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary)



class gribFile(object):
    """Wrapper for ECMWF grib files with ERA-Interim data. Provides list of
    contents and converts them into ncdf file.

    Args:
        file_grib: File name of a grib file. Derived from ECMWF, dimensions of
                   time, level, lat, lon. Multiple value variables possible.

    Example:
        file_grib = '/Users/stgruber/Desktop/wind.grib'
        file_ncdf = gribFile(file_grib).toNCDF(file_ncdf)
    """
    def __init__(self, file_grib):
        self.file_grib = file_grib
        #get list of contents
        self.jday = []
        self.date = []
        self.levs = []
        self.nams = []
        self.step = []
        grbs = pg.open(self.file_grib)
        grbs.seek(0)
        self.lats = grbs[1].latlons()[0][:,0]
        self.lons = grbs[1].latlons()[1][0,:]
        for grb in grbs:
           self.jday.append(pg.julian_to_datetime(grb.julianDay))
           self.date.append(pg.julian_to_datetime(grb.julianDay) +
                            timedelta(hours=grb.step))
           self.levs.append(grb.level)
           self.nams.append(grb.name)
           self.step.append(grb.step) #usually 0, 3/6/9/12 for accumulated
        grbs.close()
        self.nams = list(set(self.nams))
        self.step = list(set(self.step))
        self.date = list(set(self.date))
        self.jday = list(set(self.jday))
        self.levs = list(set(self.levs))
        self.lats = list(set(self.lats))
        self.lons = list(set(self.lons))
        self.date.sort()
        self.jday.sort()
        self.levs.sort()
        self.lats.sort()
        self.lons.sort()
        self.step.sort() # expect 3/6/9/12 for accumulated
        self.accumulated = ['Total precipitation',
                            'Surface thermal radiation downwards',
                            'Surface solar radiation downwards']
        self.ndate = nc.date2num(self.date, #netCDF date
                                 units = "seconds since 1970-1-1",
                                 calendar='standard') #like UNIX time

    def list(self):
        return {'Dates'       : self.date,
                'Levels'      : self.levs,
                'Variables'   : self.nams,
                'Latitudes'   : self.lats,
                'Longitudes'  : self.lons}


    def toNCDF(self, file_ncdf):
        #NCDF file name
        self.file_ncdf = file_ncdf

        #initialize new data file and create group
        ncd_root = nc.Dataset(self.file_ncdf, 'w', format='NETCDF4_CLASSIC')

        #make dimensions
        ncd_root.createDimension('level', len(self.levs))
        ncd_root.createDimension('time',  len(self.date))
        ncd_root.createDimension('lat',   len(self.lats))
        ncd_root.createDimension('lon',   len(self.lons))

        #make dimension variables
        times      = ncd_root.createVariable('time',    'd',('time',))
        levels     = ncd_root.createVariable('level',  'i4',('level',))
        latitudes  = ncd_root.createVariable('lat',    'f4',('lat',))
        longitudes = ncd_root.createVariable('lon',    'f4',('lon',))

        #assign dimensions
        times[:]      = self.ndate
        levels[:]     = self.levs
        latitudes[:]  = self.lats
        longitudes[:] = self.lons

        #make actual variables
        variables = []
        for var in self.nams:
            # isolate from [u'Geopotential']
            variables.append(ncd_root.createVariable(var,'f4',
                                                        ('time','level',
                                                         'lat','lon',)))
        #read file, get levels and times
        grbindx = pg.index(self.file_grib,'name','level','dataDate',
                                          'dataTime','step')

        levs = np.array(self.levs)
        for l in self.levs:
            for d in self.jday:
                nd = nc.date2num(d, units = "seconds since 1970-1-1",
                                 calendar = 'standard')
                var_n = 0
                for var in self.nams:
                    #distinguish forecast data to deal with accumulated fields
                    if var in self.accumulated:
                        vpre = 0 #initial for subtraction
                        for s in self.step:
                            sel = grbindx.select(name = var, level = l,
                                             dataDate = int(d.strftime("%Y%m%d")),
                                             dataTime = d.hour * 100, step = s)
                            vnow = sel[0].values - vpre
                            vpre = sel[0].values

                            #assign values to netCDF
                            nds = nd + s * 3600 #add step seconds to ncdf time
                            variables[var_n][self.ndate==nds,levs==l,::-1,:] = vnow
                    else:
                        sel = grbindx.select(name = var, level = l,
                                             dataDate = int(d.strftime("%Y%m%d")),
                                             dataTime = d.hour * 100, step = 0)

                        #assign values to netCDF
                        variables[var_n][self.ndate==nd,levs==l,::-1,:] = sel[0].values
                    var_n += 1

        #close ERA-Interim GRIB
        grbindx.close()
        #close netCDF4
        ncd_root.close()

        #return file name
        return self.file_ncdf

    def __str__(self):
        return 'Wrapper object for grib files that can make ncdf.'



class eraData(object):
    """
    Class for the manipulation of ERA-Interim data saved as ncdf.

    Example: pl = plData('/Users/stgruber/Desktop/file.nc')
    """
    def __init__(self):#, file_ncdf, dem_file
        self.g          = 9.80665 #m s-2
        self.absZero    = 273.15

    def describe(self):
        '''
        Generates and prints information on the data contained in ncdf file. The
        information is later on available as a variable.
        '''
        ncf = nc.Dataset(self.file_ncdf, 'r')
        self.lat_min = min(ncf.variables['lat'][:])
        self.lat_max = max(ncf.variables['lat'][:])
        self.lon_min = min(ncf.variables['lon'][:])
        self.lon_max = max(ncf.variables['lon'][:])
        self.lev_min = min(ncf.variables['level'][:])
        self.lev_max = max(ncf.variables['level'][:])
        self.time_min = min(ncf.variables['time'][:])
        self.time_max = max(ncf.variables['time'][:])
        self.time_min = nc.num2date(self.time_min,
                               units = "seconds since 1970-1-1",
                               calendar='standard')
        self.time_max = nc.num2date(self.time_max,
                               units = "seconds since 1970-1-1",
                               calendar='standard')

        print "Time:      "    + str(self.time_min) + " to " + str(self.time_max)
        print "Latitude:  "   + str(self.lat_min) + " to " + str(self.lat_max)
        print "Longitude: " + str(self.lon_min) + " to " + str(self.lon_max)
        print "Level:     "     + str(self.lev_min) + " to " + str(self.lev_max)
        ncf.close()


    def NCDFmerge(self, file_list, file_new):
        '''Merge multiple netCDF files with identical structure but differing
        times together into one netCDF file'''
        #TODO: sort in time dimension
        #merge netCDF files
        ncl = nc.MFDataset(file_list, aggdim='time')

        #initialize new data file and create group
        ncn = nc.Dataset(file_new, 'w', format='NETCDF4_CLASSIC')

        #make dimensions
        ncn.createDimension('level', len(ncl.variables['level'][:]))
        ncn.createDimension('time',  len(ncl.variables['time'][:]))
        ncn.createDimension('lat',   len(ncl.variables['lat'][:]))
        ncn.createDimension('lon',   len(ncl.variables['lon'][:]))

        #make dimension variables
        times      = ncn.createVariable('time',    'd',('time',))
        levels     = ncn.createVariable('level',  'i4',('level',))
        latitudes  = ncn.createVariable('lat',    'f4',('lat',))
        longitudes = ncn.createVariable('lon',    'f4',('lon',))

        #assign dimensions
        times[:]      = ncl.variables['time'][:]
        levels[:]     = ncl.variables['level'][:]
        latitudes[:]  = ncl.variables['lat'][:]
        longitudes[:] = ncl.variables['lon'][:]

        #create and assign variables
        for var in ncl.variables:
            if not var in ['time','level','lat','lon']:
                nowvar =  ncn.createVariable(var,'f4',('time','level',
                                                            'lat','lon',))
                nowvar[:,:,:,:] = ncl.variables[var][:,:,:,:]

        #cleanup
        ncl.close()
        ncn.close()

    def split_seq(self, seq, size):
        '''Split a list into chunks of defined size'''
        newseq = []
        splitsize = 1.0/size*len(seq)
        for i in range(size):
                newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
        return newseq

    def DateFile(self, filename, get='beg'):
        if get == 'beg': res = filename[-19:-13]
        if get == 'end': res = filename[ -9: -3]
        return res

    def NCDFmergeWildcard(self, files, n_to_combine):
        '''Merge multiple netCDF files with identical structure but differing
        times together into one netCDF file'''
        #get directory list and split
        all_list = self.split_seq(sorted(gl.glob(files)), n_to_combine)

        for file_list in all_list:
            sbeg = self.DateFile(file_list[ 0], get='beg')
            send = self.DateFile(file_list[-1], get='end')
            file_new = file_list[0][:-19] + 'm' + '_' + sbeg + '_' + send + '.nc'
            print file_list
            print file_new
            self.NCDFmerge(file_list, file_new)


class toposcale(object):
    """Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server, for converting grib to ncdf.

    Args:
        date: A dictionary specifying the time period desired with a begin
              and an end date given as a datetime.datetime object.

        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east.

        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels
                   needed. Unit: [m].

        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].

        directory: Directory to hold output files

    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(1999, 1, 2)}
        area  = {'north' :  40.0,
                 'south' :  41.0,
                 'west'  :  60.0,
                 'east'  :  61.0}
        elevation = {'min' :   50,
                     'max' : 2000}
        directory = '/Users/stgruber/Desktop/aaa'
        ts = toposcale(date, area, elevation, directory, 5)
        ts.retrieve()
        ts.deleteGrib()
    """
    def __init__(self, date, area, elevation, directory, increment_days):
        self.date      = date
        self.area      = area
        self.elevation = elevation
        self.directory = directory
        self.increment = increment_days
        self.nc_files  = ''
        #TODO ensure increments is smaller or equal than chosen time window

    def getFileNames(self):
        return self.nc_files
        #TODO: add methods to get file names with * in them

    def retrieve(self):
        #define variables
        """
        var_pl = ['airt', 'rh', 'wind', 'geop']
        var_sa = ['airt2', 'dewp2', 'wind10']
        var_sf = ['prec', 'swin', 'lwin']
        """
        var_pl = ['airt', 'geop']
        var_sa = ['airt2']
        var_sf = ['prec']
        #enter time loop
        date_i = {'beg' : datetime(1994, 1, 1), 'end' : datetime(1999, 1, 2)}
        slices = floor(float((self.date['end'] - self.date['beg']).days)/
                       self.increment)+1

        for ind in range (0, int(slices)):
            #prepare time slices
            date_i['beg'] = self.date['beg'] + timedelta(days =
                                                         self.increment * ind)
            date_i['end'] = self.date['beg'] + timedelta(days =
                                                         self.increment * (ind+1) - 1)
            if ind == (slices -1):
                date_i['end'] = self.date['end']

            #actual functions
            pl = ERApl(date_i, self.area, self.elevation,
                       var_pl, self.directory)
            sa = ERAsa(date_i, self.area, var_sa, self.directory)
            self.ERAli   = [pl, sa] #combine in list

            #download from ECMWF server convert to netCDF
            for era in self.ERAli:
                era.download()
                era.toNCDF()

        #topography
        top = ERAto(self.area, self.directory)
        top.download()
        top.toNCDF()

    def __str__(self):
        return "Object for toposcale data download and conversion"
