# -*- coding: utf-8 -*-
# 
# REanalysis Downscaling Cold Air Pooling Parameterization (REDCAPP)
#
# === COPYRIGHT AND LICENCE ====================================================
#
# Copyright 2013-2017 Stephan Gruber 
#           2015-2017 Bin Cao
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# === CONTRIBUTIONS ============================================================
#
# Code for ERA-Interim interaction written by Stephan Gruber, with additions and
# testing by Bin Cao. Terrain analysis, interpolation, and LSCF prediction
# written by Bin Cao with input from Stephan Gruber. 
#
# === NOTES ====================================================================
#
# Get REDCAPP at https://github.com/geocryology/REDCAPP
# Check ECMWF job status: http://apps.ecmwf.int/webmars/joblist/
#
# ==============================================================================

import numpy as np
from ecmwfapi import ECMWFDataServer
import netCDF4 as nc
import pygrib  as pg
import csv
import re

from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter,generic_filter,convolve,minimum_filter,maximum_filter
from math import radians, exp, floor
from bisect import bisect_left
from datetime import datetime, timedelta
from os import path, remove

import glob as gl



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

        print "Time:      " + str(self.time_min) + " to " + str(self.time_max)
        print "Latitude:  " + str(self.lat_min) + " to " + str(self.lat_max)
        print "Longitude: " + str(self.lon_min) + " to " + str(self.lon_max)
        print "Level:     " + str(self.lev_min) + " to " + str(self.lev_max)
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


class redcapp_get(object):
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
        rg = redcapp_get(date, area, elevation, directory, 5)
        rg.retrieve()
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
        var_pl = ['airt', 'geop']
        var_sa = ['airt2']
        
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
        return "Object for data download and conversion"




class rawData(object):
    """
    Args:
        dir_data: directory containing all raw data and output data
        
    Example:
        dir_data = 'C:/users/bincao/Desktop/data'
        dataImport = rawData(dir_data)
        sa = dataImport.saf_get()#get sa file in the given directory
        pl = dataImport.plf_get()#get pl file in the given directory
    """
    
    def __init__(self, dir_data):
        self.dir = dir_data
        
    def file_get(self, nomenclature):
        """file check"""
        
        flist = gl.glob(path.join(self.dir, nomenclature))
        if len(flist) >= 1:
            return flist[0]
        else:
            print("File " + nomenclature + " not found in directory " + 
                  self.dir + ".")
    
    def plf_get(self, nomenclature = 'era_pl_*_to_*'):
        """finds the pressure level files in the directory 
        based on 'era_pl_*_to_*'
        """
        return self.file_get(nomenclature)
            
    def saf_get(self, nomenclature = 'era_sa_*_to_*'):
        """finds the pressure level files in the directory 
        based on 'era_sa_*_to_*'
        """
        return self.file_get(nomenclature)

            
    def geopf_get(self, nomenclature = 'era_to*'):
        """finds the era-interim geopotential file in the directory 
        based on 'era_to*'
        """
        return self.file_get(nomenclature)

            
    def asciiDemHeader(self, demAsiccf):
        """Returns header information of input DEM in ASCIIGRID format"""
        header = []
        with open(demAsiccf) as f:
            reader = csv.reader(f)
            i=0
            for row in reader:
                i = i+1
                if (i <= 5):
                    header.append(row[0])
            return header
            
    def asciiDemEle(self, demAsciif):
        """
        
        Args:
           demAsciif: DEM in ASCIIGRID format and lat/lon WGS84 grid
           
         Returns:
             headerInfo: header information of input ascii file
             ele: array-like elevation of input ascii file
        """

        #read elevation
        return np.loadtxt(demAsciif, delimiter=' ', skiprows = 5)
        
    def ascii2ncdf(self, demAsciif, dem_out):
        """convert input DEM in ASCIIGRID to netcdf file
        
        Args:
            demAsciif: DEM in ASCIIGRID format and lat/lon WGS84 grid. 
                the header of ascii file should be same as the example data,
                xllcorner is the lontitude in the left low point, while
                yllcorner is the latitude in the left low point.
            dem_out: output DEM in netcdf format
            
        Returns a ASCIIGRID-based DEM in netcdf
        
        Example:
            dir_data = 'C:/users/bincao/Desktop/data'
            demAsciif = 'DEM_testArea.asc'
            dem_out  = 'C:/users/bincao/Desktop/data/DEM_fine-scale.nc'
            
            dataImport =  dataImport = rawData(dir_data)
            dataImport.ascii2ncdf(dem_file, dem_out)
            
        """
        
        #meta information
        header = self.asciiDemHeader(demAsciif)
        ncol   = int(re.findall(r'\d+', header[0])[0])
        nrow   = int(re.findall(r'\d+', header[1])[0])
        xllcorner = float(re.findall(r'\d+\.\d+', header[2])[0])
        yllcorner = float(re.findall(r'\d+\.\d+', header[3])[0])
        cellsize  = float(re.findall(r'\d+\.\d+', header[4])[0])
        
        #get variables
        ele = self.asciiDemEle(demAsciif)# elevation
        lats = np.linspace(yllcorner, yllcorner+cellsize*nrow, 
                           nrow, endpoint= True)# latitude
        lons = np.linspace(xllcorner, xllcorner+cellsize*ncol, 
                           ncol, endpoint= True)# lontitude
        
        #create nc file
        nc_root = nc.Dataset(dem_out ,'w', format = 'NETCDF4_CLASSIC')
        
        #create dimensions
        nc_root.createDimension('lat', nrow)
        nc_root.createDimension('lon', ncol)
        
        #create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes  = nc_root.createVariable('lat', 'f4', ('lat'))
        elevation  = nc_root.createVariable('elevation', 
                                           'f4', ('lat', 'lon'), zlib = True)
        
        #assign variables
        longitudes[:] = lons
        latitudes[:]  = lats[::-1]
        elevation[:]  = ele
        
        #attribute
        nc_root.description = "high-resolution topography file"
        #resolution = cellsize
        longitudes.units = 'degree_east (decimal)'
        latitudes.units  = 'degree_north (decimal)'
        elevation.units  = 'm'
        
        nc_root.close()


class downscaling(object):
    """
    Return object for downscaling that has methods for interpolationg
    upper-air temperature and surface influences at surface level
    based on disaggregating coarse-grid reanalysis and dem.
    
    Args:
        dem: A required fine-scale dem in netcdf format
    
    Example:
        dem  = 'example_alps.nc'
        geop = 'alps_geop.nc'
        sa   = 'alps_sa_79_15
        pl   = 'alps_pl_79_15.nc'
        
        downscaling = downscaling(dem, geop, sa, pl)
        
    """
    
    def __init__(self, geop, sa, pl, dem = None):
        self.g    = 9.80665 #Gravitational acceleration [m/s2]
        self.geop = nc.Dataset(geop)
        self.sa   = nc.Dataset(sa)
        self.pl   = nc.Dataset(pl)
        if not (dem is None):
            self.dem  = nc.Dataset(dem)
        
        
    def demGrid(self, stations = None):     
        """Return metadata of given stations or dem. 
        Format of stations desired [lat, lon, geop].
        
        Args: 
            stations: A list of dictionaries describing stations. If not given,
                metadata derived from given dem.
            
        Returns:
            out_xyz_dem: Metadata [lat, lon, geop] of input dem or sites
            lons: Longitude of input sites
            lats: Latitude of input sites
            shape: Shape of input dem or sites
            names: Name of input stations. Names will only be avaiable when 
                stations are inputted
        """
        
        if not (stations is None):
            lats = [s['lat'] for s in stations]
            lons = [s['lon'] for s in stations]
            names = [s['name'] for s in stations]
            siteLocation = np.asarray([[s['lat'],s['lon'],
                                        s['ele']*self.g] for s in stations])
            shape = siteLocation.shape
            return siteLocation, lats, lons, shape, names        

        #out_xyz based on dem
        lons = self.dem.variables['lon'][:]
        lats = self.dem.variables['lat'][:]
        geop = self.dem.variables['elevation'][:]*self.g
        shape = geop.shape
        
        
        lons, lats = np.meshgrid(lons, lats)
        lons = lons.reshape(lons.size)
        lats = lats.reshape(lats.size)
        geop = geop.reshape(geop.size)
        
        out_xyz_dem = np.array([lats, lons, geop]).T

        return out_xyz_dem, lats, lons, shape
        
        
    def geoGrid(self):     
        """
        Returns original (NO interpolation) coarse gird metadata
        in format of [lat, lon, geop], based on geopotential file
        """
        longitude = self.geop['lon'][:]
        latitude  = self.geop['lat'][:]
        lons, lats = np.meshgrid(longitude, latitude)
        lons = lons.reshape(lons.size)
        lats = lats.reshape(lats.size)
        geop = self.geop['Geopotential'][0,:,:]#geopotential
        geop = geop.reshape(geop.size)

        out_xyz_ori = np.array([lats, lons, geop]).T
        
        return out_xyz_ori
        
    
    def surGrid(self, lats, lons, stations):
        
        """
        Return interpolated surface geopotetial.
        
        Args: 
            lats: Latitude of intersted sites
            lons: Longitude of intersted sites
            stations: A list of dictionaries describing stations. If not given,
                metadata derived from given dem.
        
        Returns:
            out_xyz_sur: Fine-scale surface level metadata [lat, lon, geop],
            where latitude and lontitude are ontained from dem,
            while geop is interpolated from coarse geopotential file
            
        Example:
            downscaling = downscaling(dem, geop, sa, pl)
            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_ori = downscaling.geoGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, out_xyz_dem[:,:2])
        """
        
        longitude = self.geop['lon'][:]
        latitude  = self.geop['lat'][:]
        in_v = self.geop['Geopotential'][0,0,:,:] #geopotential
        fz = RegularGridInterpolator((latitude,longitude), in_v, 'linear')
        out_xy = np.array([lats, lons]).T

        if not (stations is None):
            lats = [s['lat'] for s in stations]
            lons = [s['lon'] for s in stations]
            out_xy = np.asarray([[s['lat'],s['lon']] for s in stations])
        
        z_interp = fz(out_xy)
        out_xyz_sur = np.array([lats, lons, z_interp]).T
        
        return out_xyz_sur

 

    def surTa(self, ind_time, out_xyz_sur):
        """Return interpolated 2-metre temperature.
        
        Args:
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            out_xyz_sur: 
        
        Returns:
            t_sa: interpolation fine-scale surface air temperatue based on 
            grid 2-meter temperature
            
        Example:
            dem  = 'example_alps.nc'
            geop = 'alps_geop.nc'
            sa   = 'alps_sa_79_15.nc'
            pl   = 'alps_pl_79_15.nc'
        
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)


            surTa = downscaling.surTa(0, out_xyz_sur)
        """
            
        in_v = self.sa['2 metre temperature'][ind_time,0,:,:]#geopotential
        in_v -= 273.15
        lat = self.sa.variables['lat'][:]
        lon = self.sa.variables['lon'][:]

        f_sa = RegularGridInterpolator((lat,lon), in_v, 'linear')
        t_sa = f_sa(out_xyz_sur[:,:2]) 

        return t_sa


    def gridValue(self, variable, ind_time):
        """
        Return original grid temperatures and geopotential of differnet
        pressure levels. The function are called by inLevelInterp() to
        get the input ERA-Interim values.
        
        Args: 
            variable: Given interpolated climate variable
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            
        Returns:
            gridT: Grid temperatures of different pressure levels. Retruned 
            temperature are formated in [level, lat, lon]
            gridZ: Grid geopotential of different pressure levels. Retruned 
            temperature are formated in [level, lat, lon]
            gridLon: Grid longitude of pressure level variables
            gridLat: Grid latitude of pressure level variables
        
        Example:
            gridT,gridZ,gridLat,gridLon=downscaling.gridValue('Temperature',0)
            
        """
        
        gridT = self.pl.variables[variable][ind_time,:,:,:]
        gridZ = self.pl.variables['Geopotential'][ind_time,:,:,:]
        #x and y
        
        gridLat = self.pl['lat'][:]
        gridLon = self.pl['lon'][:]
        

        return gridT,gridZ,gridLat,gridLon


    def inLevelInterp(self,gridT, gridZ, gridLat, gridLon, out_xyz):
        """
        This is a 2D interpolatation, and returns interpolated temperatures
        of different pressure levels.
        
        Args:
            gridT: Grid temperatures of different pressure levels. Retruned 
                temperature are formated in [level, lat, lon]
            gridZ: Grid geopotential of different pressure levels. Retruned 
                temperature are formated in [level, lat, lon]
            gridLat: Grid longitude of pressure level variables
            gridLon: Grid latitude of pressure level variables
            out_xyz: Given sites, which will be interpolated.
            
        Returns:
            t_interp: Interpolated temperatre of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            z_interp: Interpolated geopotential of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
        
        Examples:
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)

            #interpolate 2-meter temperature
            surTa = downscaling.surTa(0, out_xyz_sur)
            #original ERA-I values
            gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
            #interpolate temperatures and geopotential of different 
            pressure levels.

            t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
                                                           gridLat,gridLon,
                                                           out_xyz_dem)
        """
        
        shape = gridT.shape
        #create array to hold interpolation resultes
        t_interp = np.zeros([shape[0], len(out_xyz)])
        z_interp = np.zeros([shape[0], len(out_xyz)])

        #temperatue and elevation interpolation 2d
        for i in range(shape[0]):
            ft = RegularGridInterpolator((gridLat,gridLon), 
                                          gridT[i,:,:], 'linear')
            fz = RegularGridInterpolator((gridLat,gridLon), 
                                          gridZ[i,:,:], 'linear')
            t_interp[i,:] = ft(out_xyz[:,:2])#temperature
            z_interp[i,:] = fz(out_xyz[:,:2])#elevation

        t_interp -= 273.15

        return t_interp[::-1,:], z_interp[::-1,:]
        
    

        
    def fast1d(self, t_interp, z_interp, out_xyz):
        """This is a 1D interpoation. The function return interpolated 
        upper air temperature at the given sites by 
        interpolation between different pressure levels.
        
        Args:
            t_interp: Interpolated temperatre of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            z_interp: Interpolated geopotential of different pressure levels. 
                The returned values are fomrated in [level, lat, lon]
            out_xyz: Given sites with elevation, which will be interpolated.
            
        Returns:
            dG:upper-air temperature at given sites
                
        Example: 
            downscaling = downscaling(dem, geop, sa, pl)

            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)

            
            surTa = downscaling.surTa(0, out_xyz_sur)
            #original ERA-I values
            gridT,gridZ,gridLat,gridLon = downscaling.gridValue(variable,0)
            #interpolate temperatures and geopotential of different 
            pressure levels.

            t_interp, z_interp = downscaling.inLevelInterp(gridT,gridZ,
                                                           gridLat,gridLon,
                                                           out_xyz_dem)
            
            #upper air temperature at the coarse and fine scale of elevation
            pl_sa = fast1d(t_interp, z_interp, out_xyz_sur)
            pl_obs = fast1d(t_interp, z_interp, out_xyz_dem)
        """
    
        ele = out_xyz[:,2]
        size = np.arange(out_xyz.shape[0])
        n = [bisect_left(z_interp[:,i], ele[i]) for i in size]
        n = [x+1 if x == 0 else x for x in n]
        
        lowN = [l-1 for l in n]
        
        upperT = t_interp[n,size]
        upperZ = z_interp[n,size]
        dG  = upperT-t_interp[lowN,size]#<0
        dG /= upperZ-z_interp[lowN,size]#<0
        dG *= out_xyz[:,2] - upperZ#>0
        dG += upperT
             
        return dG
    
          
        
    def interpAll(self, variable, ind_time, out_xyz_sur, out_xyz_obs):
        """Returns all needed interpolated temperatures at given time.
        
        Args:
            variable: Interpolated climated variable
            ind_time: Time need to be interpolated. Time is in interger (e.g.
            0, 1, 2)
            out_xyz_sur: Interploated sites[lat, lon, geop], in which the 
            geopotential are interpolated from ERA-Interim geopotential file.
            out_xyz_obs: Interploated sites[lat, lon, geop], in which the 
            geopotential are gotten from DEM.
            
        Returns:
            pl_obs: Interpolated upper-air temperature at given sites
                with dem level
            pl_sur: Interpolated upper-air temperature at given sites
                with surface level
            
            
        Example:
            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)
            pl_obs,pl_sur,t_sa = downscaling.interpAll(variable, ind_time, 
                                                       out_xyz_sur, 
                                                       out_xyz_dem)
            
            
        """
        gridT,gridZ,gridLat,gridLon = self.gridValue(variable, ind_time)
        t_interp,z_interp = self.inLevelInterp(gridT,gridZ,gridLat,gridLon,
                                               out_xyz_obs)
        pl_sur = self.fast1d(t_interp, z_interp, out_xyz_sur)
        pl_obs = self.fast1d(t_interp, z_interp, out_xyz_obs)
        t_sa = self.surTa(ind_time, out_xyz_sur)
        dt = t_sa-pl_sur
        
        return pl_obs, dt

    def spatialMean(self, variable, daterange):
        """Return the MEAN upper-air temperature and 
        land surface influence during given date range and at given  area
        
        Args:
            variable: Climated variable to interpolate
            daterange: Date range to interpolate
            out_xyz_sur: Surface level sites to interpolate 
            out_xyz_obs: Topography sites to interpolate
                
        Returns:
            pl: Interpolated MEAN free-atmosphere during given date range and
            at given area.
            dt: Interpolated MEAN surface level land surface influences during 
            given date range and at given area.
            
        Example:
            
            dem  = 'example_alps.nc'
            geop = 'alps_geop.nc'
            sa   = 'alps_sa_79_15.nc'
        
            downscaling = downscaling(dem, geop, sa, pl)
            
            out_xyz_dem, lats, lons, shape = downscaling.demGrid()
            out_xyz_sur = downscaling.surGrid(lats, lons, None)
            pl,dt = downscaling.spatialMean(variable, daterange, 
                                            out_xyz_sur, out_xyz_dem, shape)
            
        """
        
        #obtain time range
        date_vec = nc.num2date(self.pl.variables['time'][:],
                                units = "seconds since 1970-1-1",
                                calendar='standard')


        #index of time steps to interpolate
        mask  = date_vec >= daterange.get('beg')
        mask *= date_vec <= daterange.get('end')
        ind_time_vec = np.arange(len(mask))[mask]
        out_time = date_vec[mask]

        sum_pl_obs = 0
        sum_dt = 0
        
        out_xyz_dem, lats, lons, shape = self.demGrid()
        out_xyz_sur = self.surGrid(lats, lons, None)

        print("\nConducting downscaling now, have a cup of coffee please\n")

        for ind_out, ind_time in enumerate(ind_time_vec):
            print(out_time[ind_out])
            pl_obs,dt = self.interpAll(variable,ind_time, out_xyz_sur, out_xyz_dem)
            sum_pl_obs += pl_obs
            sum_dt += dt

        dt = sum_dt/out_time.size
        pl = sum_pl_obs/out_time.size
        
        dt = dt.reshape(shape)
        pl = pl.reshape(shape)

        return pl,dt
       
    def stationTimeSeries(self, variable, daterange, stations):
        """Return upper-air temperature and land surface influence
        at given time steps
        
        Args:
            variable: Climated variable to interpolate
            daterange: Date range to interpolate
            out_xyz_sur: Surface level sites to interpolate 
            out_xyz_obs: Topography sites to interpolate
            
        Returns:
            out_valu: interpolated upper-air temperature & land-surface effects
            out_time: time series
            
        Example:
            Downscaling = downscaling(geop, sa, pl, dem_ncdf)
            stations=[{'name':'COV','lat': 46.41801198, 'lon': 9.821232448, 'ele': 3350.5},
                      {'name':'SAM','lat': 46.52639523, 'lon': 9.878944266, 'ele': 1756.2}]
            out_xyz_dem, lats, lons, shape, names = Downscaling.demGrid(stations)
            out_xyz_sur = Downscaling.surGrid(lats, lons, site)
            
            pl, dt, out_time, names = Downscaling.stationTimeSeries(variable, 
                                                                    daterange, 
                                                                    out_xyz_sur,
                                                                    out_xyz_dem)
        """
        
        #obtain time range
        date_vec = nc.num2date(self.pl.variables['time'][:],
                                units = "seconds since 1970-1-1",
                                calendar='standard')
        
        #obtain station information
        out_xyz_dem, lats, lons, shape, names = self.demGrid(stations)
        out_xyz_sur = self.surGrid(lats, lons, stations)

        #index of time steps to interpolate
        mask  = date_vec >= daterange.get('beg')
        mask *= date_vec <= daterange.get('end')
        ind_time_vec = np.arange(len(mask))[mask]
        out_time = date_vec[mask]

        out_valu = np.zeros((2, ind_time_vec.size, 
                             out_xyz_dem.shape[0]))#pl_obs, dt

        print("\nConducting downscaling now, have a cup of coffee please\n")

        for ind_out, ind_time in enumerate(ind_time_vec):
            print(out_time[ind_out])
            out_valu[:,ind_out,:] = self.interpAll(variable, ind_time,
                                                   out_xyz_sur, out_xyz_dem)

        return out_valu[0],out_valu[1], out_time, names

    
    def extractStationAirTCSV(self, daterange, variable, stations, 
                              stat_out1, stat_out2):
        """
        Extracts time series from gridded data based on a list of dictionaries
        describing stations. Time serie(s) are written into csv files,
        variables are rounded to a precision of 3 decimal places.
        
        Args:
            daterange: Date range to interpolate. Format of time range 
                desired (yy, mon, day, hh, min)
            variable: Climate variable to interpolate
            stations: Sites to interpolate. Format of stations desired
                ('name':'siteName','lat':latNumber, 'lon':lonNumber, 'ele':eleNumber)
            stat_out1: station timeseries air temperature output file name
            stat_out1: station timeseries coarse land-surface effects output
                file name
            
        Returns: 
            Returns a csv file contains time series downscaled surface air 
            temperature and land surface influences (surface level) at given 
            sites.

        Example:
            
            daterange = {'beg' : datetime(2000, 01, 11, 00, 00),
                         'end' : datetime(2000, 01, 11, 06, 00)}
                         
            vairable = 'Temperature'
            
            stations = [{'name':'TAE', 'lat':47.47986561, 'lon':8.904870734, 'ele':85.80},
                        {'name':'AAR', 'lat':47.38799123, 'lon':8.043881188, 'ele':454.87}]
            
            file_out = ['/Users/bincao/OneDrive/pl_obs.csv',
                        '/Users/bincao/OneDrive/dT.csv']
            
            downscaling.extractStationDataCSV(daterange,variable,stations,file_out) 
        """
        
        #interpolate and extract values
        
        values, time, names = self.stationTimeSeries(variable, daterange, stations)

        #write CSV
        file_out = [stat_out1, stat_out2]
        names.insert(0, 'Time_UTC')
        for i in range(2):
              with open(file_out[i], 'wb') as output_file:
                 writer = csv.writer(output_file)
                 writer.writerow(names)
                 valu = values[i].tolist()
                 for n in range(len(time)):
                     row = ['%.3f' % elem for elem in valu[n]]
                     row.insert(0,time[n])
                     writer.writerow(row)
    
                     
    
    def extractSpatialAirTNCF(self, daterange, variable, file_out):
        """
        Extracts mean of given date range from gridded data based on a 
        fine-scale DEM. Mean values are written into a netcdf file,
        variables are rounded to a precision of 3 decimal places.
        
        Args:
            daterange: Date range to interpolate. Format of time range 
                       desired (yy, mon, day, hh, min)
            variable: Climate variable to interpolate
            filout: Name of output netcdf file
        
        Returns:
            A netcdf file contains coordinate information, downscaled surface
            air temperature and land surface influence (surface level).
            

        Example:
            
            daterange = {'beg' : datetime(2000, 01, 11, 00, 00),
                         'end' : datetime(2000, 01, 11, 06, 00)}
                         
            vairable = 'Temperature'
            
            stations = [{'name':'TAE', 'lat':47.47986561, 'lon':8.904870734, 'ele':85.80},
                        {'name':'AAR', 'lat':47.38799123, 'lon':8.043881188, 'ele':454.87}]
            
            file_out  = '/Users/bincao/Desktop/subgrid.nc'
            
            downscaling.extractSpatialDataNCF(daterange, variable, file_out) 
        """
        
        #temperature and coarse land-surface effects
        pl,dt = self.spatialMean(variable, daterange)
        shape = self.dem.variables['elevation'].shape
                                  
        #create nc file
        nc_root = nc.Dataset(file_out ,'w', format = 'NETCDF4_CLASSIC')
        
        #create dimensions
        nc_root.createDimension('lat', shape[0])
        nc_root.createDimension('lon', shape[1])
        
        #create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes = nc_root.createVariable('lat', 'f4', ('lat'))
        Ta = nc_root.createVariable('surface air temperature', 
                                     'f4', ('lat', 'lon'), zlib = True)
        dT = nc_root.createVariable('coarse scale of land surface influence', 
                                     'f4', ('lat', 'lon'), zlib = True)
        
        #assign variables
        lons = self.dem.variables['lon'][:]
        lats = self.dem.variables['lat'][:]
        longitudes[:] = lons
        latitudes[:] = lats
        Ta[:] = pl
        dT[:] = dt
        
        #attribute
        nc_root.description = "Downscaled upper-air temperature and "\
                              "coarse-scale land surface influence with a"\
                              "spatial resolution of input dem"
        longitudes.units = 'degree_east (decimal)'
        latitudes.units  = 'degree_north (decimal)'
        Ta.units = 'celsius'
        dT.units = 'celsius'
        
        nc_root.close()



class topography(object):
    """
    Return object for topography that has methods for deriving topographic 
    factors based on DEM.
    
    (a) MRVBF:
        Multiresolution index of valley bottom flatness (MRVBF) is derived by 
        Gallant and Dowling (2003), detailed could be found:
        http://onlinelibrary.wiley.com/doi/10.1029/2002WR001426/abstract
        The first slope threshold is set to 50% so that the mrvbf is smoother 
        than original one and could represent the cold air pooling better.
    (b) Hypsometric Position:
        Hypsometric position is calculated as the radio of the number of cells
        with higher elevation than given site to the total number cells in the 
        surrounding region and ranges from 1 (deepest valley) to 0 
        (highest peak).
        
    Args:
         demFile: input dem file
         demResoultion: the resolution of input dem in degree
         
    Example:
        dem  = 'example_alps.nc'
        demResolution = 3./3600
        topo = topography(dem, demResoultion)
        topo.spatialTopo('/Users/bincao/Desktop/topo.nc')
    """
    

    def __init__(self, demFile, demResolution):
        self.dem = nc.Dataset(demFile, 'r')
        self.R = 6371000 # the mean radius (in meter) of Earth
        self.resolution = demResolution #units = degree
        self.lon = self.dem['lon'][:]
        self.lat = self.dem['lat'][:]
        self.shape = [len(self.lat), len(self.lon)]
        self.size = len(self.lon)*len(self.lat)
    
    def describe(self):
        """Returns summary information of DEM file
        Example:
            topo.describe()
        """
        minlon = np.min(self.lon)
        maxlon = np.max(self.lon)
        minlat = np.min(self.lat)
        maxlat = np.max(self.lat)

        print 'Ranges:  West: %s,  East: %s, South: %s, North: %s '  %(minlon, maxlon, minlat, maxlat)
        print 'Resolution: ' + str(self.resolution)
        print 'Dem shape [Lat, Lon]: ' + str(self.shape)
        print 'Dem size:' + str(self.size)

    def pixelLength(self, lat, L=1):
        """Return the cellsize for L step in meter. The function
        Args:
            lat: latitude in list
            L: interger, step number
            
        Returns:
            yL: length of cellsize in lontitude direction
            [x*cellsize for x in xL]: length of cellsize in latitude direction
            
        Example:
            yi,xi = topo.pixelLength(topo.lat, L=3)
        """

        y1 = self.resolution *np.pi*self.R/180#lon
        x1 = [np.cos(radians(elem)) *y1 for elem in lat]#lat at base scale
        if L <= 2:
            cellsize = 1
        else:
            cellsize = 3**(L-2)
        yL = y1 * cellsize
        xL = x1[(cellsize-1)/2:len(lat):cellsize]

        return [yL, [x*cellsize for x in xL]]#x1[np.arange((cellsize-1)/2,len(lat),cellsize)]
    
    def demSize(self):
        """Returns dem size in km"""
        
        yi = self.pixelLength(self.lat)[0] # cell size in m
        # dem size [lat, lon] in km
        size = [yi * dim/1000 for dim in self.shape]
        
        return size
    
    def sizeCheck(self, threshold = 0.8):
        """Check weather input dem is big enough to simulate 
        hypsometric position."""
        
        size = self.demSize()
        # valid percent (without edge effecet)
        validPer = [(size[0]-30)/size[0], (size[1]-30)/size[1]]
        
        if any([validPer[0] < threshold, validPer < threshold]):
            print ('Input DEM (width in ' + 
                               str(int(size[1])) + ' km, breadth in ' + 
                               str(int(size[0])) + ' km) ' + 
                'is too small to meaningfully ' +
                'accommodate the 30 km neighbourhood. Only one value is used ' +
                'as elevation range for the entire DEM')
            
            return True
        
        else:
            return False
        
    
    def scale(self, x, t, p):
        """
        Scale the input value onto [0,1]

        Args:
            x: Input value (> 0)
            t: Threshold parameter
            p: Shape parameter, larger values give more abrupt transitions
            between 1 (for x << t) and 0 (for x >> t).
            
        Returns:
            scaled value
         """

        return 1/(1 + (x/t)**p)


    def smoothDEM(self, ele):
        """
        Returns a smoothed DEM in 2D Gaussian kernel array format.
        The smoothing is performed with the Arc/Info focal mean function
        using an 11 × 11 Gaussian smoothing kernel with an radius of three
        cell (to corresopnd with the factor of 3 resolution change).
        Detailed introduciton could be found from equation (11) in MRVBF
        paper.
        
        Args:
            ele: Input fine-scale of elevation
            
        Returns:
            Smoothed DEM
            
        Example:
            sdemL = topo.smoothDEM(topo.dem['elevation'][:])
        """
        
        return gaussian_filter(ele,np.sqrt(4.5), truncate = 11)
        
        

    def refine(self, L, coarseValue, out_xy = None, limitSize = 80000000):
        """
        Return refined coarseValue by linear interpolation.
        
        Args:
            L: Interger, step number
            coarseValue: Values of coarse scale need be refined
            out_xy: Sites need to be refined
            limitSize: Maximum cell size input once. This could avoid RAM crash
            
            
        Returns:
            Refined value
            
        Example:
            
        """
        
        scale = 3**(L-2)
        latIndex = range((scale-1)/2,len(self.lat),scale)
        lonIndex = range((scale-1)/2,len(self.lon),scale)
        f = RegularGridInterpolator((self.lat[latIndex][::-1],self.lon[lonIndex]),
                coarseValue[::-1,:], method = 'linear',bounds_error = False)

        if not (out_xy is None):
            return f(out_xy)

        elif self.size <= limitSize:
            lon, lat = np.meshgrid(self.lon, self.lat)
            gridBase = np.array([lat.reshape(lat.size), lon.reshape(lon.size)]).T
            fineValue = f(gridBase)
            return fineValue.reshape((len(self.lat),len(self.lon)))
        else:
            fineValue = np.zeros((len(self.lat),len(self.lon)))
            chunkN = self.size/limitSize
            latLegth = len(self.lat)/chunkN

            for pieceN in range(0, chunkN-1):
                startLat = latLegth*pieceN
                endLat = latLegth*(pieceN+1)
                lon, lat = np.meshgrid(self.lon, self.lat[startLat:endLat])
                gridBase = np.array([lat.reshape(lat.size), 
                                     lon.reshape(lon.size)]).T
                fineValue[startLat:endLat,:]=f(gridBase).reshape(latLegth,len(self.lon))

            lon, lat = np.meshgrid(self.lon, self.lat[latLegth*(chunkN-1):])
            gridBase = np.array([lat.reshape(lat.size), 
                                 lon.reshape(lon.size)]).T
            fineValue[latLegth*(chunkN-1):,:] = f(gridBase).reshape(len(self.lat[latLegth*(chunkN-1):]),len(self.lon))

            return fineValue


    def flatness(self, ele, Tf, out_xy = None,  L = 1, Pf = 4):
        """
        Return flatness for step L, flatness is defined as the inversion of 
        slope and called by scale() function.
        Slope is modeled by standard finite difference techniques, and 
        described as a percentage or 100 times the tangent of slope angle.
        
        The detailed describtion of calculating slope could be found:
            http://help.arcgis.com/En/Arcgisdesktop/10.0/Help/index.html#/How_Slope_works/009z000000vz000000/

        Args:
            ele: array_like, input elevation
            Tf: threshold, transform slope to flatness by scale() function
            out_xy: array_like, the interested sites. If not given (None), 
                flatnessof all the DEM cells would be simulated. Default is 
                None.
            L: interger, step number
            Pf: Shape parameter used to scale the value to range of 0-1.

        Results:
            Return flatness value at L step
            
        Example:
            ele = topo.dem.variables['elevation'][:]
            #flatness for the finest step.
            F1 = topo.flatness(ele, out_xy = None, Tf = 50)
        """
        
        #degree to meter
        yi,xi = self.pixelLength(self.lat, L)
        xi = np.array(xi)
        ##slope
        kernelLon = np.array([[1,0,-1],[2,0,-2],[1,0,-1]])
        kernelLat = np.array([[1,2,1],[0,0,0],[-1,-2,-1]])
        dLon = convolve(ele, kernelLon)*100/(8*xi[:,None])
        dLat = convolve(ele, kernelLat)*100/(8*yi)
        slope = np.sqrt(dLon**2+dLat**2)
        del dLon, dLat

        #slope to flatness
        flatness = self.scale(slope, Tf, Pf)
        del slope
        #interpolate to the given sites
        if L <= 2:
            if not (out_xy is None):
                f = RegularGridInterpolator((self.lat[::-1], self.lon),
                    flatness[::-1,:], method = 'linear',bounds_error = False)
                flatness = f(out_xy)
            return flatness

        else:
            flatness = self.refine(L=L, coarseValue=flatness, out_xy=out_xy)
            return flatness


    def __percentile(self, x):
        """Return the rank of the element in surrouding cells.
        The function is called by lowness function"""
        return np.sum(x <= x[(x.size-1)/2])


    def lowness(self, ele, out_xy=None, L=1, Tl=0.4, Pl=3, lowRadius=13):
        """
        Returns lowness for given sites, which is measured as the radio of 
        number of points of lower elevation to the total number of points in 
        the surrounding region.

        Args:
            ele: array_like, input elevation
            out_xy: array_like, the interested sites
            L: interger, step number
            lowRadius: interger, the DEM cell size, half the number of cells
            used for the remainder of the steps
            L: interger, step number
            Tl: threshold, used to transform elevation percentile to lowness
            by scale() function
            Pl: shape parameter,used to transform elevation percentile to
            lowness by scale() function,
                  larger values give more abrupt transitions
                  
        Returns:
            Return lowness value at L step.
            
        Example:
            ele = topo.dem.variables['elevation'][:]
            #lowness for the finest step.
            L1 = topo.lowness(ele, out_xy=None, lowRadius=7)
        """
        
        size = float(lowRadius)**2
        pctl = generic_filter(ele, self.__percentile,
                                size = lowRadius)
        pctl /= size
        lowness = self.scale(pctl, Tl , Pl)
        if L <= 2:
            if not (out_xy is None):
                f = RegularGridInterpolator((self.lat[::-1], self.lon[:]),
                    lowness[::-1,:], method = 'linear',bounds_error = False)
                lowness = f(out_xy)
            return lowness
        else:
            lowness = self.refine(L = L, coarseValue=lowness, out_xy=out_xy)
            return lowness


    def finestScale(self, out_xy = None, initTf = 50.0):
        """
        Return MRVBF2 and CF2, which are combined values determined
        from the first and second-scale steps.
        
        Args:
            out_xy: Array-like, sites need to be simulated. If not 
            given (None), all the DEM cells will be simulated.
            
        Returns:
            MRVBF2: Combined MRVBF erived from finest-scale and second
                steps
            CF2: Combined flatness derived from finest-scale and second
                steps
                
        Example:
            MRVBF2, CF2 = topo.finestScale()
        """
        
        #first step
        #ele = self.dem['elevation'][:]
        F1 = self.flatness(self.dem['elevation'][:], out_xy=out_xy, Tf=initTf)
        L1 = self.lowness(self.dem['elevation'][:], out_xy=out_xy, lowRadius=7)
        PVF1 = F1*L1
        VF1 = 1 - self.scale(PVF1, 0.3, 4)
        #second step
        F2 = self.flatness(self.dem['elevation'][:], out_xy=out_xy,Tf=initTf/2)
        L2 = self.lowness(self.dem['elevation'][:], out_xy=out_xy, lowRadius=13)
        PVF2 = F2*L2
        VF2 = 1 - self.scale(PVF2, 0.3,4)
        #mrvbf2
        w2 = 1 - self.scale(VF2, 0.4, 6.68)
        MRVBF2 = w2*(1+VF2) + (1-w2)*VF1
        CF2 = F1*F2

        return MRVBF2, CF2

  
    def nmrvbf(self,out_xy = None, initTf = 50.0):
        """
        Returns smoothed mrvbf.
        
        Args:
        out_xy: Array-like, sites need to be simulated. If not 
            given (None), all the DEM cells will be simulated.
        initTf: first slope threshold
        
        Returns:
            mrvbf: smoothed and scaled mrvbf
        """
        mrvbf, cf = self.finestScale(initTf = initTf, out_xy = out_xy)
        sdemL = self.smoothDEM(self.dem['elevation'][:]) # smoothed base resolution dem
        meanKernel = np.full((3,3), 1.0/(3*3))
        for L in range(3,9):
            #print L
            shape = sdemL.shape
            Tf = initTf / (2**(L-1))
            #aggreate sdem for L step
            latIndex = range(1, shape[0], 3)
            lonIndex = range(1, shape[1], 3)
            sdemL = convolve(sdemL, meanKernel)[latIndex, :]
            sdemL = sdemL[:, lonIndex]

            #MRVBF for L step
            np.seterr(invalid = 'ignore')
            FL = self.flatness(sdemL, Tf, out_xy = out_xy,L = L)
            LL = self.lowness(sdemL, out_xy = out_xy, L= L, lowRadius = 13)
            cf = cf*FL
            PVFL = cf*LL
            VFL = 1 - self.scale(PVFL, 0.3,4)
            wL = 1 - self.scale(VFL, 0.4, np.log10((L-0.5)/0.1)/np.log10(1.5))
            mrvbf = wL*(L-1+VFL) + (1-wL)*mrvbf
			
        mrvbf/=8

        return mrvbf
        

    def aggregation(self, dem, scaleFactor = 3):
        """"
        Return a aggregated DEM, in which the cell size is is increased by
        a factor of scaleFactor.
        
        Args:
            dem: Input orginal fine-scale dem
            scaleFactor: A factor used to aggregate dem.
            
        Returns:
            aggDEM: Aggregated coarse-scale dem.
            latIndex: Index of aggregated dem in respect of origincal latitude
            lonIndex: Index of aggregated dem in respect of origincal lontitude
            
        Example:
            ele = topo.dem.variables['elevation'][:]
            aggDem, latIndex, lonIndex = topo.aggregation(ele, 5)
        """
        
        shape = dem.shape
        meanKernel = np.full((scaleFactor,scaleFactor), 
                             1.0/(scaleFactor*scaleFactor))
        latIndex = range((scaleFactor-1)/2, shape[0], scaleFactor)
        lonIndex = range((scaleFactor-1)/2, shape[1], scaleFactor)
        aggDem = convolve(dem, meanKernel)[latIndex, :]
        aggDem = aggDem[:, lonIndex]

        return [aggDem, latIndex, lonIndex]


    def aroundArea(self, centerSite, size = 30):
        """
        Return the neighbouring elevation of givens site.
        Args:
            centerSite: 
                The center cell [lat, lon] of the surrounding 
                area. All the cells with the distance <= size/2 will be clipped.
            size: 
                Diameter of surrounding area in km.
        
        Returns:
            Array-like dataframe.
            
        Example:
            centerSite = [46.749374, 9.5506067]
            eleSubset = topo.aroundArea(centerSite)
        """
        
        # area radiaus
        yi = self.pixelLength(self.lat)[0] # cell size in m
        lowRadius = int(size*1000/(2*yi)) # radius in pixel
        # nearest cell index
        latPosition = np.abs(self.lat-centerSite[0]).argmin()
        lonPosition = np.abs(self.lon-centerSite[1]).argmin()
        # subset area corner index
        latli = latPosition - lowRadius # bottom index
        latui = latPosition + lowRadius # top index
        lonli = lonPosition - lowRadius # left index
        lonri = lonPosition + lowRadius # right index
        #nearest area with bound size
        eleSubset = self.dem.variables['elevation'][latli:latui, lonli:lonri]
        return eleSubset

    
    def siteHypso(self, out_xy, bound = 30):
        """
        Returns the hypsometric position in the surrouding area.
        
        Args:
            out_xy: Site to simulate hypsometric position.
            bound: Diameter of surrounding area in km.
            
        Returns:
            lowness: Hyposmetric position in the surrouding area. Lowness
            ranges from 1 (deepest valley) to 0 (highest peak).
            
        Example:
            out_xy = np.array([[46.749374, 9.5506067, 1556.0],
                               [46.665974, 9.6340027, 939.0]])
            lowness = topo.siteHypso(out_xy)
        """
        
        lowness = []
        for i, site in enumerate(out_xy):
            ind_ele = self.aroundArea(site, size = bound)
            pctl = np.sum(ind_ele >= out_xy[i,2])/float(ind_ele.size)
            lowness.append(pctl)
        
        return lowness


    def hypso(self, bound = 30, out_xy = None):
        """
        Return hypsometric position based on fine scale of dem.

        Args:
            bound: Diameter of surrounding size
            
        Returns:
            lowness: Hypsometric position of all dem cells derived from
            fine scale of dem.
            
        Examples:
            hypso = topo.hypso(bound = 30, out_xy = None)
        """
        
        yi = self.pixelLength(self.lat)[0]
        lowRadius = int((bound*1000/yi))
        pctl = generic_filter(self.dem['elevation'][:], 
                              self.__percentile, size = lowRadius)
        pctl = pctl/(float(lowRadius)**2)
        lowness = 1-pctl
        del pctl

        if not (out_xy is None):
            lowInterp = RegularGridInterpolator((self.lat[::-1],self.lon),
                        lowness[::-1], method = 'linear', bounds_error = False)
            lowness = lowInterp(out_xy)
        return lowness


    def coarseHypso(self, bound = 30, out_xy = None):
        """
        Return hypsometric position based on coarse scale of dem.

        Args:
            bound: Diameter of surrounding size
            
        Returns:
            lowness: Hypsometric position of all dem cells derived from
            fine scale of dem.
            
        Examples:
            hypso = topo.hypso(bound = 30, out_xy = None)
        """

        yi = self.pixelLength(self.lat)[0]
        scaleFactor = int(np.ceil(500//yi) // 2 * 2 + 1)
        aggDem = self.aggregation(self.dem['elevation'][:], scaleFactor)

        #lowness
        lowRadius = 61#bound*1000/500
        pctl = generic_filter(aggDem[0], 
                              self.__percentile, 
                              size = lowRadius)/(float(lowRadius)**2)
        lowness = 1 - pctl

        #refine
        lowInterp = RegularGridInterpolator((self.lat[aggDem[1]][::-1],
                                             self.lon[aggDem[2]]),
            lowness[::-1], method = 'linear', bounds_error = False)


        if not (out_xy is None):
            lowness = lowInterp(out_xy)

        else:
            lon, lat = np.meshgrid(self.lon, self.lat)
            gridBase = np.array([lat.reshape(lat.size), 
                                 lon.reshape(lon.size)]).T
            lowness = lowInterp(gridBase).reshape((len(self.lat), len(self.lon)))

        return lowness#,pctl

    
    def eleRange(self, bound = 30, out_xy = None):
        """
        Returns elevation range for each dem cell within a area.
        If the input DEM is too small elevation range is condisered as 
        a constant value, which is derived from the center 30 km × 30 km area,
        for the entire DEM. Otherwise, elevation range is derived from the 
        neighbourhood 30 km area.
        
        Args:
            bound: Diameter of surrounding size in km
            
        Returns:
            rangeE: Array like, range of elevation range.
        
        
        Example:
            eleRange = topo.eleRange(out_xy = None, bound = 30)
        """
        
        if self.sizeCheck():
            centerlat = self.lat[len(self.lat)/2]
            centerlon = self.lon[len(self.lon)/2]
            subEle = self.aroundArea([centerlat, centerlon])
                
            rangeCon = np.max(subEle) - np.min(subEle)
                
            if not (out_xy is None):
                rangeE = np.repeat(rangeCon, len(out_xy))
            else:
                rangeE = np.ones(self.shape)*rangeCon
                
        else:
            lowRadius = bound*1000/(self.pixelLength(self.lat)[0])
            dem = gaussian_filter(self.dem['elevation'][:], np.sqrt(4.5))
            minEle = minimum_filter(dem, size = lowRadius)
            maxEle = maximum_filter(dem, size = lowRadius)
            rangeE = maxEle - minEle
            del minEle, maxEle, dem
            if not (out_xy is None):
                rangeInterp = RegularGridInterpolator((self.lat[::-1],self.lon),
                rangeE[::-1], method = 'linear', bounds_error = False)
                rangeE = rangeInterp(out_xy)
        
        return rangeE
		
    
		
class landSurCorrectionFac(object):
    """
    Returns land surface correction factor based on fine-scale DEM.
    
    Args:
        dem: fine-scale DEM in netcdf format
        demResolution: resolution of the input dem
        alpha: Adjust constant value
        beta:  A factor relating to fractional influence of surface effects on
               air temperature.
        gamma: A factor relating to cold air pooling on air temperature.
        
        The default values are derived from the Swiss Alps. Please see details
        from the REDCAPP paper.
        
    Returns:
        lscf: Land surafce correction factor
        
    Example:
        dem = 'DEM_testArea.nc'
        demRes = 3.0/3600
        lscf = landSurCorrectionFac(dem, demRes)
        LSCF = lscf.correctionFactor(hypso, mrvbf, eleRange)
    """
    
    def __init__(self, dem, demResolution, alpha= 0.61, beta= 1.56, gamma= 465):
        self.dem   = dem
        self.gamma = gamma
        self.beta  = beta
        self.alpha = alpha
        self.resolution = demResolution
        self.dem_ncdf = nc.Dataset(self.dem, 'r')
    
    def scale(self, eleRange):
        
        return(np.exp(-eleRange/self.gamma))
    
    def LSCF(self, hypso, mrvbf, eleR):
        """Returns land surface correction factor"""
        s = self.scale(eleR)
        h = hypso* (1- s) + s
        v = mrvbf*(1-s)
        
        lscf = self.alpha*h + self.beta*v
        
        return lscf
    
        
    def spatialLSCF(self, file_out):
        """"Returns and export spatialized land surface correction"""
        
        topo = topography(self.dem, self.resolution)
        mrvbf = topo.nmrvbf(out_xy = None, initTf = 50.0)
        hypso = topo.coarseHypso()
        eleR  = topo.eleRange()
        lscf = self.LSCF(hypso, mrvbf, eleR)
        
        
        # ---- export geomorphometric factors ---------------------------------
        # create nc file
        nc_root = nc.Dataset(file_out ,'w', format = 'NETCDF4_CLASSIC')
        
        # create dimensions
        nc_root.createDimension('lat', mrvbf.shape[0])
        nc_root.createDimension('lon', mrvbf.shape[1])
        
        # create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes  = nc_root.createVariable('lat', 'f4', ('lat'))
        Hypso      = nc_root.createVariable('hypso','f4',('lat','lon'),zlib=True)
        Mrvbf      = nc_root.createVariable('mrvbf','f4',('lat','lon'),zlib=True)
        RangeE     = nc_root.createVariable('eleRange','f4',('lat','lon'),zlib=True)
        Lscf       = nc_root.createVariable('lscf','f4',('lat','lon'), zlib=True)
        
        longitudes.setncatts({'long_name': u"longitude"})
        latitudes.setncatts({'long_name': u"latitude"})
        Mrvbf.setncatts({'long_name': 
                         "normalized multiresolution index of valley bottom flatness"})
        Hypso.setncatts({'long_name': u"hyposmetric position"})
        RangeE.setncatts({'long_name': u"elevation range in prescirbed neighbourhood"})
        Lscf.setncatts({'long_name': u"Land surface correction factor"})
        
        #assign variables
        longitudes[:] = self.dem_ncdf.variables['lon'][:]
        latitudes[:] = self.dem_ncdf.variables['lat'][:]
        Hypso[:] = hypso
        Mrvbf[:] = mrvbf
        RangeE[:] = eleR
        Lscf[:] = lscf
        
        #attribute
        nc_root.description = "fine-scale DEM-derived topographic factors"
        longitudes.units = 'degree_east (demical)'
        latitudes.units = 'degree_north (demical)'
        
        nc_root.close()
        
        return lscf
    
    
    def stationLSCF(self, stations, file_out):
        """Returns station land surface correction factor"""
        out_xyz = np.asarray([[s['lat'],s['lon'],s['ele']] for s in stations])
        
        topo = topography(self.dem, self.resolution)
        mrvbf = topo.nmrvbf(out_xy = out_xyz[:,:2], initTf = 50.0)
        hypso = topo.siteHypso(out_xy = out_xyz, bound = 30)
        eleR = topo.eleRange(out_xy = out_xyz[:,:2], bound = 30)
        
        lscf = self.LSCF(hypso, mrvbf, eleR)
        
        names = [s['name']for s in stations]
        with open (file_out, 'wb') as f:
            writer = csv.writer(f)
            writer.writerow(['station', 'hypso', 'mrvbf', 'eleRange', 'LSCF'])
            
            values = np.asarray([names, mrvbf, hypso, eleR, lscf]).T
            
            for row in values:
                writer.writerow(row)
        
        return lscf
    
    
class redcappTemp(object):
    """returns REDCAPP derived surface air temperature for both
    given dem area (spatialized mean air temperature) and 
    stations (air temperature time series)
    
    Args:
        geop: geopotential file, considered as coarese scale of 
            topography.
        sa: surface air temperature from renalysis
        pl: pressure level temperature from renalysis
        variable: variables to conducted here. Variable of 'Temperature' is
            used here
        daterange: the date range to be downscaled
        dem: DEM file in netcdf format, consiered as fine scale of topography
        demResolution: resolution of input dem
    
    Returns:
        returns (1) spatalized mean air temperature in netcdf format, if input
        dem file; (2) air temperature time series in csv format, if input 
        station information
        
    Examples:
        # input data
        geop = 'ecmwf_erai_to.nc'
        sa   = 'ecmwf_erai_sa_m_151201_151231.nc'
        pl   = 'ecmwf_erai_pl_m_151201_151231.nc'
        
        variable = 'Temperature'
        dem_ncf  = 'DEM_testArea.nc'
        resolution = 3.0/3600
        
        date  = {'beg' : datetime(2015,12,1,00,00),
                 'end' : datetime(2015,12,10,18,00)}
        
        # run
        Redcapp = redcappTemp(geop, sa, pl, variable, date, dem_ncdf, resolution)

        # SPATIALIZED MEAN AIR TEMPERATURE
        Redcapp.extractSpatialDataNCF(spatTemp_out)

        # AIR TEMPERATURE TIME SERIES
        Redcapp.extractStationDataCSV(statTemp_out, stations)
        
    """
    
    def __init__(self, geop, sa, pl, variable, daterange, dem, demResolution,
                 alpha= 0.61, beta= 1.56, gamma= 465):
        self.geop = geop
        self.sa = sa
        self.pl = pl
        self.variable = variable
        self.daterange = daterange
        self.dem = dem
        self.demncdf = nc.Dataset(dem, 'r')
        self.resolution = demResolution
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        
    def edgeClip(self, values):
        """
        Drops the egde of temperature without data caused by mrvbf simulation.
        """
        
        #the corner with values
        shape = values.shape
        center = [i/2 for i in shape]
        left = np.min(np.where(np.isfinite(values[center[0],:])))#left
        right = np.max(np.where(np.isfinite(values[center[0],:])))+1#right
        
        upper = np.min(np.where(np.isfinite(values[:,center[1]])))#upper
        low = np.max(np.where(np.isfinite(values[:,center[1]])))+1#low
        
        values = values[upper:low, left:right]
           
        lons = self.demncdf.variables['lon'][left:right]
        lats = self.demncdf.variables['lat'][upper:low]
        
        return values, lons, lats
        
        
    def spatialTemp(self, topo_out):
        """Returns spatialized mean air temperature."""
        #upp-air temperature and coarse land-surface effects
        Downscaling = downscaling(self.geop, self.sa, self.pl, self.dem)
        pl, dt = Downscaling.spatialMean(self.variable, self.daterange)
        
        #lscf
        print "Temperature Done!"
        print "Conducting terrain analysis..."
        LSCF = landSurCorrectionFac(self.dem, self.resolution)
        lscf = LSCF.spatialLSCF(topo_out)
        
        #redcapp temperaure
        temp = pl+lscf*dt
        
        temp, lons, lats = self.edgeClip(temp)
        
        return temp, lons, lats
    
    def stationTemp(self, stations, topo_out):
        """Returns air temperature time series."""
        
        #upp-air temperature and coarse land-surface effects
        Downscaling = downscaling(self.geop, self.sa, self.pl)
        
        pl,dt,time,names = Downscaling.stationTimeSeries(self.variable, 
                                                         self.daterange, 
                                                         stations)
        
        #lscf
        print "Temperature Done!"
        print "Conducting terrain analysis..."
        LSCF = landSurCorrectionFac(self.dem, self.resolution)
        lscf = LSCF.stationLSCF(stations, topo_out)
        
        #redcapp temperature
        temp = np.zeros((pl.shape))# hold temp

        for i in np.arange(len(names)):
            pli = pl[:,i]
            dti = pl[:,i]
            lscfi = lscf[i]
            temp[:,i] = pli + lscfi*dti
        
        return temp, time, names
        
    
    def extractSpatialDataNCF(self, topo_out, temp_out):
        """"Export spatialized mean air temperature of given dem
        in netcdf format. Please not that the area of output spatial 
        temperature will be smaller than the given dem owing to the 
        mrvbf simulation donot include the edge. Plsease also see the 
        introduction of clipEdge() fucntion"""
        
        temp, lons, lats = self.spatialTemp(topo_out)
        
        #create nc file
        nc_root = nc.Dataset(temp_out ,'w', format = 'NETCDF4_CLASSIC')
        
        #create dimensions
        nc_root.createDimension('lat', len(lats))
        nc_root.createDimension('lon', len(lons))
        
        #create variables
        longitudes = nc_root.createVariable('lon', 'f4', ('lon'))
        latitudes = nc_root.createVariable('lat', 'f4', ('lat'))
        Ta = nc_root.createVariable('surface air temperature', 
                                     'f4', ('lat', 'lon'), zlib = True)
        
        #assign variables
        longitudes[:] = lons
        latitudes[:] = lats
        Ta[:] = temp

        #attribute
        nc_root.description = "REDCAPP-derived surface air temperature"
        longitudes.units = 'degree_east (decimal)'
        latitudes.units  = 'degree_north (decimal)'
        Ta.units = 'celsius'
        
        nc_root.close()

    def extractStationDataCSV(self, stations, topo_out, temp_out):
        """
        Exports air temperature time series of the given stations 
        in csv format
        """
        
        temp,time,names = self.stationTemp(stations, topo_out)

        #write CSV
        names.insert(0, 'Time_UTC')
        with open(temp_out, 'wb') as output_file:
            writer = csv.writer(output_file)
            writer.writerow(names)
            valu = temp.tolist()
            for n in range(len(time)):
                row = ['%.3f' % elem for elem in valu[n]]
                row.insert(0,time[n])
                writer.writerow(row)
	
      