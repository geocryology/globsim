# -*- coding: utf-8 -*-
#
# Methods for downloading ERA5 data from the ECMWF server for limited
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
from __future__ import print_function

import re

import numpy   as np
import netCDF4 as nc
import cdsapi

from datetime import datetime, timedelta
from math     import exp, floor, atan2, pi
from os       import path, listdir, makedirs
from ecmwfapi.api import ECMWFDataServer
from scipy.interpolate import interp1d
from fnmatch import filter as fnfilter

import urllib3
urllib3.disable_warnings()

from globsim.generic import ParameterIO, StationListRead,  series_interpolate, variables_skip,  str_encode, GenericDownload, GenericScale, GenericInterpolate
from globsim.nc_elements import netcdf_base, new_interpolated_netcdf, new_scaled_netcdf
from globsim.meteorology import spec_hum_kgkg, pressure_from_elevation


try:
    from nco import Nco
except ImportError:
    print("*** NCO not imported, netCDF appending not possible. ***")
    pass 

try:
    import ESMF
    # Check ESMF version.  7.0.1 behaves differently than 7.1.0r 
    ESMFv = int(re.sub("[^0-9]", "", ESMF.__version__))
    ESMFnew = ESMFv > 701
    
except ImportError:
    print("*** ESMF not imported, interpolation not possible. ***")
    pass   


class ERA5generic(object):
    """
    Parent class for other ERA5 classes.
    
    """           
    CDS_DICT =  {'130.128' : 'temperature',
                     '157.128' : 'relative_humidity',
                     '131.128' : 'u_component_of_wind',
                     '132.128' : 'v_component_of_wind',
                     '129.128' : 'geopotential',
                     '167.128' : '2m_temperature',
                     '168.128' : '2m_dewpoint_temperature',
                     '206.128' : 'total_column_ozone',        # also consider surface_solar_radiation_downward_clear_sky
                     '137.128' : 'total_column_water_vapour', # also consider surface_solar_radiation_downward_clear_sky
                     '165.128' : '10m_u_component_of_wind',
                     '166.128' : '10m_v_component_of_wind',
                     '228.128' : 'total_precipitation',
                     '169.128' : 'surface_solar_radiation_downwards',
                     '175.128' : 'surface_thermal_radiation_downwards',
                     '172.128' : 'land_sea_mask'
                     }
                     
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

    def getPressureLevels(self, elevation):
        """Restrict list of ERA5 pressure levels to be downloaded"""
        Pmax = pressure_from_elevation(elevation['min']) + 55
        Pmin = pressure_from_elevation(elevation['max']) - 55 
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
          'dataset' : "era5",
          'stream'  : "oper",
          'class'   : "ea",
          'format'  : "netcdf",
          'grid'    : "0.25/0.25"}
        return dictionary_gen
 
    def getDstring(self):
        return ('_' + self.date['beg'].strftime("%Y%m%d") + "_to_" +
                      self.date['end'].strftime("%Y%m%d"))
    
    def download(self, api, storage):
        if api == 'ecmwf':
            server = ECMWFDataServer()
            print(server.trace('=== ERA5 ({}API): START ACCESS ON {} ===='.format('ECMWF', storage.upper())))
            server.retrieve(self.getDictionary())
            print(server.trace('=== ERA5 ({}API): END ACCESS ON {} ===='.format('ECMWF', storage.upper())))
        
        elif api == 'cds':
            self.download_cds(storage)
    
    def download_cds(self, storage='cds'):
        ''' 'patched' download function to allow download using cdsapi from MARS or CDS servers  '''
        server = cdsapi.Client()
        
        query = self.getDictionary()
        target = query.pop('target')
        levtype = query.pop('levtype')
        
        # select dataset
        if storage == 'ecmwf':
            dataset = 'reanalysis-era5-complete'
        elif storage == 'cds':
            dataset = {'pl' : 'reanalysis-era5-pressure-levels',
                       'sfc': 'reanalysis-era5-single-levels'
                      }[levtype]
            query = self.ECM2CDS(query)
        
        # launch download request
        print(server.info('=== ERA5 ({}API): START ACCESS ON {} ===='.format("CDS", storage.upper())))
        if path.isfile(target):
             print("WARNING: File '{}' already exists and was skipped".format(target))
        else:
            server.retrieve(dataset, query, target)
        print(server.info('=== ERA5 ({}API): END ACCESS ON {} ===='.format("CDS", storage.upper())))

    def ECM2CDS(self, query):
        ''' convert ECMWF query to CDS format '''

        # remove unnecessary keys
        for key in ['class', 'dataset', 'stream', 'type']:
            query.pop(key)
        
        # replace key('param'):val(str)  string with key('variables'):val(list)
        query['variable'] = [self.CDS_DICT[L] for L in query.pop('param').split('/')]
        
        # replace levellist if in dictionary (also string '100/200/300' to list ['100', '200', '300'])
        if 'levellist' in query:
            query['pressure_level'] = [L for L in query.pop('levellist').split('/')]
            
        # reformat date range (start/to/end --> start/end)
        query['date'] = re.sub("/to", "", query['date'])
        
        # reformat time string into list
        query['time'] = [t[0:5] for t in query['time'].split('/')]
        
        # add product type (one of Reanalysis, Ensemble members, Ensemble mean, Ensemble spread)
        query['product_type'] = 'reanalysis'
        
        return(query)
        
    def TranslateCF2ERA(self, variables, dpar):
        """
        Translate CF Standard Names into ERA5 code numbers.
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
                  "ERA5 data: {0}")
        return string.format(self.getDictionary) 


class ERA5pl(ERA5generic):
    """Returns an object for ERA5 data that has methods for querying the
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
        ERA5pl = ERA5pl(date, area, elevation, variables, directory) 
        ERA5pl.download()
    """
    def __init__(self, date, area, elevation, variables, directory):
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        outfile = 'era5_pl' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)
 
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable names. 
        dpar = {'air_temperature'   : '130.128',           # [K]
                'relative_humidity' : '157.128',           # [%]
                'wind_speed'        : '131.128/132.128'}   # [m s-1]
                
        # translate variables into those present in ERA pl data        
        self.TranslateCF2ERA(variables, dpar)
        self.param += '/129.128' # geopotential always needed [m2 s-2]
        self.time = ("00:00:00/01:00:00/02:00:00/03:00:00/"
                     "04:00:00/05:00:00/06:00:00/07:00:00/"
                     "08:00:00/09:00:00/10:00:00/11:00:00/"
                     "12:00:00/13:00:00/14:00:00/15:00:00/"
                     "16:00:00/17:00:00/18:00:00/19:00:00/"
                     "20:00:00/21:00:00/22:00:00/23:00:00")
    
    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "pl",
           'levellist': self.getPressureLevels(self.elevation),
           'time'     : self.time,
           'step'     : "0",
           'type'     : "an",
           'param'    : self.param,
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary    
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA5 data: {0}")
        return string.format(self.getDictionary) 

                                
class ERA5sa(ERA5generic):
    """
    Returns an object for ERA5 data that has methods for querying the
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
        ERA5sa = ERA5sa(date, area, variables, directory) 
        ERA5sa.download()      
    """
    def __init__(self, date, area, variables, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile = 'era5_sa' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)
        
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable names. 
        dpar = {'air_temperature'   : '167.128',  # [K] 2m values
                'relative_humidity' : '168.128',  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky' : 
                    '206.128/137.128',  # [kg m-2] Total column ozone 
                                        # [kg m-2] Total column W vapor                                                             
                'wind_speed' : '165.128/166.128'}   # [m s-1] 10m values
                
        # translate variables into those present in ERA pl data        
        self.TranslateCF2ERA(variables, dpar)
    
    def getTime(self):
        
        times = np.arange(0, 24)
        times = [str(t).zfill(2) for t in times]
        times = [t+':00:00' for t in times]
        times = '/'.join(times)
        
        return times

    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : self.getTime(),
           'step'     : "0",
           'type'     : "an",
           'param'    : self.param,
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary  
        
    def __str__(self):
        string = ("Class for ERA5 surface analysis data: {0}")
        return string.format(self.getDictionary)         


class ERA5sf(ERA5generic):
    """
    Returns an object for ERA5 data that has methods for querying the
    ECMWF server for surface forecast variables (prec, swin, lwin).
       
    Args:        
        target:    File name of the netcdf file to be created.
        
        variables:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].
        
        ERAgen:    ERA5generic() object with generic information on the area and
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
        ERA5sf = ERA5sf(date, area, variables, directory) 
        ERA5sf.download()   
    """
    def __init__(self, date, area, variables, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile = 'era5_sf' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable names. 
        # [m] total precipitation
        # [J m-2] short-wave downward
        # [J m-2] long-wave downward
        dpar = {'precipitation_amount'              : '228.128', #[kg m**-2 s**-1]
                'downwelling_shortwave_flux_in_air' : '169.128', #[J m**-2]
                'downwelling_longwave_flux_in_air'  : '175.128'} #[J m**-2]
                
        # translate variables into those present in ERA pl data        
        self.TranslateCF2ERA(variables, dpar)

    def getStep(self):
        steps = np.arange(1, 13)
        steps = [str(s) for s in steps]
        steps = '/'.join(steps)
        
        return steps

    def getDictionary(self):
    
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : ("00:00:00/01:00:00/02:00:00/03:00:00/"
                     "04:00:00/05:00:00/06:00:00/07:00:00/"
                     "08:00:00/09:00:00/10:00:00/11:00:00/"
                     "12:00:00/13:00:00/14:00:00/15:00:00/"
                     "16:00:00/17:00:00/18:00:00/19:00:00/"
                     "20:00:00/21:00:00/22:00:00/23:00:00"),
           'step'     : self.getStep(),
           'type'     : "fc",
           'param'    : self.param,
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary  
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA5 air tenperature data: {0}")
        return string.format(self.getDictionary)      


class ERA5to(ERA5generic):
    """
    Returns an object for downloading and handling ERA5
    topography (invariant).
       
    Args:      
              
    Example:
        area  = {'north' :  40.0,
                 'south' :  45.0,
                 'west'  :  60.0,
                 'east'  :  65.0}            
        directory = '/Users/stgruber/Desktop'             
        ERA5to = ERA5to(area, directory) 
        ERA5to.download()       
    """
    def __init__(self, area, directory):
        self.area       = area
        self.date       = {'beg' : datetime(2017, 1, 1),
                           'end' : datetime(2017, 1, 1)}
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'era5_to.nc')

    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "00:00:00",
           'step'     : "0",
           'type'     : "an",
           'param'    : "129.128/172.128", # geopotential and land-sea mask
           'target'   : self.file_ncdf
           } 
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary   
        
    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA5 air tenperature data: {0}")
        return string.format(self.getDictionary) 


class ERA5download(GenericDownload):
    """
    Class for ERA5 data that has methods for querying 
    the ECMWF server, returning all variables usually needed.
       
    Args:
        pfile: Full path to a Globsim Download Parameter file. 
        api  : Which API to use. Either 'cds' (default) or 'ecmwf' (deprecated'
        storage : Which server to access data from. Either 'cds'  (default) or 'ecmwf' (deprecated). Note
                 that you can use the cds api to access the ecmwf storage (MARS)
    Example:          
        ERAd = ERA5download(pfile) 
        ERAd.retrieve()
    """
        
    def __init__(self, pfile, api='cds', storage='cds'):
        super().__init__(pfile)
        par = self.par
        self._set_data_directory("era5")
        
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}

        # chunk size for downloading and storing data [days]        
        self.chunk_size = par.chunk_size   

        # set api and download server/storage
        self.api = api
        self.storage = storage

    def retrieve(self):
        """
        Retrieve all required ERA5 data from MARS/CDS server.
        """        
        # prepare time loop
        date_i = {}
        slices = floor(float((self.date['end'] - self.date['beg']).days)/
                       self.chunk_size)+1
                                                 
        # topography
        if path.isfile(path.join(self.directory,'era5_to.nc')):
            print("WARNING: File 'era5_to.nc' already exists. Skipping.")
        else: 
            top = ERA5to(self.area, self.directory)
            top.download(self.api, self.storage)
        
        for ind in range (0, int(slices)): 
            #prepare time slices   
            date_i['beg'] = self.date['beg'] + timedelta(days = 
                            self.chunk_size * ind)
            date_i['end'] = self.date['beg'] + timedelta(days = 
                            self.chunk_size * (ind+1) - 1)
            if ind == (slices-1):
                date_i['end'] = self.date['end']
            
            #actual functions                                                                           
            pl = ERA5pl(date_i, self.area, self.elevation, 
                       self.variables, self.directory) 
            sa = ERA5sa(date_i, self.area, self.variables, self.directory) 
            sf = ERA5sf(date_i, self.area, self.variables, self.directory) 
        
            #download from ECMWF server convert to netCDF  
            ERAli = [pl, sa, sf]
            for era in ERAli:
                era.download(self.api, self.storage)          
        
        # report inventory
        self.inventory()  
        
    def inventory(self):
        """
        Report on data avaialbe in directory: time slice, variables, area 
        """
        print("\n\n\n")
        print("=== INVENTORY FOR GLOBSIM ERA5 DATA ===\n")
        print("Download parameter file: \n" + self.pfile + "\n")
        # loop over filetypes, read, report
        file_type = ['era5_pl_*.nc', 'era5_sa_*.nc', 
                     'era5_sf_*.nc', 'era5_t*.nc']
        for ft in file_type:
            infile = path.join(self.directory, ft)
            nf = len(fnfilter(listdir(self.directory), ft))
            print(str(nf) + " FILE(S): " + infile)
            
            if nf > 0:
                # open dataset
                ncf = nc.MFDataset(infile, 'r', aggdim='time')
                
                # list variables
                keylist = [str_encode(x) for x in ncf.variables.keys()]
                print("    VARIABLES:")
                print("        " + str(len(keylist)) + 
                      " variables, including dimensions")
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
        return "Object for ERA5 data download and conversion"


class ERA5interpolate(GenericInterpolate):
    """
    Collection of methods to interpolate ERA5 netCDF files to station
    coordinates. All variables retain their original units and time stepping.
    """
        
    def __init__(self, ifile):
        super().__init__(ifile)
        par = self.par        
       
        self.dir_inp = path.join(par.project_directory,'era5') 
        
        #convert longitude to ERA notation if using negative numbers  
        self.stations['longitude_dd'] = self.stations['longitude_dd'] % 360             

    def ERA2station(self, ncfile_in, ncfile_out, points,
                    variables = None, date = None):
        
        """
        Biliner interpolation from fields on regular grid (latitude, longitude) 
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all ERA5 files). The type 
        of variable and file structure are determined from the input.
        
        This function creates an empty of netCDF file to hold the interpolated 
        results, by calling new_interpolated_netcdf(). Then, data is 
        interpolated in temporal chunks and appended. The temporal chunking can 
        be set in the interpolation parameter file.
        
        Args:
        ncfile_in: Full path to an ERA5 derived netCDF file. This can
                   contain wildcards to point to multiple files if temporal
                  chunking was used.
            
        ncfile_out: Full path to the output netCDF file to write.     
        
        points: A dictionary of locations. See method StationListRead in
                generic.py for more details.
    
        variables:  List of variable(s) to interpolate such as 
                    ['r', 't', 'u','v', 't2m', 'u10', 'v10', 'ssrd', 'strd', 'tp'].
                    Defaults to using all variables available.
    
        date: Directory to specify begin and end time for the derived time 
                series. Defaluts to using all times available in ncfile_in.
        
        cs: chunk size, i.e. how many time steps to interpolate at once. This 
            helps to manage overall memory usage (small cs is slower but less
            memory intense).          
        """
                
        # read in one type of mutiple netcdf files       
        ncf_in = nc.MFDataset(ncfile_in, 'r', aggdim ='time')
        
        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # build the output of empty netCDF file
        rootgrp = new_interpolated_netcdf(ncfile_out, self.stations, ncf_in, 
                                      time_units='hours since 1900-01-01 00:00:0.0') 
        rootgrp.source = 'ERA5, interpolated bilinearly to stations'
        rootgrp.close()   
                                  
        # open the output netCDF file, set it to be appendable ('a')
        ncf_out = nc.Dataset(ncfile_out, 'a')

        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]
        #"hours since 1900-01-01 00:00:0.0"
        t_unit = ncf_in.variables['time'].units 
        try :
            t_cal = ncf_in.variables['time'].calendar
        except AttributeError :  # attribute doesn't exist
            t_cal = u"gregorian" # standard
        time = nc.num2date(nctime, units = t_unit, calendar = t_cal)
        
        # detect invariant files (topography etc.)
        invariant = True if len(time) == 1 else False                                                                         
    
        # restrict to date/time range if given
        if date is None:
            tmask = time < datetime(3000, 1, 1)
        else:
            tmask = (time < date['end']) * (time >= date['beg'])
                          
        # get time vector for output
        time_in = nctime[tmask]     

        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter  = len(time_in)//self.cs
        niter += ((len(time_in) % self.cs) > 0)

        # loop over chunks
        for n in range(niter):
            # indices (relative to index of the output file)
            beg = n * self.cs
            # restrict last chunk to lenght of tmask plus one (to get last time)
            end = min(n*self.cs + self.cs, len(time_in))-1
            
            # time to make tmask for chunk 
            beg_time = nc.num2date(time_in[beg], units=t_unit, calendar=t_cal)
            if invariant:
                # allow topography to work in same code, len(nctime) = 1
                end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                #end = 1
            else:
                end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)
                
            #'<= end_time', would damage appending
            tmask_chunk = (time <= end_time) * (time >= beg_time)
            if invariant:
                # allow topography to work in same code
                tmask_chunk = [True]
                 
            # get the interpolated variables
            dfield, variables = self.interp2D(ncfile_in, ncf_in, 
                                                 self.stations, tmask_chunk,
                                                 variables=None, date=None) 
            # append time
            ncf_out.variables['time'][:] = np.append(ncf_out.variables['time'][:], 
                                                     time_in[beg:end+1])
                                  
            #append variables
            for i, var in enumerate(variables):
                if variables_skip(var):
                    continue
                                                              
                if pl:
                    if ESMFnew:
                        # dfield has dimensions (station, variables, time, pressure levels
                        ncf_out.variables[var][beg:end+1,:,:] = dfield.data[:,i,:,:].transpose((1,2,0))
                    else:
                        # dimension: time, level, station (pressure level files)
                        ncf_out.variables[var][beg:end+1,:,:] = dfield.data[i,:,:,:]        
                else:   
                    if ESMFnew:
                        # dfield has dimensions (station, variables, time)
                        ncf_out.variables[var][beg:end+1,:] = dfield.data[:,i,:].transpose((1,0))
                    else:
                        # dfield has dimensions time, station (2D files)
                        ncf_out.variables[var][beg:end+1,:] = dfield.data[i,:,:]

                                     
        #close the file
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
        #TODO check the aggdim does not work
        ncf = nc.MFDataset(ncfile_in, 'r', aggdim='time')
        height = ncf.variables['height'][:]
        nt = len(ncf.variables['time'][:])
        nl = len(ncf.variables['level'][:])
        
        # list variables
        varlist = [str_encode(x) for x in ncf.variables.keys()]
        for V in ['time', 'station', 'latitude', 'longitude', 'level','height','z']:
            varlist.remove(V)

        # === open and prepare output netCDF file ==============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = netcdf_base(ncfile_out, len(height), nt, 'hours since 1900-01-01 00:00:0.0')
        rootgrp.source = 'ERA5, interpolated (bi)linearly to stations'
        
        # assign base variables
        time[:]      = ncf.variables['time'][:]
        station[:]   = ncf.variables['station'][:]
        latitude[:]  = ncf.variables['latitude'][:]
        longitude[:] = ncf.variables['longitude'][:]
        height[:]    = ncf.variables['height'][:]
        
        # create and assign variables from input file
        for var in varlist:
            tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))    
            tmp.long_name = str_encode(ncf.variables[var].long_name)
            tmp.units     = str_encode(ncf.variables[var].units)  

        # add air pressure as new variable
        var = 'air_pressure'
        varlist.append(var)
        tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))    
        tmp.long_name = var.encode('UTF8')
        tmp.units     = 'hPa'.encode('UTF8') 
        # end file prepation ===================================================
                                                                                                
        # loop over stations
        for n, h in enumerate(height): 
            # convert geopotential [mbar] to height [m], shape: (time, level)
            elevation = ncf.variables['z'][:,:,n] / 9.80665
            # TODO: check if height of stations in data range
            
            # difference in elevation, level directly above will be >= 0
            elev_diff = elevation - h
            # vector of level indices that fall directly above station. 
            # Apply after ravel() of data.
            va = np.argmin(elev_diff + (elev_diff < 0) * 100000, axis=1) 
            # mask for situations where station is below lowest level
            mask = va < (nl-1)
            va += np.arange(elevation.shape[0]) * elevation.shape[1]
            
            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask # +1 when OK, +0 when below lowest level
            
            wa, wb = self.calculate_weights(elev_diff, va, vb)
            
            #loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                if var == 'air_pressure':
                    # pressure [Pa] variable from levels, shape: (time, level)
                    data = np.repeat([ncf.variables['level'][:]],
                                      len(time),axis=0).ravel() 
                else:    
                    #read data from netCDF
                    data = ncf.variables[var][:,:,n].ravel()
                    
                ipol = data[va]*wa + data[vb]*wb   # interpolated value                    
                rootgrp.variables[var][:,n] = ipol # assign to file   
    
        rootgrp.close()
        ncf.close()
        # closed file ==========================================================    

    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to 
        the more generically ERA-like interpolation functions.
        """                       

        # 2D Interpolation for Invariant Data      
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.            
        dummy_date  = {'beg' : datetime(1979, 1, 1, 12, 0),
                       'end' : datetime(1979, 1, 1, 12, 0)}     
        self.ERA2station(path.join(self.dir_inp,'era5_to.nc'), 
                         path.join(self.dir_out,'era5_to_' + 
                                   self.list_name + '.nc'), self.stations,
                                   ['z', 'lsm'], date = dummy_date)  
                                  
        # === 2D Interpolation for Surface Analysis Data ===    
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['t2m'],  # [K] 2m values
                'relative_humidity' : ['d2m'],  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky' : 
                    ['tco3', 'tcwv'],   # [kg m-2] Total column ozone 
                                        # [kg m-2] Total column W vapor                                                             
                'wind_speed' : ['u10', 'v10']}   # [m s-1] 10m values   
        varlist = self.TranslateCF2short(dpar)                      
        self.ERA2station(path.join(self.dir_inp,'era5_sa_*.nc'), 
                         path.join(self.dir_out,'era5_sa_' + 
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)          
        
        # 2D Interpolation for Surface Forecast Data    'tp', 'strd', 'ssrd' 
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.       
        # [m] total precipitation
        # [J m-2] short-wave downward
        # [J m-2] long-wave downward
        dpar = {'precipitation_amount'              : ['tp'],   
                'downwelling_shortwave_flux_in_air' : ['ssrd'], 
                'downwelling_longwave_flux_in_air'  : ['strd']} 
        varlist = self.TranslateCF2short(dpar)                           
        self.ERA2station(path.join(self.dir_inp,'era5_sf_*.nc'), 
                         path.join(self.dir_out,'era5_sf_' + 
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)          
                         
        
        # === 2D Interpolation for Pressure Level Data ===
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['t'],           # [K]
                'relative_humidity' : ['r'],           # [%]
                'wind_speed'        : ['u', 'v']}      # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('z')
        self.ERA2station(path.join(self.dir_inp,'era5_pl_*.nc'), 
                         path.join(self.dir_out,'era5_pl_' + 
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)
                                   
        # 1D Interpolation for Pressure Level Data
        self.levels2elevation(path.join(self.dir_out,'era5_pl_' + 
                                        self.list_name + '.nc'), 
                              path.join(self.dir_out,'era5_pl_' + 
                                        self.list_name + '_surface.nc'))      
                                                   

class ERA5scale(GenericScale):
    """
    Class for ERA5 data that has methods for scaling station data to
    better resemble near-surface fluxes.
    
    Processing kernels have names in UPPER CASE.
       
    Args:
        sfile: Full path to a Globsim Scaling Parameter file. 
              
    Example:          
        ERAd = ERA5scale(sfile) 
        ERAd.process()
    """
    NAME = "ERA5"
       
    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par
        
        # input file handles
        self.nc_pl = nc.Dataset(path.join(self.intpdir, 'era5_pl_' + 
                                self.list_name + '_surface.nc'), 'r')
        self.nc_sa = nc.Dataset(path.join(self.intpdir, 'era5_sa_' + 
                                self.list_name + '.nc'), 'r')
        self.nc_sf = nc.Dataset(path.join(self.intpdir, 'era5_sf_' + 
                                self.list_name + '.nc'), 'r')
        self.nc_to = nc.Dataset(path.join(self.intpdir, 'era5_to_' + 
                                self.list_name + '.nc'), 'r')
        self.nstation = len(self.nc_to.variables['station'][:])                     
                              
        # output file 
        self.output_file = self.getOutNCF(par, 'era5')
        
        # time vector for output data
        # get time and convert to datetime object
        nctime = self.nc_pl.variables['time'][:]
        self.t_unit = self.nc_pl.variables['time'].units #"hours since 1900-01-01 00:00:0.0"
        self.t_cal  = self.nc_pl.variables['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)
        
        # interpolation scale factor
        self.time_step = par.time_step * 3600    # [s] scaled file
        self.interval_in = (time[1]-time[0]).seconds #interval in seconds
        self.interpN = floor(self.interval_in/self.time_step)
        
        #number of time steps for output, include last value
        self.nt = int(floor((max(time) - min(time)).total_seconds()
                            / 3600 / par.time_step)) + 1 

        # vector of output time steps as datetime object
        # 'seconds since 1900-01-01 00:00:0.0'
        mt = min(time)
        self.times_out = [mt + timedelta(seconds = (x*self.time_step)) 
                          for x in range(0, self.nt)]                                                                   
                                      
        # vector of output time steps as written in ncdf file [s]
        self.scaled_t_units = 'seconds since 1900-01-01 00:00:0.0'
        self.times_out_nc = nc.date2num(self.times_out, 
                                        units = self.scaled_t_units, 
                                        calendar = self.t_cal) 
        
    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one 
        variable and adds it to the netCDF file.
        """    
        self.rg = new_scaled_netcdf(self.output_file, 
                                 self.nc_pl, self.times_out_nc,
                                 t_unit = self.scaled_t_units, 
                                 station_names = self.stations['station_name'])
        
        
        # iterate thorugh kernels and start process
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                print(kernel_name)
                getattr(self, kernel_name)()
            
        # close netCDF files   
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()
        
    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """        
        # add variable to ncdf file
        vn = 'PRESS_pl' # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))    
        var.long_name = 'air_pressure {} pressure levels only'.format(self.NAME)
        var.units     = 'Pa' 
        var.standard_name = 'surface_air_pressure'
        
        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)  
        values  = self.nc_pl.variables['air_pressure'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            #scale from hPa to Pa 
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc, 
                                        time_in*3600, values[:, n]) * 100          

    def AIRT_C_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """        
        # add variable to ncdf file
        vn = 'AIRT_pl' # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))    
        var.long_name = 'air_temperature {} pressure levels only'.format(self.NAME)
        var.units     = 'degrees_C' 
        var.standard_name = 'air_temperature'
        
        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)  
        values  = self.nc_pl.variables['t'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc, 
                                        time_in*3600, values[:, n]-273.15)          

    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """
        # add variable to ncdf file
        vn = 'AIRT_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '2_metre_temperature {} surface only'.format(self.NAME)
        var.units     = 'degrees_C' 
        var.standard_name = 'air_temperature'
        
        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)      
        values  = self.nc_sa.variables['t2m'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc, 
                                                    time_in*3600, 
                                                    values[:, n]-273.15)           
        
    def AIRT_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP.
        """       
        print("AIRT_redcapp")            

    def PREC_mm_sur(self):
        """
        Precipitation sum in mm for the time step given.
        """   
        # add variable to ncdf file
        vn  = 'PREC_sur' # variable name
        var = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Total precipitation {} surface only'.format(self.NAME)
        var.units     = 'kg m-2 s-1'
        var.standard_name = 'precipitation_amount'
        
        # interpolation scale factor
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)
        
        # total prec [mm] in 1 second  
        values  = self.nc_sf.variables['tp'][:]*1000/self.interval_in
        
        # interpolate station by station
        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            f = interp1d(time_in*3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)*self.time_step          
            
    def RH_per_sur(self):
        """
        Relative humdity derived from surface data, exclusively. Clipped to
        range [0.1,99.9]. Kernel AIRT_sur must be run before.
        """         
        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)  
        values  = self.nc_sa.variables['d2m'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            dewp[:, n] = series_interpolate(self.times_out_nc, 
                                            time_in*3600, values[:, n]-273.15) 
                                                    
        # add variable to ncdf file
        vn = 'RH_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Relative humidity {} surface only'.format(self.NAME)
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'
        
        # simple: https://doi.org/10.1175/BAMS-86-2-225
        RH = 100 - 5 * (self.rg.variables['AIRT_sur'][:, :]-dewp[:, :])
        self.rg.variables[vn][:, :] = RH.clip(min=0.1, max=99.9)    

        
    def WIND_sur(self):
        """
        Wind speed and direction temperature derived from surface data, 
        exclusively.
        """    
        # temporary variable, interpolate station by station
        U = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)  
        values  = self.nc_sa.variables['u10'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            U[:, n] = series_interpolate(self.times_out_nc, 
                                         time_in*3600, values[:, n]) 
        
        # temporary variable, interpolate station by station
        V = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)  
        values  = self.nc_sa.variables['v10'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            V[:, n] = series_interpolate(self.times_out_nc, 
                                         time_in*3600, values[:, n]) 

        # wind speed, add variable to ncdf file, convert
        vna = 'WSPD_sur' # variable name
        var           = self.rg.createVariable(vna,'f4',('time', 'station'))    
        var.long_name = '10 wind speed {} surface only'.format(self.NAME)
        var.units     = 'm s-1'  
        var.standard_name = 'wind_speed'
                
        # wind direction, add variable to ncdf file, convert, relative to North 
        vnb = 'WDIR_sur' # variable name
        var           = self.rg.createVariable(vnb,'f4',('time', 'station'))    
        var.long_name = '10 wind direction {} surface only'.format(self.NAME)
        var.units     = 'degree'
        var.standard_name = 'wind_from_direction'

        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            WS = np.sqrt(np.power(V,2) + np.power(U,2))
            WD = [atan2(V[i, n], U[i, n])*(180/pi) + 
                  180 for i in np.arange(V.shape[0])]
            self.rg.variables[vna][:, n] = WS
            self.rg.variables[vnb][:,n] = WD                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    def SW_Wm2_sur(self):
        """
        Short-wave downwelling radiation derived from surface data, exclusively.  
        This kernel only interpolates in time.
        """   
        
        # add variable to ncdf file
        vn = 'SW_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Surface solar radiation downwards {} surface only'.format(self.NAME)
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_shortwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)  
        values  = self.nc_sf.variables['ssrd'][:]/3600/self.interval_in#[w/m2/s]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            f = interp1d(time_in*3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)*self.time_step

    def LW_Wm2_sur(self):
        """
        Long-wave downwelling radiation derived from surface data, exclusively.  
        This kernel only interpolates in time.
        """   
        
        # add variable to ncdf file
        vn = 'LW_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Surface thermal radiation downwards {} surface only'.format(self.NAME)
        var.standard_name = 'surface_downwelling_longwave_flux'
        var.units     = 'W m-2'
        
        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:].astype(np.int64)  
        values  = self.nc_sf.variables['strd'][:]/3600/self.interval_in #[w m-2 s-1]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            f = interp1d(time_in*3600, values[:, n], kind='linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc)*self.time_step   

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg]
        https://crudata.uea.ac.uk/cru/pubs/thesis/2007-willett/2INTRO.pdf
        '''
        
        # temporary variable,  interpolate station by station
        dewp = np.zeros((self.nt, self.nstation), dtype=np.float32)
        time_in = self.nc_sa.variables['time'][:].astype(np.int64)  
        values  = self.nc_sa.variables['d2m'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            dewp[:, n] = series_interpolate(self.times_out_nc, 
                                            time_in*3600, values[:, n]-273.15) 

        # compute
        SH = spec_hum_kgkg(dewp[:, :], 
                           self.rg.variables['PRESS_pl'][:, :])  
        
        # add variable to ncdf file
        vn = 'SH_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Specific humidity {} surface only'.format(self.NAME)
        var.units     = '1'  
        var.standard_name = 'specific_humidity'
        self.rg.variables[vn][:, :] = SH                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

