#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Copyright Xiaojing Quan & Stephan Gruber
#==============================================================================
# a script for checking the continuity of Time with given the chunk for downloaded netCDF files from ERA_Interim and MERRA-2,
# 
# -- Step 1: read in downloaded netCDF files with a wildcard expression;
# -- Step 2: check the time coverage of varialbes to be continous; 
# -- Step 3: if there was any missing one, to list the missing series of time 
#
#===============================================================================
from datetime          import datetime, timedelta, date
from os                import path
from netCDF4           import Dataset, MFDataset
from generic           import ParameterIO, StationListRead
from __future__        import print_function

import numpy as np
import csv
import netCDF4 as nc
import math
import time as tc
import sys


class MERRA2DataCheck(object):
    """
        To check  the continuity of the time coverage of from all the downloaded 
        netCDF files with a wildward expression
          
        Args:
        ncfile_in: Full path to an MERRA-2 derived netCDF file. This can
        contain wildcards to point to multiple files if temporal
      
    """

    def __init__(self, ifile):
        #read parameter file
        self.ifile = ifile
        par = ParameterIO(self.ifile)
        self.dir_inp = path.join(par.project_directory,'merra2') 
        self.dir_out = path.join(par.project_directory,'station')
        self.variables = par.variables
        self.list_name = par.list_name
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        
        #read station points 
        self.stations = StationListRead(self.stations_csv)  
        
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}
                      
                      
    def DataReadin(self, file_in):    
        """
        To read in all the downloaded netCDF files with a wildward expression
          
        Args:
            ncfile_in: Full path to an MERRA-2 derived netCDF file. This can
                        contain wildcards to point to multiple files if temporal
                        chunking was used.
              
            ncfile_out: Full path to the output netCDF file to write.  
          """

        # open netcdf file handle, can be one file of several with wildcards
        ncf = nc.MFDataset(file_in, 'r', aggdim ='time')
                                                                 
        #get variable 'time' 
        nctime = ncf.variables['time'][:]
        t_unit = ncf.variables['time'].units 
        try :
            t_cal = ncf.variables['time'].calendar
        except AttributeError : # Attribute doesn't exist
            t_cal = u"gregorian" # or standard
 
        # check the continourity of the time coverage
        if file_in == path.join(self.dir_inp,'merra_pl_*.nc'):
             missingIndex = list(set(xrange(nctime[0], nctime[-1] + 1, 6)) - set(nctime))       
             if missingIndex == []:
                print("NO MISSING TIME FOR", file_in)
             else:    
                 missingIndex.sort()
                 missingTime = nc.num2date(missingIndex, units = t_unit, calendar = t_cal)
                 
                 print('FOR', file_in)
                 
                 print("Number of Missing Indices of Time : ", len(missingIndex))
                 
                 print("Unit of Time: ", t_unit)
                 
                 print('Missing Time List:', list(missingTime))
        
            
        else:
             missingIndex = list(set(xrange(nctime[0], nctime[-1] + 1)) - set(nctime))      
             if missingIndex == []:
                print("NO MISSING TIME FOR", file_in)
             else:    
                 missingIndex.sort()
                 missingTime = nc.num2date(missingIndex, units = t_unit, calendar = t_cal)

                 print('FOR', file_in)

                 print("Number of Missing Indices of Time : ", len(missingIndex))
                 
                 print("Unit of Time: ", t_unit)

                 print('Missing Time List:', list(missingTime)              )

    def process(self):
        """
        combine the specific type of single or mutiple netCDF files
        """                       

        # === Surface  Data ===    
        #                
        self.DataReadin(path.join(self.dir_inp,'merra_sa_*.nc'))          

        # === Pressure-level Data ===    
        # 
        #                
        self.DataReadin(path.join(self.dir_inp,'merra_pl_*.nc'))          

        # === radiation Data ===    
        # 
        #                
        self.DataReadin(path.join(self.dir_inp,'merra_sr_*.nc'))          


class ERADataCheck(object):
    """
        To check  the continourity of the time coverage of from all the downloaded 
        netCDF files with a wildward expression
          
        Args:
        ncfile_in: Full path to an ERA_Interim derived netCDF file. This can
        contain wildcards to point to multiple files if temporal
      
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
        
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}
                      
                      
    def DataReadin(self, file_in):    
        """
        To read in all the downloaded netCDF files with a wildward expression
          
        Args:
            ncfile_in: Full path to an MERRA-2 derived netCDF file. This can
                        contain wildcards to point to multiple files if temporal
                        chunking was used.
              
            ncfile_out: Full path to the output netCDF file to write.  
          """
	  
	# open netcdf file handle, can be one file of several with wildcards
	ncf = nc.MFDataset(file_in, 'r', aggdim ='time') 
                  
        #get variable 'time' 
        nctime = ncf.variables['time'][:]
        t_unit = ncf.variables['time'].units 
        try :
            t_cal = ncf.variables['time'].calendar
        except AttributeError : # Attribute doesn't exist
            t_cal = u"gregorian" # or standard
 
        # check the continourity of the time coverage
        if file_in == path.join(self.dir_inp,'era_sf_*.nc'):
             missingIndex = list(set(xrange(nctime[0], nctime[-1] + 1, 3)) - set(nctime))       
             if missingIndex == []:
                print("NO MISSING TIME FOR", file_in)
             else:    
                 missingIndex.sort()
                 missingTime = nc.num2date(missingIndex, units = t_unit, calendar = t_cal)
                 
                 print('FOR', file_in)
                 
                 print("Number of Missing Indices of Time : ", len(missingIndex))
                 
                 print("Unit of Time: ", t_unit)
                 
                 print('Missing Time List:', list(missingTime)              )
        
            
        else:
             missingIndex = list(set(xrange(nctime[0], nctime[-1] + 1, 6)) - set(nctime))      
             if missingIndex == []:
                print("NO MISSING TIME FOR", file_in)
             else:    
                 missingIndex.sort()
                 missingTime = nc.num2date(missingIndex, units = t_unit, calendar = t_cal)

                 print('FOR', file_in)

                 print("Number of Missing Indices of Time : ", len(missingIndex))
                 
                 print("Unit of Time: ", t_unit)

                 print('Missing Time List:', list(missingTime)              )

    def process(self):
        """
        combine the specific type of single or mutiple netCDF files
        """                       

        # === Surface  Data ===    
        #                
        self.DataReadin(path.join(self.dir_inp,'era_sa_*.nc'))          

        # === Pressure-level Data ===    
        # 
        #                
        self.DataReadin(path.join(self.dir_inp,'era_pl_*.nc'))          

        # === radiation Data ===    
        # 
        #                
        self.DataReadin(path.join(self.dir_inp,'era_sf_*.nc'))          

                
# ===========For Run DataCheck ================================================ 
# Args:
# ifile : Full path to a Globsim Interpolate Paramter file under operating latap or virtual machine
  
ifile = '/home/xquan/src/globsim/examples/par/examples.globsim_interpolate'

# to check merra2 data 
MERRA2check = MERRA2DataCheck(ifile)

MERRA2check.process()

# # to check eraint data
# ERAcheck = ERADataCheck(ifile)
# 
# ERAcheck.process()

