# -*- coding: utf-8 -*-
# 
# Globsim
#
# === COPYRIGHT AND LICENCE ====================================================
#
# Copyright 2017-2018 Stephan Gruber & Xiaojing Quan
#           
#
# 
#
# === CONTRIBUTIONS ============================================================
#
#
#
# === NOTES ====================================================================
#
# Get globsim at https://github.com/geocryology/globsim
# Check ECMWF job status: http://apps.ecmwf.int/webmars/joblist/
# TODO -multiprocessing
# TODO -Downloading climate reanalysis data (ERA_Interim, MERRA, etc.)
# TODO -interpolate per ERA-gridcell nearest to center of points given
#      -stratify list of points to do that  
# TODO -Plotting and presenting the interpolation results
#
# ==============================================================================

from datetime import datetime, timedelta
from ecmwfapi import ECMWFDataServer
from math     import exp, floor, radians
from os       import path, remove
from scipy.interpolate import griddata,RegularGridInterpolator, NearestNDInterpolator,LinearNDInterpolator
from scipy.ndimage import gaussian_filter,generic_filter,convolve,minimum_filter,maximum_filter
from itertools import izip
from bisect import bisect_left
from netCDF4  import Dataset as NetCDFFile 

import numpy   as np
import pygrib  as pg
import netCDF4 as nc
import glob    as gl
import csv
import re



class Interp2d(object):
    """
    This is a 2D interpolatation, and returns interpolated variables in 
    specific time and pressure level at required single and multipil locations 
        
    Args:
        sa: ERA_Interim raw metadata in surface level   
        pl: ERA_Interim raw metadata in pressure levels 
        dem: A required fine-scale dem in netcdf format
            
    Returns:
        Interpolated variables in the format of [lats, lons, values] at gives locations
        
    Example:
            
        dem  = 'example_Ottawa.nc'
        sa   = 'era_sa_20160101_to_20160105.nc'
        pl   = 'era_pl_20160101_to_20160105.nc'
        
        Interp2d = Interp2d(dem, pl)
        
    """             
    
    def __init__(self, sa, pl, dem =None):      
#        self.sa     = NetCDFFile(sa)
        self.pl     = NetCDFFile(pl)
        self.g       = 9.80665 # Gravitational acceleration [m/s2]
        self.absZero = 273.15  
        if not (dem is None):
            self.dem = NetCDFFile(dem) 
     
    def gridVariable(self,variable, ind_time, ind_lev):
        """
        Return original grid variables and geopotential of differnet
        pressure levels. 
        
        Args: 
            variable: Given interpolated climate variable
            ind_time: Time need to be interpolated. The format is in interger 
            (e.g. 0, 1, 2,...)
            ind_lev: Pressure level need to be interpolated.The format is in 
            interger (e.g. 1, 2, ...)
            
        Returns: 
          
           Temp: Grid temperatures of at one specific pressure level. Retruned 
             are formated in [lat, lon]
            lat: Grid latitude of pressure level variables
            lon: Grid longitude of pressure level variables
        
        Example:
            
            dem  = 'example_Ottawa.nc'
            geop = 'era_to.nc'
            sa   = '/Users/xquan/data/era_sa_20160101_to_20160105.nc'
            pl   = '/Users/xquan/data/era_pl_20160101_to_20160105.nc'
           
            Interp2d = Interp2d(dem, pl)  
            
            ind_time=1
            ind_lev=1
            
            Temp = Interp2d.gridVariable('Temperature',ind_time,ind_lev)
            
        """
        
      #  ind_time=0
      #  ind_lev=0
        
        va= self.pl.variables[variable][ind_time,ind_lev,:,:]
        
        lat=self.pl.variables['lat'][:]
        lon=self.pl.variables['lon'][:]
        lev=self.pl.variables['level'][:]
        
        return va, lat, lon, lev
        

    def interVariable(self, variable, lats, lons):
        """
        Return interpolated variable in one specific pressure level at given 
        pixel
        
        Args: 
            variable: Given interpolated climate variable
            lats: The latitude at given locations
            lons: The longitude at given locations
        
        Returns: 
            interp_va: interpolated variable at given locations in the format of
            [lats, lons, values]           
        
        Example:
            
            dem  = 'example_Ottawa.nc'
            geop = 'era_to.nc'
            sa   = '/Users/xquan/data/era_sa_20160101_to_20160105.nc'
            pl   = '/Users/xquan/data/era_pl_20160101_to_20160105.nc'
           
            Interp2d = Interp2d(dem, pl)  
            
            ind_time=1
            ind_lev=1
            
            temp = Interp2d.gridVariable('Temperature',ind_time, ind_lev)
            
            lats= [45.4, 45.45]
            lons= [284.30,284.35] # convert from: +360
            
            intep_temp=Interp2d.interVariable(temp,lats,lons)
          
            
        """

        lat=self.pl.variables['lat'][:]
        lon=self.pl.variables['lon'][:]
     
        f_va=RegularGridInterpolator((lat, lon), variable, 'linear')        
        out_xy = np.array([lats, lons]).T
        interp_va = f_va(out_xy)
        out_variable= np.array([lats,lons,interp_va]).T


        print out_variable
     
        
        return out_variable