#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Copyright Xiaojing Quan & Stephan Gruber
# =============================================================================    
# REVISION HISTORY 
# 20170309 -- Initial version created 
#
#==============================================================================
# An example for processing and analysing ERA_Interim data:
# --Air Temperature at pressure levels 
# --Relative Humidity 
# --Wind Speed/Direction 
# --Air Temperature at 2 meter 
# --Suface Solar Radiation Downwards 
# --Suface Thermal radiation Downwards 
# --Total Precipitation
# --Others
#
# Files Format: NetCDF
#==============================================================================
#PURPOSE 
#To demonstrate how to read and write data with NetCDF files using a NetCDF file 
# from ERA_Interim(ECMWF)
# Plotting using Matplotlib and Basemap
#==============================================================================
# REFERENCE
#
#==============================================================================

from datetime import datetime # Python standard library datetime module 
from netCDF4  import Dataset # http://code.google.com/p/netcdf4-python/
from os       import path, remove
from mpl_toolkits.basemap import Basemap,shiftgrid,shiftgrid
from scipy.interpolate import griddata

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

# Reading in variables

#nc_f='User/xquan/data/era_pl_20160101_to_20160105.nc' #file name 
#nc_fid=Dataset(nc_f,'r') # open the file

#root_grp= Dataset ('User/xquan/data/era_pl_20160101_to_20160105.nc')






#Creat the basemap object






        


#Populate a gridded array for plotting 








# Creat and customize plot






