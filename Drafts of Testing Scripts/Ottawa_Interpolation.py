#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Copyright Xiaojing Quan & Stephan Gruber
# =============================================================================    
# REVISION HISTORY 
# 20170316 -- Initial version created 
#
#==============================================================================
# An example for interpolating ERA_Interim data:
# --Air Temperature at pressure levels [time*level*lat*lon]
# --Relative Humidity [time*level*lat*lon]
# --Wind Speed/Direction [time*level*lat*lon]
# --Air Temperature at 2 meter [time*level*lat*lon]
# --Surface Solar Radiation Downwards [time*level*lat*lon] (level=1)
# --Surface Thermal Radiation Downwards [time*level*lat*lon] (level=1)
# --Total Precipitation [time*level*lat*lon] (level=1)
# --Others
#
# Files Format: NetCDF
#
#==============================================================================
# REFERENCE
#
# -Source Code: 
# https://github.com/geocryology/REDCAPP/redcapp.py (COPYRIGHT: Stephan Gruber & Bin Can)
#
#
#==============================================================================
#PURPOSE 
#-Step1: Interpolating ERA parameters at one single pixel in 2D demission (X=lon, Y=lat) 

#-Step2: Interpolating ERA parameters at mutiple pixels in 2D demission 

#-Step3: Interpolating ERA paremeters at mutiple pixels in 2D demisison at one specific pressure level

#==============================================================================
#NOTES
# Looping Conception: 2D(longitude,latitude) x 1D(Pressre Levels)
# Saved as netcdf file (by time seris)

#==============================================================================

from datetime import datetime # Python standard library datetime module 
from os       import path, remove
from scipy.interpolate import griddata

import netCDF4 as nc
import datetime as dt
import numpy as np


dir_data= '/Users/xquan/data'
dir_src= '/Users/xquan/src/globsim'

execfile(path.join(dir_src, 'globsim.py'))


# Get in ERA_Interim raw data at pressure levels

dem  = 'User/xquan/data/example_Ottawa.nc'
pl = '/Users/xquan/data/era_pl_20160101_to_20160105.nc'

num  = nc.Dataset('/Users/xquan/data/era_pl_20160101_to_20160105.nc')

Interp2d = Interp2d(dem, pl)     

# Read in interpoalted variable at specific time and pressue level indexs

ind_time=0
ind_lev=0


temp=num.variables['Temperature'][ind_time,ind_lev,:,:]


#temp=Interp2d.gridVariable('Temperature',ind_time, ind_lev)
#print temp

rh=Interp2d.gridVariable('Relative humidity',ind_time, ind_lev)
u=Interp2d.gridVariable('U component of wind',ind_time, ind_lev)
v=Interp2d.gridVariable('V component of wind',ind_time, ind_lev)
gp=Interp2d.gridVariable('Geopotential',ind_time, ind_lev)





#==========================Step 1==============================================
#Interpolated Tempreture at given pixels in one given time index and pressure level

lats= [45.35, 45.40, 45.45]
lons= [284.25,284.30,284.35]       # convert from: +360


interpTemp=Interp2d.interVariable(temp,lats,lons)



#-----------------------Step 2-------------------------------------------------





#-----------------------Step 3-------------------------------------------------





#-------------------------------------------------------------------------------

# Saved result as netCDF file

#f=NetCDFFile('Sample_Result_1.nc', 'w', format='NETCDF4')    # creat a netCDF file

#tempgrp=f.creatGroup('Temp_data')                            # creat a data group for temperature result output

#tempgrp.creatDimension ('lon',len(lon))                      # Specify the dimension of data                   
#tempgrp.creatDimension('lat', len(lat))        
#tempgrp.creatDimension('time',None)                          

#longitude=tempgrp.creatVariable('Longitude', 'f4', 'lon')    # Building output variables    f4:32 bit float
#latitude=tempgrp.creatVariable('Latitude','f4', 'lat')       #                              i4: 32 bit integer
#levels=tempgrp.creatVariable('Temperature', 'f4', ('time', 'lon', 'lat'))
#time= tempgrp.creatVariable('Time','i4', 'time')

#longitude[:]= lon                                            # Pass the values of interpolated results to the output variables
#latitude[:]= lat
#temp[:,:,:]= temp_data

#longitude.units = 'degrees east'                             # Add local attributes to variable instances
#latitude.units = 'degrees north'
#time.units = 'days since Jan 01, 0001'
#temp.units = 'Kelvin'

#f.close()                                                    # Close the dataset                                                  

