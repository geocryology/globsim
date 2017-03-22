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
# --Suface Solar Radiation Downwards [time*level*lat*lon] (level=1)
# --Suface Thermal radiation Downwards [time*level*lat*lon] (level=1)
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

from datetime import datetime, timedelta
from netCDF4  import Dataset as NetCDFFile 

from os import path
from scipy.ndimage import gaussian_filter,generic_filter,convolve,minimum_filter,maximum_filter

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

# Reading in variables


nc=NetCDFFile('/Users/xquan/data/era_pl_20160101_to_20160105.nc') #open the file

# Read the variables from the netCDF file and assign them to Python variables

lat=nc.variables['lat'][:]
lon=nc.variables['lon'][:]
time=nc.variables['time'][:]
mslp=nc.variables['level'][:]
rh=nc.variables['Relative humidity'][:]
temp=nc.variables['Temperature'][:]
u=nc.variables['U component of wind'][:]
v=nc.variables['V component of wind'][:]




#==========================Step 1==============================================

dir_data= '/Users/xquan/data'
dir_src= '/Users/xquan/src/globsim'

execfile(path.join(dir_src, 'redcapp_XQ.py'))

#-----Option 1(Utilizing Classes from redcapp)---------------------------------Q: Option 1 or 2, which one is a right way to conduct the 2D interpolation?

#get the raw data from directory containing all raw data and output data
dataImport=rawData(dir_data)
sa=dataImport.saf_get() #get da file in the given directory 
pl=dataImport.plf_get() #get pl file in the given directory
geop=dataImport.geopf_get()  # geopotential file


# 2D interpolation 

dem='example_Ottawa.nc'                                                        #Q: what is the conception to generate the file? Input or output variable?
geop ='era_to.nc'                                                              #Q: geop: era_to.nc? if it wasn't, how to generate this data file?
a='era_sa_20160101_to_20160105.nc'                                          
pl='era_pl_20160101_to_20160105.nc'

downscaling=downscaling(dem,geop,sa,pl)

out_xyz_dem, lats, lons, shape= downscaling.demGrid()
out_xyz_sur = downscaling.surGrid(lats, lons, None)

# interpolate 2-meter temperature
surTa=downscaling.surTa(0, out_xyz_sur)

# original ERA-I values
gridT, gridZ, gridLat, gridLon= downscaling. gridValue(variable,0)

#interpolate temperatures and geopotential of different pressue levels
t_interp, z_interp=downscaling.inLevelInterp(gridT, gridZ, gridLat, gridLon, out_xyz_dem) 


#---Option 2(interploting variables at one single pixel at given time index----


time_idx=1



















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


