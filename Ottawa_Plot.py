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
from netCDF4  import Dataset as NetCDFFile # http://code.google.com/p/netcdf4-python/
from os       import path, remove
from mpl_toolkits.basemap import Basemap,shiftgrid,shiftgrid,addcyclic
from scipy.interpolate import griddata

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

# Reading in variables

nc=NetCDFFile('/Users/xquan/data/era_pl_20160101_to_20160105.nc') #open the file

# Read the variables from the netCDF file and assign them to Python variables

lat=nc.variables['lat'][:]
lon=nc.variables['lon'][:]
time=nc.variables['time'][:]
mslp=nc.variables['level'][:]
RH=nc.variables['Relative humidity'][:]
t=nc.variables['Temperature'][:]
u=nc.variables['U component of wind'][:]
v=nc.variables['V component of wind'][:]

time_idx=20
#Fix the problem of Python and the renalaysis are slightly off in time
#offset=dt.timedelta(hour=48)
#List of all times in the file as datetime objects
#dt_time = [dt.date(1, 1, 1) + dt.timedelta(hours=t) - offset\
#           for t in time]
#cur_time = dt_time[time_idx]

#Plot of Georeferenced Temperature 
#Setup the map. (See http://matplotlib.org/basemap/users/mapsetup.html)
map=Basemap(projection='moll', llcrnrlat=-90, urcrnrlat=90,\
            llcrnrlon=0, urcrnrlon=360, resolution='c', lon_0=0)# projection, lat/lon extents and resolution of polygons to draw
                                                               # resolutions: c - crude, l - low, i - intermediate, h - high, f - full
#Others
map.drawcoastlines()
map.drawmapboundary()
map.drawcountries()

# Make the plot continuous
#air_cyclic, lons_cyclic = addcyclic(t[time_idx, :, :], lon)

# Shift the grid so lons go from -180 to 180 instead of 0 to 360.
#air_cyclic, lons_cyclic = shiftgrid(180., air_cyclic, lons_cyclic, start=False)

#Make the grid so lons go from -180 to 180 instead of 0 to 360.
#air_cyclic, lons_cyclic=shiftgrid(180., air_cyclic, lons_cyclic, start=False)

#Creat 2D lat/lon/ arrays for Basemap
lons, lats=np.meshgrid(lon-180, lat)

# Transofrm lat/lon into plotting coordinates for projection
x,y=map(lons, lats)

# Plot of air temperature with 11 contour intervals

#cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
#cbar.set_label("%s (%s)" % (nc.variables['Temperature'].var_desc,\
#                            nc.variables['Temperature'].units))
#plt.title("%s on %s" % (nc.variables['Temeprature'].var_desc, cur_time))

temp=map.contourf(x,y,t[1,1,:,:])
cb = map.colorbar(temp,"bottom", size="5%", pad="2%")
plt.title('Temperature')
cb.set_label('Temperature (K)')



#Save figure
plt.show()
plt.savefig('temp.png')






