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
#==============================================================================
#PURPOSE 
#-Step1: Interpolating ERA parameters at one single pixel in 2D demission (X=lon, Y=lat) 

#-Step2: Interpolating ERA parameters at mutiple pixels in 2D demission 

#-Step3: Interpolating ERA paremeters at mutiple pixels in 2D demisison at one specific pressure level

#Looping Conception: 2D(longitude,latitude) x 1D(Pressre Levels)
#==============================================================================
# REFERENCE
#
# -Source Code: 
# https://github.com/geocryology/REDCAPP/redcapp.py (COPYRIGHT: Stephan Gruber & Bin Can)
#
#==============================================================================
#
