#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 08:21:28 2017

@author: xquan
"""

#==============================================================================
# An example for downloading and manipulating ERA-Interim data
# Location: OTTAWA AREA (Latitude:45.4, Longitude:-75.7)
# Saved as ncdf
#==============================================================================

from datetime import datetime
from os import path

##### HOW TO RUN THIS #########################################################
#
# (1) Register to ECMWF (free) https://apps.ecmwf.int/registration/
#
# (2) Follow the instructions for "Installing your API key" on
# https://software.ecmwf.int/wiki/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch
#
# (3) Adapt the script below (settings and location) and run it
#
###############################################################################

import ERA_Interim

# set parameters
date      = {'beg' : datetime(2016, 1, 1), 'end' : datetime(2016, 2, 1)}
area      = {'north':  40.0, 'south': 50.0, 'west': -70.0, 'east': -80.0}
elevation = {'min': 0, 'max': 1000}           
directory = '/home/xquan/data'             

# run
ts = ERAbatch(date, area, elevation, directory, 15) 
ts.retrieve()
ts.deleteGrib()
