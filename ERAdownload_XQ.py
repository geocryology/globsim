#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# (C) Copyright Xiaojing Quan & Stephan Gruber
#
# For variable codes and units of ERA-Interim, see: 
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
#==============================================================================
# An example for downloading and manipulating ERA-Interim data at one single location
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

#settings
dir_data = '/Users/xquan/data'
dir_src  = '/Users/xquan/src/globsim'

execfile(path.join(dir_src, 'ERA_Interim_Download.py'))

# Location: OTTAWA AREA (Latitude:45.4, Longitude:-75.7)
date  = {'beg' : datetime(2016,1,1),
         'end' : datetime(2016,2,1)}
         
area  = {'north' : 40.00,
         'south' : 50.00,
         'west'  : -70.00, 
         'east'  : -80.00}
         
elevation = {'min' : 0, 
             'max' : 2000}
#run                                    
ts = ERAbatch(date, area, elevation, dir_data, 5)  
ts.retrieve()

ERAdownload_XQ = eraData()
ERAdownload_XQ.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_sa_*'),1)
ERAdownload_XQ.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_pl_*'),1)
ERAdownload_XQ.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_sf_*'),1)



