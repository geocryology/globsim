#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# (C) Copyright Bin Cao & Stephan Gruber
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# For variable codes and units of ERA-Interim, see: 
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
#==============================================================================
#A example of for the manipulation, especially download of ERA-Interim data 
#saved as ncdf.
#==============================================================================

from datetime import datetime
from os import path

##### HOW TO RUN THIS ##########################################################
#
# (1) Register to ECMWF (free) https://apps.ecmwf.int/registration/
#
# (2) Follow the instructions for "Installing your API key" on
# https://software.ecmwf.int/wiki/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch
#
# (3) Adapt the script below (settings and location) and run it
#
################################################################################

#settings
dir_data = '/Users/xquan/data'
dir_src  = '/Users/xquan/src/globsim'

execfile(path.join(dir_src, 'download.py'))


#location: alps
date  = {'beg' : datetime(2015,1,1),
         'end' : datetime(2015,1,10)}
         
area  = {'north' : 46.65,
         'south' : 46.35,
         'west'  : 9.70, #positive is westwards of Greenwich
         'east'  : 9.95} #positive is westwards of Greenwich
         
elevation = {'min' : 0, 
             'max' : 4500}
#run                                    
ts = toposcale(date, area, elevation, dir_data, 5) 
ts.retrieve()


eraDownload = eraData()
eraDownload.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_sa_*'),1)
eraDownload.NCDFmergeWildcard(path.join(dir_data, 'ecmwf_erai_pl_*'),1)
