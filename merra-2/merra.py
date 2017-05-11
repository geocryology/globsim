#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Copyright Xiaojing Quan & Stephan Gruber
# =============================================================================    
# REVISION HISTORY 
# 20170510 -- Initial version created 
#
#==============================================================================
# A scripts for downloading MERRA-2 reanalysis data:
# --Air Temperature at pressure levels [time*level*lat*lon]
# --Relative Humidity [time*level*lat*lon]
# --Wind Speed/Direction [time*level*lat*lon]
# --Air Temperature at 2 meter [time*level*lat*lon]
# --Surface Solar Radiation Downwards [time*level*lat*lon] (level=1)
# --Surface Thermal Radiation Downwards [time*level*lat*lon] (level=1)
# --Total Precipitation [time*level*lat*lon] (level=1)
# --Others
# Saved as netCDF 4
#====================HOW TO RUN THIS ==========================================
#
# (1) Register a New User in Earthdata Login: 
#  https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login
#
# (2) Authorize NASA GESDISC DATA ARCHIVE in Earthdata Login
# https://disc.gsfc.nasa.gov/registration/authorizing-gesdisc-data-access-in-earthdata_login
#
# (3) Adapt the script below with authrized USERNAME and PASSWROD
# 
# (4) Obtaining the URL address of the objected single dataset at:
#     https://disc.sci.gsfc.nasa.gov/daac-bin/FTPSubset2.pl
# 
# (5) Obtianing the mutiple dataset with spefici spacial and temporal)
#==============================================================================
# IMPORTANT Notes: 
# Selected URLs list:
# url1 = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4'
#       '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')                     # 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics

# url2 = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4'
#       '/2016/02/MERRA2_400.inst3_3d_asm_Np.20160201.nc4')                     # 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
                                                                               
                                                                               
# url3 = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXRAD.5.12.4'
#        '/2016/02/MERRA2_400.tavg1_2d_rad_Nx.20160202.nc4')                    # 2d, 1-Hourly, Time-Averaged, Single-Level, Assimilation, Radiation Diagnostics 
                                                                               

#==============================================================================

from pydap.client import open_url
from pydap.cas.urs import setup_session
from datetime import datetime, timedelta
from os import path

import numpy as np
import csv
import netCDF4 as nc

username = "quanxj17"
password = "Qxj17carleton"

url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4'
       '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')                     # 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics

#== Read in dataset============================================================
session = setup_session(username, password, check_url=url)
ds = open_url(url, session=session)

#===========get variable keys==================================================
print ds.keys

#===========get latitudes,longitude,Level======================================
lat = ds.lat[:]
lon = ds.lon[:]
#lev = ds.lev[:]
#===========find shape of variables============================================
print lat.shape
print lon.shape
#==============================================================================
class MERRAgeneric(object):
    """
    Parent class for other merra classes.
    """
     
    def areaget(self,lat, lon):
        """Gets the specific area with given latitude and longitude"""
        # get the index of lat, lon 
        id_lat = lat.size
        id_lon = lon.size
         
        for i in range(0,(id_lat-1)):
            if lat[i] > 40.0: 
               id_south = i
               print id_south  
            elif lat[i] < 41.0: 
               id_north = i
               print id_north
            else: 
               print "Got the range of latitude" 
        for j in range(0,(id_lon-1)):
            if lon[j] > 60.0:
               id_west = j
               print id_west
            elif lon[j] < 61.0:
               id_east = j 
            else:
               print "Got the range of longitude"
                 
        return(id_south, id_north,id_west,i)
 






