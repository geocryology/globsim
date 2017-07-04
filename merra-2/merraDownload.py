#==============================================================================
# CODE CONCEPTION 
# 1.Commands for reading in urls list of the files with required variables (5 files types mentioned above in total)
# 2.Commands for looping the urls list with specific temporal range (each file by each file)
# 3.Function for defining the spatial subset indices
# 4.Commands for writing and saving defined dataset with specific spatial and temporal ranges, save the dataset as NetCDF file in a local directory.
#
#==============================================================================
#NOTE
#1. how to loop multiple urls to read in all required variables at one time?
#
#==============================================================================

from pydap.client import open_url
from pydap.cas.urs import setup_session
from datetime import datetime, timedelta
from os import path

import numpy as np
import csv
import netCDF4 as nc



#Account for Database Access
username = "quanxj17"
password = "Qxj17carleton"

#settings directory 
dir_data = '/Users/xquan/data'
dir_src  = '/Users/xquan/src/globsim/merra-2'

#Given wanted datatype, mydate, area, elevation

date_start = "2016-01-01"
date_end   = "2016-02-01" 
data_type = '3dm_ana'                                                          # options: '2dm', '3dm_ana','3dm_asm' and '2dr'
area = [ 40.0, 45.0, 60.0, 65.0]
 
elevation = {'min' : 50, 'max' : 2000}  # needed to be test furtherly 

# test the functions in merra.py

#==============================================================================

#Downloading test sample

# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4'
#         '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')                     # 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics

# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4'
#        '/2016/01/MERRA2_400.inst3_3d_asm_Np.20160201.nc4')                     # 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 

# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4'
#        '/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4')                      # 3d,3-hourly,Instantaneous,Pressure-Level, Analyzed Meteorological Fields                                                                                                                                                                                                                                                               
                                                     
url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXFLX.5.12.4'
       '/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4')                      # 2d, 1-hourly, 2D, single-level, full horizontal resolution, Surface Flux Diagnostics

#url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T3NPRAD.5.12.4'
#      '/2016/01/MERRA2_400.tavg3_3d_rad_Np.20160102.nc4')                      # 3d, 3-Hourly, Time-Averaged, Pressure-Level, Assimilation, Radiation Diagnostics


#Read in dataset
session = setup_session(username, password, check_url=url)
ds1 = open_url(url, session=session)

#get variable keys
print ds1.keys

#get latitudes,longitude,Level
# lat = ds.lat[:]
# lon = ds.lon[:]
# lev = ds.lev[:]

# print lat.size
# print lon.size

#get the indices of selected range of area

#execfile(path.join(dir_src, 'merra.py'))

area = [40.0, 45.0, 60.0, 65.0]

#tmp = MERRAgeneric()
#lat_id, lon_id = tmp.areaGET(lat, lon, area)

#lat_id = MERRAgeneric(lat, area)
#lon_id = MERRAgeneric(lon, area)
# id_lat = np.where((lat[:] > area[0]) & (lat[:] < area[1])) 
# id_lon = np.where((lon[:] > area[2]) & (lon[:] < area[3])) 

# print id_lat
# print id_lon

# #find shape of variables
# ds.T2M.shape
# #get variable from dataset
# 
#t2m = ds.T2M[:,:,:]                                                            
# 
# 
# #get subset of variables
# 
#t2m = t2m[:,id_lat,:]
#t2m = t2m[:,:,id_lon]
# 
# print t2m.shape




