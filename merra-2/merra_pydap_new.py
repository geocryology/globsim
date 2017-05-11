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

#=====Account for Database Access==============================================
username = "quanxj17"
password = "Qxj17carleton"

#=====Set up urls==============================================================
#url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/'
#       'M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4')          # test url 

#url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4'
#       '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')                     # 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics

#url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4'
#       '/2016/02/MERRA2_400.inst3_3d_asm_Np.20160201.nc4')                     # 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
                                                                               
                                                                               
#url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXRAD.5.12.4'
#        '/2016/02/MERRA2_400.tavg1_2d_rad_Nx.20160202.nc4')                    # 2d, 1-Hourly, Time-Averaged, Single-Level, Assimilation, Radiation Diagnostics 
                                                                               

url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T3NPRAD.5.12.4'
      '/2016/01/MERRA2_400.tavg3_3d_rad_Np.20160102.nc4')                      # 3d, 3-Hourly, Time-Averaged, Pressure-Level, Assimilation, Radiation Diagnostics


#== Read in dataset============================================================
session = setup_session(username, password, check_url=url)
ds = open_url(url, session=session)

#===========get variable keys==================================================
print ds.keys

#===========get latitudes,longitude,Level======================================
# lat = ds.lat[:]
# lon = ds.lon[:]
# lev = ds.lev[:]


#===========find shape of variables============================================
# ds.T2M.shape
# ds.V10M.shape
# ds.U10M.shape
# ds.T.shape
# ds.RH.shape
# ds.V.shape
# ds.U.shape
# ds.PHIS.shape
# ds.LWGAB.shape
# ds.LWGNT.shape
# ds.SWGDN.shape
ds.DTDTSWR.shape
ds.DTDTLWR.shape

#===========get subset of one variable=========================================
# need to set up the spatial subset indices
#
# t2m = ds.T2M[:,0:2,0:2]                                                      # 2m Surface Temperature (Unit:K) [time*lat*lon]
# v10m = ds.V10M[:,:,:]                                                        # 10-meter_northward_wind 
# u10m = ds.U10M[:,:,:]                                                        # 10-meter_eastward_wind      
# temp = ds.T[:,:,0:2,0:2]                                                     # Air Temperatuere at Pressure Levels (42) (Unit:K) [time*lev*lat*lon]
# rh  = ds.RH[1,1,0:2,0:2]                                                     # Relative Humity at Pressure Levels (42) (Unit:1) [time*lev*lat*lon]
# v   = ds.V[1,1,0:2,0:2]                                                      # Northward_wind (Unit:m-s) [time*lev*lat*lon]
# u   = ds.U[1,1,0:2,0:2]                                                      # Eastward_wind (Unit:m-s) [time*lev*lat*lon]
# gp  = ds.PHIS[:,:,:]                                                         # Surface geopotential height (units: m+2 s-2) [time*lat*lon]
#lwgab = ds.LWGAB[:,:,:]                                                       # Surface absorbed longwave radiation (unit:W m-2) [time*lat*lon] 
#lwgnt = ds.LWGNT[:,:,:]                                                       # Surface net downward longwave flux (units: W m-2) [time*lat*lon]
#swgdn = ds.SWGDN[:,:,:]                                                       # Surface incoming shortwave flux (units: W m-2) [time*lat*lon]
dtdtswr = ds.DTDTSWR[:,:,:,:]                                                  # Air temperature tendency due to shortwave (units:Ks-1) [time*lev*lat*lon]
dtdtlwr = ds.DTDTLWR[:,:,:,:]                                                  # Air temperature tendency due to longwave (units:Ks-1 ) [time*lev*lat*lon]


#==================write out abstracted variables==============================





#==================Save output as netCDf file==================================

#f=NetCDFFile('Sample_Result_1.nc', 'w', format='NETCDF4')                     # creat a netCDF file

#tempgrp=f.creatGroup('Temp_data')                                             # creat a data group for temperature result output

#tempgrp.creatDimension('time',None)   
#tempgrp.creatDimension('lat', len(lat))  
#tempgrp.creatDimension ('lon',len(lon))                                       # Specify the dimension of data                   
      
                       

#longitude=tempgrp.creatVariable('Longitude', 'f4', 'lon')                     # Building output variables    f4:32 bit float
#latitude=tempgrp.creatVariable('Latitude','f4', 'lat')                        #                              i4: 32 bit integer
#levels=tempgrp.creatVariable('Temperature', 'f4', ('time', 'lon', 'lat'))
#time= tempgrp.creatVariable('Time','i4', 'time')

#longitude[:]= lon                                                             # Pass the values of interpolated results to the output variables
#latitude[:]= lat
#temp[:,:,:]= temp

#longitude.units = 'degrees east'                                              # Add local attributes to variable instances
#latitude.units = 'degrees north'
#time.units = 'days since Jan 01, 0001'
#temp.units = 'Kelvin'

#f.close()                                                                     # Close the dataset                                                  




#==============================================================================