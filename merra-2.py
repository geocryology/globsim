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

# Samples of Selected URLs list:

# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4'
#        '/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4')                     # 3d,6-hourly,Instantaneous,Pressure-Level, Analyzed Meteorological Fields                                                                                                                                                                                                                                                               

# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4'
#        '/2016/01/MERRA2_400.inst3_3d_asm_Np.20160201.nc4')                     # 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
                                                     
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4'
#         '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')                    # 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics,meteorological Diagnostics data 

# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXFLX.5.12.4'
#       '/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4')                      # 2d,1-hourly, single-level, full horizontal resolution, Surface Flux Diagnostics

# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T3NPRAD.5.12.4'
#      '/2016/01/MERRA2_400.tavg3_3d_rad_Np.20160102.nc4')                       # 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics

#==============================================================================

from pydap.client import open_url
from pydap.cas.urs import setup_session
from datetime import datetime, timedelta
from os import path
from netCDF4 import Dataset

import pydap.lib
import numpy as np
import csv
import netCDF4
import itertools
import pandas



class MERRAgeneric():
    """
    Parent class for other merra classes.
    """
           
    def getURLs(self, beg, end):                                                                                                                                          
        """ Set up urls by given range of date and type of data
            to getobjected url address (2d, 3d meterological fields and radiation datasets)
            url_2dm: 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
            url_3dm_ana: 3d,6-hourly,Instantaneous,Pressure-Level,Analyzed Meteorological Fields 
            url_3dm_asm: 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
            url_2dr: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
            
            Args:
            beg = "2016/01/01"
            end = "2016/12/31"

            urls_3dmana: get type of url for 3d Analyzed Meteorological Fields data
            urls_3dmasm: get type of url for 3d Assimilated Meteorological Fields data
            urls_2dm:    get type of url for 2d meterological data
            urls_2ds:    get type of url for 2d surface flux data
            urls_2dr:    get type of url fro 2d radiation diagnostics data
            
        """
        #setup the based url strings    
        baseurl_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/') # baseurl for 2d dataset
        baseurl_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/') # baseurl for 3d dataset
 
        baseurl_3dn = ('M2I6NPANA.5.12.4/','/MERRA2_400.inst6_3d_ana_Np.')            # sub url of 3d Analyzed Meteorological Fields data           
        baseurl_3da = ('M2I3NPASM.5.12.4/','/MERRA2_400.inst3_3d_asm_Np.')            # sub url of 3d Assimilated Meteorological Fields data
        baseurl_2dm = ('M2I1NXASM.5.12.4/','/MERRA2_400.inst1_2d_asm_Nx.')            # sub url of 2d meteorological Diagnostics data     
        baseurl_2dr = ('M2T1NXRAD.5.12.4/','/MERRA2_400.tavg1_2d_rad_Nx.')            # sub url of 2d radiation Diagnostics data
        baseurl_2ds = ('M2T1NXFLX.5.12.4/','/MERRA2_400.tavg1_2d_flx_Nx.')            # sub url of 2d suface flux Diagnostics data                                                              
        format = ('.nc4')

        #Setup the start and end of dates
        begin = datetime.strptime(beg, '%Y/%m/%d')
        end  = datetime.strptime(end, "%Y/%m/%d")

        #Setup the based string of dates for urls 
        res1 = [d.strftime("%Y/%m") for d in pandas.date_range(begin,end)]
        res2 = [d.strftime("%Y%m%d") for d in pandas.date_range(begin,end)]        
        
       # get the urls list
       
        urls_3dmana = []
        urls_3dmasm = []
        urls_2dm = []
        urls_2dr = []
        urls_2ds = []
        for i in range(0,len(res1)): 
                 urls_3dmana.append(baseurl_3d + baseurl_3dn[0] + res1[i] + baseurl_3dn[1] + res2[i] + format)  # urls of 3d Analyzed Meteorological Fields datasets with temporal subset             
                 urls_3dmasm.append(baseurl_3d + baseurl_3da[0] + res1[i] + baseurl_3da[1] + res2[i] + format)  # urls of 3d Assimilated Meteorological Fields datasets with temporal subset                  
                 urls_2dm.append(baseurl_2d + baseurl_2dm[0] + res1[i] + baseurl_2dm[1] + res2[i] + format)     # urls of 2d meteorological Diagnostics datasets temporal subset 
                 urls_2ds.append(baseurl_2d + baseurl_2ds[0] + res1[i] + baseurl_2ds[1] + res2[i] + format)     # urls of 2d suface flux Diagnostics datasets with temporal subset  
                 urls_2dr.append(baseurl_2d + baseurl_2dr[0] + res1[i] + baseurl_2dr[1] + res2[i] + format)     # urls of 2d radiation Diagnostics datasets with temporal subset 
 
        return urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr
 
    def download(self, username, password, urls, chunk_size):
        """ Access the MERRA server by account information and defiend urls
            Args:
            username = "xxxxxx"
            password = "xxxxxx"
            urls = urls_3dmana,urls_3dmasm,urls_2dm,urls_2ds, urls_2dr: a full list of urls by specific date range of wanted types of dataset
            size: the wanted size of urls list for each chunk
                              
        """
        chunk_size = 5
        urls_chunks = [urls[x:x+chunk_size] for x in xrange(0, len(urls), chunk_size)]      

        print ('===== MERRA: START ======')
        print ('TIME TO GET A COFFEE')        
        ds = {}
        for i in range(len(urls_chunks)):
            ds[i] = {}
            ####ds[0] = {}
            url = urls_chunks[i]
            for j in range(len(url)): 
                session = setup_session(username, password, check_url=url[j])        
                ds[i][j] = open_url(url[j], session=session) 
                ###ds[0][j] = open_url(url[j], session=session)  #chris
                print ('------COMPLETED------','CHUNK', i, 'URL', j )
                print url[j]
            print ds[i][j].keys
            ###return ds #Chris     
        print ('===== MERRA: COMPLETED =======')
        print ('Type:', type(ds), 'Length:', len(urls))      
        
        return ds        
    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 datasets        
           variable = ['T','V','U','H','lat','lon','lev','time']                         # for extrating from 3d Analyzed Meteorological Fields datasets 
                    = ['RH','lat','lon','lev','time']                                    # for extracting from 3d Assimilated Meteorological Fields datasets
                    = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']  # for extracting from 2d meteorological Diagnostics datasets
                    = ['PRECTOT','lat','lon','time']                                     # for extracing from 2d suface flux Diagnostics datasets  
                    = ['SWGNT','LWGNT','lat','lon','time']                               # for extracing from 2d radiation Diagnostics datasets 
                    
           ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, chunk_size)
        """
       
        out_variable = {}
        for i in range(0, len(ds)):
            out_variable[i] = {}
            for j in range(len(ds[i])):
                print "run", i, j
                outputVar = []
                for x in range(0,len(variable)):
                    outputVar.append(variable[x])

                var = ds[i][j].keys()
                for l in range(len(outputVar)):
                    foundVariable = False
                    if outputVar[l] in var:
                        for k in range(len(var)):
                            if foundVariable != True:
                                if var[k] == outputVar[l]:
                                    temp = "" + var[k]
                                    outputVar[l] = ds[i][j][temp]
                                    foundVariable = True
                out_variable[i][j] = outputVar

        print "Length of Out_Variable", len(out_variable[0][0])
        
        return out_variable
    
    def getArea(self, area, ds): 
        """Gets the specific area with given latitude and longitude
           For example: 
           area = {'bbS':60.0, 'bbN': 65.0, 'bbW':-115.0, 'bbE': -110.0}
           ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm, urls_2dm, urls_2dr, urls_2ds, chunk_size)
        """
             
        # pass the value of individual row lat and lon to Lat and Lon for the area subset
        Lat = ds[0][0].lat[:]
        Lon = ds[0][0].lon[:]
                        
        # get the indices of selected range of Lat,Lon
        id_lon = np.where((Lon[:] > area['bbW']) & (Lon[:] < area['bbE'])) 
        id_lat = np.where((Lat[:] > area['bbS']) & (Lat[:] < area['bbN'])) 
       
        # convert id_lat, id_lon from tuples to string
        id_lat = list(itertools.chain(*id_lat))
        id_lon = list(itertools.chain(*id_lon))   
        
        return id_lat, id_lon

    def latLon_3d(self, out_variable, p1, p2, p3, p4, id_lat, id_lon): 
        """
        Get Latitude, Longitude, Levels, and Time for pressure levels datasets
        args:
        out_variable = MERRAgeneric().getVariables(variable, ds) 
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
        p1 = id_Latitude
        p2 = id_Longitude
        p3 = id_Level
        p4 = id_Time  
        """
   
        Lat  = {}
        Lon  = {}
        lev  = {}
        time = {}
        for i in range(0, len(out_variable)):
            Lat[i] = {}
            Lon[i] = {}
            lev[i] = {}
            time[i] = {}
            for j in range(0, len(out_variable[i])):
                print "run", i, j
                Lat[i][j]   = out_variable[i][j][p1][:]
                Lon[i][j]   = out_variable[i][j][p2][:]
                lev[i][j]   = out_variable[i][j][p3][:]
                time[i][j]  = out_variable[i][j][p4][:]   
        
        #For Latitude and Longitude 
       
        lat = {}
        lon = {}               
        for i in range(0, len(Lat)):
            lat[i] = {}
            lon[i] = {}       
            for j in range(len(Lat[i])):
                lat[i][j] = Lat[i][j][id_lat]                
                lon[i][j] = Lon[i][j][id_lon]
        
         
        return lat, lon, lev, time    

    def latLon_2d(self, out_variable, p1, p2, p3, id_lat, id_lon): 
        """
        Get Latitude, Longitude, Levels, and Time for surface datasets
        args:
        out_variable = MERRAgeneric().getVariables(variable, ds) 
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
        p1 = id_Latitude
        p2 = id_Longitude
        p3 = id_Time 
        """
   
        Lat  = {}
        Lon  = {}
        time = {}
        for i in range(0, len(out_variable)):
            Lat[i] = {}
            Lon[i] = {}
            time[i] = {}
            for j in range(0, len(out_variable[i])):
                print "run", i, j
                Lat[i][j]   = out_variable[i][j][p1][:]
                Lon[i][j]   = out_variable[i][j][p2][:]
                time[i][j]  = out_variable[i][j][p3][:]   
        
        #For Latitude and Longitude 
       
        lat = {}
        lon = {}               
        for i in range(0, len(Lat)):
            lat[i] = {}
            lon[i] = {}       
            for j in range(len(Lat[i])):
                lat[i][j] = Lat[i][j][id_lat]                
                lon[i][j] = Lon[i][j][id_lon]
        
        return lat, lon, time


    def dataStuff_3d(self, position, lat, lon, lev, time, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
         Args:
           lat, lon, lev, time = MERRAgeneric().latLon_3d(out_variable, p1, p2, p3, p4, id_lat, id_lon)  
           id_lat, id_lon =  MERRAgeneric().getArea(area, ds) 
           
           3d Analyzed Meteorological Fields datasets:
           temp = MERRAgeneric().dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v = MERRAgeneric().dataStuff(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u = MERRAgeneric().dataStuff(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
           h = MERRAgeneric().dataStuff(3, lat, lon, lev, time, id_lat, id_lon, out_variable)
           
           3d Assimilated Meteorological Fields datasets:
           rh = MERRAgeneric().dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)    
                  
           position: [t=0,v=1,u=2, h=3]
                     [rh=0]
        """
        
        data = {}
        for i in range(0, len(out_variable)):
            data[i] = {}
            for j in range(0, len(out_variable[i])):
                print "run", i, j 
                data[i][j] = out_variable[i][j][position][:]

        # Restrict the area for data set
        data_area = {}
        for i in range(0, len(data)): 
            data_area[i] = {}
            for j in range(0, len(data[i])):
                print "run", i, j
                data_area[i][j] = data[i][j][:,:,id_lat,:]
            for j in range(0, len(data_area[i])):
                data_area[i][j] = data_area[i][j][:,:,:,id_lon]
            
        del data 

        return data_area
    
    def dataStuff_2d(self, position, lat, lon, time, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
         Args:
           lat, lon, time = MERRAgeneric().latLon_2d(out_variable, p1, p2, p3, id_lat, id_lon)  
           id_lat, id_lon =  MERRAgeneric().getArea(area, ds) 
                      
           2d meteorological Diagnostics datasets:
           t2m = MERRAgeneric().dataStuff_2d(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u2m = MERRAgeneric().dataStuff_2d(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v2m = MERRAgeneric().dataStuff_2d(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
           u10m = MERRAgeneric().dataStuff_2d(3, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v10m = MERRAgeneric().dataStuff_2d(4, lat, lon, lev, time, id_lat, id_lon, out_variable)
           
           2d suface flux Diagnostics datasets:
           prectot = MERRAgeneric().dataStuff_2d(0, lat, lon, lev, time, id_lat, id_lon, out_variable
           
           2d radiation Diagnostics datasets:
           t2m = MERRAgeneric().dataStuff_2d(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u2m = MERRAgeneric().dataStuff_2d(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
       
           position: [t2m=0, u2m=1, v2m=2, u10m=3, v10m=4]
                     [prectot = 0]
                     [swgnt=0, lwgnt=1]           
        """
        data = {}
        for i in range(0, len(out_variable)):
            data[i] = {}
            for j in range(0, len(out_variable[i])):
                print "run", i, j 
                data[i][j] = out_variable[i][j][position][:]

        # Restrict the area for data set
        data_area = {}
        for i in range(0, len(data)): 
            data_area[i] = {}
            for j in range(0, len(data[i])):
                print "run", i, j
                data_area[i][j] = data[i][j][:,id_lat,:]
            for j in range(0, len(data_area[i])):
                data_area[i][j] = data_area[i][j][:,:,id_lon]
            
        del data 

        return data_area

        
class MERRApl_ana(MERRAgeneric):
    """Returns variables from downloaded MERRA 3d Analyzed Meteorological Fields datasets  
       which are abstracted with specific temporal and spatial range  
       
    Args:
        beg, end: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
                      
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['T','V','U','H','lat','lon','lev','time']
        
              
    """
    
    def getdDs(username, password, urls, size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_3dmana
           size = 5
           ds = MERRAgeneric().download(username, password, urls, size)
        """    
     
        return ds
     
    
    def getVariables(variable, ds):
        """Return the objected variables from the specific MERRA-2 3D Analyzed Meteorological Fields datasets        
           Args:
           variable = ['T','V','U','H','lat','lon','lev','time']
           ds = MERRAgeneric().download( username, password, urls_3dmana, size)
           
        """
        out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)

        return out_variable_3dmana 
        
    def getlatLon (out_variable_3dmana, p1, p2, p3, p4, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)
            p1 = 4 (id_Latitude)
            p2 = 5 (id_Longitude)
            p3 = 6 (id_Level)
            p4 = 7 (id_Time)  

        """       
        p1 = 4
        p2 = 5
        p3 = 6
        p4 = 7

        lat, lon, lev, time1 = MERRAgeneric().latLon_3d(out_variable_3dmana, id_lat, id_lon)
        
        return lat, lon, lev, time1
    

class MERRApl_asm(MERRAgeneric):
    """Returns variables from downloaded MERRA 3d Assimilated Meteorological Fields data, 
       which are abstracted with specific temporal and spatial range        

    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['RH','lat','lon','lev','time']
        

    """

    def getDs(username, password, urls, size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_3dmasm
           size = 5
           ds = MERRAgeneric().download(username, password, urls, size)
        """    
     
        return ds



    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 3D datasets        
            variable = ['RH','lat','lon','lev','time']
            ds = MERRAgeneric.download( username, password, urls_3dmasm, size)
        """
        out_variable_3dmasm = MERRAgeneric.getVariables(variable, ds)

        return out_variable_3dmasm
        
    def getlatLon (out_variable_3dmasm, p1, p2, p3, p4, id_lat, id_lon):       # Do I need to rename the out_variable or not??
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmasm = MERRAgeneric().getVariables(variable, ds)
            p1 = 1 (id_Latitude)
            p2 = 2 (id_Longitude)
            p3 = 3 (id_Level)
            p4 = 4 (id_Time)  

        """       
        p1 = 1
        p2 = 2
        p3 = 3
        p4 = 4

        time2 = MERRAgeneric().getlatLon_3d(p1, p2, p3, p4, out_variable_3dmasm, id_lat, id_lon)
        
        return time2
        

class saveNCDF_pl(temp, v, u, h, rh, time1, time2, lev, lat, lon):                                   # for saving abstracted pressure-levels variables
        """ write output netCDF file for abstracted variables from original meteorological data 
            at pressure levels
            demension: time, level, lat, lon
            variables: time1: array([   0,  360,  720, 1080], dtype=int32, Unit:minute)
                       time2: array([   0,  180,  360,  540,  720,  900, 1080, 1260], dtype=int32, Unit: minute)
                       temperature(time1,lev,lat,lon), 
                       U component of wind(time1,lev,lat,lon), 
                       V component of wind(time1,lev,lat,lon),
                       geopotential heights(time1,lev,lat,lon)
                       relative humidity(time,lev,lat,lon), time, level, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_pl.nc') # edit the name of saving nc file with specific date later, to save the data by each chunk
            temp = dataStuff(0, lat, lon, lev, time1, id_lat, id_lon, out_variable_3dmana) 
            v = dataStuff(1, lat, lon, lev, time1, id_lat, id_lon, out_variable_3dmana)
            u = dataStuff(2, lat, lon, lev, time1, id_lat, id_lon, out_variable_3dmana)       
            h = dataStuff(3, lat, lon, lev, time1, id_lat, id_lon, out_variable_3dmana)        
            rh = dataStuff(0, lat, lon, lev, time2, id_lat, id_lon, out_variable_3dmasm)  
            lat, lon, lev, time1 = MERRAgeneric().latLon(out_variable_3dmana, id_lat, id_lon)
            lat, lon, lev, time2 = MERRAgeneric().latLon(out_variable_3dmasm, id_lat, id_lon) 
                     
        """
  
        
        
        def getTime(self, start_date, end_date, time_step):                                                                                                                                          
            """set up time range 
             beg = "2016/01/01"
             end = "2016/12/31"
             """
            time_step = '6H'
            start = datetime.strptime(beg, '%Y/%m/%d')
            end  = datetime.strptime(end, "%Y/%m/%d")
            ind = pandas.date_range(start, end, freq = time_step) 
            
            time_diff = end - start
            print time_diff.days
            # get list of wanted date
            date = [start + timedelta(days=x) for x in range(time_diff.days + 1)]

            return date, ind
   

        def saveData(temp, v, u, h, rh, time1, time2, lev, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_pl.nc')                     # edit the name of saving nc file with specific date later, to save the data by each chunk

            """
            
            file_ncdf  = path.join(dir_data,'merra_pl.nc') 
            rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
            print(rootgrp.file_format)
            rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at pressure levels'
            rootgrp.featureType = "3_Dimension"
        
            #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
            TIME1 = time1[0][0]/60           # unit conversion minutues to hours
            TIME2 = time2[0][0]/60
            LEV = lev[0][0]
            LAT = lat[0][0]
            LON = lon[0][0]

            chunk_size = 5
            #dimensions
 #          times1  = rootgrp.createDimension('times1', len(TIME1)*chunk_size)
            times2  = rootgrp.createDimension('times', len(TIME2)*chunk_size)
    #       times = rootgrp.createDimension('times', None)
            levels = rootgrp.createDimension('levels', len(LEV))
            lats   = rootgrp.createDimension('lats', len(LAT))
            lons   = rootgrp.createDimension('lons', len(LON))
        
            #Temperature at pressure levels                      
            Temp            = rootgrp.createVariable('Temperature', 'f4', ('times','levels','lats', 'lons'),fill_value=9.9999999E14)
            Temp.standard_name = "air_temperature"
            Temp.units         = "K" 
            Temp.missing_value = 9.9999999E14
            Temp.fmissing_value = (9.9999999E14, 'f')
            Temp.vmax = (9.9999999E14, 'f')
            Temp.vmin = (-9.9999999E14, 'f')
            # get the data with subset of area
            temp = MERRAgeneric().dataStuff_3d(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
            #restructing the shape 
            temp_temp = []
            for i in range(0, len(temp)):
                for j in range(0,len(temp[i])):
                    for k in range(0,len(temp[i][j])):
                        temp_temp.append(temp[i][j][k][:])
            temp_temp = np.asarray(temp_temp, dtype = float)
            # pass the values
            Temp[:,:,:,:] = temp_temp[0:20,:,:,:]             
           
            #Northward of wind at pressure levels 
            V               = rootgrp.createVariable('V_component_of_wind','f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
            V.standard_name = "northward_wind"
            V.unit          = "m s-1" 
            V.missing_value = 9.9999999E14
            V.fmissing_value = 9.9999999E14, 'f'
            V.vmax = 9.9999999E14, 'f'
            V.vmin = -9.9999999E14, 'f'   
            # get the data wtih subset of area
            v = MERRAgeneric().dataStuff_3d(1, lat, lon, lev, time, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            v_temp = []
            for i in range(0, len(v)):
                for j in range(0,len(v[i])):
                    for k in range(0,len(v[i][j])):
                        v_temp.append(v[i][j][k][:])
            v_temp = np.asarray(v_temp, dtype = float)
            # pass the values
            V[:,:,:,:] = v_temp[0:20,:,:,:]             

            #Eastward of wind at pressure levels
            U               = rootgrp.createVariable('U_component_of_wind', 'f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
            U.standard_name = "eastward_wind"
            U.unit          = "m s-1" 
            U.missing_value = 9.9999999E14
            U.fmissing_value = (9.9999999E14, 'f')
            U.vmax = (9.9999999E14, 'f')
            U.vmin = (-9.9999999E14, 'f')
            #get the data with subset of area
            u = MERRAgeneric().dataStuff_3d(2, lat, lon, lev, time, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            u_temp = []
            for i in range(0, len(u)):
                for j in range(0,len(u[i])):
                    for k in range(0,len(u[i][j])):
                        u_temp.append(u[i][j][k][:])
            u_temp = np.asarray(u_temp, dtype = float)
            # pass the values
            U[:,:,:,:] = u_temp[0:20,:,:,:]             
            
            # Geopotential Height at pressure levels
            H = rootgrp.createVariable('H', 'f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
            H.standard_name = "geopotential height"
            H.units         = "m+2 s-2" 
            H.missing_value = 9.9999999E14
            H.fmissing_value = (9.9999999E14,'f')
            H.vmax = (9.9999999E14, 'f')
            H.vmin = (-9.9999999E14, 'f')
            #get the data with subset of area
            h = MERRAgeneric().dataStuff_3d(3, lat, lon, lev, time, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            h_temp = []
            for i in range(0, len(h)):
                for j in range(0,len(h[i])):
                    for k in range(0,len(h[i][j])):
                        h_temp.append(h[i][j][k][:])                       
            h_temp = np.asarray(h_temp, dtype = float)
            # pass the values
            H[:,:,:,:] = h_temp[0:20,:,:,:]             

            # Relative Humidity after moist at pressure levels
            RH               = rootgrp.createVariable('Relative_humidity', 'f4', ('times2','levels', 'lats', 'lons'),fill_value=9.9999999E14)
            RH.standard_name = "relative humidity after moist"
            RH.units       = "1"  
            RH.missing_value = 9.9999999E14
            RH.fmissing_value = (9.9999999E14, 'f')
            RH.vmax = (9.9999999E14, 'f')
            RH.vmin = (-9.9999999E14, 'f')
            #get the data with subset of area
            rh = MERRAgeneric().dataStuff_3d(0, lat, lon, lev, time2, id_lat, id_lon, out_variable_3dmasm)
            #restructing the shape 
            rh_temp = []
            for i in range(0, len(rh)):
                for j in range(0,len(rh[i])):
                    for k in range(0,len(rh[i][j])):
                        rh_temp.append(rh[i][j][k][:])
            rh_temp = np.asarray(rh_temp, dtype = float)
            # pass the values
            RH[:,:,:,:] = rh_temp[0:40,:,:,:]             
        
            Time1               = rootgrp.createVariable('time1', 'i4', ('times'))
            Time1.standard_name = "time"
            Time1.units         = "hour since" + beg.strftime("%Y-%m-%d 00:00:00") # needed to add the beging date into it
            Time1.calendar      = "standard"
            # pass the values
            Time1[:] = time1[0][0][:]                                           # 6-hourly time step (for Temp, U, V, H)                                           
            # Time[:] = time2[0][0][:]                                         # 3-hourly time step (for RH)

            Time2               = rootgrp.createVariable('time2', 'i4', ('times'))
            Time2.standard_name = "time"
            Time2.units         = "hour since" + beg.strftime("%Y-%m-%d 00:00:00") # needed to add the beging date into it
            Time2.calendar      = "standard"
            # pass the values
            Time2[:] = time2[0][0][:]                                           # 6-hourly time step (for Temp, U, V, H)                                           
                               

            Level               = rootgrp.createVariable('level','i4', ('levels'))
            Level.standard_name = "air_pressure"
            Level.units         = "hPa"
            # pass the values
            Level[:]      = lev[0][0][:]                    

            Latitudes               = rootgrp.createVariable('latitudes', 'f4',('lats'))
            Latitudes.standard_name = "latitude"
            Latitudes.units         = "degrees_north"
            Latitudes.axis          = 'Y'
            Latitudes[:]  = lat[0][0][:]                    # pass the values of latitude

            Longitudes               = rootgrp.createVariable('longitudes', 'f4',('lons'))
            Longitudes.standard_name = "longitude"
            Longitudes.units         = "degrees_east"
            Longitudes.axis          = 'X'
            Longitudes[:] = lon[0][0][:]                    # pass the values of longitudes
        
        
            #close the root group

            rootgrp.close()          


class MERRAsm(MERRAgeneric):
    """Returns variables from downloaded MERRA 2d meteorological Diagnostics data, 
       which are abstracted with specific temporal and spatial range        
       
    Args:
        beg, end: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['T2M','U2M','V2M','U10M','V10M','lat','lon','time'].
                   T2M: 2-meter_air_temperature
                   U2M: 2-meter_eastward_wind
                   V2M: 2-meter_northward_wind
                   U10M: 10-meter_eastward_wind
                   V10M: 10-meter_northward_wind          
                      
     """
    def getDs(username, password, urls, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_2dm
           chunk_size = 5
           ds = MERRAgeneric().download(username, password, urls, chunk_size)
        """    
     
        return ds

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 3D datasets        
            variable = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']
            ds = MERRAgeneric.download( username, password, urls_2dm, chunk_size)
            
        """
        out_variable_2dm = MERRAgeneric().getVariables(variable, ds)

        return out_variable_2dm
         
    def getlatLon (out_variable_2dm, p1, p2, p4, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)
            p1 = 7 (id_Latitude)
            p2 = 8 (id_Longitude)
            p4 = 9 (id_Time) 

        """       
        p1 = 7
        p2 = 8
        p3 = 9

        lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dm, id_lat, id_lon)
        
        return lat, lon, time
        


class MERRAsf(MERRAgeneric):
    """Returns variables from downloaded MERRA 2d suface flux Diagnostics data, 
       which are abstracted with specific temporal and spatial range        
       
    Args:
        beg, end: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
                      
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['PRECTOT','lat','lon','time'].

              
    """
    
    def getDs(username, password, urls, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_2dm
           chunk_size = 5
           ds = MERRAgeneric().download(username, password, urls, chunk_size)
        """    
     
        return ds
    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D      
           variable = ['PRECTOT','lat','lon','time']
           ds = MERRAgeneric.download( username, password, urls_2ds, chunk_size)
        """        
        out_variable_2ds = MERRAgeneric.getVariables(variable, ds)

    return out_variable_2ds

    def getlatLon_2d (out_variable_2ds, p1, p2, p3, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)
            p1 = 1 (id_Latitude)
            p2 = 2 (id_Longitude)
            p4 = 3 (id_Time) 

        """       
        p1 = 1
        p2 = 2
        p3 = 3

        lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2ds, id_lat, id_lon)
        
        return lat, lon, time


class saveNCDF_sa(t2m, u2m, v2m, u10m, v10m, prectot, time, lev, lat, lon):                                  
        """ write output netCDF file for abstracted variables from original  2d meteorological Diagnostics datasetsc and suface flux Diagnostics datasets
            demension: time, lat, lon
            variables: time: array([   0,   60,  120,  180,  240,  300,  360,  420,  480,  540,  600,
                             660,  720,  780,  840,  900,  960, 1020, 1080, 1140, 1200, 1260,
                             1320, 1380], dtype=int32, Unit: minute)
                       t2m(time,llat,lon), 
                       u2m(time,lat,lon), 
                       v2m(time,lat,lon),
                       u10m(time,lat,lon)
                       v10m(time,lat,lon)
                       prectot(time,lat,lon), 
                       time, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_pl.nc') # edit the name of saving nc file with specific date later, to save the data by each chunk
            t2m = dataStuff_2d(0, lat, lon,time, id_lat, id_lon, out_variable_2dm) 
            u2m = dataStuff_2d(1, lat, lon, time, id_lat, id_lon, out_variable_2dm)
            v2m = dataStuff_2d(2, lat, lon, time, id_lat, id_lon, out_variable_2dm)       
            u10m = dataStuff_2d(3, lat, lon, time, id_lat, id_lon, out_variable_2dm)        
            v10m = dataStuff_2d(4, lat, lon, time, id_lat, id_lon, out_variable_2dm) 
            prectot = dataStuff_2d(0, lat, lon, time, id_lat, id_lon, out_variable_2ds) 
            lat, lon,time = MERRAgeneric().latLon_2d(out_variable_2dm, id_lat, id_lon)
                     
        """
     

        def saveData(t2m, u2m, v2m, u10m, v10m, time, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_sa.nc')                     # edit the name of saving nc file with specific date later, to save the data by each chunk

            """
            
            file_ncdf  = path.join(dir_data,'merra_sa.nc') 
            rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
            print(rootgrp.file_format)
            rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at surface level'
            rootgrp.featureType = "2_Dimension"
        
            #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
            TIME = time[0][0]/60           # unit conversion minutues to hours
            LAT = lat[0][0]
            LON = lon[0][0]

            chunk_size = 5
            #dimensions
            times  = rootgrp.createDimension('times', len(TIME)*chunk_size)
    #       times = rootgrp.createDimension('times', None)
            lats   = rootgrp.createDimension('lats', len(LAT))
            lons   = rootgrp.createDimension('lons', len(LON))
        
            #Temperature at 2 meter                      
            T2M            = rootgrp.createVariable('2_meter_air_temperature', 'f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            T2M.standard_name = "2_meter_air_temperature"
            T2M.units         = "K" 
            T2M.missing_value = 9.9999999E14
            T2M.fmissing_value = (9.9999999E14, 'f')
            T2M.vmax = (9.9999999E14, 'f')
            T2M.vmin = (-9.9999999E14, 'f')
            # get the data with subset of area
            t2m = MERRAgeneric().dataStuff(0, lat, lon, time, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            t2m_temp = []
            for i in range(0, len(t2m)):
                for j in range(0,len(t2m[i])):
                    for k in range(0,len(t2m[i][j])):
                        t2m_temp.append(t2m[i][j][k][:])
            t2m_temp = np.asarray(t2m_temp, dtype = float)
            del t2m
            # pass the values
            T2M[:,:,:] = t2m_temp[0:120,:,:]   # len(TIME)*(chunk_size)          
            
            # Eastward Wind at 2 meter surface level
            U2M               = rootgrp.createVariable('2_meter_eastward_wind','f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            U2M.standard_name = "2_meter_eastward_wind"
            U2M.unit          = "m s-1" 
            U2M.missing_value = 9.9999999E14
            U2M.fmissing_value = 9.9999999E14, 'f'
            U2M.vmax = 9.9999999E14, 'f'
            U2M.vmin = -9.9999999E14, 'f'   
            # get the data wtih subset of area
            u2m = MERRAgeneric().dataStuff_2d(1, lat, lon, time, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            u2m_temp = []
            for i in range(0, len(u2m)):
                for j in range(0,len(u2m[i])):
                    for k in range(0,len(u2m[i][j])):
                        u2m_temp.append(u2m[i][j][k][:])
            u2m_temp = np.asarray(u2m_temp, dtype = float)
            del u2m
            # pass the values
            U2M[:,:,:] = u2m_temp[0:120,:,:]             

            #Northward Wind at 2 meter surface level
            V2M               = rootgrp.createVariable('2_meter_northward_wind','f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            V2M.standard_name = "2_meter_eastward_wind"
            V2M.unit          = "m s-1" 
            V2M.missing_value = 9.9999999E14
            V2M.fmissing_value = 9.9999999E14, 'f'
            V2M.vmax = 9.9999999E14, 'f'
            V2M.vmin = -9.9999999E14, 'f'   
            # get the data wtih subset of area
            v2m = MERRAgeneric().dataStuff_2d(2, lat, lon, time, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            v2m_temp = []
            for i in range(0, len(v2m)):
                for j in range(0,len(v2m[i])):
                    for k in range(0,len(v2m[i][j])):
                        v2m_temp.append(v2m[i][j][k][:])
            v2m_temp = np.asarray(v2m_temp, dtype = float)
            del v2m
            # pass the values
            V2M[:,:,:] = v2m_temp[0:120,:,:]             
            
            # Eastward Wind at 10 meter surface level
            U10M               = rootgrp.createVariable('10_meter_eastward_wind','f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            U10M.standard_name = "10_meter_eastward_wind"
            U10M.unit          = "m s-1" 
            U10M.missing_value = 9.9999999E14
            U10M.fmissing_value = 9.9999999E14, 'f'
            U10M.vmax = 9.9999999E14, 'f'
            U10M.vmin = -9.9999999E14, 'f'   
            # get the data wtih subset of area
            u10m = MERRAgeneric().dataStuff_2d(3, lat, lon, time, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            u10m_temp = []
            for i in range(0, len(u10m)):
                for j in range(0,len(u10m[i])):
                    for k in range(0,len(u10m[i][j])):
                        u10m_temp.append(u10m[i][j][k][:])
            u10m_temp = np.asarray(u10m_temp, dtype = float)
            del u10m
            # pass the values
            U10M[:,:,:] = u10m_temp[0:120,:,:]             

            #Northward Wind at 10 meter surface level
            V10M               = rootgrp.createVariable('10_meter_northward_wind','f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            V10M.standard_name = "10_meter_eastward_wind"
            V10M.unit          = "m s-1" 
            V10M.missing_value = 9.9999999E14
            V10M.fmissing_value = 9.9999999E14, 'f'
            V10M.vmax = 9.9999999E14, 'f'
            V10M.vmin = -9.9999999E14, 'f'   
            # get the data wtih subset of area
            v10m = MERRAgeneric().dataStuff_2d(4, lat, lon, time, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            v10m_temp = []
            for i in range(0, len(v10m)):
                for j in range(0,len(v10m[i])):
                    for k in range(0,len(v10m[i][j])):
                        v10m_temp.append(v10m[i][j][k][:])
            v10m_temp = np.asarray(v10m_temp, dtype = float)
            del v10m
            # pass the values
            V10M[:,:,:] = v10m_temp[0:120,:,:]   
            
            # Total Precipitation          
            PRECTOT               = rootgrp.createVariable('total_precipitation','f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            PRECTOT.standard_name = "total_precipitation"
            PRECTOT.unit          = "kg m-2 s-1" 
            PRECTOT.missing_value = 9.9999999E14
            PRECTOT.fmissing_value = 9.9999999E14, 'f'
            PRECTOT.vmax = 9.9999999E14, 'f'
            PRECTOT.vmin = -9.9999999E14, 'f'   
            # get the data wtih subset of area
            prectot = MERRAgeneric().dataStuff_2d(0, lat, lon, time, id_lat, id_lon, out_variable_2ds)
            #restructing the shape 
            prectot_temp = []
            for i in range(0, len(prectot)):
                for j in range(0,len(prectot[i])):
                    for k in range(0,len(prectot[i][j])):
                        prectot_temp.append(prectot[i][j][k][:])
            prectot_temp = np.asarray(prectot_temp, dtype = float)
            del prectot
            # pass the values
            PRECTOT[:,:,:] = prectot_temp[0:120,:,:]   
       
            Time               = rootgrp.createVariable('time', 'i4', ('times'))
            Time.standard_name = "time"
            Time.units         = "hour since" + beg.strftime("%Y/%m/%d 00:00:00") # needed to add the beging date into it
            Time.calendar      = "standard"
            # pass the values
            Time[:] = time[0][0][:]                                               # 1-hourly time step                                            
  
            Latitudes               = rootgrp.createVariable('latitudes', 'f4',('lats'))
            Latitudes.standard_name = "latitude"
            Latitudes.units         = "degrees_north"
            Latitudes.axis          = 'Y'
            Latitudes[:]  = lat[0][0][:]                    # pass the values of latitude

            Longitudes               = rootgrp.createVariable('longitudes', 'f4',('lons'))
            Longitudes.standard_name = "longitude"
            Longitudes.units         = "degrees_east"
            Longitudes.axis          = 'X'
            Longitudes[:] = lon[0][0][:]                    # pass the values of longitudes
        
        
            #close the root group

            rootgrp.close()          


class MERRAsr(MERRAgeneric):
    """Returns variables from downloaded MERRA 2d radiation Diagnostics datasets, 
       which are abstracted with specific temporal and spatial range        
       
    Args:
        beg, end: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['SWGNT','LWGNT','lat','lon','time'].
                   SWGNT:surface net downward shortwave flux(time*lat*lon)
                   LWGNT:surface net downward longwave flux(time*lat*lon)
        
        directory: Directory to hold output files
              
    Example:
        from datetime import datetime
        beg   = "2016/01/01"
        end   = "2016/02/01" 
        area = {'bbS':60.0, 'bbN': 65.0, 'bbW':-115.0, 'bbE': -110.0}
        dir_data = '/Users/xquan/data'             
        MERRAsr = MERRAsr(beg, end, area, elevation, variable, dir_data) 
        MERRAsr.download()
    """
    
    def getDs(username, password, urls, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 2d radiation Diagnostics datasets 
           Args:
           username = ******
           password = ******
           urls = urls_2dr
           chunk_size = 5
           ds = MERRAgeneric().download(username, password, urls, chunk_size)
        """    
     
        return ds

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D radiation Diagnostics datasets        
            variable = ['SWGNT','LWGNT', 'lat','lon','time']
            urls = urls_2dr
            ds = MERRAgeneric.download( username, password, urls, chunk_size)
            
            
        """
        out_variable_2dr = MERRAgeneric().getVariables(variable, ds)

        return out_variable_2dr
         
    def getlatLon (out_variable_2dr, p1, p2, p3, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)
            p1 = 2 (id_Latitude)
            p2 = 3 (id_Longitude)
            p4 = 4 (id_Time) 

        """       
        p1 = 2
        p2 = 3
        p3 = 4

        lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dr, id_lat, id_lon)
        
        return lat, lon, time
        
class saveNCDF_sr(swgnt, lwgnt, time, lat, lon):                                  
        """ write output netCDF file for abstracted variables from original 2d radiation Diagnostics datasets  datasets 
            demension: time, lat, lon
            variables: time: array([   0,   60,  120,  180,  240,  300,  360,  420,  480,  540,  600,
                             660,  720,  780,  840,  900,  960, 1020, 1080, 1140, 1200, 1260,
                             1320, 1380], dtype=int32, Unit: minute)
                       swgnt(time,llat,lon), 
                       lwgnt(time,lat,lon), 
                       time, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_sr.nc') # edit the name of saving nc file with specific date later, to save the data by each chunk
            swgnt = dataStuff_2d(0, lat, lon,time, id_lat, id_lon, out_variable_2dr) 
            lwgnt = dataStuff_2d(1, lat, lon, time, id_lat, id_lon, out_variable_2dr)
            lat, lon,time = MERRAgeneric().latLon_2d(out_variable_2dr, id_lat, id_lon)
                     
        """
     

        def saveData(swgnt, lwgnt, time, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_sr.nc')                     # edit the name of saving nc file with specific date later, to save the data by each chunk

            """
            
            file_ncdf  = path.join(dir_data,'merra_sr.nc') 
            rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
            print(rootgrp.file_format)
            rootgrp.source      = 'Merra, abstrated radiation variables from metadata at surface level'
            rootgrp.featureType = "2_Dimension"
        
            #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
            TIME = time[0][0]/60           # unit conversion minutues to hours
            LAT = lat[0][0]
            LON = lon[0][0]

            chunk_size = 5
            #dimensions
            times  = rootgrp.createDimension('times', len(TIME)*chunk_size)
    #       times = rootgrp.createDimension('times', None)
            lats   = rootgrp.createDimension('lats', len(LAT))
            lons   = rootgrp.createDimension('lons', len(LON))
        
            #surface net downward shortwave flux                      
            SWGNT            = rootgrp.createVariable('surface_net_downward_shortwave_flux', 'f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            SWGNT.standard_name = "surface_net_downward_shortwave_flux"
            SWGNT.units         = "W m-2" 
            SWGNT.missing_value = 9.9999999E14
            SWGNT.fmissing_value = (9.9999999E14, 'f')
            SWGNT.vmax = (9.9999999E14, 'f')
            SWGNT.vmin = (-9.9999999E14, 'f')
            # get the data with subset of area
            swgnt = MERRAgeneric().dataStuff(0, lat, lon, time, id_lat, id_lon, out_variable_2dr)
            #restructing the shape 
            swgnt_temp = []
            for i in range(0, len(swgnt)):
                for j in range(0,len(swgnt[i])):
                    for k in range(0,len(swgnt[i][j])):
                        swgnt_temp.append(swgnt[i][j][k][:])
            swgnt_temp = np.asarray(swgnt_temp, dtype = float)
            del swgnt
            # pass the values
            SWGNT[:,:,:] = swgnt_temp[0:120,:,:]   # len(TIME)*(chunk_size)          

            #surface net downward longwave flux                      
            LWGNT            = rootgrp.createVariable('surface_net_downward_longwave_flux', 'f4', ('times','lats', 'lons'),fill_value=9.9999999E14)
            LWGNT.standard_name = "surface_net_downward_longwave_flux"
            LWGNT.units         = "W m-2" 
            LWGNT.missing_value = 9.9999999E14
            LWGNT.fmissing_value = (9.9999999E14, 'f')
            LWGNT.vmax = (9.9999999E14, 'f')
            LWGNT.vmin = (-9.9999999E14, 'f')
            # get the data with subset of area
            lwgnt = MERRAgeneric().dataStuff(0, lat, lon, time, id_lat, id_lon, out_variable_2dr)
            #restructing the shape 
            lwgnt_temp = []
            for i in range(0, len(lwgnt)):
                for j in range(0,len(lwgnt[i])):
                    for k in range(0,len(lwgnt[i][j])):
                        lwgnt_temp.append(lwgnt[i][j][k][:])
            lwgnt_temp = np.asarray(lwgnt_temp, dtype = float)
            del lwgnt
            # pass the values
            LWGNT[:,:,:] = lwgnt_temp[0:120,:,:]   # len(TIME)*(chunk_size)          
                   
            Time               = rootgrp.createVariable('time', 'i4', ('times'))
            Time.standard_name = "time"
            Time.units         = "hour since" + beg.strftime("%Y-%m-%d 00:00:00") # needed to add the beging date into it
            Time.calendar      = "standard"
            # pass the values
            Time[:] = time[0][0][:]                                               # 1-hourly time step                                            
  
            Latitudes               = rootgrp.createVariable('latitudes', 'f4',('lats'))
            Latitudes.standard_name = "latitude"
            Latitudes.units         = "degrees_north"
            Latitudes.axis          = 'Y'
            Latitudes[:]  = lat[0][0][:]                    # pass the values of latitude

            Longitudes               = rootgrp.createVariable('longitudes', 'f4',('lons'))
            Longitudes.standard_name = "longitude"
            Longitudes.units         = "degrees_east"
            Longitudes.axis          = 'X'
            Longitudes[:] = lon[0][0][:]                    # pass the values of longitudes
        
        
            #close the root group

            rootgrp.close()          
               
#=========================For Run MERRA-2======================================
# Output Datasets Types:

# - merra_pl_20160101_to_20160105: The extracted meteorological variables with defined chunk at pressure levels from original Merra-2 data products with specific temporal and spatial setups  
# - merra_sa_20160101_to_20160105: The extracted meteorological variables with defined chunk at surface level from original Merra-2 data products with specific temporal and spatial setups
# - merra_sr_20160101_to_20160105: The extracted radiation variables with with defined chunk at surface level from original Merra-2 data products with specific temporal and spatial setups  
#     
#==============================================================================    


from pydap.client import open_url
from pydap.cas.urs import setup_session
from datetime import datetime, timedelta
from os import path
from netCDF4 import Dataset

import pydap.lib
import numpy as np
import csv
import netCDF4 as nc
import itertools
import pandas


#Account for Database Access
username = "quanxj17"
password = "Qxj17carleton"

#settings directory 
dir_data = '/Users/xquan/data'
# dir_data = '/Users/xquan/src/globsim/merra2'

dir_src  = '/Users/xquan/src/globsim/'

#Given wanted datatype, mydate, area, elevation

beg   = "2016/01/01"
end   = "2016/02/01" 
                                                                               
area = {'bbS':60.0, 'bbN': 65.0, 'bbW':-115.0, 'bbE': -110.0}
 
elevation = {'min' : 50, 'max' : 2000}  


#==============================================================================
