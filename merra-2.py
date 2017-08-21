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
#         '/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')                    # 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics,meteorological Diagnostics  

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
            end = "2016/02/01"
            
            Return:
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
        Begin = datetime.strptime(beg, '%Y/%m/%d')
        End  = datetime.strptime(end, "%Y/%m/%d")

        #Setup the based string of dates for urls 
        res1 = [d.strftime("%Y/%m") for d in pandas.date_range(Begin,End)]
        res2 = [d.strftime("%Y%m%d") for d in pandas.date_range(Begin,End)]        
        
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
            chunk_size: the wanted size of urls list for each chunk
            Return:
            ds: the each individual original dataset (The dictionary of basetype with children['variables']) 
                                                           of 3d Analyzed Meteorological Fields, 
                                                              3d Assimilated Meteorological Fields datasets,
                                                              2d meteorological Diagnostics datasets,
                                                              2d surface flux Diagnostics datasets,
                                                              2d radiation Diagnostics datasets
            ds_structure: [lengths of total_chunks * chunk_size]
                                                              
                              
        """
        chunk_size = 5
        urls_chunks = [urls[x:x+chunk_size] for x in xrange(0, len(urls), chunk_size)]      

        print ('================ MERRA-2: START ===============')
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
        print ('=============== MERRA-2: COMPLETED ================')
        infor = urls[0].split('/')
        print 'Dataset:', infor[2], infor[3],infor[4], infor[8]
        print 'Type:', type(ds)
        print 'Days:', len(urls)    
        
        return ds        
    
    def Variables(self, variable, ds):
        """Get the objected variables from the specific MERRA-2 datasets        
           variable = ['T','V','U','H','lat','lon','lev','time']                         # for extrating from 3d Analyzed Meteorological Fields datasets 
                    = ['RH','lat','lon','lev','time']                                    # for extracting from 3d Assimilated Meteorological Fields datasets
                    = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']  # for extracting from 2d meteorological Diagnostics datasets
                    = ['PRECTOT','lat','lon','time']                                     # for extracing from 2d suface flux Diagnostics datasets  
                    = ['SWGNT','LWGNT','lat','lon','time']                               # for extracing from 2d radiation Diagnostics datasets 
                    
           ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, chunk_size)
           Return:
           out_variable: The extracted variables in the order of given list 
                         Type: BaseType with the given data baseProxy
                         
           out_variable structure: [lenghs of total chunks * chunk_size * lenghs of variables]
           
        """
       
        out_variable = {}
        for i in range(0, len(ds)):
            out_variable[i] = {}
            for j in range(len(ds[i])):
                print "Run", "Chunk:", i, "NO.:", j
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

        print "Length of Out_Variable:", len(out_variable[0][0])
        
        return out_variable
    
    def getArea(self, area, ds): 
        """Gets the specific indexs  of the latitudes and longitudes of given area
           For example: 
           area = {'bbS':60.0, 'bbN': 65.0, 'bbW':-115.0, 'bbE': -110.0}
           ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm, urls_2dm, urls_2dr, urls_2ds, chunk_size)
           Return:
           id_lat: wanted indexs of latitudes from the original merra global dataset
           id_lon: wanted indexs of longitudes from the original merra global dataset
        """
             
        # pass the value of individual row lat and lon to Lat and Lon for the area subset
        Lat = ds[0][0].lat[:]
        Lon = ds[0][0].lon[:]
                        
        # get the indices of selected range of Lat,Lon
        id_lon = np.where((Lon[:] > area['bbW']) & (Lon[:] < area['bbE'])) 
        id_lat = np.where((Lat[:] > area['bbS']) & (Lat[:] < area['bbN'])) 
       
        # convert id_lat, id_lon from tuples to string
        id_lon = list(itertools.chain(*id_lon))   
        id_lat = list(itertools.chain(*id_lat))

        
        return id_lat, id_lon

    def latLon_3d(self, out_variable, p1, p2, p3, p4, id_lat, id_lon): 
        """
        Get Latitude, Longitude, Levels, and Time for datasets at the pressure levels
        Args:
        out_variable = MERRAgeneric().getVariables(variable, ds) 
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
        p1 = id_Latitude
        p2 = id_Longitude
        p3 = id_Level
        p4 = id_Time 
        Return: 
        lat, lon, lev, time: the extracted latitudes, longitudes in the given area from 3D original dateset
        lat_structure: lengths of total_chunks * chunk_size * lengths of extracted latitudes
        lon_structure: lengths of total_chunks * chunk_size * lengths of extracted longitudes
        lev_structure: lengths of total_chunks * chunk_size * total number of pressure levels from original dataset 
        time_structure: lengths of total_chunks * chunk_size * total number of hours per day from original dataset                                        
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
                print "run", "Chunk:", i, "NO.:", j
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
        Get Latitude, Longitude, Levels, and Time for datasets at surface level
        Args:
        out_variable = MERRAgeneric().getVariables(variable, ds) 
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
        p1 = id_Latitude
        p2 = id_Longitude
        p3 = id_Time 
        Return:
        lat, lon, lev, time: the extracted latitudes, longitudes in the given area from 2D original dateset
        lat_structure: length of total_chunks * chunk_size * lengths of extracted latitudes
        lon_structure: length of total_chunks * chunk_size * lengths of extracted longitudes
        time_structure: length of total_chunks * chunk_size * total number of hours per day from original dataset                               
        """
   
        Lat  = {}
        Lon  = {}
        time = {}
        for i in range(0, len(out_variable)):
            Lat[i] = {}
            Lon[i] = {}
            time[i] = {}
            for j in range(0, len(out_variable[i])):
                print "Run", "Chunk:",i, "NO.:", j
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


    def dataStuff_3d(self, position, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
           Args:
           id_lat, id_lon =  MERRAgeneric().getArea(area, ds) 
           position: [t=0,v=1,u=2, h=3]
                     [rh=0]

           3d Analyzed Meteorological Fields datasets:
           t = MERRAgeneric().dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v = MERRAgeneric().dataStuff(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u = MERRAgeneric().dataStuff(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
           h = MERRAgeneric().dataStuff(3, lat, lon, lev, time, id_lat, id_lon, out_variable)
           
           3d Assimilated Meteorological Fields datasets:
           rh = MERRAgeneric().dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)    
                  
           Return:
           data_area: the extracted individual variable at pressure levels at the given area
           data_area_Structure: length of total_chunks * chunk_size * [time*level*lat*lon]
        """
        print "get data"
        data = {}
        for i in range(0, len(out_variable)):
            data[i] = {}
            for j in range(0, len(out_variable[i])):
                print "Run","Chunk", i, "NO.:", j 
                data[i][j] = out_variable[i][j][position][:]

        # Restrict the area for data set
        print "restrict area"
        data_area = {}
        for i in range(0, len(data)): 
            data_area[i] = {}
            for j in range(0, len(data[i])):
                print "Run", "Chunk", i, "NO.:", j
                data_area[i][j] = data[i][j][:,:,id_lat,:]
            for j in range(0, len(data_area[i])):
                data_area[i][j] = data_area[i][j][:,:,:,id_lon]
            
        del data 

        return data_area
    
    def dataStuff_2d(self, position, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
         Args:            
           id_lat, id_lon =  MERRAgeneric().getArea(area, ds) 
           position: [t2m=0, u2m=1, v2m=2, u10m=3, v10m=4]
                     [prectot = 0]
                     [swgnt=0, lwgnt=1]
          
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
           Return:
           data_area: the extracted individual variable at surface levels at the given area
           data_area_structure: length of total_chunks * chunk_size * [time*lat*lon]
           
        """
        print "get data"
        data = {}
        for i in range(0, len(out_variable)):
            data[i] = {}
            for j in range(0, len(out_variable[i])):
                print "Run","Chunk", i, "NO.:", j
                data[i][j] = out_variable[i][j][position][:]

        # Restrict the area for data set
        print "restrict area"
        data_area = {}
        for i in range(0, len(data)): 
            data_area[i] = {}
            for j in range(0, len(data[i])):
                print "Run","Chunk", i, "NO.:", j
                data_area[i][j] = data[i][j][:,id_lat,:]
            for j in range(0, len(data_area[i])):
                data_area[i][j] = data_area[i][j][:,:,id_lon]
            
        del data 

        return data_area
        
    def getTime(self, beg, end):                                                                                                                                          
        """set up date and time series for netCDF output results 
            beg = "2016/01/01"
            end = "2016/02/01"
            Return: 
            Date: a string list of date in the range of begin and end (format: yearmonthday)
            time_ind1: the time series for 6-hours step in the range of begin and end dates  
            time_ind2: the time series for 3-hours step in the range of begin and end dates
            time_ind3: the time series for 1-hour step in the range of begin and end dates
        """ 
            
            
        Start = datetime.strptime(beg, '%Y/%m/%d')
        End  = datetime.strptime(end, '%Y/%m/%d')
            
        # Set up the wanted time step
        time_step1 = '6H'
        time_step2 = '3H'
        time_step3 = '1H'
            
        #get extra one one more day for getting the full range of time series
        End1 = End + timedelta(days=1)           
            
        #get the Datetimeindex with time_step 
        time_ind1 = pandas.date_range(Start, End1, freq = time_step1)       
        time_ind2 = pandas.date_range(Start, End1, freq = time_step2) 
        time_ind3 = pandas.date_range(Start, End1, freq = time_step3)
        
        # time_ind1 = Start + np.arange(24) * timedelta(hours=6)
        # time_ind1 = Start + np.arange(24) * timedelta(hours=3)
        # time_ind1 = Start + np.arange(24) * timedelta(hours=1)        

        #Convert the Datetimeindex to numpy and get rid of one extra time in the end
        time_ind1 = (np.array(time_ind1))[0:-1]
        time_ind2 = (np.array(time_ind2))[0:-1]
        time_ind3 = (np.array(time_ind3))[0:-1]
                         
        # get list of wanted date series
        date_diff = End - Start
        date = [Start + timedelta(days=x) for x in range(date_diff.days + 1)]
        date = [d.strftime('%Y%m%d') for d in date]

        return date, time_ind1, time_ind2, time_ind3
 
    def restruDatastuff(self, data_area):                                                                                                                                          
        """ Restructuring the dimension of abstracted data stuff for preparing to save netCDF output results furtherly
            return: 
            data_total: for 3d: (len(date)*len(time/day), level, lat, lon)
                        for 2d: (len(date)*len(time/day), lat, lon)
        """ 
        data_total = []                                                   
        for i in range(0, len(data_area)):                                           
            for j in range(0,len(data_area[i])):
                for k in range(0,len(data_area[i][j])):
                    data_total.append(data_area[i][j][k][:])
        
        data_total = np.asarray(data_total, dtype = float)

        return data_total
        
        
class MERRApl_ana():
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
    
    def getDs(self, beg, end, username, password, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           beg = "2016/01/01"
           end = "2016/02/01"
           chunk_size = 5
           ds = MERRAgeneric().download(username, password, urls, chunk_size)
        """    

        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr = MERRAgeneric().getURLs(beg, end)
        print urls_3dmana
        urls = urls_3dmana
        ds = MERRAgeneric().download(username, password, urls, chunk_size)
        
        return ds
     
    
    def getVariables(self, ds):                                    
        """Return the objected variables from the specific MERRA-2 3D Analyzed Meteorological Fields datasets        
           Args:
           ds = MERRAgeneric().download( username, password, urls_3dmana, size)
           
        """
        variable = ['T','V','U','H','lat','lon','lev','time']
        out_variable_3dmana = MERRAgeneric().Variables(variable, ds)

        return out_variable_3dmana 
        
    def getlatLon_3d (self, area, ds, out_variable_3dmana, id_lat, id_lon):
        # old: def getlatLon (self, area, ds, out_variable_3dmana, p1, p2, p3, p4, id_lat, id_lon):
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
        
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)

        lat, lon, lev, time = MERRAgeneric().latLon_3d(out_variable_3dmana, p1, p2, p3, p4, id_lat, id_lon)
        
        return lat, lon, lev, time
    

class MERRApl_asm():
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

    def getDs(self, beg, end, username, password, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_3dmasm
           size = 5
           ds = MERRAgeneric().download(username, password, urls, size)
        """    
        
        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr = MERRAgeneric().getURLs(beg, end)
        urls = urls_3dmasm
        ds = MERRAgeneric().download(username, password, urls, chunk_size)
        return ds
        

    def getVariables(self, ds):
        """Return the objected variables from the specific MERRA-2 3D datasets        
            variable = ['RH','lat','lon','lev','time']
            ds = MERRAgeneric.download( username, password, urls_3dmasm, chunk_size)
        """
        variable = ['RH','lat','lon','lev','time']
        out_variable_3dmasm = MERRAgeneric().Variables(variable, ds)

        return out_variable_3dmasm
        
    def getlatLon_3d (self, area, ds, out_variable_3dmasm, id_lat, id_lon):       # Do I need to rename the out_variable or not??
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
        
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
        
        lat, lon, lev, time = MERRAgeneric().latLon_3d(out_variable_3dmasm, p1, p2, p3, p4, id_lat, id_lon)
        
        return lat, lon, lev, time
        

class SaveNCDF_pl_3dmana():                                                         # for saving abstracted pressure-levels variables
        """ write output netCDF file for abstracted variables from original meteorological data 
            at pressure levels
            demension: time, level, lat, lon
            variables: time: array([   0,  360,  720, 1080], dtype=int32, Unit:minute)
                       temperature(time,lev,lat,lon), 
                       U component of wind(time,lev,lat,lon), 
                       V component of wind(time,lev,lat,lon),
                       geopotential heights(time,lev,lat,lon)
                       time, level, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            date, time_ind1 = MERRAgeneric().getTime(beg, end)
            t = MERRAgeneric().dataStuff_3d(0, id_lat, id_lon, out_variable_3dmana) 
            v = MERRAgeneric().dataStuff_3d(1, id_lat, id_lon, out_variable_3dmana)
            u = MERRAgeneric().dataStuff_3d(2,  id_lat, id_lon, out_variable_3dmana)       
            h = MERRAgeneric().dataStuff_3d(3, id_lat, id_lon, out_variable_3dmana)        
            lat, lon, lev, time1 = MERRAgeneric().latLon(out_variable_3dmana, id_lat, id_lon)
                     
        """
  
           

        def saveData(self, out_variable_3dmana, time, lev, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
                               
            """
            date, time_ind1,time_ind2, time_ind3 = MERRAgeneric().getTime(beg, end)
            
            #Setup size of saving file
            chunk_size = 5
            date_size = len(date)
            # for t,v, u, h, rh (double hour_size)
            hour_size = len(time[0][0])
            int_size = date_size//chunk_size
            res_type = (date_size*hour_size)%(chunk_size*hour_size)
            
            if (res_type > 0):
                size_type = [chunk_size*hour_size]*int_size + [res_type]
            
            else:
                size_type = [chunk_size*hour_size]*int_size           
            
            # get the data with subset of area for t
            t = MERRAgeneric().dataStuff_3d(0, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            t_total = MERRAgeneric().restruDatastuff(t)
            del t

            # get the data wtih subset of area for V
            v = MERRAgeneric().dataStuff_3d(1, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            v_total = MERRAgeneric().restruDatastuff(v)
            del v
            
            #get the data with subset of area for U
            u = MERRAgeneric().dataStuff_3d(2, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            u_total = MERRAgeneric().restruDatastuff(u)
            del u
            
            #get the data with subset of area
            h = MERRAgeneric().dataStuff_3d(3, id_lat, id_lon, out_variable_3dmana)
            #restructing the shape 
            h_total = MERRAgeneric().restruDatastuff(h)
            del h
            
            #Set up the list of output variables
            var_list = [['Temperature',"air_temperature", "K", t_total],
                        ['V_component_of_wind', "northward_wind","m s-1", v_total],
                        ['U_component_of_wind',"eastward_wind","m s-1", u_total],
                        ['H',"geopotential height","m+2 s-2", h_total]]
            
            # save nc file 
            var_low = 0
            var_up = 0
            for i in range(0, 1):
            #for i in range(0, len(size_type)):
                var = size_type[i]
                var_low = var_up
                var_up = var_low + var
                
                #set up file path and names 
                file_ncdf  = path.join(dir_data,("merra_pl-1" + "_" + (date[var_low/len(time[0][0])]) + "_" + "to" + "_" +(date[var_up/len(time[0][0]) - 1]) + ".nc"))      
                rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
                print("File Type:", rootgrp.file_format)
                rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at pressure levels'
                rootgrp.featureType = "3_Dimension"
    
                #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
                LEV = lev[0][0]
                LAT = lat[0][0]
                LON = lon[0][0]
                #dimensions
                times  = rootgrp.createDimension('times', var)
                levels = rootgrp.createDimension('levels', len(LEV))
                lats   = rootgrp.createDimension('lats', len(LAT))
                lons   = rootgrp.createDimension('lons', len(LON))
                
                #Output the results of output variables
                for x in range(0,len(var_list)):
                    out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('times', 'levels','lats', 'lons'),fill_value=9.9999999E14)       
                    out_var.standard_name = var_list[x][1]
                    out_var.units         = var_list[x][2] 
                    out_var.missing_value = 9.9999999E14
                    out_var.fmissing_value = (9.9999999E14, 'f')
                    out_var.vmax = (9.9999999E14, 'f')
                    out_var.vmin = (-9.9999999E14, 'f')   
                    out_var[:,:,:,:] = var_list[x][3][var_low:var_up,:,:,:]    #data generic name with data stored in it
                
    
                Time = rootgrp.createVariable('times', 'i4', ('times'))
                Time.standard_name = "time"
                Time.units  = "hour since " + str(datetime.strptime(beg, '%Y/%m/%d'))                 
                Time.calendar = "standard"   
                # pass the values
                Time[:] = time_ind1[var_low:var_up]                                  # 6-hourly time step (for Temp, U, V, H)                                           
                # Time[:] = nc.date2num(time_ind1[var_low:var_up],unit = Time.units,calendar = Time.calendar)
                                                        
                Level = rootgrp.createVariable('level','i4', ('levels'))
                Level.standard_name = "air_pressure"
                Level.units = "hPa"
                # pass the values
                Level[:] = lev[0][0][:]                    

                Latitudes               = rootgrp.createVariable('latitudes', 'f4',('lats'))
                Latitudes.standard_name = "latitude"
                Latitudes.units         = "degrees_north"
                Latitudes.axis          = 'Y'
                Latitudes[:]  = lat[0][0][:]                                   # pass the values of latitude

                Longitudes               = rootgrp.createVariable('longitudes', 'f4',('lons'))
                Longitudes.standard_name = "longitude"
                Longitudes.units         = "degrees_east"
                Longitudes.axis          = 'X'
                Longitudes[:] = lon[0][0][:]                                   # pass the values of longitudes
    
                #close the root group
                rootgrp.close()          

class SaveNCDF_pl_3dmasm():                                                         # for saving abstracted pressure-levels variables
        """ write output netCDF file for abstracted variables from original meteorological data 
            at pressure levels
            demension: time, level, lat, lon
            variables: time: array([   0,  180,  360,  540,  720,  900, 1080, 1260], dtype=int32, Unit: minute)
                       relative humidity(time,lev,lat,lon), 
                       time, level, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            date, time_ind2 = MERRAgeneric().getTime(beg, end)
            rh = MERRAgeneric().dataStuff_3d(0, id_lat, id_lon, out_variable_3dmasm)  
            lat, lon, lev, time = MERRAgeneric().latLon(out_variable_3dmasm, id_lat, id_lon) 
                     
        """
    
        def saveData(self, out_variable_3dmasm, time, lev, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
                               
            """
            date,time_ind1,time_ind2, time_ind3 = MERRAgeneric().getTime(beg, end)
            
            #Setup size of saving file
            chunk_size = 5
            date_size = len(date)
            # for rh
            hour_size = len(time[0][0])
            int_size = date_size//chunk_size
            res_type = (date_size*hour_size)%(chunk_size*hour_size)
            
            if (res_type > 0):
                size_type = [chunk_size*hour_size]*int_size + [res_type]
            
            else:
                size_type = [chunk_size*hour_size]*int_size           

            #get the data with subset of area for rh
            rh = MERRAgeneric().dataStuff_3d(0, id_lat, id_lon, out_variable_3dmasm)
            #restructing the shape 
            rh_total = MERRAgeneric().restruDatastuff(rh)
            del rh

            #Set up the list of output variables
            var_list = [['Relative_humidity',"relative humidity after moist", "1", rh_total]]
            
            # save nc file 
            var_low = 0
            var_up = 0
            for i in range(0, len(size_type)):
                var = size_type[i]
                var_low = var_up
                var_up = var_low + var
                
                #set up file path and names 
                file_ncdf  = path.join(dir_data,("merra_pl-2" + "_" + (date[var_low/len(time[0][0])]) + "_" + "to" + "_" +(date[var_up/len(time[0][0]) - 1]) + ".nc"))      
                rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
                print("File Type:",rootgrp.file_format)
                rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at pressure levels'
                rootgrp.featureType = "3_Dimension"
    
                #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
                LEV = lev[0][0]
                LAT = lat[0][0]
                LON = lon[0][0]
                #dimensions
                times  = rootgrp.createDimension('times', var)
                levels = rootgrp.createDimension('levels', len(LEV))
                lats   = rootgrp.createDimension('lats', len(LAT))
                lons   = rootgrp.createDimension('lons', len(LON))
                
                #Output the results of output variables
                for x in range(0,len(var_list)):
                    out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('times', 'levels','lats', 'lons'),fill_value=9.9999999E14)       
                    out_var.standard_name = var_list[x][1]
                    out_var.units         = var_list[x][2] 
                    out_var.missing_value = 9.9999999E14
                    out_var.fmissing_value = (9.9999999E14, 'f')
                    out_var.vmax = (9.9999999E14, 'f')
                    out_var.vmin = (-9.9999999E14, 'f')   
                    out_var[:,:,:,:] = var_list[x][3][var_low:var_up,:,:,:]       #data generic name with data stored in it
                
                Time = rootgrp.createVariable('times', 'i4', ('times'))
                Time.standard_name = "time"
                Time.units = "hour since" + str(datetime.strptime(beg, '%Y/%m/%d'))                  
                Time.calendar = "standard"
                # pass the values
                Time[:] = time_ind2[var_low:var_up]                                  # 3-hourly time step (for RH)                                           
                            
                Level = rootgrp.createVariable('level','i4', ('levels'))
                Level.standard_name = "air_pressure"
                Level.units = "hPa"
                # pass the values
                Level[:] = lev[0][0][:]                    

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


class MERRAsm():
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
    def getDs(self, beg, end, username, password, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_2dm
           chunk_size = 5
           ds = MERRAgeneric().download(username, password, urls, chunk_size)
        """    
        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr = MERRAgeneric().getURLs(beg, end)
        urls = urls_2dm
        ds = MERRAgeneric().download(username, password, urls, chunk_size)
        
        return ds

    
    def getVariables(self, ds):
        """Return the objected variables from the specific MERRA-2 3D datasets        
            variable = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']
            ds = MERRAgeneric.download( username, password, urls_2dm, chunk_size)
            
        """
        variable = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']
        out_variable_2dm = MERRAgeneric().Variables(variable, ds)

        return out_variable_2dm
         
    def getlatLon_2d(self, area, ds, out_variable_2dm, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)
            p1 = 7 (id_Latitude)
            p2 = 8 (id_Longitude)
            p3 = 9 (id_Time) 

        """       
        
                
        p1 = 7
        p2 = 8
        p3 = 9

        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)

        lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dm, p1, p2, p3, id_lat, id_lon)
        
        return lat, lon, time
        


class MERRAsf():
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
    
    def getDs(self, beg, end, username, password, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 3d Analyzed Meteotological 
           Fields data products
           Args:
           username = ******
           password = ******
           urls = urls_2ds
           chunk_size = 5
           ds = MERRAgeneric().download(username, password, urls, chunk_size)
        """    
        
        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr = MERRAgeneric().getURLs(beg, end)
        urls = urls_2ds
        ds = MERRAgeneric().download(username, password, urls, chunk_size)
        
        return ds
    
    def getVariables(self, ds):
        """Return the objected variables from the specific MERRA-2 2D suface flux Diagnostics data    
           
           ds = MERRAgeneric.download( username, password, urls_2ds, chunk_size)
        """        
        variable = ['PRECTOT','lat','lon','time']
        out_variable_2ds = MERRAgeneric().Variables(variable, ds)

        return out_variable_2ds

    def getlatLon_2d(self, area, ds, out_variable_2ds, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 2D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_3dmana = MERRAgeneric().getVariables(variable, ds)
            p1 = 1 (id_Latitude)
            p2 = 2 (id_Longitude)
            p3 = 3 (id_Time) 

        """       
        p1 = 1
        p2 = 2
        p3 = 3

        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
        
        lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2ds, p1, p2, p3, id_lat, id_lon)
        
        return lat, lon, time


class SaveNCDF_sa():                                  
        """ write output netCDF file for abstracted variables from original 2d meteorological Diagnostics datasetsc and suface flux Diagnostics datasets
            demension: time, lat, lon
            variables: time: array([   0,   60,  120,  180,  240,  300,  360,  420,  480,  540,  600,
                             660,  720,  780,  840,  900,  960, 1020, 1080, 1140, 1200, 1260,
                             1320, 1380], dtype=int32, Unit: minute)
                       t2m(time*lat*lon), 
                       u2m(time*lat*lon), 
                       v2m(time*lat*lon),
                       u10m(time*lat*lon)
                       v10m(time*lat*lon)
                       prectot(time*lat*lon), 
                       time, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            date, time_ind3 = MERRAgeneric().getTime(beg, end)
            t2m = dataStuff_2d(0, id_lat, id_lon, out_variable_2dm) 
            u2m = dataStuff_2d(1, id_lat, id_lon, out_variable_2dm)
            v2m = dataStuff_2d(2, id_lat, id_lon, out_variable_2dm)       
            u10m = dataStuff_2d(3,id_lat, id_lon, out_variable_2dm)        
            v10m = dataStuff_2d(4, id_lat, id_lon, out_variable_2dm) 
            prectot = dataStuff_2d(0, out_variable_2ds) 
            lat, lon,time = MERRAgeneric().latLon_2d(out_variable_2dm, id_lat, id_lon)
                     
        """
     

        def saveData(self, out_variable_2dm, out_variable_2ds, time, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
            """
            
            date, time_ind1, time_ind2, time_ind3 = MERRAgeneric().getTime(beg, end)
            
            #Setup size of saving file
            chunk_size = 5
            date_size = len(date)
            hour_size = len(time[0][0])
            int_size = date_size//chunk_size
            res_type = (date_size*hour_size)%(chunk_size*hour_size)
                        
            if (res_type > 0):
                size_type = [chunk_size*hour_size]*int_size + [res_type]
            
            else:
                size_type = [chunk_size*hour_size]*int_size           
    
            # get the data with subset of area for t2m
            t2m = MERRAgeneric().dataStuff_2d(0, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            t2m_total = MERRAgeneric().restruDatastuff(t2m)
            del t2m

            # get the data wtih subset of area for u2m
            u2m = MERRAgeneric().dataStuff_2d(1, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            u2m_total = MERRAgeneric().restruDatastuff(u2m)
            del u2m

            # get the data wtih subset of area
            v2m = MERRAgeneric().dataStuff_2d(2, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            v2m_total = MERRAgeneric().restruDatastuff(v2m)
            del v2m

            # get the data wtih subset of area
            u10m = MERRAgeneric().dataStuff_2d(3, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            u10m_total = MERRAgeneric().restruDatastuff(u10m)
            del u10m
            
            # get the data wtih subset of area
            v10m = MERRAgeneric().dataStuff_2d(4, id_lat, id_lon, out_variable_2dm)
            #restructing the shape 
            v10m_total = MERRAgeneric().restruDatastuff(v10m)
            del v10m
             
            # get the data with subset of area
            prectot = MERRAgeneric().dataStuff_2d(0, id_lat, id_lon, out_variable_2ds)
            #restructing the shape 
            prectot_total = MERRAgeneric().restruDatastuff(prectot)
            del prectot
  
            var_list = [['2_meter_air_temperature',"2_meter_air_temperature", "K", t2m_total],
                        ['2_meter_eastward_wind',"2_meter_eastward_wind","m s-1", u2m_total],
                        ['2_meter_northward_wind',"2_meter_northward_wind","m s-1", v2m_total],
                        ['10_meter_eastward_wind',"10_meter_eastward_wind","m s-1", u10m_total],
                        ['10_meter_northward_wind',"10_meter_northward_wind","m s-1", v10m_total],
                        ['total_precipitation',"total_precipitation","kg m-2 s-1", prectot_total]]
                    
            
            #save nc file
            var_low = 0
            var_up = 0
            for i in range(0, 1):
            # for i in range(0, len(size_type)):
                var = size_type[i]
                var_low = var_up
                var_up = var_low + var
    
                #set up file path and names 
                file_ncdf  = path.join(dir_data,("merra_sa" + "_" + (date[var_low/len(time[0][0])]) + "_" + "to" + "_" +(date[var_up/len(time[0][0]) - 1]) + ".nc"))
                rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
                print("File Type:", rootgrp.file_format)
                rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at surface level'
                rootgrp.featureType = "2_Dimension"
            
                #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
                TIME = time[0][0]/60           # unit conversion minutues to hours
                LAT = lat[0][0]
                LON = lon[0][0]
    
                chunk_size = 5
                #dimensions
                times  = rootgrp.createDimension('times', var)
                lats   = rootgrp.createDimension('lats', len(LAT))
                lons   = rootgrp.createDimension('lons', len(LON))
                
                #Output the results of extracted variables
                for x in range(0,len(var_list)):
                    out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('times','lats', 'lons'),fill_value=9.9999999E14)       
                    out_var.standard_name = var_list[x][1]
                    out_var.units         = var_list[x][2] 
                    out_var.missing_value = 9.9999999E14
                    out_var.fmissing_value = (9.9999999E14, 'f')
                    out_var.vmax = (9.9999999E14, 'f')
                    out_var.vmin = (-9.9999999E14, 'f')   
                    out_var[:,:,:] = var_list[x][3][var_low:var_up,:,:]        #data generic name with data stored in it
        
                Time  = rootgrp.createVariable('time', 'i4', ('times'))
                Time.standard_name = "time"
                Time.units         = "hour since" + str(datetime.strptime(beg, '%Y/%m/%d'))
                Time.calendar      = "standard"
                # pass the values
                Time[:] = time_ind3[var_low:var_up]                                               # 1-hourly time step                                            
    
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


class MERRAsr():
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
        
    """
    
    def getDs(self, beg, end, username, password, chunk_size):
        """Return the orginal datasets structured with defined chuncks form the specific MERRA-2 2d radiation Diagnostics datasets 
           Args:
           username = ******
           password = ******
        """    
        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr = MERRAgeneric().getURLs(beg, end)
        urls = urls_2dr
        ds = MERRAgeneric().download(username, password, urls, chunk_size)
        
        return ds

    
    def getVariables(self, ds):
        """Return the objected variables from the specific MERRA-2 2D radiation Diagnostics datasets        
            variable = ['SWGNT','LWGNT', 'lat','lon','time']
            urls = urls_2dr
            ds = MERRAgeneric.download( username, password, urls, chunk_size)   
        """
        variable = ['SWGNT','LWGNT', 'lat','lon','time']
        out_variable_2dr = MERRAgeneric().Variables(variable, ds)

        return out_variable_2dr
         
    def getlatLon_2d(self, area, ds, out_variable_2dr, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Time from specific MERRA-2 2D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable_2dr = MERRAgeneric().getVariables(variable, ds)
            p1 = 2 (id_Latitude)
            p2 = 3 (id_Longitude)
            p3 = 4 (id_Time) 

        """       
        id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
       
        p1 = 2
        p2 = 3
        p3 = 4

        lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dr, p1, p2, p3, id_lat, id_lon)
        
        return lat, lon, time
        
class SaveNCDF_sr():                                  
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
            date, time_ind3 = MERRAgeneric().getTime(beg, end)
            swgnt = MERRAgeneric().dataStuff_2d(0, id_lat, id_lon, out_variable_2dr) 
            lwgnt = MERRAgeneric().dataStuff_2d(1, id_lat, id_lon, out_variable_2dr)
            lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dr, id_lat, id_lon)
                     
        """

        def saveData(self, out_variable_2dr, time, lat, lon):
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
            """
            Args: 
            dir_data  = '/Users/xquan/data'  
            
            """
            date, time_ind1,time_ind2, time_ind3 = MERRAgeneric().getTime(beg, end)

            #Setup size of saving file
            chunk_size = 5
            date_size = len(date)
            hour_size = len(time[0][0])
            int_size = date_size//chunk_size
            res_type = (date_size*hour_size)%(chunk_size*hour_size)
            if (res_type > 0):
                size_type = [chunk_size*hour_size]*int_size + [res_type]
            else:
                size_type = [chunk_size*hour_size]*int_size           


            # get the data with subset of area
            swgnt = MERRAgeneric().dataStuff_2d(0, id_lat, id_lon, out_variable_2dr)
            #restructing the shape 
            swgnt_total = MERRAgeneric().restruDatastuff(swgnt)
            del swgnt
            
            # get the data with subset of area
            lwgnt = MERRAgeneric().dataStuff_2d(0, id_lat, id_lon, out_variable_2dr)
            #restructing the shape 
            lwgnt_total = MERRAgeneric().restruDatastuff(lwgnt)
            del lwgnt

            #Set up the list of output variables
            var_list = [['surface_net_downward_shortwave_flux',"surface_net_downward_shortwave_flux", "W m-2", swgnt_total],
                        ['surface_net_downward_longwave_flux',"surface_net_downward_longwave_flux","W m-2", lwgnt_total]]

            var_low = 0
            var_up = 0
            # for i in range(0, 1):
            for i in range(0, len(size_type)):
                var = size_type[i]
                var_low = var_up
                var_up = var_low + var
    
                # set up file path and names  
                file_ncdf  = path.join(dir_data,("merra_sr" + "_" + (date[var_low/len(time[0][0])]) + "_" + "to" + "_" +(date[var_up/len(time[0][0]) - 1]) + ".nc"))
                rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
                print("File Type:", rootgrp.file_format)
                rootgrp.source      = 'Merra, abstrated radiation variables from metadata at surface level'
                rootgrp.featureType = "2_Dimension"
            
                #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup 
                TIME = time[0][0]/60           # unit conversion minutues to hours
                LAT = lat[0][0]
                LON = lon[0][0]
    
                #dimensions
                times  = rootgrp.createDimension('times', var)
                lats   = rootgrp.createDimension('lats', len(LAT))
                lons   = rootgrp.createDimension('lons', len(LON))
            
                #Output the results of extracted variables
                for x in range(0,len(var_list)):
                    out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('times','lats', 'lons'),fill_value=9.9999999E14)       
                    out_var.standard_name = var_list[x][1]
                    out_var.units         = var_list[x][2] 
                    out_var.missing_value = 9.9999999E14
                    out_var.fmissing_value = (9.9999999E14, 'f')
                    out_var.vmax = (9.9999999E14, 'f')
                    out_var.vmin = (-9.9999999E14, 'f')   
                    out_var[:,:,:] = var_list[x][3][var_low:var_up,:,:]         
                                    
                Time               = rootgrp.createVariable('time', 'i4', ('times'))
                Time.standard_name = "time"
                Time.units         = "hour since" + str(datetime.strptime(beg, '%Y/%m/%d')) # needed to add the beging date into it
                Time.calendar      = "standard"
                # pass the values
                Time[:] = time_ind3[var_low:var_up]                                               # 1-hourly time step                                            
    
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
from datetime import datetime, timedelta, date
from os import path
from netCDF4 import Dataset
from dateutil.rrule import rrule, DAILY

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

# execfile(path.join(dir_src, 'merra-2.py'))

#Given wanted datatype, mydate, area, elevation

beg   = "2016/01/01"
end   = "2016/02/01" 
                                                                               
area = {'bbS':50.0, 'bbN': 70.0, 'bbW':-120.0, 'bbE': -100.0}
 
elevation = {'min' : 50, 'max' : 2000}  

chunk_size = 5

# Get merra-2 3d meteorological analysis variables at pressure levels
beg   = "2016/01/01"
end   = "2016/01/12" 

# startDay = date(2016, 01, 01)
# endDay   = date(2016, 01, 10)
startDay = datetime.strptime(beg, '%Y/%m/%d')
endDay   = datetime.strptime(end, '%Y/%m/%d')

x = 0
for dt in rrule(DAILY, dtstart = startDay, until = endDay):
        currentDay = (str(dt.strftime("%Y")) + "/" + str(dt.strftime("%m")) + "/" + str(dt.strftime("%d")))
        x += 1
        if (x == 1):
            beg = currentDay
    
        if (x == chunk_size or currentDay == endDay):
            x = 0
            end = currentDay
            
            #get merra-2 meterological varaibles at pressure levels
            print ("-----Get wanted variables from Merra-2 3d, 6-hourly Pressure-Level, Analyzed Meteorological Fields-----")
            ds_ana = MERRApl_ana().getDs(beg, end, username, password, chunk_size)
            
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds_ana)
            
            out_variable_3dmana = MERRApl_ana().getVariables(ds_ana)
            
            lat, lon, lev, time = MERRApl_ana().getlatLon_3d(area, ds_ana, out_variable_3dmana, id_lat, id_lon)
           
            # Output merra-2 meteorological analysis variable at pressure levels
            #For T, V, U, H
            SaveNCDF_pl_3dmana().saveData(out_variable_3dmana, time, lev, lat, lon)
            print ("----------------------------------------Result NO.1: Completed----------------------------------------")

            # Get merra-2 3d meteorological assimilated variables at pressure levels
            print ("-----Get wanted variables from Merra-2 3d,3-hourly,Pressure-Level,Assimilated Meteorological Fields-----")
            ds_asm = MERRApl_asm().getDs(beg, end, username, password, chunk_size)                       
            
            out_variable_3dmasm = MERRApl_asm().getVariables(ds_asm)
            
            lat, lon, lev, time = MERRApl_asm().getlatLon_3d(area, ds_asm, out_variable_3dmasm, id_lat, id_lon)
            
            # Output meteorological assimilated variable at pressure levels
            # For RH
            SaveNCDF_pl_3dmasm().saveData(out_variable_3dmasm, time, lev, lat, lon)
            print ("----------------------------------------Result NO.2: Completed----------------------------------------")
            
            #get merra-2 meterological varaibles at surface level
            # Get merra-2 2d meteorological Diagnostics variables at surface level
            print ("-----Get wanted variables from Merra-2 1-hourly,Single-level,meteorological Diagnostics -----")
            ds_2dm = MERRAsm().getDs(beg, end, username, password, chunk_size)

            out_variable_2dm = MERRAsm().getVariables(ds_2dm)

            lat, lon, time = MERRAsm().getlatLon_2d(area, ds_2dm, out_variable_2dm, id_lat, id_lon)
            
            
            # Get merra-2 2d suface flux Diagnostics variables at surface level
            print ("-----Get wanted variables from Merra-2 1-hourly, single-level, Surface Flux Diagnostics -----")
            ds_2ds = MERRAsf().getDs(beg, end, username, password, chunk_size)

            out_variable_2ds = MERRAsf().getVariables(ds_2ds)
            
            lat, lon, time = MERRAsf().getlatLon_2d(area, ds_2ds, out_variable_2ds, id_lat, id_lon)
            
            # Output marra-2 variable at surface level 
            SaveNCDF_sa().saveData(out_variable_2dm, out_variable_2ds, time, lat, lon)
            print ("----------------------------------------Result NO.3: Completed----------------------------------------")

            #get merra-2 radiation varaibles
            # Get merra-2 2d radiation variables
            print ("-----Get wanted variables from Merra-2 2d,1-Hourly,Time-Averaged,Single-Level, Radiation Diagnostics -----")
            ds_2dr = MERRAsr().getDs(beg, end, username, password, chunk_size)

            out_variable_2dr = MERRAsr().getVariables(ds_2dr)

            lat, lon, time = MERRAsr().getlatLon_2d(area, ds_2dr, out_variable_2dr, id_lat, id_lon)

            #Output merra-2 radiation variables 
            SaveNCDF_sr().saveData(out_variable_2dr, time, lat, lon)
            print ("----------------------------------------Result NO.4: Completed----------------------------------------")

    