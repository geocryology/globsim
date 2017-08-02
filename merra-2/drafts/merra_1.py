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
        
    
    def getDate(self, start_date, end_date):                                                                                                                                          
        """set up time range 
        date_start = "2016-01-01"
        date_end   = "2016-12-31"
        """
        
        start = datetime.strptime(date_start, '%Y-%m-%d')
        end  = datetime.strptime(date_end, "%Y-%m-%d")
        time_diff = end - start
        print time_diff.days
        # get list of wanted date
        date = [start + timedelta(days=x) for x in range(time_diff.days + 1)]

        return date
   
    def getURLs(self, date_start, date_end):                                                                                                                                          
        """ Set up urls by given range of date and type of data
            to getobjected url address (2d, 3d meterological fields and radiation datasets)
            url_2dm: 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
            url_3dm_ana: 3d,6-hourly,Instantaneous,Pressure-Level,Analyzed Meteorological Fields 
            url_3dm_asm: 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
            url_2dr: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
            
            Args:
            date_start = "2016-01-01"
            date_end   = "2016-12-31"

            urls_3dmana: get type of url for 3d Analyzed Meteorological Fields data
            urls_3dmasm: get type of url for 3d Assimilated Meteorological Fields data
            urls_2dm:    get type of url for 2d meterological data
            urls_2ds:    get type of url for 2d surface flux data

            
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
        start = datetime.strptime(date_start, '%Y-%m-%d')
        end  = datetime.strptime(date_end, "%Y-%m-%d")

        #Setup the based string of dates for urls 
        res1 = [d.strftime("%Y/%m") for d in pandas.date_range(start,end)]
        res2 = [d.strftime("%Y%m%d") for d in pandas.date_range(start,end)]        
        
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
                 urls_2dr.append(baseurl_2d + baseurl_2dr[0] + res1[i] + baseurl_2dr[1] + res2[i] + format)     # urls of 2d radiation Diagnostics datasets with temporal subset 
                 urls_2ds.append(baseurl_2d + baseurl_2ds[0] + res1[i] + baseurl_2ds[1] + res2[i] + format)     # urls of 2d suface flux Diagnostics datasets with temporal subset  

        return urls_3dmana, urls_3dmasm, urls_2dm, urls_2dr, urls_2ds 
 
    def download(self, username, password, urls, size):
        """ Access the MERRA server by account information and defiend urls
            Args:
            username = "xxxxxx"
            password = "xxxxxx"
            urls = urls_3dmana,urls_3dmasm,urls_2dm,urls_2dr,urls_2ds: a full list of urls by specific date range of wanted types of dataset
            size: the wanted size of urls list for each chunk
                              
        """
        size = 5
        urls_chunks = [urls[x:x+size] for x in xrange(0, len(urls), size)]      

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
        print ('Type:', type(ds), 'Length:', len(ds))      
        
        return ds        
    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 datasets        
           variable = ['T','V','U','H','lat','lon','lev','time']                         # for extrating from 3d Analyzed Meteorological Fields datasets 
                    = ['RH','lat','lon','lev','time']                                    # for extracting from 3d Assimilated Meteorological Fields datasets
                    = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']  # for extracting from 2d meteorological Diagnostics datasets
                    = ['PRECTOT','lat','lon','time']                                     # for extracing from radiation Diagnostics datasets 
                    = ['SWGNT','LWGNT','lat','lon','time']                               # for extracing from 2d suface flux Diagnostics datasets
                    
           ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm, urls_2dm, urls_2dr, urls_2ds, size)
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

        print "Length of Out_Variable", len(out_variable)
        
        return out_variable
    
    def getArea(self, area, ds): 
        """Gets the specific area with given latitude and longitude
           For example: 
           area = {'lat':[40.0, 45.0], 'lon':[60.0, 65.0]}
           ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm, urls_2dm, urls_2dr, urls_2ds, size)
        """
             
        # pass the value of individual row lat and lon to Lat and Lon for the area subset
        Lat = ds[0][0].lat[:]
        Lon = ds[0][0].lon[:]
                        
        # get the indices of selected range of Lat,Lon
        id_lon = np.where((Lon[:] > area['lon'][0]) & (Lon[:] < area['lon'][1])) 
        id_lat = np.where((Lat[:] > area['lat'][0]) & (Lat[:] < area['lat'][1])) 
       
        # convert id_lat, id_lon from tuples to string
        id_lat = list(itertools.chain(*id_lat))
        id_lon = list(itertools.chain(*id_lon))   
        
        return id_lat, id_lon

    def latLon(self, out_variable, p1, p2, p3, p4, id_lat, id_lon): 
        """
        Get Latitude, Longitude, Levels, and Time
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


    def dataStuff(self, position, lat, lon, lev, time, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
         Args:
           lat, lon, lev, time = latLon()
           id_lat, id_lon =  MERRAgeneric().getArea(area, ds) 
           temp = dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v = dataStuff(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u = dataStuff(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
           h = dataStuff(3, lat, lon, lev, time, id_lat, id_lon, out_variable)        
           position: [t=0,v=1,u=2, h=3]
                     [rh=0]
                     [t2m=0, u2m=1, v2m=2, u10m=3, v10m=4]
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
                data_area[i][j] = data[i][j][:,:,id_lat,:]
            for j in range(0, len(data_area[i])):
                data_area[i][j] = data_area[i][j][:,:,:,id_lon]
            del data 

        return data_area
        
class MERRApl_ana(MERRAgeneric):
    """Returns variables from downloaded MERRA 3d Analyzed Meteorological Fields datasets  
       which are abstracted with specific temporal and spatial range  
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
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
        out_variable = MERRAgeneric().getVariables(variable, ds)

        return out_variable 
        
    def getlatLon (out_variable, p1, p2, p3, p4, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable = MERRAgeneric().getVariables(variable, ds)
            p1 = 4 (id_Latitude)
            p2 = 5 (id_Longitude)
            p3 = 6 (id_Level)
            p4 = 7 (id_Time)  

        """       
        p1 = 4
        p2 = 5
        p3 = 6
        p4 = 7

        lat, lon, lev, time = MERRAgeneric().latLon(out_variable_3dmana, id_lat, id_lon)
        
        return lat, lon, lev, time

    def getdataStuff(position, lat, lon, lev, time, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
        Args:
           lat, lon, lev, time = latLon()
           temp = dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v = dataStuff(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u = dataStuff(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
           h = dataStuff(3, lat, lon, lev, time, id_lat, id_lon, out_variable)        
           position: [t=0,v=1,u=2,h=3]
        """     
        temp = dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
        v = dataStuff(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
        u = dataStuff(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
        h = dataStuff(3, lat, lon, lev, time, id_lat, id_lon, out_variable)        

        return temp, v, u, h
    


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
            ds = MERRAgeneric.download( username, password, urls_3dasm, size)
        """
        out_variable = MERRAgeneric.getVariables(variable, ds)

        return out_variable
        
    def getlatLon (out_variable, p1, p2, p3, p4, id_lat, id_lon):
        """
        Return the objected Latitude, Longitude, Levels, Time from specific MERRA-2 3D datasets
        Args:
            id_lat, id_lon =  MERRAgeneric().getArea(area, ds)
            out_variable = MERRAgeneric().getVariables(variable, ds)
            p1 = 1 (id_Latitude)
            p2 = 2 (id_Longitude)
            p3 = 3 (id_Level)
            p4 = 4 (id_Time)  

        """       
        p1 = 1
        p2 = 2
        p3 = 3
        p4 = 4

        lat, lon, lev, time = MERRAgeneric().latLon(out_variable_3dmasm, id_lat, id_lon)
        
        return lat, lon, lev, time

    def getdataStuff(position, lat, lon, lev, time, id_lat, id_lon, out_variable):  
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area 
        Args:
           lat, lon, lev, time = latLon()
           temp = dataStuff(0, lat, lon, lev, time, id_lat, id_lon, out_variable)
           v = dataStuff(1, lat, lon, lev, time, id_lat, id_lon, out_variable)
           u = dataStuff(2, lat, lon, lev, time, id_lat, id_lon, out_variable)       
           h = dataStuff(3, lat, lon, lev, time, id_lat, id_lon, out_variable)        
           position: [t=0,v=1,u=2,h=3]
        """     
        
        return temp, v, u, h
        


class saveNCDF_pl(self, file_ncdf, t, v, u, h, rh, time, lev, lat, lon):        # for saving abstracted pressure-levels variables
        """ write output netCDF file for abstracted variables from original meteorological data 
            at pressure levels
            demension: time, level, lat, lon
            variables: temperature(time,lev,lat,lon), U component of wind(time,lev,lat,lon), 
                       V component of wind(time,lev,lat,lon),
                       geopotential heights(time,lev,lat,lon)
                       relative humidity(time,lev,lat,lon), time, level, lat, lon.
            Args: 
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_pl.nc')
            T = restrictArea.T
            V = restrictArea.V
            U = restrictArea.U
            H = restrictArea.H
            Time = restrictArea.Time
            Lev = restrictArea.Lev
            Lat = restrictArea.Lat
            Lon = restrictArea.Lon                    
        """   
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4')
        print(rootgrp.file_format)
        rootgrp.source      = 'Merra, abstrated variables from metadata at pressure levels'
        rootgrp.featureType = "3_Dimension"
        
        #Arrange the format of dimensions for time, levels, latitude and longitude
        
        # For time : convert the time (len(time)*8) to series (one demension)
        # Time = 
        Lev = lev[0]
        Lat = lat[0]
        Lon = lon[0]

        #dimensions
        times  = rootgrp.createDimension('times', len(Time))
#       times = rootgrp.createDimension('times', None)
        levels = rootgrp.createDimension('levels', len(Lev))
        lats   = rootgrp.createDimension('lats', len(Lat))
        lons   = rootgrp.createDimension('lons', len(Lon))
        
        #variables             
        t               = rootgrp.createVariable('Temperature', 'f4', ('times','levels','lats', 'lons'),fill_value=9.9999999E14)
        t.standard_name = "air_temperature"
        t.units         = "K" 
        t.missing_value = 9.9999999E14
        t.fmissing_value = (9.9999999E14, 'f')
        t.vmax = (9.9999999E14, 'f')
        t.vmin = (-9.9999999E14, 'f')
    
        v               = rootgrp.createVariable('V_component_of_wind','f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
        v.standard_name = "northward_wind"
        v.unit          = "m s-1" 
        v.missing_value = 9.9999999E14
        v.fmissing_value = 9.9999999E14, 'f'
        v.vmax = 9.9999999E14, 'f'
        v.vmin = -9.9999999E14, 'f'
                            
        u               = rootgrp.createVariable('U_component_of_wind', 'f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
        u.standard_name = "eastward_wind"
        u.unit          = "m s-1" 
        u.missing_value = 9.9999999E14
        u.fmissing_value = (9.9999999E14, 'f')
        u.vmax = (9.9999999E14, 'f')
        u.vmin = (-9.9999999E14, 'f')
        
        h = rootgrp.createVariable('H', 'f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
        h.standard_name = "geopotential height"
        h.units         = "m+2 s-2" 
        h.missing_value = 9.9999999E14
        h.fmissing_value = (9.9999999E14,'f')
        h.vmax = (9.9999999E14, 'f')
        h.vmin = (-9.9999999E14, 'f')

        rh               = rootgrp.createVariable('Relative_humidity', 'f4', ('times','levels', 'lats', 'lons'),fill_value=9.9999999E14)
        rh.standard_name = "relative humidity after moist"
        rh.units       = "1"  
        rh.missing_value = 9.9999999E14
        rh.fmissing_value = (9.9999999E14, 'f')
        rh.vmax = (9.9999999E14, 'f')
        rh.vmin = (-9.9999999E14, 'f')
        
        time               = rootgrp.createVariable('time', 'i4', ('times'))
        time.standard_name = "time"
        time.units         = "minutes since" + date.strftime("%Y-%m-%d 00:00:00")
        time.calendar      = "standard"
 
        level               = rootgrp.createVariable('level','i4', ('levels'))
        level.standard_name = "air_pressure"
        level.units         = "hPa"

        latitudes               = rootgrp.createVariable('latitudes', 'f4',('lats'))
        latitudes.standard_name = "latitude"
        latitudes.units         = "degrees_north"
        latitudes.axis          = 'Y'

        longitudes               = rootgrp.createVariable('longitudes', 'f4',('lons'))
        longitudes.standard_name = "longitude"
        longitudes.units         = "degrees_east"
        longitudes.axis          = 'X'
        
        # assign values
        t[:,:,:,:]    = T[:,:,:,:]                # pass the values of temperature
        v[:,:,:,:]    = V[:,:,:,:]                # pass the values of northward wind
        u[:,:,:,:]    = U[:,:,:,:]                # pass the values of eastward wind
        h[:,:,:,:]    = H[:,:,:,:]               # pass the values of geopotential height
        rh[:,:,:,:]   = RH[:,:,:,:]               # pass the values of relative humidity 
      

        time[:]       = Time[:]
        level[:]      = Lev[:]                    # pass the values of pressure level
        latitudes[:]  = Lat[:]                    # pass the values of latitude
        longitudes[:] = Lon[:]                    # pass the values of longitudes
        
        #close the root group
        rootgrp.close()          


class MERRAsm(MERRAgeneric):
    """Returns variables from downloaded MERRA 2d meteorological Diagnostics data, 
       which are abstracted with specific temporal and spatial range        
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
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

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 3D datasets        
            variable = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']
            ds = MERRAgeneric.download( username, password, urls_2dm, size)
        """
        out_variable = MERRAgeneric.getVariables(variable, ds)

        return out_variable 


class MERRAsf(MERRAgeneric):
    """Returns variables from downloaded MERRA 2d suface flux Diagnostics data, 
       which are abstracted with specific temporal and spatial range        
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
                      
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['PRECTOT','lat','lon','time'].

              
    """
    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D      
           variable = ['PRECTOT','lat','lon','time']
           ds = MERRAgeneric.download( username, password, urls_2ds, size)
        """        
        out_variable = MERRAgeneric.getVariables(variable, ds)

    return out_variable



# class saveNCDF_sa():   # for saving abstracted surface variables


class MERRAsr(MERRAgeneric):
    """Returns variables from downloaded MERRA 2d radiation Diagnostics data, 
       which are abstracted with specific temporal and spatial range        
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['SWGNT','LWGNT''lat','lon','time'].
                   SWGNT:surface net downward shortwave flux(time*lat*lon)
                   LWGNT:surface net downward longwave flux(time*lat*lon)
        
        directory: Directory to hold output files
              
    Example:
        from datetime import datetime
        date = datetime(2016, 01, 01)   
        area = {'lat':[40.0,45.0], 'lon':[60.0,65.0]}
        elevation = {'min' :    0, 
                     'max' : 8850}
        dir_data = '/Users/xquan/data'             
        MERRAsr = MERRAsr(mydate, area, elevation, variable, dir_data) 
        MERRAsr.download()
    """

    def __init__(self, date, area, elevation, dir_data):
        self.date     = date
        self.area       = area
        self.elevation  = elevation
        self.dir_data  = dir_data
        self.file_ncdf  = path.join(self.dir_data,'merra_sa.nc') # needed to change later

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D and 3D datasets        
           variable = ['SWGNT','LWGNT','lat','lon','time']
           ds =  MERRAgeneric.download( username, password, urls_2dr, size)
        """
     
        out_variable = MERRAgeneric.getVariables(variable, ds)
       
    return out_variable 


# class saveNCDF_sr():  # for saving abstracted radiation variables
 



# class MERRAinterp(object)        # for interpolating the abstracted variables at required pressure levels at each station
 

 
 
       
                   
    
        
        
        
        
        
        