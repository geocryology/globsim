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
                                                                               
url3 = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4'
       '/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4')                      # 3d,6-hourly,Instantaneous,Pressure-Level,Analyzed Meteorological Fields                                                                                                                                                                                                                                         
                                                                               
# url4 = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXRAD.5.12.4'
#        '/2016/02/MERRA2_400.tavg1_2d_rad_Nx.20160202.nc4')                    # 2d, 1-Hourly, Time-Averaged, Single-Level, Assimilation, Radiation Diagnostics 
                                                                               

#==============================================================================

from pydap.client import open_url
from pydap.cas.urs import setup_session
from datetime import datetime, timedelta
from os import path
from netCDF4 import Dataset

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
   
    def getURLs(self, date_start, date_end, data_type):                                                                                                                                          
        """ Set up urls by given range of date and type of data
            to getobjected url address (2d, 3d meterological fields and radiation datasets)
            url_2dm: 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
            url_3dm_ana: 3d,6-hourly,Instantaneous,Pressure-Level,Analyzed Meteorological Fields 
            url_3dm_asm: 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
            url_2dr: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
            
            Args:
            date_start = "2016-01-01"
            date_end   = "2016-12-31"

            data_type = '3dm_ana' # get type of url for 3d Analyzed Meteorological Fields data
            data_type = '3dm_asm' # get type of url for 3d Assimilated Meteorological Fields dara
            data_type = '2dm' # get type of url for 2d meterological data
            data_type = '2dr' # get type of url for 2d radiation data
            data_type = '2ds' # get type of url for 2d surface flux data

            
        """
        #setup the based url strings    
        baseurl_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/') # baseurl for 2d dataset
        baseurl_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/') # baseurl for 3d dataset
 
        baseurl_3dn = ('M2I6NPANA.5.12.4/','/MERRA2_400.inst6_3d_ana_Np.')            # sub url of 3d Analyzed Meteorological Fields data           
        baseurl_3da = ('M2I3NPASM.5.12.4/','/MERRA2_400.inst3_3d_asm_Np.')            # sub url of 3d Assimilated Meteorological Fields data
        baseurl_2dm = ('M2I1NXASM.5.12.4/','/MERRA2_400.inst1_2d_asm_Nx.')            # sub url of 2d meteorological Diagnostics data     
        baseurl_2dr = ('M2T1NXRAD.5.12.4/','/MERRA2_400.tavg1_2d_rad_Nx.')            # sub url of 2d radiation Diagnostics data
        baseurl_2ds = ('M2T1NXFLX.5.12.4/','/MERRA2_400.tavg1_2d_flx_Nx.')            # sub url of 2d suface flux Diagnostics data                                                              # sub url of 2d surface flux Diagnostics data
        format = ('.nc4')

        #Setup the start and end of dates
        start = datetime.strptime(date_start, '%Y-%m-%d')
        end  = datetime.strptime(date_end, "%Y-%m-%d")

        #Setup the based string of dates for urls 
        res1 = [d.strftime("%Y/%m") for d in pandas.date_range(start,end)]
        res2 = [d.strftime("%Y%m%d") for d in pandas.date_range(start,end)]        
        
       # get the urls list
        urls = []
        for i in range(0,len(res1)):
              if data_type == '3dm_ana':
                 urls.append(baseurl_3d + baseurl_3dn[0] + res1[i] + baseurl_3dn[1] + res2[i] + format)              
              elif data_type == '3dm_asm':
                 urls.append(baseurl_3d + baseurl_3da[0] + res1[i] + baseurl_3da[1] + res2[i] + format)              
              elif data_type == '2dm':        
                 urls.append(baseurl_2d + baseurl_2dm[0] + res1[i] + baseurl_2dm[1] + res2[i] + format)
              elif data_type == '2dr':
                 urls.append(baseurl_2d + baseurl_2dr[0] + res1[i] + baseurl_2dr[1] + res2[i] + format)
              elif data_type == '2ds':
                 urls.append(baseurl_2d + baseurl_2ds[0] + res1[i] + baseurl_2ds[1] + res2[i] + format)      
        return urls
 
    def download(self, username, password, urls):
        """ Access the MERRA server by account information and defiend urls
            Args:
            username = "xxxxxx"
            password = "xxxxxx"                  
        """
        
        print ('===== MERRA: START ======')
        print ('TIME TO GET A COFFEE')        
        ds = {}
        for i in range(0, len(urls)): 
            session = setup_session(username, password, check_url=urls[i])        
            ds[i] = open_url(urls[i], session=session) 
        print ('===== MERRA: STOP =======')
        print type(ds)
        print ds[0].keys
        
        return ds        

    def getArea(self, area, ds): 
        """Gets the specific area with given latitude and longitude
           For example: 
           area = [ 40.0, 45.0, 60.0, 65.0]
           ds: original datasets from def download()
        """
        
        # get the row lat and lon      
        lat = {}
        lon = {}      
        for i in range(len(ds)):
            lat[i] = ds[i].lat[:]
            lon[i] = ds[i].lon[:]
             
        Lat = lat[0]
        Lon = lon[0]
                        
        # get the indices of selected range of lat,lon
        id_lat = np.where((Lat[:] > area[0]) & (Lat[:] < area[1])) 
        id_lon = np.where((Lon[:] > area[2]) & (Lon[:] < area[3])) 
        
        # convert id_lat, id_lon from tuples to string
        id_lat = list(itertools.chain(*id_lat))
        id_lon = list(itertools.chain(*id_lon))   
        return id_lat, id_lon

    def getPressure(self, elevation):                                          #
        """Convert elevation into air pressure using barometric formula"""
        g  = 9.80665   #Gravitational acceleration [m/s2]
        R  = 8.31432   #Universal gas constant for air [N·m /(mol·K)]    
        M  = 0.0289644 #Molar mass of Earth's air [kg/mol]
        P0 = 101325    #Pressure at sea level [Pa]
        T0 = 288.15    #Temperature at sea level [K]
        #http://en.wikipedia.org/wiki/Barometric_formula
        return P0 * exp((-g * M * elevation) / (R * T0)) / 100 #[hPa] or [bar]
    
    def getPressureLevels(self, elevation):                                    #
        """Restrict list of MERRA pressure levels to be downloaded"""
        Pmax = self.getPressure(elevation['min']) + 55
        Pmin = self.getPressure(elevation['max']) - 55 
        levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775,
                         800, 825, 850, 875, 900, 925, 950, 975, 1000])     
        mask = (levs > Pmin) * (levs <= Pmax) #select
        levs = '/'.join(map(str, levs[mask]))
        return levs


    def getNCDF(self):                                                         #
        return self.file_ncdf   
        
class MERRApl_ana(MERRAgeneric):
    """Returns an object for MERRA data that has methods for querying the server  
       of Goddard Earth Sciences Data and Information Services Center.
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['T','V','U','H','lat','lon','lev','time']
        
        directory: Directory to hold output files
              
    Example:
        from datetime import datetime
        date = datetime(2016, 01, 01)   
        area  = [40.0, 45.0, 60.0, 65.0]
        elevation = {'min' :    0, 
                     'max' : 8850}
        dir_data = '/Users/xquan/data'             
        MERRApl = MERRApl(mydate, area, elevation, variable, dir_data) 
        MERRApl.download()
    """

    def __init__(self, date, area, elevation, dir_data):
        self.date     = date
        self.area       = area
        self.elevation  = elevation
        self.dir_data  = dir_data
        self.file_ncdf  = path.join(self.dir_data,'merra_pl_ana.nc')               # Need to set up the date string into the output .nc file with range of date

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D and 3D datasets        
           variable = ['T','V','U','H','lat','lon','lev','time']
        """
        
        out_variable = {}
        for i in range(0, len(ds)):
            print "run", i
            outputVar = []
            for x in range(0,len(variable)):
                outputVar.append(variable[x])
            
            var = ds[i].keys()
            for j in range(len(outputVar)):
                foundVariable = False
                if outputVar[j] in var:
                    for l in range(len(var)):
                        if foundVariable != True:
                            if var[l] == outputVar[j]:
                                temp = "" + var[l]
                                outputVar[j] = ds[i][temp]
                                foundVariable = True
            out_variable[i] = outputVar
        return out_variable 

    def restrictArea(self, id_lat, id_lon, out_variable):
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area          
        """

        # pass the values
        # For Air Temperature   
        T = {}
        for i in range(0, len(out_variable)):
            print "run", i
            T[i] = out_variable[i][0][:]
       # For Wind Component 
        V = {}
        U = {}
        for i in range(0, len(out_variable)):
            print "run", i    
            V[i] = out_variable[i][1][:]
            U[i] = out_variable[i][2][:]
       # For Geopotential Height
        for i in range(0, len(out_variable)):
            print "run", i
        
        # Latitude, Longitude, Levels, and Time
        Lat  = {}
        Lon  = {}
        lev  = {}
        time = {}
        for i in range(0, len(out_variable)):
            print "run", i
            Lat[i]   = out_variable[i][4][:]
            Lon[i]   = out_variable[i][5][:]
            lev[i]   = out_variable[i][6][:]
            time[i]  = out_variable[i][7][:]

        # Restrict the area
        # For Air Temperature
        t = {}
        for i in range(0, len(T)):
            print "run", i
            t[i] = T[i][:,:,id_lat,:]
        for i in range(0, len(t)):
            t[i] = t[i][:,:,:,id_lon]

        # For Wind Component V
        v = {}
        for i in range(0, len(V)):
            print "run", i
            v[i] = V[i][:,:,id_lat,:]
        for i in range(0, len(v)):
            v[i] = v[i][:,:,:,id_lon]
        
        # For Wind Componnet U
        u = {}
        for i in range(0, len(U)):
            print "run", i
            u[i] = U[i][:,:,id_lat,:]
        for i in range(0, len(u)):
            u[i] =u[i][:,:,:,id_lon]

        # For Geopotential Height
        h = {}
        for i in range(0, len(H)):
            print "run", i
            h[i] = H[i][:,:,id_lat,:]
        for i in range(0, len(h)):
            h[i] = h[i][:,:,:,id_lon]

        #For Latitude and Longitude 
        lat = {}
        lon = {}        
        for i in range(len(Lat)):
            lat[i] = Lat[i][id_lat]
            lon[i] = Lon[i][id_lon]
        return t, v, u, h, lat, lon, lev, time

class MERRApl_asm(MERRAgeneric):
    """Returns an object for MERRA data that has methods for querying the server  
       of Goddard Earth Sciences Data and Information Services Center.
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
        
        variable:  List of variable(s) to download that can include one, several
        
        directory: Directory to hold output files
              
    Example:
    """

    def __init__(self, date, area, elevation, dir_data):
        self.date     = date
        self.area       = area
        self.elevation  = elevation
        self.dir_data  = dir_data
        self.file_ncdf  = path.join(self.dir_data,'merra_pl_asm.nc')               # Need to set up the date string into the output .nc file with range of date

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D and 3D datasets        
        """
        
        out_variable = {}
        for i in range(0, len(ds)):
            print "run", i
            outputVar = []
            for x in range(0,len(variable)):
                outputVar.append(variable[x])
            
            var = ds[i].keys()
            for j in range(len(outputVar)):
                foundVariable = False
                if outputVar[j] in var:
                    for l in range(len(var)):
                        if foundVariable != True:
                            if var[l] == outputVar[j]:
                                temp = "" + var[l]
                                outputVar[j] = ds[i][temp]
                                foundVariable = True
            out_variable[i] = outputVar
        return out_variable 

    def restrictArea(self, id_lat, id_lon, out_variable):
        """Define the outputs ones &
           pass the values of abstrated variables to the output ones and 
           restrict the area          
        """
        # pass the values
        RH   = {}
        for i in range(0, len(out_variable)):
            print "run", i    

        # Restrict the area
        # For Relative Humidity
        rh = {}
        for i in range(0, len(RH)):
            print "run", i
            rh[i] = RH[i][:,:,id_lat,:]
        for i in range(0, len(rh)):
            rh[i] = rh[i][:,:,:,id_lon]

        #For Latitude and Longitude 
        lat = {}
        lon = {}        
        for i in range(len(Lat)):
            lat[i] = Lat[i][id_lat]
            lon[i] = Lon[i][id_lon]
        
        """ write output netCDF file for abstracted variables from original meteorological data 
            at pressure levels
            demension: time, level, lat, lon
            dir_data  = '/Users/xquan/data'  
            file_ncdf  = path.join(dir_data,'merra_pl.nc')
            T = restrictArea.T
            V = restrictArea.V
            U = restrictArea.U
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
        rh[:,:,:,:]   = RH[:,:,:,:]               # pass the values of relative humidity 
        time[:]       = Time[:]
        level[:]      = Lev[:]                    # pass the values of pressure level
        latitudes[:]  = Lat[:]                    # pass the values of latitude
        longitudes[:] = Lon[:]                    # pass the values of longitudes
        
        #close the root group
        rootgrp.close()          

class MERRAsa(MERRAgeneric):
    """Returns an object for MERRA data that has methods for querying the server  
       of Goddard Earth Sciences Data and Information Services Center.
       
    Args:
        date: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg],to get 
              the indies of defined latitudes and longitudes.  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['T2M','U2M','V2M','U10M','V10M','lon','lev','time'].
                   T2M: 2-meter_air_temperature
                   U2M: 2-meter_eastward_wind
                   V2M: 2-meter_northward_wind
                   U10M: 10-meter_eastward_wind
                   V10M: 10-meter_northward_wind          
        
        directory: Directory to hold output files
              
    Example:
        from datetime import datetime
        date = datetime(2016, 01, 01)   
        area  = [40.0, 45.0, 60.0, 65.0]
        elevation = {'min' :    0, 
                     'max' : 8850}
        dir_data = '/Users/xquan/data'             
        MERRAsa = MERRApl(mydate, area, elevation, variable, dir_data) 
        MERRAsa.download()
    """

    def __init__(self, date, area, elevation, dir_data):
        self.date     = date
        self.area       = area
        self.elevation  = elevation
        self.dir_data  = dir_data
        self.file_ncdf  = path.join(self.dir_data,'merra_sa.nc')

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D and 3D datasets        
           variable = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lon','lev','time']
        """
        
 
# class MERRAsr(MERRAgeneric)      # for downloading 2d radiation variables

# class MERRAsf(MERRAgeneric)      # for downloading 2d surface flux variables
 
# class MERRAinterp(object)        # for interpolating the abstracted variables at required pressure levels at each station
 

 
 
       
                   
    
        
        
        
        
        
        