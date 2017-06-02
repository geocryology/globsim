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
from netCDF4 import Dataset

import numpy as np
import csv
import netCDF4
import itertools



class MERRAgeneric():
    """
    Parent class for other merra classes.
    """
    
    def getURL(self, mydate, data_type):                                                                                                                                          
        """Converts datetime objects into string &
            get objected url address (2d, 3d meterological fields and radiation datasets)
            url_2dm: 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
            url_3dm: 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields 
            url_2dr: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
            
            Example:
            mydate = datetime(2016, 1, 1) 
            data_type = '2dm' # get type of url for 2d meterological data
            data_type = '3dm' # get type of url for 3d meterological data
            data_type = '2dr' # get type of url for 2d radiation data
        """

        format = ('.nc4')
        res1 = mydate.strftime("%Y/%m") 
        res2 = mydate.strftime("%Y%m%d")  
        baseurl_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/') # baseurl for 2d dataset
        baseurl_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/') # baseurl for 3d dataset
        baseurl1 = ('M2I1NXASM.5.12.4/','/MERRA2_400.inst1_2d_asm_Nx.')            # sub url of 2d meteorological data     
        baseurl2 = ('M2I3NPASM.5.12.4/','/MERRA2_400.inst3_3d_asm_Np.')            # sub url of 3d meteorological data           
        baseurl3 = ('M2T1NXRAD.5.12.4/','/MERRA2_400.tavg1_2d_rad_Nx.')            # sub url of 2d radiation data

        if data_type == '2dm':        
           url = baseurl_2d + baseurl1[0] + res1 + baseurl1[1] + res2 + format
        elif data_type == '3dm':
           url = baseurl_3d + baseurl2[0] + res1 + baseurl2[1] + res2 + format
        elif data_type == '2dr':
           url = baseurl_2d + baseurl3[0] + res1 + baseurl3[1] + res2 + format        
           # return (url)
 
    def download(self, username, password, url):
        """ Access the MERRA server by account information and defiend urls"""
        session = setup_session(username, password, check_url=url)        
        print ('=== MERRA: START ====')
        ds = open_url(url, session=session)
        print ('=== MERRA: STOP =====')
        print type(ds)
        print ds.keys()

    def getArea(self, area): 
        """Gets the specific area with given latitude and longitude
           For example: 
           area = [ 40.0, 45.0, 60.0, 65.0] 
        """
        
        # get the row lat and lon
        lat = ds.lat[:]
        lon = ds.lon[:]
        
        # get the indices of selected range of lat,lon
        id_lat = np.where((lat[:] > area[0]) & (lat[:] < area[1])) 
        id_lon = np.where((lon[:] > area[2]) & (lon[:] < area[3])) 
        
        # convert id_lat, id_lon from tuples to string
        id_lat = list(itertools.chain(*id_lat))
        id_lon = list(itertools.chain(*id_lon))    

    #?? needed to revised furtherly (updated on 2rd June)
    def getPressure(self, elevation):                     
        """Convert elevation into air pressure using barometric formula"""
        g  = 9.80665   #Gravitational acceleration [m/s2]
        R  = 8.31432   #Universal gas constant for air [N·m /(mol·K)]    
        M  = 0.0289644 #Molar mass of Earth's air [kg/mol]
        P0 = 101325    #Pressure at sea level [Pa]
        T0 = 288.15    #Temperature at sea level [K]
        #http://en.wikipedia.org/wiki/Barometric_formula
        return P0 * exp((-g * M * elevation) / (R * T0)) / 100 #[hPa] or [bar]
    
    def getPressureLevels(self, elevation):
        """Restrict list of MERRA pressure levels to be downloaded"""
        Pmax = self.getPressure(elevation['min']) + 55
        Pmin = self.getPressure(elevation['max']) - 55 
        levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 
                         800, 825, 850, 875, 900, 925, 950, 975, 1000])        # replace by MERRA pressure levels
        mask = (levs > Pmin) * (levs <= Pmax) #select
        levs = '/'.join(map(str, levs[mask]))
        return levs


    def getNCDF(self):
        return self.file_ncdf   
        

class MERRApl(MERRAgeneric):
    """Returns an object for MERRA data that has methods for querying the server  
       of Goddard Earth Sciences Data and Information Services Center.
       
    Args:
        mydate: A dictionary specifying the specific date desired as a datetime.datetime object.
              
        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg].  
              
        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels 
                   needed. Unit: [m].
        
        variable:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].
        
        directory: Directory to hold output files
              
    Example:
        from datetime import datetime
        mydate = datetime(2016, 01, 01)   
        area  = [40.0, 45.0, 60.0, 65.0]
        elevation = {'min' :    0, 
                     'max' : 8850}
        variable  = ['PHIS','RH','V','T','U']             
        directory = '/Users/xquan/data'             
        MERRApl = MERRApl(mydate, area, elevation, variable, directory) 
        MERRApl.download()
    """

    def __init__(self, date, area, elevation, variable, directory):
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'merra_pl.nc')

    
    def getVariables(self, variable, ds):
        """Return the objected variables from the specific MERRA-2 2D and 3D datasets        
           variable = ['PHIS','RH','T','V','U','lat','lon','lev','time']
        """
        var = ds.keys()
        for i in range(len(variable)):
            foundVariable = False
            if variable[i] in var:
               for j in range(len(var)):
                   if foundVariable != True:
                        if var[j] == variable[i]:
                            temp = "" + var[j]        #? why ""
                            variable[i] = ds[temp]
                            foundVariable = True
     
        # # Make sure it works loop              
        # for i in range(len(variable)):
        #    print variable[i]
        sgp  = variable[0]
        RH   = variable[1]
        airt = variable[2]
        V    = variable[3]
        U    = variable[4]
        lat  = variable[5]
        lon  = variable[6]
        lev  = variable[7]
        Time = variable[8]
        
    def saveNCDF(self, variable):
        """ write output netCDF file for abstracted variables from original meteorological data 
            at pressure levels
            demension: time, level, lat, lon
            variables: surface geopotential(time,lat, lon), relative humidity(time,lev,lat,lon), 
                       temperature(time,lev,lat,lon), U component of wind(time,lev,lat,lon), 
                       V component of wind(time,lev,lat,lon), time, level, lat, lon.
        """   
        # creat a NetCDF file for saving output variables (Dataset object, also the root group).
        rootgrp = Dataset('/Users/xquan/data/merra_pl.nc', 'w', format='NETCDF4')
        print(rootgrp.file_format)
        rootgrp.source      = 'Merra, abstrated variables from metadata at pressure levels'
        rootgrp.featureType = "3_Dimension"
        
        #dimensions
        times  = rootgrp.createDimension('times', len(Time))
        levels = rootgrp.createDimension('levels', len(lev))
        lats   = rootgrp.createDimension('lats', len(lat))
        lons   = rootgrp.createDimension('lons', len(lon))
        
        #variables
        phis = rootgrp.createVariable('PHIS', 'f4', ('times', 'lats', 'lons'))
        phis.standard_name = "surface geopotential height"
        phis.units         = "m+2 s-2" 
        rh               = rootgrp.createVariable('Relative_humidity', 'f4', ('times','levels', 'lats', 'lons'))
        rh.standard_name = "relative humidity after moist"
        rh.units       = "1"     
        t               = rootgrp.createVariable('Temperature', 'f4', ('times','levels','lats', 'lons'))
        t.standard_name = "air_temperature"
        t.units         = "K"      
        v               = rootgrp.createVariable('V_component_of_wind','f4', ('times','levels', 'lats', 'lons'))
        v.standard_name = "northward_wind"
        v.unit          = "m s-1"                      
        u               = rootgrp.createVariable('U_component_of_wind', 'f4', ('times','levels', 'lats', 'lons'))
        u.standard_name = "eastward_wind"
        u.unit          = "m s-1" 
        time               = rootgrp.createVariable('time', 'i4', ('times'))
        time.standard_name = "time"
        time.units         = "minutes since 1980-01-01 00:00:00"
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
        i = 0
        for i in range(len(variable[8])):
            phis[i,:,:] = variable[0][i,:,:]        # pass the values of surface geopotential height      
            rh[i,:,:,:] = variable[1][i,:,:,:]      # pass the values of relative humidity 
            t[i,:,:,:]  = variable[2][i,:,:,:]      # pass the values of temperature
            v[i,:,:,:]  = variable[3][i,:,:,:]      # pass the values of northward wind
            u[i,:,:,:]  = variable[4][i,:,:,:]      # pass the values of eastward wind
            time[i]     = variable[8][i]            # pass the values of time
      
        level[:] = lev[:]             # pass the values of pressure level
        latitudes[:] = lat[:]         # pass the values of latitude
        longitudes[:] = lon[:]        # pass the values of longitudes
        
        #close the root group
        rootgrp.close()
                
                  
                         

 
 
 
 
 
 
 
 
 
 
 
 
 
# class MERRAsa(MERRAgeneric)      # for downloading 2d surface variables
 
# class MERRAsf(MERRAgeneric)      # for downloading 2d radiation variables
 
# class MERRAinterp(object)        # for interpolating the abstracted variables at required pressure levels at each station
 
# etc. 
 
 
       
                   
    
        
        
        
        
        
        