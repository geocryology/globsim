#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright Xiaojing Quan & Stephan Gruber
# =============================================================================
# REVISION HISTORY
# 20170510 -- Initial Version Created
# 20171208 -- First Draft Completed
#
# ==============================================================================
# A scripts for downloading MERRA-2 reanalysis data:
# -- Geoppotential Height at pressure levels [time*level*lat*lon] (Unit:m) (6-hourly/day)
# -- Air Temperature at pressure levels [time*level*lat*lon] (Unit:K) (6-hourly/day)
# -- Relative Humidity at Pressure Levels[time*level*lat*lon] (Unit:1) (3-hourly/day)
# -- Easteward Wind at Pressure Levels [time*level*lat*lon] (Unit:m/s) (6-hourly/day)
# -- Northward Wind at Pressure Levels [time*level*lat*lon] (Unit:m/s) (6-hourly/day)
# -- Air Temperature at 2 Meter [time*lat*lon] (Unit:K) (1-hourly/day)
# -- Eastward Wind at 2 Meter [time*lat*lon] (Unit:K) (1-hourly/day)
# -- Northward Wind at 2 Meter [time*lat*lon] (Unit: m/s) (1-hourly/day)
# -- Eastward Wind at 10 Meter  [time*lat*lon] (Unit: m/s) (1-hourly/day)
# -- Northward Wind at 10 Meter [time*lat*lon] (Unit: m/s) (1-hourly/day)
# -- Precipitation Flux [time*lat*lon] (Unit: kg/m2/s) (1-hourly/day)
# -- Surface Incoming Shoertwave Flux [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Incoming Shortwave Flux Assuming Clear Sky [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Net Downward Longwave Flux [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Net Downward Longwave Flux Assuming Clear Sky [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Longwave Flux Emitted from Surface [time*lat*lon] (Unit:W/m2) (1-hourly/day)
#
# -- Surface Absorbed Longwave Flux [time*lat*lon] (Unit:W/m2) (1-hourly/day)
# -- Surface Absorbed Longwave Flux Assuming Clear Sky [time*lat*lon] (Unit:W/m2) (1-hourly/day)
#
# Saved as netCDF 4
# ====================HOW TO RUN THIS ==========================================
#
# (1) Register a New User in Earthdata Login:
#  https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login
#
# (2) Authorize NASA GESDISC DATA ARCHIVE in Earthdata Login
# https://disc.gsfc.nasa.gov/registration/authorizing-gesdisc-data-access-in-earthdata_login
#
# (3) Adapt the script below with:
#    - Authrized Username and Password (setup in .merrarc file),
#    - Input parameters: Date, Area, Elevation, Chunk_size, Variables, etc.
#      (setup in Globsim download parameter file )
#
# (4) Obtaining the URL addresses of the objected datasets at:
#     https://disc.sci.gsfc.nasa.gov/daac-bin/FTPSubset2.pl
#
# (5) Obtianing the mutiple datasets with spefici spacial and temporal)
#
# (6) Get all varialbes which are needed, and saved in NetCDF files
#
# ==============================================================================
# IMPORTANT Notes:

# 1. Samples of Selected URLs list:

# 3d,6-hourly,Instantaneous,Pressure-Level, Analyzed Meteorological Fields
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4')

# 3d,3-hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NPASM.5.12.4/2016/01/MERRA2_400.inst3_3d_asm_Np.20160201.nc4')

# 2d,1-hourly,Instantaneous,Single-level,Assimilation,Single-Level Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4')

# 2d,1-hourly, single-level, full horizontal resolution, Surface Flux Diagnostics
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXFLX.5.12.4/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4')

# 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics
# url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T3NPRAD.5.12.4/2016/01/MERRA2_400.tavg3_3d_rad_Np.20160102.nc4')

# 2d 1-Hourly,Time-Averaged,Single-Level,Assimilation,Single-Level Diagnostics V5.12.4
# url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXSLV.5.12.4/2016/01/MERRA2_400.tavg1_2d_slv_Nx.20160101.nc4')

# 2. Radiation Variables Processing:
# Considering the variables below as equal utilization:
#
# downwelling_shortwave_flux_in_air_assuming_clear_sky = Surface_Incoming_Shoertwave_Flux
#
# downwelling_shortwave_flux_in_air_assuming_clear_sky = Surface_Incoming_Shortwave_Flux_Assuming_Clear_Sky
#
# downwelling_longwave_flux_in_air = Longwave Flux Emitted from Surface + Surface Net Downward Longwave Flux
#
# downwelling_longwave_flux_in_air_assuming_clear_sky = Longwave Flux Emitted from Surface + Surface Net Downward Longwave Flux Assuming Clear Sky

# ==============================================================================

from __future__        import print_function

import numpy as np
import netCDF4 as nc
import itertools
import pandas
import time as tc
import sys

from pydap.client      import open_url
from pydap.cas.urs     import setup_session
from datetime          import datetime, timedelta
from os                import path, makedirs
from netCDF4           import Dataset
from dateutil.rrule    import rrule, DAILY
from math              import floor, atan2, pi

from scipy.interpolate import interp1d, griddata
from time              import sleep

from globsim.common_utils import str_encode, series_interpolate, variables_skip, GenericDownload, GenericScale, GenericInterpolate
from globsim.meteorology import LW_downward, pressure_from_elevation
from globsim.nc_elements import netcdf_base, new_scaled_netcdf


class MERRAgeneric():
    """
    Parent class for other merra classes.
    """

    def getURLs(self, date):
        """ Set up urls by given range of date and type of data to get objected
        url address

        Args:
        baseurl_2d: baseurl for 2D dataset
        baseurl_3d: baseurl for 3D dataset
        baseurl_3dn_*: sub url of 3D Pressure Levels Analyzed Meteorological Fields data
        baseurl_3da_*: sub url of 3D Pressure Levels Assimilated Meteorological Fields data
        baseurl_2dm_*: sub url of 2D Single-Level Diagnostics
        baseurl_2dr_*: sub url of 2D radiation Diagnostics data
        baseurl_2ds_*: sub url of 2D suface flux Diagnostics data
        baseurl_2dv_*: sub url of 2D Single-Level,Assimilation Single-Level Diagnostics data

        Return:
        urls_3dmana: urls of 3D Analyzed Meteorological Fields datasets
        urls_3dmasm: urls of 3D Assimilated Meteorological Fields datasets
        urls_2dm: urls of 2D meteorological Diagnostics datasets
        urls_2ds: urls of 2D suface flux Diagnostics datasets
        urls_2dr: urls of 2D radiation Diagnostics datasets
        urls_2dv: urls of 2D Assimilation Single-Level Diagnostics datasets
        url_2dc: urls of 2D single-level constant model parameters datasets
        """
        # Setup the based url strings
        baseurl_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/')
        baseurl_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/')

        baseurl_3dn = (baseurl_3d + 'M2I6NPANA.5.12.4/{YM}/MERRA2_{FN}.inst6_3d_ana_Np.{YMD}.nc4')
        baseurl_3da = (baseurl_3d + 'M2I3NPASM.5.12.4/{YM}/MERRA2_{FN}.inst3_3d_asm_Np.{YMD}.nc4')
        baseurl_2dm = (baseurl_2d + 'M2I1NXASM.5.12.4/{YM}/MERRA2_{FN}.inst1_2d_asm_Nx.{YMD}.nc4')
        baseurl_2dr = (baseurl_2d + 'M2T1NXRAD.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_rad_Nx.{YMD}.nc4')
        baseurl_2ds = (baseurl_2d + 'M2T1NXFLX.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_flx_Nx.{YMD}.nc4')
        baseurl_2dv = (baseurl_2d + 'M2T1NXSLV.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_slv_Nx.{YMD}.nc4')
        
        # build the urls list
        urls_3dmana = []
        urls_3dmasm = []
        urls_2dm = []
        urls_2dr = []
        urls_2ds = []
        urls_2dv = []

        for d in pandas.date_range(date['beg'], date['end']):
            ym = d.strftime("%Y/%m")
            ymd = d.strftime("%Y%m%d")
            fn = self.get_file_number(d.year)

            urls_3dmana.append(baseurl_3dn.format(YM=ym, FN=fn, YMD=ymd))
            urls_3dmasm.append(baseurl_3da.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dm.append(baseurl_2dm.format(YM=ym, FN=fn, YMD=ymd))
            urls_2ds.append(baseurl_2ds.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dr.append(baseurl_2dr.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dv.append(baseurl_2dv.format(YM=ym, FN=fn, YMD=ymd))

        # Setup URL for getting constant model parameters (2D, single-level, full horizontal resolution)
        url_2dc = ['https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4']

        return urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv
    
    @staticmethod
    def get_file_number(year):
        """
        The file names consist of a number and a meta data string. 
        The number changes over the years. 1980 until 1991 it is 100, 
        1992 until 2000 it is 200, 2001 until 2010 it is  300 
        and from 2011 until now it is 400.
        """
        file_number = ''
        
        if year >= 1980 and year < 1992:
            file_number = '100'
        elif year >= 1992 and year < 2001:
            file_number = '200'
        elif year >= 2001 and year < 2011:
            file_number = '300'
        elif year >= 2011:
            file_number = '400'
        else:
            raise Exception('The specified year is out of range.')
        return file_number

    def download(self, username, password, urls):  # TODO: rename to 'get_opendap_connections'
        """ Access the MERRA server by account information and defined urls
            Args:
                 username = "xxxxxx"
                 password = "xxxxxx"
                 urls = urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr
            Return:
                ds:  The dictionary of each individual original dataset with
                     children['variables']
        """

        print('================ MERRA-2 SERVER ACCESS: START ================')
        print('TIME TO GET A COFFEE')
        ds = {}
        for i in range(len(urls)):
            for delay in range(0,60):
                try:  # try to download the file
                    session = setup_session(username, password, check_url=urls[i])
                    ds[i] = open_url(urls[i], session=session)
                    break
                except Exception as e:
                    if delay < 59:
                        print("Error downloading file: {}. Trying again ({})".format(urls[i], str(delay)))
                        print("Error message : \n{}".format(e))
                        sleep(delay)
                        pass
                    else:
                        print("Error downloading file: " + urls[i] + ". Giving up.")
                        print("Error message : \n{}".format(e))
                        raise RuntimeError("==> Unsuccessful after 60 attempts.")
            # print('------COMPLETED------','URL NO.:', i+1)
            print(urls[i])
        print(ds[0].keys())
        print('================ MERRA-2 SERVER ACCESS: COMPLETED ================')
        infor = urls[0].split('/')
        print('Dataset:', infor[2], infor[3],infor[4])
        print('Type:', type(ds))
        print('Days:', len(urls))

        return ds

    def Variables(self, variable, ds):
        """Get the objected variables from the specific MERRA-2 datasets
        Args:
        variable = ['T','U','V','H','lat','lon','lev','time']
        variable = ['RH','H','lat','lon','lev','time']
        variable = ['U2M','T2M','TQL','V2M','V10M','U10M','QV2M','lat','lon','time']
        variable = ['PRECTOT','PRECTOTCOR', 'lat','lon','time']
        variable = ['SWGDN','LWGNT', 'SWGDNCLR', 'LWGNTCLR', 'LWGEM', 'lat','lon','time']
        variable = ['PHIS', 'FRLAND', 'FRLANDICE', 'lat', 'lon','time']
        ds = MERRAgeneric().download (username, password, urls_3dmana, urls_3dmasm,
             urls_2dm, urls_2dr, urls_2ds, chunk_size)
        Return:
        out_variable: The extracted variables in the order of given list
        Structure of out_variable: [chunk_size * lenghs of variables]
        """

        out_variable = {}
        for i in range(0, len(ds)):
            # print("Run", "NO.:", i+1)
            outputVar = []
            for x in range(0,len(variable)):
                outputVar.append(variable[x])

            var = list(ds[i].keys())
            for j in range(len(outputVar)):
                foundVariable = False
                if outputVar[j] in var:
                    for k in range(len(var)):
                        if not foundVariable:
                            if var[k] == outputVar[j]:
                                temp = "" + var[k]
                                outputVar[j] = ds[i][temp]
                                foundVariable = True
            out_variable[i] = outputVar

        return out_variable

    def getArea(self, area, ds):
        """Get the specific indices  of the latitudes and longitudes of given area
        An example:
            area = {'north':65.0, 'south': 60.0, 'west': -115.0, 'east': -110.0}
            ds = MERRAgeneric().download(username, password, urls_3dmana, urls_3dmasm,
                 urls_2dm, urls_2dr, urls_2ds, chunk_size)
        Return:
            id_lat: indices of latitudes
            id_lon: indices of longitudes
        """

        # pass the value of individual row lat and lon to Lat and Lon for the area subset
        Lat = ds[0].lat[:]
        Lon = ds[0].lon[:]

        # get the indices of selected range of Lat,Lon
        id_lon = np.where((Lon[:] >= area['west']) & (Lon[:] <= area['east']))
        id_lat = np.where((Lat[:] >= area['south']) & (Lat[:] <= area['north']))

        # convert id_lat, id_lon from tuples to string
        id_lon = list(itertools.chain(*id_lon))
        id_lat = list(itertools.chain(*id_lat))

        return id_lat, id_lon

    def getPressureLevels(self, elevation):
        """Restrict list of MERRA-2 pressure levels to be download
        
        elevation : dict
            dictionary with keys 'max' and 'min' from configuration
        """
        Pmax = pressure_from_elevation(elevation['min']) + 55
        Pmin = pressure_from_elevation(elevation['max']) - 55

        levs = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750,
                         725, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250,
                         200, 150, 100, 70, 50, 40, 30, 20, 10, 7.0, 5.0, 4.0,
                         3.0, 2.0, 1.0, 0.7, 0.5, 0.4, 0.3, 0.1])

        # Get the indices of selected range of elevation
        id_lev = np.where((levs >= Pmin) & (levs <= Pmax))
        id_lev = list(itertools.chain(*id_lev))

        return id_lev

    def latLon_3d(self, out_variable, id_lat, id_lon, id_lev):
        """Get Latitude, Longitude, Levels, and Time for datasets at the
        pressure levels

           Args:
           out_variable = MERRAgeneric().getVariables(variable, ds)
           id_lat, id_lon = MERRAgeneric().getArea(area, ds)
           -4 = id_Latitude
           -3 = id_Longitude
           -2 = id_Level
           -1 = id_Time

           Return:
           lat, lon, lev, time: the extracted latitudes, longitudes, levels
           Structure of lat: chunk_size * lengths of extracted latitudes
           Structure of lon: chunk_size * lengths of extracted longitudes
           Structure of lev: chunk_size * lengths of extracted levels
           Structure of time: chunk_size * total number of hours per day
        """

        Lat  = {}
        Lon  = {}
        Lev  = {}
        time = {}
        for i in range(0, len(out_variable)):
            for delay in range(0,60):
                try:  # try to obtain the data of Lat, Lon, Lev, time each by each
                    # print("run", "NO.:", i+1)
                    Lat[i]  = out_variable[i][-4][:]
                    Lon[i]  = out_variable[i][-3][:]
                    Lev[i]  = out_variable[i][-2][:]
                    time[i] = out_variable[i][-1][:]
                    break
                except Exception:
                    if delay < 59:
                        print("Error downloading data: " + ". Trying again (" + str(delay) + ")")
                        sleep(delay)
                        pass
                    else:
                        print("Error downloading data: " + ". Giving up.")
                        raise RuntimeError("==> Unsuccessful after 60 attempts.")

        # For Latitude and Longitude
        lat = {}
        lon = {}
        for i in range(0, len(Lat)):
            lat[i] = Lat[i][id_lat]
            lon[i] = Lon[i][id_lon]

        # For elevation
        lev = {}
        for i in range(0, len(Lev)):
            lev[i] = Lev[i][id_lev]

        return lat, lon, lev, time

    def latLon_2d(self, out_variable, id_lat, id_lon):
        """Get Latitude, Longitude, Levels, and Time for datasets at surface
        level
           Args:
           out_variable = MERRAgeneric().getVariables(variable, ds)
           id_lat, id_lon = MERRAgeneric().getArea(area, ds)
           -3 = id_Latitude
           -2 = id_Longitude
           -1 = id_Time
           Return:
           lat, lon, lev, time: the extracted latitudes, longitudes (2D)
           Structure of lat: [chunk_size * lengths of extracted latitudes
           Structure of lon: [chunk_size * lengths of extracted longitudes]
           Structure of time: [lchunk_size * total number of hours per day]
        """
        Lat  = {}
        Lon  = {}
        time = {}
        for i in range(0, len(out_variable)):
            for delay in range(0,60):
                try:  # try to obtain the data of Lat, Lon, time each by each
                    # print("run", "NO.:", i+1)
                    Lat[i]  = out_variable[i][-3][:]
                    Lon[i]  = out_variable[i][-2][:]
                    time[i] = out_variable[i][-1][:]
                    break
                except Exception:
                    if delay < 59:
                        print("Error downloading data: " + ". Trying again (" + str(delay) + ")")
                        sleep(delay)
                        pass
                    else:
                        print("Error downloading data: " + ". Giving up.")
                        raise RuntimeError("==> Unsuccessful after 60 attempts.")

        # For Latitude and Longitude
        lat = {}
        lon = {}
        for i in range(0, len(Lat)):
            lat[i] = Lat[i][id_lat]
            lon[i] = Lon[i][id_lon]

        return lat, lon, time

    def dataStuff_3d(self, position, id_lat, id_lon, id_lev, out_variable):
        """Define the outputs ones and pass the values of extracted variables
        to the output ones and restrict the area (3D data product).
        Args:
        id_lat, id_lon = MERRAgeneric().getArea(area, ds)
        position: index of the variable in get_variables (setup automatically)
        And examples:
        For 3d Analyzed Meteorological Fields datasets:
        t = MERRAgeneric().dataStuff (position, lat, lon, lev, time,
        id_lat, id_lon, out_variable)
        u = MERRAgeneric().dataStuff(position, lat, lon, lev, time,
        id_lat, id_lon, out_variable)
        v = MERRAgeneric().dataStuff(position, lat, lon, lev, time,
        id_lat, id_lon, out_variable)
        h = MERRAgeneric().dataStuff(position, lat, lon, lev, time,
        id_lat, id_lon, out_variable)
       Return:
       data_area: the extracted individual variable at pressure levels at the
       given area
       Structure of data_area: [chunk_size *[time*level*lat*lon]]
        """
        print("Get Data")
        data = {}
        for i in range(0, len(out_variable)):
            for delay in range(0,60):
                try:  # try to obtain the data of variable each by each
                    print(f"Run Day NO.: {i + 1}")
                    data[i] = out_variable[i][position][:]
                    break
                except Exception:
                    if delay < 59:
                        print("Error downloading data:  Trying again (" + str(delay) + ")")
                        sleep(delay)
                        pass
                    else:
                        print("Error downloading data: Giving up.")
                        raise RuntimeError("==> Unsuccessful after 60 attempts.")

        # Restrict the area for data set
        print("Restrict Area and Elevation")
        data_area = {}
        for i in range(0, len(data)):
            print("Run", "Day NO.:", i + 1)
            data_area[i] = data[i][:,id_lev,:,:]

        for j in range(0, len(data_area)):
            data_area[j] = data_area[j][:,:,id_lat,:]
        for j in range(0, len(data_area)):
            data_area[j] = data_area[j][:,:,:,id_lon]

        del data

        return data_area

    def dataStuff_2d(self, position, id_lat, id_lon, out_variable):
        """Define the outputs ones and pass the values of extracted variables
        to the output ones and restrict the area (2D data products).
        Args:
        id_lat, id_lon = MERRAgeneric().getArea(area, ds)
        position: index of the variable in get_variables (setup automatically)
        An examples:
        For 2D meteorological Diagnostics datasets:
        t2m = MERRAgeneric().dataStuff_2d(position, lat, lon, time, id_lat,
              id_lon, out_variable_2dm)
        u2m = MERRAgeneric().dataStuff_2d(position, lat, lon, time, id_lat,
              id_lon, out_variable_2dm)
        v2m = MERRAgeneric().dataStuff_2d(position, lat, lon, time, id_lat,
              id_lon, out_variable_2dm)
        u10m = MERRAgeneric().dataStuff_2d(position, lat, lon, time, id_lat,
              id_lon, out_variable_2dm)
        v10m = MERRAgeneric().dataStuff_2d(position, lat, lon, time, id_lat,
               id_lon, out_variable_2dm)
        Return:
        data_area: the extracted individual variable at surface levels
        Structure of data_area: chunk_size * [time*lat*lon]
        """
        print("Get Data")
        data = {}
        for i in range(0, len(out_variable)):
            for delay in range(0,60):
                try:  # try to obtain the data of variable each by each
                    print(f"Run Day NO.: {i+1}")
                    data[i] = out_variable[i][position][:]
                    break
                except Exception:
                    if delay < 59:
                        print("Error downloading data: " + ". Trying again (" + str(delay) + ")")
                        sleep(delay)
                        pass
                    else:
                        print("Error downloading data: " + ". Giving up.")
                        raise RuntimeError("==> Unsuccessful after 60 attempts.")

        # Restrict the area for data set
        print("Restrict Area")
        data_area = {}
        for i in range(0, len(data)):
            print(f"Run Day NO.: {i+1}")
            data_area[i] = data[i][:,id_lat,:]
        for j in range(0, len(data_area)):
            data_area[j] = data_area[j][:,:,id_lon]

        del data

        return data_area

    def getTime(self, date):
        """set up date and time series for netCDF output results
        Return:
        date_ind: a string list of date in the range of begin and end
        time_ind1: the time series for 6-hours step in the range of date
        time_ind2: the time series for 3-hours step in the range of date
        time_ind3: the time series for 1-hour step in the range of date
        """
        Start = date['beg']
        End   = date['end']

        # Set up the time step
        time_step1 = '6H'
        time_step2 = '3H'
        time_step3 = '1H'

        # get extra one one more day for getting the full range of time series
        End1 = End + timedelta(days=1)

        # get the Datetime index associated with time_step
        time_ind1 = (pandas.date_range(Start, End1, freq=time_step1))[0:-1]
        time_ind2 = (pandas.date_range(Start, End1, freq=time_step2))[0:-1]
        time_ind3 = (pandas.date_range(Start, End1, freq=time_step3))[0:-1]

        # Convert datetime
        try:                         # for pandas 0.19
            time_ind1.to_datetime()
            time_ind2.to_datetime()
            time_ind3.to_datetime()
        except Exception:                       # for pandas 0.23
            time_ind1 = pandas.to_datetime(time_ind1)
            time_ind2 = pandas.to_datetime(time_ind2)
            time_ind3 = pandas.to_datetime(time_ind3)

        # get list of date series
        date_diff = End - Start
        date_ind = [Start + timedelta(days=x) for x in range(date_diff.days + 1)]
        date_ind = [d.strftime('%Y%m%d') for d in date_ind]

        return date_ind, time_ind1, time_ind2, time_ind3

    def restruDatastuff(self, data_area):
        """ Restructuring the dimension of abstracted data stuff for preparing
        to save netCDF output results.
        return:
        data_total: [len(date)*len(time/day), level, lat, lon] (3D)
                    [len(date)*len(time/day), lat, lon] (2D)
        """
        data_total = []
        for i in range(0, len(data_area)):
            for j in range(0,len(data_area[i])):
                data_total.append(data_area[i][j][:])

        data_total = np.asarray(data_total, dtype=float)

        return data_total

    def tempExtrapolate(self, t_total, h_total, elevation):
        """ Processing 1D vertical extrapolation for Air Temperature, at where
            the values are lacking  (marked by 9.9999999E14) from merra-2 3D
            Analyzed Meteorological Fields datasets.
            IMPORTANT TIP:
            To set up 'ele_max = 2500' (meter) or higher
            Reason: to make sure. get enough levels of geopotential height for
            conducting 1dinterp (linear) (2 points of values needed at least)
        """

        # restructure t_total [time*lev*lat*lon] to [lat*lon*time*lev]
        t_total = t_total[:,:,:,:].transpose((2,3,0,1))
        h_total = h_total[:,:,:,:].transpose((2,3,0,1))

        # find the value gap and conduct 1d extrapolation
        for i in range(0, len(t_total)):
            for j in range(0, len(t_total[0])):
                t_time = t_total[i][j][:]
                h_time = h_total[i][j][:]
                for k in range(0, len(t_time)):
                    t_lev = t_time[k][:]
                    h_lev = h_time[k][:]
                    id_interp = []
                    for z in range(0, len(t_lev)):
                        # find the indices of levels with missing values
                        if t_lev[z] > 99999:
                            id_interp.append(z)

                            if id_interp != []:
                                # get the levels of geopotential heights with missing values
                                lev_interp = h_lev[id_interp]
                                # pass the index of first found level with existing value to z_top
                                z_top = id_interp[-1] + 1
                                # get values at the lowest 3 levels of geopotential heights with existed values
                                lev_3p = h_lev[z_top:z_top + 3]
                                # get values at the lowest 3 levels of air temperature with existed values
                                t_3p = t_lev[z_top:z_top + 3]
                                # Using spicy.interpolate.interp1d function-------------------------
                                # Require >= 2 points of levs and t in minimum
                                if len(lev_3p) >= 2:
                                    # build linear function based on given values at lowest 3 levels of air temperature and geopotential heights
                                    f = interp1d(lev_3p, t_3p, kind='linear', fill_value='extrapolate')
                                    # use built function to calculate the values of air temperature at the found missing-values levels
                                    t_interp = f(lev_interp)
                                    # fill the calculated values into missing-values levels
                                    t_lev[id_interp] = t_interp
                                else:
                                    print('Numbers of points for extrapolation are too low (less then 2):', len(lev_3p))
                                    print('Failed to conduct extrapolation at some points in the output')
                                    print('Current ele_max =', elevation['max'])
                                    print('Higher Value of "ele_max" is needed to reset: > 2500')
                                    sys.exit(0)

                        else:
                            t_lev[z] = t_lev[z]
                        h_lev[z] = h_lev[z]

                    #assign back
                    t_time[k][:] = t_lev
                    h_time[k][:] = h_lev

                #replace the extrapolated value [time * level] to each individual cell
                t_total[i][j][:] = t_time
                h_total[i][j][:] = h_time

        #restructure back
        t_total = t_total[:,:,:,:].transpose((2,3,0,1))
        h_total = h_total[:,:,:,:].transpose((2,3,0,1))

        return t_total

    def wind_rh_Extrapolate(self, data_total):
        """Processing 1D vertical extrapolation for UV and RH, at where the
           values are lacking (marked by 9.9999999E14) from merra-2 3d Analyzed
           Meteorological Fields datasets.
           Wind (U,V) and Relative Humidity (RH) are utilized the value of at
           lowest pressure levels to the ones with value gaps
        """
        #restructure u_total,v_total [time*lev*lat*lon] to [lat*lon*time*lev]
        data_total = data_total[:,:,:,:].transpose((2,3,0,1))

        #find and fill the value gap
        for i in range(0, len(data_total)):
            for j in range(0,len(data_total[0])):
                data_time = data_total[i][j][:]
                for k in range(0,len(data_time)):
                    data_lev = data_time[k][:]
                    id_interp = []
                    for z in range(0, len(data_lev)):
                        if data_lev[z] > 99999:
                           id_interp.append(z)

                        if id_interp != []:
                            z_top = id_interp[-1] + 1
                            data_lev[id_interp] = data_lev[z_top]
                        else:
                            data_lev[z] = data_lev[z]

                    data_time[k][:] = data_lev

                #replace the interpolation value to each single pixel
                data_total[i][j][:] = data_time

        #restructure back
        data_total = data_total[:,:,:,:].transpose((2,3,0,1))

        return data_total

    def MERRA_skip(self, merralist):
        ''' To remove the extra variables from downloaded MERRA2 data'''
        for var in ['LWGEM', 'LWGNT', 'LWGNTCLR']:
            if var in merralist :
               merralist.remove(var)

    def netCDF_empty(self, ncfile_out, stations, nc_in):
        #TODO change date type from f4 to f8 for lat and lon
        '''
        Creates an empty station file to hold interpolated reults. The number of
        stations is defined by the variable stations, variables are determined by
        the variable list passed from the gridded original netCDF.

        ncfile_out: full name of the file to be created
        stations:   station list read with common_utils.StationListRead()
        variables:  variables read from netCDF handle
        lev:        list of pressure levels, empty is [] (default)
        '''
        rootgrp = netcdf_base(nc_in, ncfile_out, len(stations), None,
                              'hours since 1980-01-01 00:00:00')

        station = rootgrp["station"]
        latitude = rootgrp["latitude"]
        longitude = rootgrp["longitude"]
        height = rootgrp["height"]

        # assign station characteristics
        station[:]   = list(stations['station_number'])
        latitude[:]  = list(stations['latitude_dd'])
        longitude[:] = list(stations['longitude_dd'])
        height[:]    = list(stations['elevation_m'])

        # extra treatment for pressure level files
        try:
            lev = nc_in.variables['level'][:]
            print("== 3D: file has pressure levels")
            level           = rootgrp.createDimension('level', len(lev))
            level           = rootgrp.createVariable('level','i4',('level'))
            level.long_name = 'pressure_level'
            level.units     = 'hPa'
            level[:] = lev
        except:
            print("== 2D: file without pressure levels")
            lev = []

        #remove extra variables
        varlist_merra = [str_encode(x) for x in nc_in.variables.keys()]
        self.MERRA_skip(varlist_merra)

        # create and assign variables based on input file
        for n, var in enumerate(varlist_merra):
            if variables_skip(var):
                continue
            print("VAR: ", str_encode(var))
            # extra treatment for pressure level files
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4', ('time', 'level', 'station'))
            else:
                tmp = rootgrp.createVariable(var,'f4', ('time', 'station'))
            tmp.long_name = str_encode(nc_in.variables[var].long_name)  # for merra2
            tmp.units     = str_encode(nc_in.variables[var].units)

        #close the file
        rootgrp.close()

class SaveNCDF_pl_3dm():
        """ write output netCDF file for abstracted variables from original
            meteorological data at pressure levels
        """

        def varList(self, date, get_variables_3dmasm, get_variables_3dmana, id_lat,
                    id_lon, id_lev, out_variable_3dmasm, out_variable_3dmana, chunk_size,
                    time, lev, lat, lon, dir_data, elevation):
            """
            conduct 1D Extrapolation for T,U,V, and build the variable list for
            output netcdf file
            """
            date_ind, time_ind1,time_ind2, time_ind3 = MERRAgeneric().getTime(date)

            # Setup size of saving file
            date_size = len(date_ind)

            # for t,v, u, h
            hour_size = len(time[0])
            int_size = date_size//chunk_size
            res_type = (date_size*hour_size)%(chunk_size*hour_size)

            if (res_type > 0):
                size_type = [chunk_size*hour_size]*int_size + [res_type]

            else:
                size_type = [chunk_size*hour_size]*int_size

            # Get the variables and Set up the list of output variables
            rh_total = []
            h_total = []
            t_total = []
            u_total = []
            v_total = []

            var_out = {'RH':['relative_humidity','relative humidity','1', rh_total],
                       'H':['geopotential_height','geopotential_height', 'm', h_total],
                       'T':['air_temperature', 'air_temperature','K', t_total],
                       'U':['eastward_wind','eastward_wind_component','m/s', u_total],
                       'V':['northward_wind','northward_wind_component', 'm/s', v_total]}

            var_list = []
            # get RH
            for i in range(0, len(get_variables_3dmasm[0:-4])):
                for x in var_out.keys():
                    if x == get_variables_3dmasm[i]:
                        print("------Get Subset of Variable at Pressure Levels------", get_variables_3dmasm[i])
                        var = MERRAgeneric().dataStuff_3d(i, id_lat, id_lon, id_lev, out_variable_3dmasm)
                        # restructuring the shape
                        var_total = MERRAgeneric().restruDatastuff(var)
                        if x == 'RH':
                           rh_total = var_total
                           print("Conduct 1D Extrapolation for 'RH'")
                           var_total = MERRAgeneric().wind_rh_Extrapolate(rh_total)   #1D Extrapolation for Relative Humidity
                           rh_total = var_total

                           # Extracting RH at same time indice as geopotential height
                           # rh_total[double_time, level, lat, lon] to rh_total[even position of time, level, lat, lon]
                           rh_total = rh_total[::2,:,:,:]
                        del var
                        var_out[x][3] = rh_total
                        var_list.append([get_variables_3dmasm[i], var_out[x][0], var_out[x][1], var_out[x][2], var_out[x][3]])
            # get H,T,U,V
            for i in range(0, len(get_variables_3dmana[0:-4])):
                for x in var_out.keys():
                    if x == get_variables_3dmana[i]:
                        print("------Get Subset of Variable at Pressure Levels------", get_variables_3dmana[i])
                        var = MERRAgeneric().dataStuff_3d(i, id_lat, id_lon, id_lev, out_variable_3dmana)
                        # restructing the shape
                        var_total = MERRAgeneric().restruDatastuff(var)
                        del var
                        if x == 'H':
                           h_total = var_total
                        elif x == 'T':
                           t_total = var_total
                           print("Conduct 1D Extrapolation for 'T'")
                           var_total = MERRAgeneric().tempExtrapolate(t_total, h_total, elevation) # 1D Extrapolation for Air Temperature
                        elif x == 'U':
                           u_total = var_total
                           print("Conduct 1D Extrapolation for 'U'"                                       )
                           var_total = MERRAgeneric().wind_rh_Extrapolate(u_total)   # 1D Extrapolation for Eastward Wind
                        elif x == 'V':
                           v_total = var_total
                           print("Conduct 1D Extrapolation for 'V'")
                           var_total = MERRAgeneric().wind_rh_Extrapolate(v_total)   # 1D Extrapolation for Northward Wind

                        var_out[x][3] = var_total

                        var_list.append([get_variables_3dmana[i],var_out[x][0], var_out[x][1], var_out[x][2], var_out[x][3]])


            return var_list, time_ind1, date_ind, size_type

        def saveData(self, date, get_variables_3dmasm, get_variables_3dmana,
                     id_lat, id_lon, id_lev, out_variable_3dmasm, out_variable_3dmana,
                     chunk_size, time, lev, lat, lon, dir_data, elevation):
            """
            build NetCDF file for saving output variables (T, U, V, H, RH, lat, lon, levels, time)
            """

            # get var_list, time indices
            var_list, time_ind1, date_ind, size_type = self.varList(
                    date, get_variables_3dmasm, get_variables_3dmana, id_lat,
                    id_lon, id_lev, out_variable_3dmasm, out_variable_3dmana,
                    chunk_size, time, lev, lat, lon, dir_data, elevation)

            # save nc file
            var_low = 0
            var_up = 0
            for i in range(0, 1):
            # for i in range(0, len(size_type)):
                var = size_type[i]
                var_low = var_up
                var_up = var_low + var

                # set up file path and names
                file_ncdf  = path.join(dir_data,("merra_pl" + "_" + (date_ind[var_low // len(time[0])]) + "_" + "to" + "_" +(date_ind[var_up // len(time[0]) - 1]) + ".nc"))
                rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
                print("Saved File Type:", rootgrp.file_format)
                rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at pressure levels'
                rootgrp.featureType = "3_Dimension"

                # Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup
                LEV = lev[0]
                LAT = lat[0]
                LON = lon[0]

                # dimensions
                time  = rootgrp.createDimension('time',  None)
                level = rootgrp.createDimension('level', len(LEV))
                lats  = rootgrp.createDimension('lats',  len(LAT))
                lons  = rootgrp.createDimension('lons',  len(LON))

                # Output the results of output variables
                for x in range(0, len(var_list)):
                    out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('time', 'level', 'lats', 'lons'),fill_value=9.9999999E14)
                    out_var.standard_name = var_list[x][1]
                    out_var.long_name = var_list[x][2]
                    out_var.units         = var_list[x][3]
                    out_var.missing_value = 9.9999999E14

                    # out_var.fmissing_value = (9.9999999E14, 'f')
                    # out_var.vmin = (-9.9999999E14, 'f')
                    # out_var.vmax = (9.9999999E14, 'f')

                    out_var[:,:,:,:] = var_list[x][4][var_low:var_up,:,:,:]    # data generic name with data stored in it

                Time = rootgrp.createVariable('time', 'i4', ('time'))
                Time.standard_name = "time"
                Time.units  = "hours since 1980-1-1 00:00:00"
                Time.calendar = "gregorian"

                # pass the values
                netCDFTime = []
                for x in range(0, len(time_ind1)):
                    netCDFTime.append(nc.date2num(time_ind1[x], units = Time.units, calendar = Time.calendar))
                Time[:] = netCDFTime[var_low:var_up]

                Level = rootgrp.createVariable('level','i4', ('level'))
                Level.standard_name = "air_pressure"
                Level.long_name = "vertical level"
                Level.units = "hPa"
                Level.positive = "down"
                Level.axis = "Z"
                Level[:] = [L for L in lev[0]]

                Latitudes               = rootgrp.createVariable('latitude', 'f4',('lats'))
                Latitudes.standard_name = "latitude"
                Latitudes.units         = "degrees_north"
                Latitudes.axis          = "Y"
                Latitudes[:]  = [L for L in lat[0]]

                Longitudes               = rootgrp.createVariable('longitude', 'f4',('lons'))
                Longitudes.standard_name = "longitude"
                Longitudes.units         = "degrees_east"
                Longitudes.axis          = "X"
                Longitudes[:] = [L for L in lon[0]]

                # close the root group
                rootgrp.close()


class SaveNCDF_sa():
        """ write output netCDF file for abstracted variables from original 2D
            meteorological Diagnostics dataset and suface flux Diagnostics datasets
        """

        def varList(self, date, get_variables_2dm, id_lat, id_lon, out_variable_2dm, chunk_size, time, lat, lon, dir_data):
            """ build the variables list for output netcdf file in the format
            of (T2M,U2M, V2M, U10M,V10M,QV2M, lat, lon, time)"""

            date_ind, time_ind1, time_ind2, time_ind3 = MERRAgeneric().getTime(date)

            # Setup size of saving file
            date_size = len(date_ind)
            hour_size = len(time[0])
            int_size = date_size // chunk_size
            res_type = (date_size * hour_size) % (chunk_size * hour_size)

            if (res_type > 0):
                size_type = [chunk_size * hour_size] * int_size + [res_type]

            else:
                size_type = [chunk_size * hour_size] * int_size

            # Get the variables and set up the list for saving in netCDF file
            t2m_total = []
            u2m_total = []
            v2m_total = []
            u10m_total = []
            v10m_total = []
            sh2m_total = []

            var_out = {'T2M': ['2-metre_air_temperature',
                               'temperature_at_2m_above_the_displacement_height',
                               'K',     t2m_total],
                       'U2M': ['2-metre_eastward_wind',
                               'eastward_wind_at _2m_above_the_displacement_height',
                               'm/s',   u2m_total],
                       'V2M': ['2-metre_northward_wind',
                               'northward_wind_at_2m_above_the_displacement_height',
                               'm/s',   v2m_total],
                       'U10M':['10-metre_eastward_wind',
                               'eastward_wind_at_10m_above_displacement_height',
                               'm/s',   u10m_total],
                       'V10M':['10-metre_northward_wind',
                               'northward_wind_at_10m_above_the_displacement_height',
                               'm/s',   v10m_total],
                       'QV2M':['2-metre_specific_humidity',
                               'specific_humidity_at_2m',
                               'kg/kg', sh2m_total]}

            var_list = []
            for i in range(0, len(get_variables_2dm[0: -3])):
                for x in var_out.keys():
                    if x == get_variables_2dm[i]:
                        print("------Get Subset of Variable at Surface Level------", get_variables_2dm[i])
                        # the position of T2M, U2M, V2M, U10M, V10M in out_variable_2ds is the position in the get_variables
                        var = MERRAgeneric().dataStuff_2d(i, id_lat, id_lon, out_variable_2dm)

                        # restructure the shape
                        var_total = MERRAgeneric().restruDatastuff(var)
                        del var
                        var_out[x][3] = var_total
                        del var_total
                        var_list.append([get_variables_2dm[i], var_out[x][0], var_out[x][1], var_out[x][2], var_out[x][3]])

            return var_list, time_ind3, date_ind, size_type

        def saveData(self, date, get_variables_2dm, id_lat, id_lon, out_variable_2dm, chunk_size, time, lat, lon, dir_data):
            """build a NetCDF file for saving output variables
            (T2M,U2M, V2M, U10M,V10M, PRECTOT, PRECTOTCORR,T2MDEW,QV2M,
            lat, lon, time)"""

            # get the varable list and time indices
            var_list, time_ind3, date_ind, size_type = self.varList(date, get_variables_2dm, id_lat, id_lon, out_variable_2dm, chunk_size, time, lat, lon, dir_data)

            # save nc file
            var_low = 0
            var_up = 0

            for i in range(0, 1):
            # for i in range(0, len(size_type)):
                var = size_type[i]
                var_low = var_up
                var_up = var_low + var

                # set up file path and names
                file_ncdf  = path.join(dir_data,("merra_sa" + "_" + (date_ind[var_low // len(time[0])]) + "_" + "to" + "_" +(date_ind[var_up // len(time[0]) - 1]) + ".nc"))
                rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')

                print("Saved File Type:", rootgrp.file_format)

                rootgrp.source      = 'Merra, abstrated meteorological variables from metadata at surface level'
                rootgrp.featureType = "2_Dimension"

                # Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup
                LAT = lat[0]
                LON = lon[0]

                # dimensions
                time  = rootgrp.createDimension('time', None)
                lats  = rootgrp.createDimension('lats', len(LAT))
                lons  = rootgrp.createDimension('lons', len(LON))

                # Output the results of extracted variables
                for x in range(0, len(var_list)):
                    out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('time','lats','lons'), fill_value=9.9999999E14)
                    out_var.standard_name = var_list[x][1]
                    out_var.long_name = var_list[x][2]
                    out_var.units         = var_list[x][3]
                    out_var.missing_value = 9.9999999E14
                    # out_var.fmissing_value = (9.9999999E14, 'f')
                    # out_var.vmax = (9.9999999E14, 'f')
                    # out_var.vmin = (-9.9999999E14, 'f')
                    out_var[:,:,:] = var_list[x][4][var_low:var_up,:,:]     # data generic name with data stored in it

                Time  = rootgrp.createVariable('time', 'i4', ('time'))
                Time.standard_name = "time"
                # Time.units         = "hour since " + str(datetime.strptime(beg, '%Y/%m/%d'))
                Time.units  = "hours since 1980-1-1 00:00:00"
                Time.calendar      = "gregorian"
                # pass the values
                netCDFTime = []
                for x in range(0, len(time_ind3)):
                    netCDFTime.append(nc.date2num(time_ind3[x], units = Time.units, calendar = Time.calendar))
                Time[:] = netCDFTime[var_low:var_up]

                Latitudes               = rootgrp.createVariable('latitude', 'f4',('lats'))
                Latitudes.standard_name = "latitude"
                Latitudes.units         = "degrees_north"
                Latitudes.axis          = "Y"
                Latitudes[:]  = [L for L in lat[0]]    # pass the values of latitude

                Longitudes               = rootgrp.createVariable('longitude', 'f4',('lons'))
                Longitudes.standard_name = "longitude"
                Longitudes.units         = "degrees_east"
                Longitudes.axis          = "X"
                Longitudes[:] = [L for L in lon[0]]   # pass the values of longitudes


                #close the root group
                rootgrp.close()


class SaveNCDF_sf():
    """ write output netCDF file for abstracted variables from original 2D
        radiation Diagnostics datasets  datasets
    """
    def varList(self, date, get_variables_2dr, get_variables_2ds, get_variables_2dv,
               id_lat, id_lon, out_variable_2dr,out_variable_2ds, out_variable_2dv,
               chunk_size, time, lat, lon, dir_data):

        """build the variables list for output netcdf file
        (PRECTOT, PRECTOTCORR,T2MDEW,SWGDN, SWGDNCLR, LWRNT, LWGEM, LWGNTCLR, lat,lon,time)"""
        date_ind, time_ind1, time_ind2, time_ind3 = MERRAgeneric().getTime(date)

        # Set up time_ind3 with the begin at year-mm-dd 00:30:00
        time_ind3 = time_ind3 + timedelta(minutes=30)

        # Setup size of saving file
        date_size = len(date_ind)
        hour_size = len(time[0])
        int_size = date_size // chunk_size
        res_type = (date_size * hour_size) % (chunk_size*hour_size)

        if (res_type > 0):
            size_type = [chunk_size * hour_size] * int_size + [res_type]
        else:
            size_type = [chunk_size * hour_size] * int_size

        # Get the variables and Set up the list of output variables

        prectot_total = []
        prectotcorr_total = []
        t2mdew_total = []
        swgdn_total = []
        swgdnclr_total = []
        lwgnt_total = []
        lwgem_total = []
        lwgntclr_total = []

        var_out = {'PRECTOT':    ['total_precipitation',
                                  'total_surface_precipitation_flux',
                                  'kg/m2/s', prectot_total],
                   'PRECTOTCORR':['total_precipitation_corrected',
                                  'total_surface_precipitation_flux_corrected',
                                  'kg/m2/s', prectotcorr_total],
                   'T2MDEW':     ['2-metre_dew_point_temperature',
                                  'dew_point_temperature_at_2m',
                                  'K',       t2mdew_total],
                   'SWGDN':      ['surface_incoming_shortwave_flux',
                                  'surface_incoming_shortwave_flux',
                                  'W/m2',    swgdn_total],
                   'SWGDNCLR':   ['surface_incoming_shortwave_flux_assuming_clear_sky',
                                  'surface_incoming_shortwave_flux_assuming_clear_sky',
                                  'W/m2', swgdnclr_total],
                   'LWGNT':      ['surface_net_downward_longwave_flux',
                                  'surface_net_downward_longwave_flux',
                                  'W/m2', lwgnt_total],
                   'LWGEM':      ['longwave_flux_emitted_from_surface',
                                  'longwave_flux_emitted_from_surface',
                                  'W/m2', lwgem_total],
                   'LWGNTCLR':   ['surface_net_downward_longwave_flux_assuming_clear_sky',
                                  'surface_net_downward_longwave_flux_assuming_clear_sky',
                                  'W/m2', lwgntclr_total]}

        var_list = []
        for i in range(0, len(get_variables_2ds[0:-3])):
            for x in var_out.keys():
                if x == get_variables_2ds[i]:
                    print("------Get Subset of Variable at Surface Level------", get_variables_2ds[i])

                    # the position of PRECTOT in out_variable_2ds is 0
                    var = MERRAgeneric().dataStuff_2d(i, id_lat, id_lon, out_variable_2ds)

                    # restructure the shape
                    var_total = MERRAgeneric().restruDatastuff(var)
                    del var
                    var_out[x][3] = var_total
                    del var_total
                    var_list.append([get_variables_2ds[i], var_out[x][0], var_out[x][1], var_out[x][2], var_out[x][3]])

        for i in range(0, len(get_variables_2dv[0:-3])):
            for x in var_out.keys():
                if x == get_variables_2dv[i]:
                    print("------Get Subset of Variable at Surface Level------", get_variables_2dv[i])
                    # the position of PRECTOT in out_variable_2ds is 0
                    var = MERRAgeneric().dataStuff_2d(i, id_lat, id_lon, out_variable_2dv)
                    # restructing the shape
                    var_total = MERRAgeneric().restruDatastuff(var)
                    del var
                    var_out[x][3] = var_total
                    del var_total
                    var_list.append([get_variables_2dv[i], var_out[x][0], var_out[x][1], var_out[x][2], var_out[x][3]])

        for i in range(0, len(get_variables_2dr[0:-3])):
            for x in var_out.keys():
                if x == get_variables_2dr[i]:
                    print("------Get Subset of Variable at Surface Level------", get_variables_2dr[i])
                    var = MERRAgeneric().dataStuff_2d(i, id_lat, id_lon, out_variable_2dr)

                    # restructure the shape
                    var_total = MERRAgeneric().restruDatastuff(var)
                    del var
                    var_out[x][3] = var_total

                    if x == 'LWGNT':
                        lwgnt_total = var_total
                    elif x == 'LWGNTCLR':
                        lwgntclr_total = var_total
                    elif x == 'LWGEM':
                        lwgem_total = var_total

                    del var_total
                    var_list.append([get_variables_2dr[i],
                                     var_out[x][0],
                                     var_out[x][1],
                                     var_out[x][2],
                                     var_out[x][3]])

        # Getting downwelling longwave radiation flux conversed by the function below :
        # - downwelling longwave flux in air =  Upwelling longwave flux from surface + surface net downward longwave flux
        # - downwelling longwave flux in air assuming clear sky =  Upwelling longwave flux from surface + surface net downward longwave flux assuming clear sky

        lwgdn_total = lwgnt_total + lwgem_total
        lwgdnclr_total = lwgntclr_total + lwgem_total

        #append LWGDN, LWGDNCLR
        var_list.append(['LWGDN', 'downwelling_longwave_flux_in_air','downwelling_longwave_flux_in_air','W/m2', lwgdn_total])
        var_list.append(['LWGDNCLR','downwelling_longwave_flux_in_air_assuming_clear_sky','downwelling_longwave_flux_in_air_assuming_clear_sky','W/m2', lwgdnclr_total])

        return var_list, time_ind3, date_ind, size_type

    def saveData(self, date, get_variables_2dr, get_variables_2ds, get_variables_2dv,
                 id_lat, id_lon, out_variable_2dr, out_variable_2ds, out_variable_2dv,
                 chunk_size, time, lat, lon, dir_data):
        """Build a NetCDF file for saving output variables (SWGDN, SWGDNCLR,
           LWRNT, LWGEM, LWGNTCLR, LWGAB, LWGABCLR, lat, lon, time)
           """

        # get variable list and time indices
        var_list, time_ind3, date_ind, size_type = self.varList(
                date, get_variables_2dr, get_variables_2ds,
                get_variables_2dv, id_lat, id_lon, out_variable_2dr,
                out_variable_2ds, out_variable_2dv, chunk_size, time, lat,
                lon, dir_data)

        var_low = 0
        var_up = 0
        for i in range(0, 1):
        # for i in range(0, len(size_type)):
            var = size_type[i]
            var_low = var_up
            var_up = var_low + var

            # set up file path and names
            beg = date_ind[var_low // len(time[0])]
            end = date_ind[var_up // len(time[0]) - 1]
            fname = "merra_sf" + "_" + beg + "_" + "to" + "_" + end + ".nc"
            file_ncdf  = path.join(dir_data, fname)
            rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')

            print("Saved File Type:", rootgrp.file_format)

            rootgrp.source      = 'Merra, abstrated radiation variables from metadata at surface level'
            rootgrp.featureType = "2_Dimension"

            #Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup
            LAT = lat[0]
            LON = lon[0]

            #dimensions
            time  = rootgrp.createDimension('time', None)
            lats  = rootgrp.createDimension('lats', len(LAT))
            lons  = rootgrp.createDimension('lons', len(LON))

            #Output the results of extracted variables
            for x in range(0, len(var_list)):
                out_var = rootgrp.createVariable(var_list[x][0], 'f4', ('time', 'lats', 'lons'), fill_value=9.9999999E14)
                out_var.standard_name = var_list[x][1]
                out_var.long_name = var_list[x][2]
                out_var.units         = var_list[x][3]
                out_var.missing_value = 9.9999999E14
                # out_var.fmissing_value = (9.9999999E14, 'f')
                # out_var.vmax = (9.9999999E14, 'f')
                # out_var.vmin = (-9.9999999E14, 'f')
                out_var[:,:,:] = var_list[x][4][var_low:var_up,:,:]

            Time               = rootgrp.createVariable('time', 'i4', ('time'))
            Time.standard_name = "time"
            Time.units  = "hours since 1980-1-1 00:30:00"
            Time.calendar      = "gregorian"
            # pass the values
            netCDFTime = []
            for x in range(0, len(time_ind3)):
                netCDFTime.append(nc.date2num(time_ind3[x], units = Time.units, calendar = Time.calendar))
            Time[:] = netCDFTime[var_low:var_up]

            Latitudes               = rootgrp.createVariable('latitude', 'f4', ('lats'))
            Latitudes.standard_name = "latitude"
            Latitudes.units         = "degrees_north"
            Latitudes.axis          = 'Y'
            Latitudes[:]  = [L for L in lat[0]]                    # pass the values of latitude

            Longitudes               = rootgrp.createVariable('longitude', 'f4', ('lons'))
            Longitudes.standard_name = "longitude"
            Longitudes.units         = "degrees_east"
            Longitudes.axis          = 'X'
            Longitudes[:] = [L for L in lon[0]]                    # pass the values of longitudes

            #close the root group
            rootgrp.close()


class SaveNCDF_sc():
        """ write output netCDF file for abstracted variables from original 2D
            Constant Model Parameters
        """
        def varList(self, get_variables_2dc, id_lat, id_lon, out_variable_2dc,
                    chunk_size, time, lat, lon, dir_data):
            """Create a NetCDF file for saving output variables
            (Dataset object, also the root group)"""

            # Get the variables and set up the list for saving in netCDF file
            phis_total = []
            frlake_total = []
            frland_total = []
            frlandice_total = []
            frocean_total = []
            sgh_total = []

            var_out = {'PHIS':['surface_geopotential_height',
                               'surface_geopotential_height','m2/s2', phis_total],
                       'FRLAKE':['fraction_of_lake',
                                 'fraction_of_lake','1',frlake_total],
                       'FRLAND':['fraction_of_land',
                                 'fraction_of_land', '1', frland_total],
                       'FRLANDICE':['fraction_of_land_ice',
                                    'fraction_of_land_ice', '1', frlandice_total],
                       'FROCEAN':['fraction_of_ocean',
                                  'fraction_of_ocean', '1',frocean_total],
                       'SGH':['isotropic_stdv_of_GWD_topography',
                              'isotropic_stdv_of_GWD_topography', 'm', sgh_total]}

            var_list = []
            for i in range(0, len(get_variables_2dc[0:-3])):
                for x in var_out.keys():
                    if x == get_variables_2dc[i]:
                        print("------Get Subset of Constant Model Parameters------", get_variables_2dc[i])

                        # the position of T2M, U2M, V2M, U10M, V10M in out_variable_2ds is the position in the get_variables
                        var = MERRAgeneric().dataStuff_2d(i, id_lat, id_lon, out_variable_2dc)

                        # restructure the shape
                        var_total = MERRAgeneric().restruDatastuff(var)
                        del var
                        var_out[x][3] = var_total
                        del var_total
                        var_list.append([get_variables_2dc[i], var_out[x][0], var_out[x][1], var_out[x][2], var_out[x][3]])

            return var_list

        def saveData(self, get_variables_2dc, id_lat, id_lon, out_variable_2dc, chunk_size, time, lat, lon, dir_data):
            """create a NetCDF file for saving output variables (Dataset object, also the root group)"""

            # get var_list
            var_list = self.varList(get_variables_2dc, id_lat, id_lon, out_variable_2dc, chunk_size, time, lat, lon, dir_data)

            #set up file path and names
            file_ncdf  = path.join(dir_data,("merra_sc" + ".nc"))
            rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
            print("Saved File Type:", rootgrp.file_format)
            rootgrp.source      = 'Merra, abstrated constant model parameters'
            rootgrp.featureType = "2_Dimension"

            # Arrange the format of dimensions for time, levels, latitude and longitude for dimension setup

            # dimensions
            time  = rootgrp.createDimension('time', None)
            lats  = rootgrp.createDimension('lats', len(lat[0]))
            lons  = rootgrp.createDimension('lons', len(lon[0]))

            # Output the results of extracted variables
            for x in range(0, len(var_list)):
                out_var = rootgrp.createVariable(var_list[x][0], 'f4',
                                                 ('time','lats','lons'),
                                                 fill_value=9.9999999E14)
                out_var.standard_name = var_list[x][1]
                out_var.long_name = var_list[x][2]
                out_var.units         = var_list[x][3]
                out_var.missing_value = 9.9999999E14
                # out_var.fmissing_value = (9.9999999E14, 'f')
                # out_var.vmax = (9.9999999E14, 'f')
                # out_var.vmin = (-9.9999999E14, 'f')
                out_var[:,:,:] = var_list[x][4][:,:,:] # data generic name with data stored in it

            Time  = rootgrp.createVariable('time', 'i4', ('time'))
            Time.standard_name = "time"
            Time.units  = "hours since 1992-01-02 03:00:00"
            Time.calendar      = "gregorian"

            # Set up the value of time (one single value)
            time_ind4 = datetime.combine(
                    datetime.strptime("1992-01-02", "%Y-%m-%d"),
                    datetime.strptime("0300","%H%M").time())
            time_ind4 = (pandas.date_range(time_ind4, time_ind4, freq = '1H'))

            # pass the values
            netCDFTime = []
            for x in range(0, len(time_ind4)):
                 netCDFTime.append(nc.date2num(time_ind4[x],
                                               units = Time.units,
                                               calendar = Time.calendar))
            Time[:] = netCDFTime[:]

            Latitudes               = rootgrp.createVariable('latitude',
                                                             'f4', ('lats'))
            Latitudes.standard_name = "latitude"
            Latitudes.units         = "degrees_north"
            Latitudes.axis          = "Y"
            Latitudes[:]  = [L for L in lat[0]]  # pass the values of latitude

            Longitudes               = rootgrp.createVariable('longitude',
                                                              'f4', ('lons'))
            Longitudes.standard_name = "longitude"
            Longitudes.units         = "degrees_east"
            Longitudes.axis          = "X"
            Longitudes[:] = [L for L in lon[0]] # pass the values of longitudes

            # close the root group
            rootgrp.close()


"""
Referenced from era_interim.py (Dr.Stephan Gruber): Class ERAdownload()

Class for accessing the parameter file for downloading Merra-2 specified
variables,latitude and longitude coordinates, start, end date, minimum and
maximum elevations.

Args:
    pfile: Full path to a Globsim Download Parameter file.


"""
# TODO: update functions of checking available files
class MERRAdownload(GenericDownload):
        # build full dictionary between variable names from input parameter
        # file and original merra2 data products
        full_variables_dic = {
                'air_temperature': ['air_temperature',
                                    '2-meter_air_temperature'],
                'relative_humidity' : ['relative_humidity',
                                       '2-metre_dewpoint_temperature',
                                       '2-metre_specific_humidity'],
                'precipitation_amount': ['total_precipitation',
                                         'total_precipitation_corrected'],
                'wind_from_direction': ['eastward_wind',
                                        'northward_wind',
                                        '2-meter_eastward_wind',
                                        '2-meter_northward_wind',
                                        '10-meter_eastward_wind',
                                        '10-meter_northward_wind'],
                'wind_speed':['eastward_wind',
                              'northward_wind',
                              '2-meter_eastward_wind',
                              '2-meter_northward_wind',
                              '10-meter_eastward_wind',
                              '10-meter_northward_wind'],
                'downwelling_shortwave_flux_in_air': ['surface_incoming_shortwave_flux'],
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['surface_incoming_shortwave_flux_assuming_clear_sky'],
                'downwelling_longwave_flux_in_air': ['surface_net_downward_longwave_flux',
                                                     'longwave_flux_emitted_from_surface'],
                'downwelling_longwave_flux_in_air_assuming_clear_sky':['surface_net_downward_longwave_flux_assuming_clear_sky',
                                                                       'longwave_flux_emitted_from_surface']}

        # build variables Standards Names and referenced Names for downloading
        # from orginal MERRA-2 datasets
        # 3D Analyzed Meteorological fields data
        full_variables_pl_ana = {'geopotential_height':'H',
                                      'air_temperature':'T',
                                      'eastward_wind':'U',
                                      'northward_wind': 'V'}
        # 3D Assimilated Meteorological fields data
        full_variables_pl_asm = {'relative_humidity': 'RH'}
        # 2D Single-Level Diagnostics data (instantaneous)
        full_variables_sm = {'2-meter_air_temperature': 'T2M',
                                  '2-meter_eastward_wind': 'U2M',
                                  '2-meter_northward_wind':'V2M',
                                  '10-meter_eastward_wind':'U10M',
                                  '10-meter_northward_wind':'V10M',
                                  '2-metre_specific_humidity':'QV2M'}
        # 2D surface flux diagnostics data
        full_variables_sf = {'total_precipitation': 'PRECTOT',
                                  'total_precipitation_corrected': 'PRECTOTCORR'}
        # 2D single-level diagnostics (time-averageed)
        full_variables_sv = {'2-metre_dewpoint_temperature': 'T2MDEW'}
        # 2D radiation diagnostics data
        full_variables_sr = {'surface_incoming_shortwave_flux': 'SWGDN',
                                  'surface_incoming_shortwave_flux_assuming_clear_sky': 'SWGDNCLR',
                                  'surface_net_downward_longwave_flux':'LWGNT',
                                  'longwave_flux_emitted_from_surface': 'LWGEM',
                                  'surface_net_downward_longwave_flux_assuming_clear_sky': 'LWGNTCLR'}

    """
    Class for MERRA-2 data that has methods for querying NASA GES DISC server,
    and returning all variables usually wanted.
    """
    def __init__(self, pfile):
        super().__init__(pfile)
        par = self.par

        self._set_data_directory("merra2")


        # time bounds
        self.date  = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),#get_begin_date(par, 'merra2', ["merra_pl*", "merra_sa*", "merra_sf*"]),
                      'end': datetime.strptime(par['end'], '%Y/%m/%d')}

        # credential
        self.credential = path.join(par['credentials_directory'], ".merrarc")
        self.account = open(self.credential, "r")
        self.inf = self.account.readlines()
        # pass the first line to username  (covert list to str)
        self.username = ''.join(self.inf[0].split())
         # pass the second line to passworrd (covert list to str)
        self.password = ''.join(self.inf[1].split())

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size']

    def getVariables(self, full_variables_dic, full_variables_type):
        """
        build the major variable list for retrieving between lists from
        download parameter file and the product
        """
        get_variables = []

        for i in range(0, len(self.variables)):
            for var in full_variables_dic.keys():
                if  var == self.variables[i]:
                    var_names = full_variables_dic[var]
                    #  Set up the variables list for accassing original type
                    # of MERRA-2 datasets (3D and 2D)
                    for j in range(0, len(var_names)):
                        for var1 in full_variables_type.keys():
                            if var1 == var_names[j]:
                                get_variables.append(full_variables_type[var1])
        # set the list
        get_variables = list(set(get_variables))

        #add the basic variables into the list
        if 'relative_humidity' in list(full_variables_type.keys()) :
            get_variables.extend(['lat','lon','lev','time'])
        elif 'air_temperature' in list(full_variables_type.keys()):
            # !ADD Geopotential Height in the first element of downloading
            # list. Must be the first one
            get_variables.insert(0,'H')
            # add the variables names of latitude, longitude, levels and time
            get_variables.extend(['lat','lon','lev','time'])
        else:
            # add the variables names of latitude, longitude and time
            get_variables.extend(['lat','lon','time'])

        return get_variables

    def retrieve(self):
        """
        Retrive all required MERRA-2 data from NASA Goddard Earth Sciences Data
        and Information Services Center
        """

        # Get merra-2 3d meteorological analysis variables at pressure levels
        t_start = tc.time()

        #Chunk size for spliting files and download [days], Format:Integer
        chunk_size = int(self.chunk_size)

        # download all variables by looping over the date with chunk_size
        startDay = self.date['beg']
        endDay   = self.date['end']

        # Get merra-2 2d Constant Model Parameters (outside of time & date looping!)
        print('''========== Get Variables From MERRA2 2d Time-Invariant ==========''')
        print('''========== Single-level, Constant Model Parameters ==========''')

        self.download_merra_to()

        print("----------- Result NO.0: Completed -----------")

        # Begin downloading time-dependent variables
        x = 0
        for dt in rrule(DAILY, dtstart = startDay, until = endDay):
            currentDay = (str(dt.strftime("%Y")) + "/" + str(dt.strftime("%m")) + "/" + str(dt.strftime("%d")))
            x += 1
            if (x == 1):

                self.date['beg'] = currentDay

                print('DOWNLOADING BEGINS ON:', self.date['beg'])

                # convert date['beg'] from string back to datetime object
                self.date['beg'] = datetime.strptime(self.date['beg'],'%Y/%m/%d')

            if (x == chunk_size or dt == endDay):

                x = 0
                self.date['end'] = currentDay

                print('DOWNLOADING ENDS ON: {}'.format(self.date['end']))

                # convert date['beg'] from string back to datetime object
                self.date['end'] = datetime.strptime(self.date['end'],'%Y/%m/%d')

                # get the urls for all types of data products
                urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = MERRAgeneric().getURLs(self.date)

                # Get merra-2 3d meteorological assimilated variables at pressure levels
                print('''========== Get Variables From Merra-2 3d ==========
                         ========== 3-hourly, Pressure-Level, Assimilated Meteorological Fields ==========''')

                # get the shared variables dictionaries and pass the information to the build-in dictionaries
                get_variables_3dmasm = self.getVariables(self.full_variables_dic,self.full_variables_pl_asm)

                print(get_variables_3dmasm)

                ds_asm = MERRAgeneric().download(self.username, self.password, urls_3dmasm)

                id_lat, id_lon =  MERRAgeneric().getArea(self.area, ds_asm)

                id_lev = MERRAgeneric().getPressureLevels(self.elevation)

                out_variable_3dmasm = MERRAgeneric().Variables(get_variables_3dmasm, ds_asm)

                # get merra-2 meteorological varaibles at pressure levels
                print('''========== Get Variables From Merra-2 3d ==========
                         ========== 6-hourly, Pressure-Level, Analyzed Meteorological Fields ==========''')

                # get the shared variables dictionaries and pass the information to the build-in dictionaries
                get_variables_3dmana = self.getVariables(self.full_variables_dic, self.full_variables_pl_ana)
                print(get_variables_3dmana)

                ds_ana = MERRAgeneric().download(self.username, self.password, urls_3dmana)

                id_lat, id_lon =  MERRAgeneric().getArea(self.area, ds_ana)

                out_variable_3dmana = MERRAgeneric().Variables(get_variables_3dmana, ds_ana)

                # print(out_variable_3dmana)
                lat, lon, lev, time = MERRAgeneric().latLon_3d(out_variable_3dmana, id_lat, id_lon, id_lev)

                # Output merra-2 meteorological analysis variable at pressure levels
                # For T, V, U, H
                SaveNCDF_pl_3dm().saveData(self.date, get_variables_3dmasm, get_variables_3dmana, id_lat, id_lon, id_lev,
                                           out_variable_3dmasm, out_variable_3dmana, chunk_size, time, lev, lat, lon, self.directory, self.elevation)

                print("----------- Result NO.1: Completed -----------")

                # Get merra-2 2d meteorological Diagnostics variables at surface level
                print('''========== Get Variables From Merra-2 ==========
                      ========== 2d, 1-hourly, Single-level, Meteorological Diagnostics ==========''')

                # get the shared variables dictionaries and pass the information to the build-in dictionaries
                get_variables_2dm = self.getVariables(self.full_variables_dic, self.full_variables_sm)

                print(get_variables_2dm)

                ds_2dm = MERRAgeneric().download(self.username, self.password, urls_2dm)

                out_variable_2dm = MERRAgeneric().Variables(get_variables_2dm, ds_2dm)

                lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dm, id_lat, id_lon)

                # Output merra-2 variable at surface level
                SaveNCDF_sa().saveData(self.date, get_variables_2dm, id_lat, id_lon, out_variable_2dm, chunk_size, time, lat, lon, self.directory)

                print("----------- Result NO.2: Completed -----------")

                # Get merra-2 2d surface flux Diagnostics variables at surface level
                print('''========== Get Variables From Merra-2 2d ==========
                      ========== 1-hourly, Single-level, Surface Flux Diagnostics ==========''')

                # get the shared variables dictionaries and pass the information to the build-in dictionaries
                get_variables_2ds = self.getVariables(self.full_variables_dic, self.full_variables_sf)

                print(get_variables_2ds)

                ds_2ds = MERRAgeneric().download(self.username, self.password, urls_2ds)

                out_variable_2ds = MERRAgeneric().Variables(get_variables_2ds, ds_2ds)

                print('''========== Get Variables From Merra-2 2d ==========
                      ========== 1-hourly, Single-Level, Assimilation,Single-Level Diagnostics ==========''')

                # get the shared variables dictionaries and pass the information to the build-in dictionaries
                get_variables_2dv = self.getVariables(self.full_variables_dic, self.full_variables_sv)

                print(get_variables_2dv)

                ds_2dv = MERRAgeneric().download(self.username, self.password, urls_2dv)

                out_variable_2dv = MERRAgeneric().Variables(get_variables_2dv, ds_2dv)

                # Get merra-2 2d radiation variables
                print('''========== Get Variables From Merra-2 2d ==========
                      ========== 1-Hourly, Single-Level, Radiation Diagnostics ==========''')

                # get the shared variables dictionaries and pass the information to the build-in dictionaries
                get_variables_2dr = self.getVariables(self.full_variables_dic, self.full_variables_sr)

                print(get_variables_2dr)

                ds_2dr = MERRAgeneric().download(self.username, self.password, urls_2dr)

                out_variable_2dr = MERRAgeneric().Variables(get_variables_2dr, ds_2dr)

                lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dr, id_lat, id_lon)

                # Output merra-2 radiation variables
                SaveNCDF_sf().saveData(self.date, get_variables_2dr,
                           get_variables_2ds, get_variables_2dv,
                           id_lat, id_lon, out_variable_2dr,
                           out_variable_2ds, out_variable_2dv,
                           chunk_size, time, lat, lon, self.directory)

                print("----------- Result NO.3: Completed -----------")


        t_end = tc.time()
        t_total = int((t_end - t_start) // 60)
        print("Total Time (Minutes):", t_total)

    def download_merra_to(self):
            # check if it already exists
            if path.isfile(path.join(self.directory, ("merra_to" + ".nc"))):
                print("WARNING:  file 'merra_to.nc' already exists and is being skipped")
                return()

            urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = MERRAgeneric().getURLs(self.date)

            # get the shared variables dictionaries and pass the information to the build-in dictionaries
            get_variables_2dc = ['PHIS','FRLAKE','FRLAND','FRLANDICE',
                                 'FROCEAN','SGH','lat','lon','time']

            print(get_variables_2dc)

            ds_2dc = MERRAgeneric().download(self.username, self.password,
                                 url_2dc)

            out_variable_2dc = MERRAgeneric().Variables(get_variables_2dc,
                                           ds_2dc)

            id_lat, id_lon =  MERRAgeneric().getArea(self.area, ds_2dc)

            lat, lon, time = MERRAgeneric().latLon_2d(out_variable_2dc,
                                         id_lat, id_lon)

            # Output merra-2 variable at surface level
            SaveNCDF_sc().saveData(get_variables_2dc, id_lat, id_lon,
                       out_variable_2dc, int(self.chunk_size),
                       time, lat, lon, self.directory)




