#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic classes, methods, functions used for more than one reanalysis.
#
#
# (C) Copyright Stephan Gruber (2017)
#         
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
from __future__  import print_function

from datetime    import datetime, timedelta
from csv         import QUOTE_NONE
from os          import mkdir, path, makedirs
from math        import floor

import pandas  as pd
import netCDF4 as nc
import numpy as np
import glob
import os
# handle python 3 string types
try:
    basestring
except:
    basestring = str

class ParameterIO(object):
    """
    Reads generic parameter files and makes values available as dictionary.
    
    # read file
    par = ParameterIO('examples/par/examples.globsim_download')
    
    # access first reanalysis variable
    par.variables[0]
    
    # access data_directory
    par.data_directory
    
    # access north end of bounding box
    par.bbN
    """

    def __init__(self, pfile):
        """
        Instantiate a new object and set conventions.
        """
        self.fmt_date = "%Y/%m/%d"
        self.pfile    = pfile
        self.comment  = "#"
        self.assign   = "="
        self.file_read()
        
    def __getitem__(self, item):
        attr = getattr(self, item)
        return(attr)

    def file_read(self):
        """
        Read parameter file into a list of strings (each line is one entry) and
        parse content into a dictionary.
        """
        # read file
        with open(self.pfile, "r") as myfile:
            inpts_str = myfile.readlines()

        # parse content
        for line in inpts_str:
            d = self.line2dict(line)
            if d is not None:
                self.__dict__[list(d.keys())[0]] = list(d.values())[0]

    def __is_only_comment(self, lin):
        # checks whether line contains nothing but comment
        for c in lin:
            if c != " ":
                if c == self.comment:
                    return True
                else:
                    return False

    def __string2datetime(self, valu):
        # checks if value is a date string. If true, a datetime object is
        # returned. If false, value is returned unchanged.
        if not isinstance(valu, basestring):
            return valu

        # see if time conversion is possible
        try:
            valu = datetime.strptime(valu, self.fmt_date)
        except ValueError:
            pass
        return valu
        
    def __string2logical(self, valu):
        # checks if value is a true / false. If true, a logical object is
        # returned. If false, value is returned unchanged.
        if not isinstance(valu, basestring):
            return valu

        # see if time conversion is possible
        if valu.lower() == "false":
            valu = False
        if valu.lower() == "true":
            valu = True
        return valu

    def __string2datetime_list(self, dates):
        # convert list of date strings to datetime
        return [self.__string2datetime(date) for date in dates]

    def line2dict(self, lin):
        """
        Converts one line of a parameter file into a dictionary. Comments
        are recognised and ignored, float vectors are preserved, datetime
        converted from string.
        """
        # Check if this is valid
        if self.__is_only_comment(lin):
            return None

        # Remove possible trailing comment form line
        lin = lin.split(self.comment)[0]

        # Discard lines without value assignment
        if len(lin.split(self.assign)) != 2:
            return None

        # Extract name and value, strip of leading/trailling blanks
        name = lin.split(self.assign)[0].strip()
        valu = lin.split(self.assign)[1].strip()

        # Make a vector is commas are found
        if valu.find(",") > 0:
            # Convert to float or list of float if numeric
            try:
                valu = list(map(float, valu.split(",")))
            except ValueError:
                valu = list(valu.split(","))
                valu = [v.strip() for v in valu]
        else:
            try:
                valu = float(valu)
            except ValueError:
                pass
                    
        # Convert to datetime if it is datetime
        valu = self.__string2datetime(valu)
        
        # Convert to logical if logical
        valu = self.__string2logical(valu)
        
        # Make dictionary and return
        return {name: valu}
    
class GenericDownload(object):
    """
    Generic functionality for download classes
    """
        
    def __init__(self, pfile):
        # read parameter file
        self.pfile = pfile
        self.par = par = ParameterIO(self.pfile)
        
        self._set_elevation(par)        
        self._set_area(par)
        self._check_area(par)
            
        self.variables = par.variables
    
    @staticmethod
    def _check_area(par):
        if (par.bbN < par.bbS) or (par.bbE < par.bbW):
            raise ValueError("Bounding box is invalid: {}".format(self.area))
            
        if (np.abs(par.bbN - par.bbS) < 1.5) or (np.abs(par.bbE - par.bbW) < 1.5):
            raise ValueError("Download area is too small to conduct interpolation.")
    
    @staticmethod        
    def _set_area(par):
        self.area  = {'north':  par.bbN,
                      'south':  par.bbS,
                      'west' :  par.bbW,
                      'east' :  par.bbE}
    
    @staticmethod
    def _set_elevation(par):
        self.elevation = {'min' : par.ele_min, 
                          'max' : par.ele_max}
                          
    def _set_data_directory(self, name):
        self.directory = path.join(self.par.project_directory, name)  
        if path.isdir(self.directory) == False:
            makedirs(self.directory) 

class GenericInterpolate:
    
    def __init__(self, ifile):
        #read parameter file
        self.ifile = ifile
        self.par = par = ParameterIO(self.ifile)
        self.dir_out = self.makeOutDir(par)
        self.variables = par.variables
        self.list_name = par.station_list.split(path.extsep)[0]
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        
        #read station points 
        self.stations = StationListRead(self.stations_csv)  
        
        # time bounds, add one day to par.end to include entire last day
        self.date  = {'beg' : par.beg,
                      'end' : par.end + timedelta(days=1)}
                      
        # chunk size: how many time steps to interpolate at the same time?
        # A small chunk size keeps memory usage down but is slow.
        self.cs  = int(par.chunk_size)
                      
    def makeOutDir(self, par):
        '''make directory to hold outputs'''
        
        dirIntp = path.join(par.project_directory, 'interpolated')
        
        if not (path.isdir(dirIntp)):
            makedirs(dirIntp)
            
        return dirIntp
        
class GenericScale:
    
    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        self.par = par = ParameterIO(self.sfile)
        self.intpdir = path.join(par.project_directory, 'interpolated')
        self.scdir = self.makeOutDir(par)
        self.list_name = par.station_list.split(path.extsep)[0]
        
        # get the station file
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        #read station points 
        self.stations = StationListRead(self.stations_csv)  
        
        # read kernels
        self.kernels = par.kernels
        if not isinstance(self.kernels, list):
            self.kernels = [self.kernels]
    
    def getOutNCF(self, par, src, scaleDir = 'scale'):
        '''make out file name'''
        
        timestep = str(par.time_step) + 'h'
        src = '_'.join(['scaled', src, timestep])
        src = src + '.nc'
        fname = path.join(self.scdir, src)
        
        return fname
    
    def makeOutDir(self, par):
        '''make directory to hold outputs'''
        
        dirSC = path.join(par.project_directory, 'scaled')
        
        if not (path.isdir(dirSC)):
            makedirs(dirSC)
            
        return dirSC
        
def variables_skip(variable_name):
        '''
        Which variable names to use? Drop the ones that are dimensions.  
        '''
        skip = 0
        dims = ('time', 'level', 'latitude', 'longitude', 'station', 'height')
        if variable_name in dims:
            skip = 1      
        return skip 

def StationListRead(sfile):  
    '''
    Reads ASCII station list and returns a pandas dataframe.
    
    # read station list
    stations = StationListRead('examples/par/examples_list1.globsim_interpolate')
    print(stations['station_number'])
    '''
    # read file
    raw = pd.read_csv(sfile)    
    raw = raw.rename(columns=lambda x: x.strip())
    return(raw)


def ScaledFileOpen(ncfile_out, nc_interpol, times_out, t_unit, station_names=None):
    '''
    Open netCDF file for scaled results (same for all reanalyses) or create it 
    if it does not exist. Returns the file object so that kernel functions can 
    successively write variables to it.
    
    '''
    
    try:
        # read file if it exists
        rootgrp = nc.Dataset(ncfile_out, 'a')
        
        #TODO: make sure the times and stations match if file exists
    
    except: 
        # make netCDF outfile, variables are written in kernels
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'Reanalysis data interpolated and scaled to stations'
        rootgrp.featureType = "timeSeries"
        
        name = ncfile_out[-9:-3]

        # dimensions
        station = rootgrp.createDimension('station', 
                     len(nc_interpol.variables['station'][:]))
        time    = rootgrp.createDimension('time', len(times_out))

        # base variables
        time           = rootgrp.createVariable('time', 'i8',('time'))
        time.long_name = 'time'
        time.axis      = 'T'
        
        time.units = t_unit
        time.calendar  = 'gregorian'
        
        station             = rootgrp.createVariable('station', 'i4',('station'))
        station.long_name   = 'station index for time series data'
        station.units       = ''
        station.standard_name = 'platform_id'
        
        latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
        latitude.standard_name  = 'latitude'
        latitude.units      = 'degrees_N'
        latitude.long_name  = 'WGS84 latitude'
        latitude.axis       = 'Y'
        
        longitude           = rootgrp.createVariable('longitude','f4',('station'))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'WGS84 longitude'
        longitude.units     = 'degrees_E' 
        longitude.axis      = 'X'
        
        height           = rootgrp.createVariable('height','f4',('station'))
        height.standard_name = 'height_above_reference_ellipsoid'
        height.long_name = 'height_above_reference_ellipsoid'
        height.units     = 'm' 
        height.axis      = 'Z'
        height.positive  = 'up'
        
        crs           = rootgrp.createVariable('crs','i4')
        crs.long_name = 'coordinate system'
        crs.grid_mapping_name = 'latitude_longitude'
        crs.longitude_of_prime_meridian = 0.0
        crs.semi_major_axis = 6378137
        crs.inverse_flattening = 298.2572236
        crs.wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
        
        # assign base variables
        time[:]      = times_out
        station[:]   = nc_interpol.variables['station'][:]
        latitude[:]  = nc_interpol.variables['latitude'][:]
        longitude[:] = nc_interpol.variables['longitude'][:]
        height[:]    = nc_interpol.variables['height'][:]
        
        # add station names to netcdf
        if station_names is not None:
            # first convert to character array
            names_out = nc.stringtochar(np.array(station_names, 'S32'))
            
            # create space in the netcdf
            nchar        = rootgrp.createDimension('name_strlen', 32) 
            st           = rootgrp.createVariable('station_name', "S1", ('station', 'name_strlen'))
            st.standard_name = 'platform_name'
            st.units     = ''
            
            # add data
            st[:] = names_out
            
    return rootgrp
    

def convert_cummulative(data):
    """
    Convert values that are serially cummulative, such as precipitation or 
    radiation, into a cummulative series from start to finish that can be 
    interpolated on for sacling. 
    data: 1-dimensional time series 
    """                       
    # get increment per time step
    diff = np.diff(data)
    diff = np.concatenate(([data[0]], diff))
    
    # where new forecast starts, the increment will be smaller than 0
    # and the actual value is used
    mask = diff < 0
    diff[mask] = data[mask]
    
    #get full cummulative sum
    return np.cumsum(diff, dtype=np.float64)

def cummulative2total(data, time):
    """
    Convert values that are serially cummulative, such as precipitation or 
    radiation, into a cummulative series from start to finish that can be 
    interpolated on for sacling. 
    data: 1-dimensional time series 
    """                       
    # get increment per time step
    diff = np.diff(data)
    diff = np.concatenate(([data[0]], diff))
    
    # where new forecast starts, the increment will be smaller than 0
    # and the actual value is used
    
    mask = [timei.hour in [3,15] for timei in time]
    diff[mask] = data[mask]
    
    mask = diff < 0
    diff[mask] = 0
    
    return diff
    
def get_begin_date(par, data_folder, match_strings):
    """ Get the date to begin downloading when some files already exist
    
    Parameters
    ----------
    par : ParameterIO object
    data_folder : str
        name of subdirectory containing data files. Examples: merra2, era5
    match_strings : list
        list of glob-style strings to check. Examples ["merra_pl*", "merra_sa*","merra_sf*"]
    
    Returns
    -------
    datetime
        datetime object corresponding to the desired begin date (replaces par['beg'])
        
    This makes an inventory of all the files that have been downloaded so far and
    returns the next date to begin downloading.  If all match_strings are downloaded up to the same
    day, then the following day is returned. Otherwise, the 
    """
    directory = par['project_directory']
    print("Searching for existing files in directory")
    
    if not all([len(glob.glob(os.path.join(directory, data_folder, s))) > 0 for s in match_strings]):
        print("No existing files found. Starting download from {}".format(par['beg'].strftime("%Y-%m-%d")))
        return par['beg']
        
    datasets = [nc.MFDataset(os.path.join(directory, data_folder, s)) for s in match_strings]
    dates = [nc.num2date(x['time'][:], x['time'].units, x['time'].calendar) for x in datasets]
    
    latest = [max(d) for d in dates]
    latest = [dt.replace(hour=0, minute=0, second=0, microsecond=0) for dt in latest]
    latest_complete = min(latest)

    begin_date = latest_complete + timedelta(days=1)
    
    print("Found some files in directory. Beginning download on {}".format(begin_date.strftime("%Y-%m-%d"))    )
    return(begin_date)
    
def series_interpolate(time_out, time_in, value_in, cum=False):
    '''
    Interpolate single time series. Convenience function for usage in scaling 
    kernels.
    time_out: Array of times [s] for which output is desired. Integer. 
    time_in:  Array of times [s] for which value_in is given. Integer. 
    value_in: Value time series. Must have same length as time_in.
    cum:      Is valiable serially cummulative like LWin? Default: False.
    '''
    time_step_sec = time_out[1]-time_out[0]
    
    # convert to continuous cummulative, if values are serially cummulative
    if cum:
        value_in = convert_cummulative(value_in)

    # interpolate            
    vi = np.interp(time_out, time_in, value_in)
 
    # convert from cummulative to normal time series if needed
    if cum:
        vi = np.diff(vi) / time_step_sec
        vi = np.float32(np.concatenate(([vi[0]], vi)))
            
    return vi        

def globsimScaled2Pandas(ncdf_in, station_nr):
    '''
    Read a scaled (or interpolated) globsim netCDF file and return all values
    for one station as a Pandas data frame.
    
    ncdf_in: full path to a globsim netCDF (by station)
    
    station_nr: station_number, as given in the stations .csv file to identify 
                the station.
    
    '''
    # open file                                     
    ncf = nc.Dataset(ncdf_in, 'r')
    
    # station mask
    sm = ncf.variables['station'][:] == int(station_nr)                                                               
    # list variables
    varlist = [x.encode('UTF8') for x in ncf.variables.keys()]                                                                                                                                                                                             
    
    # get and convert time
    time = ncf.variables['time'][:]
    t_unit = ncf.variables['time'].units 
    t_cal = ncf.variables['time'].calendar   
    time = nc.num2date(time, units = t_unit, calendar = t_cal)
                
    # make data frame with time   
    df = pd.DataFrame(data=time,columns=['time'])    
    # add variables
    for var in varlist:
        if variables_skip(var):
            continue
        data = ncf.variables[var][:,sm]                                                                                                                                                                                                                                                                                                                                                                                               
        df = pd.concat([df,pd.DataFrame(data=data,columns=[var])],axis=1) 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    return df 
    
def globsim2GEOtop(ncdf_globsim, txt_geotop):
    """
    Convert globsim scaled netCDF to GEOtop meteo file.
    """
        
    outfile = '/Users/stgruber/Supervision/MSc/Mary_Pascale_Laurentian/reanalysis/station/Meteo_0001.txt'
    
    #read time object        
    time = self.rg.variables['time']        
    date = self.rg.variables['time'][:]
    
    #read all other values
    columns = ['Date','AIRT_ERA_C_pl','AIRT_ERA_C_sur','PREC_ERA_mm_sur','RH_ERA_per_sur','SW_ERA_Wm2_sur','LW_ERA_Wm2_sur','WSPD_ERA_ms_sur','WDIR_ERA_deg_sur']
    metdata = np.zeros((len(date),len(columns)))
    metdata[:,0] = date
    for n, vn in enumerate(columns[1:]):
        metdata[:,n+1] = self.rg.variables[vn][:, 0]
    
    #make data frame
    data = pd.DataFrame(metdata, columns=columns)
    data[['Date']] = nc.num2date(date, time.units, calendar=time.calendar)
    
    # round
    decimals = pd.Series([2,1,1,1,1,1,1,1], index=columns[1:])
    data.round(decimals)

    #export to file
    fmt_date = "%d/%m/%Y %H:%M"
    data.to_csv(outfile, date_format=fmt_date, index=False, float_format='%.2f')
        
                   
def globsim2CLASS(ncdf_globsim, met_class, station_nr):
    """
    Convert globsim scaled netCDF to CLASS-CTEM .met file.
    
    ncdf_globsim: full path to a globsim scaled netCDF (by station)
    
    met_class: full path to the CLASS-CTEM met file to write.
    
    station_nr: station_number, as given in the stations .csv file to identify 
                the station.
    
    The columns in CLASS-CTEM MET files are:            
    1)  Hour
    2)  Minute
    3)  Day of year
    4)  Year YYYY
    5)  Shortwave Radiation (W/m2)
    6)  Longwave Radiation (W/m2)
    7)  Precip (mm/s)
    8)  Temp.(°C)
    9)  Specific Humidity (Kg/Kg)
    10) Wind Speed (m/s)
    11) Pressure (Pa)
    
    """
    
    # columns to export
    columns = ['time','SW_ERA_Wm2_sur', 'LW_ERA_Wm2_sur', 'PREC_ERA_mmsec_sur',
               'AIRT_ERA_C_sur','SH_ERA_kgkg_sur', 'WSPD_ERA_ms_sur',
               'AIRT_PRESS_Pa_pl']
    
    # output ASCII formatting
    formatters={"time":             "  {:%H %M  %j  %Y}".format,
                "SW_ERA_Wm2_sur":   "{:8.2f}".format,
                "LW_ERA_Wm2_sur":   "{:8.2f}".format,
                "PREC_ERA_mmsec_sur":  "{:13.4E}".format,
                "AIRT_ERA_C_sur":   "{:8.2f}".format,
                "SH_ERA_kgkg_sur":  "{:11.3E}".format,
                "WSPD_ERA_ms_sur":  "{:7.2f}".format,
                "AIRT_PRESS_Pa_pl": "{:11.2f}".format}
    
    # get data
    df = globsimScaled2Pandas(ncdf_globsim, station_nr)                   
    
    # convery precipitation
    df["PREC_ERA_mmsec_sur"] = df["PREC_ERA_mm_sur"]/1800.
    
    # write FORTRAN formatted ASCII file
    with open(met_class, 'w') as f:
        f.write(' ')
        f.write(df.to_string(columns = columns,
                formatters=formatters, 
                header=False, index=False))
    f.close()    
     
def str_encode(value, encoding = "UTF8"):
    '''
    handles encoding to allow compatibility between python 2 and 3
    specifically with regards to netCDF variables.   Python 2 imports 
    variable names as unicode, whereas python 3 imports them as str.
    '''
    if type(value) == str:
        return(value)
    else:
        return(value.encode(encoding))
        
def create_globsim_directory(target_dir, name):
    """
    creates globsim directory
    """
    # create top-level
    TL = path.join(target_dir, name)
    mkdir(TL)
    
    # create subdirectories
    mkdir(path.join(TL, "eraint"))
    mkdir(path.join(TL, "Grib"))
    mkdir(path.join(TL, "jra55"))
    mkdir(path.join(TL, "merra2"))
    mkdir(path.join(TL, "par"))
    mkdir(path.join(TL, "scale"))
    mkdir(path.join(TL, "station"))
    mkdir(path.join(TL, "era5"))
    
    return(True)    
    
def nc_to_clsmet(ncd, out_dir, src, start=None, end=None):
    """
    @args
    ncd: netcdf dataset
    src: data source ("ERAI","ERA5" "MERRA2", "JRA55")
    """
    #open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)
    
    # find number of stations
    nstn = len(n['station'][:])
    
    # get date / time columns
    time = nc.num2date(n['time'][:], 
                        units=n['time'].units, 
                        calendar=n['time'].calendar)
    time_step = (time[1]-time[0]).seconds # in second
    
    
    HH =   [x.timetuple().tm_hour for x in time]
    MM =   [x.timetuple().tm_min  for x in time]
    DDD =  [x.timetuple().tm_yday for x in time]
    YYYY = [x.timetuple().tm_year for x in time]
    
    TIME = np.stack((HH, MM, DDD, YYYY))
    
    # get shortwave
    SW = "SW_{}_Wm2_sur".format(src)
    SW = n[SW][:]
    
    # get longwave
    LW = "LW_{}_Wm2_sur".format(src)
    LW = n[LW][:]
    
    # get precip
    PREC = "PREC_{}_mm_sur".format(src)
    PREC = n[PREC][:] 
    PREC /= time_step # [mm/0.5h] to [mm/s]
    
    # get temp
    AIRT = "AIRT_{}_C_sur".format(src)
    AIRT = n[AIRT][:]
    
    # get specific humidity
    SH = "SH_{}_kgkg_sur".format(src)
    SH = n[SH][:]
    
    # get wind speed
    WSPD = "WSPD_{}_ms_sur".format(src)
    WSPD = n[WSPD][:]
    
    # get pressure
    PRESS = "PRESS_{}_Pa_pl".format(src)
    PRESS = n[PRESS][:]
    
    # get site names
    NAMES = nc.chartostring(n['station_name'][:])
    
    data = np.stack((SW, LW, PREC, AIRT, SH, WSPD, PRESS))
    
    # write output files
    for i in range(nstn):
        
        # massage data into the right shape
        out_array = data[:,:,i]
        out_array = np.concatenate((TIME, out_array), 0)
        out_array = np.transpose(out_array)
        
        #get station name 
        st_name = NAMES[i]
        
        filename = "{}_{}_{}.MET".format(i, st_name, src)
        savepath = path.join(out_dir, filename)
        
        # create file
        np.savetxt(savepath, out_array, 
                fmt = [" %02u", "%02u", " %03u", " %2u", "%8.2f", "%8.2f", 
                        "%13.4E", "%8.2f", "%11.3E", "%7.2f", "%11.2f" ])
                        
                        
def nc_to_gtmet(ncd, out_dir, src, start=None, end=None):
    """
    @args
    ncd: netcdf dataset
    src: data source ("ERAI","ERA5" "MERRA2", "JRA55")
    """
    #open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)
    
    # find number of stations
    nstn = len(n['station'][:])
    
    # get date / time column
    time = nc.num2date(n['time'][:], 
                        units=n['time'].units, 
                        calendar=n['time'].calendar)
    
    time = [x.strftime('%d/%m/%Y %H:%M') for x in time]
    time = pd.DataFrame(time)
    
    # get precip
    PREC = "PREC_{}_mm_sur".format(src)
    PREC = n[PREC][:]
    
    # get wind velocity
    WSPD = "WSPD_{}_ms_sur".format(src)
    WSPD = n[WSPD][:]
    
    # get wind direction
    WDIR = "WDIR_{}_deg_sur".format(src)
    WDIR = n[WDIR][:]
    
    # get windx and windy
    
    # get RH 
    RH = "RH_{}_per_sur".format(src)
    RH = n[RH][:]
        
    # get air temp
    AIRT = "AIRT_{}_C_sur".format(src)
    AIRT = n[AIRT][:]
    
    # get dew temp (missing?)
    
    # get air pressure
    PRESS = "PRESS_{}_Pa_pl".format(src)
    PRESS = n[PRESS][:]
    PRESS *= 1e-5      # convert to bar for geotop
   
    # get shortwave solar global (direct / diffuse missing?)
    SW = "SW_{}_Wm2_sur".format(src)
    SW = n[SW][:]
    
    # get longwave incoming
    LW = "LW_{}_Wm2_sur".format(src)
    LW = n[LW][:]

    # get site names
    NAMES = nc.chartostring(n['station_name'][:])

    # combine data variables into array
    data = np.stack((PREC, WSPD, WDIR, RH, AIRT, PRESS, SW, LW))
    
    # write output files
    for i in range(nstn):
        
        # massage data into the right shape
        out_df = pd.DataFrame(np.transpose(data[:,:,i]))
        out_df = pd.concat([time, out_df], axis=1)
        out_df.columns = ["Date", "IPrec", "WindVelocity", "WindDirection", "RH",
                            "AirTemp", "AirPress", "SWglobal", "LWin"]
        

        #get station name 
        st_name = NAMES[i]
        
        # prepare paths
        filename = "{}-{}_{}_Forcing_0001.txt".format(i, st_name, src)
        savepath = path.join(out_dir, filename)
        
        # create file
        out_df.to_csv(savepath, index=False, float_format="%10.5f")
        
def create_empty_netcdf(ncfile_out, stations, nc_in, time_units):
        #TODO change date type from f4 to f8 for lat and lon
        '''
        Creates an empty station file to hold interpolated reults. The number of 
        stations is defined by the variable stations, variables are determined by 
        the variable list passed from the gridded original netCDF.
        
        ncfile_out: full name of the file to be created
        stations:   station list read with generic.StationListRead() 
        variables:  variables read from netCDF handle
        lev:        list of pressure levels, empty is [] (default)
        '''
        
        #Build the netCDF file
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4_CLASSIC')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.featureType = "timeSeries"
                                                
        # dimensions
        station = rootgrp.createDimension('station', len(stations))
        time    = rootgrp.createDimension('time', None)
                
        # base variables
        time           = rootgrp.createVariable('time', 'i4',('time'))
        time.long_name = 'time'
        time.units     = time_units
        time.calendar  = 'gregorian'
        station             = rootgrp.createVariable('station', 'i4',('station'))
        station.long_name   = 'station for time series data'
        station.units       = '1'
        latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
        latitude.long_name  = 'latitude'
        latitude.units      = 'degrees_north'    
        longitude           = rootgrp.createVariable('longitude','f4',('station'))
        longitude.long_name = 'longitude'
        longitude.units     = 'degrees_east' 
        height           = rootgrp.createVariable('height','f4',('station'))
        height.long_name = 'height_above_reference_ellipsoid'
        height.units     = 'm'  
        
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
                    
        # create and assign variables based on input file
        for n, var in enumerate(nc_in.variables):
            if variables_skip(var):
                continue                 
            print("VAR: ", str_encode(var))
            # extra treatment for pressure level files           
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4',('time', 'level', 'station'))
            else:
                tmp = rootgrp.createVariable(var,'f4',('time', 'station'))     
            tmp.long_name = nc_in.variables[var].long_name.encode('UTF8') 
            tmp.units     = nc_in.variables[var].units.encode('UTF8')  
                    
        return rootgrp
        