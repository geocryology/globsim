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
from datetime import datetime
import pandas  as pd
import netCDF4 as nc
import numpy as np

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
                self.__dict__[d.keys()[0]] = d.values()[0]

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

        # Make dictionary and return
        return {name: valu}

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


def ScaledFileOpen(ncfile_out, nc_interpol, times_out):
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
        
        if name == 'eraint':
	   time.units = 'seconds since 1900-01-01 00:00:0.0' #! For Era_Interim Scaling
	elif name == 'merra2' :
	   time.units = 'hours since 1980-01-01 00:00:0.0'  #! For MERRA2 Scaling
        else: 
	   time.units = 'hours since 1900-01-01 00:00:0.0' #! For JRA55 Scaling

        time.units = 'seconds since 1900-01-01 00:00:0.0' #! For Era_Interim Scaling

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
        
        # assign base variables
        time[:]      = times_out
        station[:]   = nc_interpol.variables['station'][:]
        latitude[:]  = nc_interpol.variables['latitude'][:]
        longitude[:] = nc_interpol.variables['longitude'][:]
        height[:]    = nc_interpol.variables['height'][:]

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
 
    # convert from cummulative to noremal time series if needed
    if cum:
        vi = np.diff(vi) / time_step_sec
        vi = np.float32(np.concatenate(([vi[0]], vi)))
            
    return vi        
                                                  
                
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
        
                   
def globsim2CLASS(ncdf_globsim, met_class):
    """
    Convert globsim scaled netCDF to CLASS-CTEM .met file.
    """
                       
                                           