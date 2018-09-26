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

from datetime    import datetime
from csv         import QUOTE_NONE


import pandas  as pd
import netCDF4 as nc
import numpy as np

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
            time.units = 'seconds since 1980-01-01 00:00:0.0'  #! For MERRA2 Scaling
        else: 
            time.units = 'seconds since 1900-01-01 00:00:0.0' #! For JRA55 Scaling

#        time.units = 'seconds since 1900-01-01 00:00:0.0' #! For Era_Interim Scaling

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
     
     
def satvapp_kPa_fT(T):
    '''
    Saturation water vapour pressure [kPa] following the Tetens formula, Eq 4.2
    in Stull, Practical Meteorology.
    
    T: Temperature [C]
    '''                          
    e0 = 0.6113  # [kPa]
    b  = 17.2694 # fitting constant
    T1 = 273.15  # [K]
    T2 = 35.86   # [K]
    T += T1
    return e0 * np.exp((b*(T-T1))/(T-T2))
                                                                                                                       

def vapp_kPa_fTd(Td): 
    '''
    Water vapour pressure [hPa] derived from dewpoint temperature [C]. Taken 
    from www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    where it is attributed to (Bolton 1980) 
    https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
    
    Td: Dew point temperature [C] 
    
    '''
    #(Bolton 1980) 
    #https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
    return 6.112 * np.exp((17.67 * Td)/(Td + 243.5)) 
    
    #https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
    #return 6.112 * np.power(10,(7.5 * Td)/(237.3 + Td))
                    
def spec_hum_kgkg(Td, Pr):
    '''
    Specific humidity [Kg/Kg]. Eq 4.7 in Stull, Practical Meteorology.
    Td: Dewpoint temperature [C]
    Pr:  Air pressure [Pa]
    '''
    E = 0.622 # density of vater vapour / density of dry air
    e  = vapp_kPa_fTd(Td) / 10
    P = Pr/1000. # from Pa to kPa 
    spec_hum = E * e / (P - e * (1 - E)) 
    return spec_hum   
   
def water_vap_pressure(RH,T):
    '''
    water vapour pressure [unit:1], Eq C9,C10 in Fiddes and Gruber (2014)
    RH: relative humidity (%)
    Tair: air temperature (kelvin)
    '''   
    es0 = 6.11 # reference saturation vapour pressure at 0ºC 
    T0 = 273.15 # Kelvin     
    lv = 2.5 * 1000000 #latent heat of vaporization of water 
    Rv = 461.5 # gas constant for water vapour 
    es = es0 * np.exp((lv)/Rv * (1/T0 - 1/T)) 
    pv = (RH * es)/100  
    return pv
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
def emissivity_clear_sky(RH,T): 
    '''
    clear sky emissivity, Eq(1) in Fiddes and Gruber (2014)
    pv: water vapour pressure (1)
    T: air temperature (kelvin)
    '''  
    pv = water_vap_pressure(RH, T)
    x1 = 0.43 
    x2 = 5.7 
    e_clear = 0.23 + x1*(pv/T)**(1/x2)  
    return e_clear


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
        
def LW_downward(RH,T,N):
    '''
    incoming longware radiation [W/m2], Eq(14) in Fiddes and Gruber (2014)
    e_clear: clear sky emissivity
    N: cloud cover 
    T: air temperature 
    ''' 
    e_clear = emissivity_clear_sky(RH,T)
    p1 = 6
    p2 = 4
    e_as = 0.979
    con = 5.67 * 10**(-8) # J/s/m/K4 Stefan-Boltzmann constant 
    lw = e_clear*(1-N**p1)+(e_as*(N**p2))*con*T**4
    return lw
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  