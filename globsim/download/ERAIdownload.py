from __future__ import print_function

import glob
import numpy as np
import netCDF4 as nc

from datetime import datetime, timedelta
from ecmwfapi.api import ECMWFDataServer
from math import floor
from os import path, listdir, remove
from fnmatch import filter

from globsim.download.GenericDownload import GenericDownload

from globsim.common_utils import  str_encode

try:
    from nco import Nco
except ImportError:
    print("*** NCO not imported, netCDF appending not possible. ***")
    pass


class ERAIgeneric(object):
    """
    Parent class for other ERA-Interim classes.
    """

    def areaString(self, area):
        """Converts numerical coordinates into string: North/West/South/East"""
        res = str(round(area['north'],2)) + "/"
        res += str(round(area['west'], 2)) + "/"
        res += str(round(area['south'],2)) + "/"
        res += str(round(area['east'], 2))
        return(res)

    def dateString(self, date):
        """Converts datetime objects into string"""
        res = (date['beg'].strftime("%Y-%m-%d") + "/to/" +
               date['end'].strftime("%Y-%m-%d"))
        return(res)

    def getPressureLevels(self, elevation):
        """Restrict list of ERA-interim pressure levels to be downloaded"""
        Pmax = pressure_from_elevation(elevation['min']) + 55
        Pmin = pressure_from_elevation(elevation['max']) - 55
        levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775,
                         800, 825, 850, 875, 900, 925, 950, 975, 1000])
        mask = (levs >= Pmin) * (levs <= Pmax)  # select
        levs = '/'.join(map(str, levs[mask]))
        return levs

    def getDictionaryGen(self, area, date):
        """
        Makes dictionary of generic variables for a server call
        """
        dictionary_gen = {
          'area'    : self.areaString(area),
          'date'    : self.dateString(date),
          'dataset' : "interim",
          'stream'  : "oper",
          'class'   : "ei",
          'format'  : "netcdf",
          'grid'    : "0.75/0.75"}
        return dictionary_gen

    def getDstring(self):
        return ('_' + self.date['beg'].strftime("%Y%m%d") + "_to_" +
                self.date['end'].strftime("%Y%m%d"))

    def download(self):
        query = self.getDictionary()
        target = query.pop('target')

        if path.isfile(target):
            print("WARNING: File '{}' already exists and was skipped".format(target))
        else:
            server = ECMWFDataServer()
            print(server.trace('=== ERA Interim: START ===='))
            server.retrieve(self.getDictionary())
            print(server.trace('=== ERA Interim: STOP ====='))

    def TranslateCF2ERA(self, variables, dpar):
        """
        Translate CF Standard Names into ERA-Interim code numbers.
        """
        self.param = ''
        for var in variables:
            try:
                self.param += dpar.get(var) + '/'
            except TypeError:
                pass
        self.param = self.param.rstrip('/')  # fix last

    def __str__(self):
        string = ("List of generic variables to query ECMWF server for "
                  "ERA-Interim data: {0}")
        return string.format(self.getDictionary)

    def netCDF_merge(self, directory):
        """
        To combine mutiple downloaded erai netCDF files into a large file with
        specified chunk_size(e.g. 500),
        -- give the full name of merged file to the output = outfile
        -- pass all data from the first input netfile to the merged file name
        -- loop over the files_list, append file one by one into the merge file
        -- pass the mergae netcdf file to interpolation module to process
        (to use nc.MFDataset by reading it)

        Args:
            ncfile_in: the full name of downloaded files (file directory + files names)
        e.g.:
              '/home/xquan/src/globsim/examples/erai/era_sa_*.nc'
              '/home/xquan/src/globsim/examples/erai/era_pl_*.nc'
              '/home/xquan/src/globsim/examples/erai/era_sf_*.nc'

        Output: merged netCDF files
        era_all_0.nc, era_all_1.nc, ...,

        """
        # set up nco operator
        nco = Nco()

        # loop over filetypes, read, report
        file_type = ['erai_sa_*.nc', 'erai_sf_*.nc', 'erai_pl_*.nc']
        for ft in file_type:
            ncfile_in = path.join(directory, ft)

            # get the file list
            files_list = glob.glob(ncfile_in)
            files_list.sort()
            num = len(files_list)

            # set up the name of merged file
            if ncfile_in[-7:-5] == 'sa':
                merged_file = path.join(ncfile_in[:-11],'erai_sa_all_'+ files_list[0][-23:-15] + "_" + files_list[num-1][-11:-3] +'.nc')
            elif ncfile_in[-7:-5] == 'sf':
                merged_file = path.join(ncfile_in[:-11],'erai_sf_all_' + files_list[0][-23:-15] + '_' + files_list[num-1][-11:-3] + '.nc')
            elif ncfile_in[-7:-5] == 'pl':
                merged_file = path.join(ncfile_in[:-11],'erai_pl_all_'+ files_list[0][-23:-15] + '_' + files_list[num-1][-11:-3] +'.nc')
            else:
                print('There is not such type of file')

            # combined files into merged files
            nco.ncrcat(input=files_list,output=merged_file, append=True)

            print('The Merged File below is saved:')
            print(merged_file)

            # clear up the data
            for fl in files_list:
                remove(fl)

    def mergeFiles(self, ncfile_in):
        """
        To combine mutiple downloaded erai netCDF files into a large file.

        Args:
            ncfile_in: the full name of downloaded files (file directory + files names)
        e.g.:
              '/home/xquan/src/globsim/examples/erai/era_sa_*.nc'
              '/home/xquan/src/globsim/examples/erai/era_pl_*.nc'
              '/home/xquan/src/globsim/examples/erai/era_sf_*.nc'

        Output: merged netCDF files
        erai_sa_all.nc, erai_sf_all.nc, erai_pl_all.nc

        """

        # read in one type of mutiple netcdf files
        ncf_in = nc.MFDataset(ncfile_in, 'r', aggdim='time')

        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # get spatial dimensions
        lat  = ncf_in.variables['latitude'][:]
        lon  = ncf_in.variables['longitude'][:]
        if pl: # only for pressure level files
            lev  = ncf_in.variables['level'][:]
            nlev = len(lev)

        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]

        #set up the name of merged file
        if ncfile_in[-7:-5] == 'sa':
            ncfile_out = path.join(ncfile_in[:-11],'erai_sa_all' + '.nc')
        elif ncfile_in[-7:-5] == 'sf':
            ncfile_out = path.join(ncfile_in[:-11],'erai_sf_all' + '.nc')
        elif ncfile_in[-7:-5] == 'pl':
            ncfile_out = path.join(ncfile_in[:-11],'erai_pl_all' + '.nc')
        else:
            print('There is not such type of file'    )

        # get variables
        varlist = [str_encode(x) for x in ncf_in.variables.keys()]
        varlist.remove('time')
        varlist.remove('latitude')
        varlist.remove('longitude')
        if pl: #only for pressure level files
            varlist.remove('level')

        #Build the netCDF file
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4_CLASSIC')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'ERA_Interim, merged downloaded original files'
        rootgrp.featureType = "timeSeries"

        # dimensions
        latitude = rootgrp.createDimension('latitude', len(lat))
        longitude = rootgrp.createDimension('longitude', len(lon))
        time    = rootgrp.createDimension('time', None)

        # base variables
        time           = rootgrp.createVariable('time', 'i4',('time'))
        time.long_name = 'time'
        time.units     = 'hours since 1900-01-01 00:00:0.0'
        time.calendar  = 'gregorian'
        latitude            = rootgrp.createVariable('latitude', 'f4',('latitude'))
        latitude.long_name  = 'latitude'
        latitude.units      = 'degrees_north'
        longitude           = rootgrp.createVariable('longitude','f4',('longitude'))
        longitude.long_name = 'longitude'
        longitude.units     = 'degrees_east'

        # assign station characteristics
        latitude[:]  = lat[:]
        longitude[:] = lon[:]
        time[:]    = nctime[:]

        # extra treatment for pressure level files
        try:
            lev = ncf_in.variables['level'][:]
            print("== 3D: file has pressure levels")
            level = rootgrp.createDimension('level', len(lev))
            level           = rootgrp.createVariable('level','i4',('level'))
            level.long_name = 'pressure_level'
            level.units     = 'hPa'
            level[:] = lev
        except:
            print("== 2D: file without pressure levels")
            lev = []

        # create and assign variables based on input file
        for n, var in enumerate(varlist):
            print("VAR: ", var)
            # extra treatment for pressure level files
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4',('time', 'level', 'latitude', 'longitude'))
            else:
                tmp = rootgrp.createVariable(var,'f4',('time', 'latitude', 'longitude'))
            tmp.long_name = ncf_in.variables[var].long_name.encode('UTF8') # for erai
            tmp.units     = ncf_in.variables[var].units.encode('UTF8')

            # assign values
            if pl: # only for pressure level files
                tmp[:] = ncf_in.variables[var][:,:,:,:]
            else:
                tmp[:] = ncf_in.variables[var][:,:,:]


        #close the file
        rootgrp.close()
        ncf_in.close()

        #get the file list
        files_list = glob.glob(ncfile_in)
        files_list.sort()

        #clear up the data
        for fl in files_list:
            remove(fl)


class ERAIpl(ERAIgeneric):
    """Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server.

    Args:
        date: A dictionary specifying the time period desired with a begin
              and an end date given as a datetime.datetime object.

        area: A dictionary delimiting the area to be queried with the latitudes
              north and south, and the longitudes west and east [decimal deg].

        elevation: A dictionary specifying the min/max elevation of the area of
                   interest. This is used to determine the pressure levels
                   needed. Unit: [m].

        variables:  List of variable(s) to download based on CF Standard Names.

        directory: Directory to hold output files

    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(2013, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  15.0,
                 'west'  :  60.0,
                 'east'  : 105.0}
        elevation = {'min' :    0,
                     'max' : 8850}
        variables  = ['air_temperature', 'relative_humidity']
        directory = '/Users/stgruber/Desktop'
        ERAIpl = ERAIpl(date, area, elevation, variables, directory)
        ERAIpl.download()
    """
    def __init__(self, date, area, elevation, variables, directory):
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        outfile = 'erai_pl' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable names.
        dpar = {'air_temperature'   : '130.128',           # [K]
                'relative_humidity' : '157.128',           # [%]
                'wind_speed'        : '131.128/132.128'}   # [m s-1]

        # translate variables into those present in ERA pl data
        self.TranslateCF2ERA(variables, dpar)
        self.param += '/129.128' # geopotential always needed [m2 s-2]

    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "pl",
           'levellist': self.getPressureLevels(self.elevation),
           'time'     : "00/06/12/18",
           'step'     : "0",
           'type'     : "an",
           'param'    : self.param,
           'target'   : self.file_ncdf
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim data: {0}")
        return string.format(self.getDictionary)


class ERAIsa(ERAIgeneric):
    """
    Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server for surface analysis variables (airt2,dewp2, ozone, vapor,
    wind10).

    Args:

        target:    File name of the netcdf file to be created.

        variable:  List of variable(s) to download based on CF Standard Names

    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(2013, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  15.0,
                 'west'  :  60.0,
                 'east'  : 105.0}
        variables  = ['air_temperature', 'relative_humidity']
        directory = '/Users/stgruber/Desktop'
        ERAIsa = ERAIsa(date, area, variables, directory)
        ERAIsa.download()
    """
    def __init__(self, date, area, variables, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile = 'erai_sa' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable names.
        dpar = {'air_temperature'   : '167.128',  # [K] 2m values
                'relative_humidity' : '168.128',  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky' :
                    '206.128/137.128',  # [kg m-2] Total column ozone
                                        # [kg m-2] Total column W vapor
                'wind_speed' : '165.128/166.128'}   # [m s-1] 10m values

        # translate variables into those present in ERA pl data
        self.TranslateCF2ERA(variables, dpar)
        print(self.param)


    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "00/06/12/18",
           'step'     : "0",
           'type'     : "an",
           'param'    : self.param,
           'target'   : self.file_ncdf
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("Class for ERA-Interim surface analysis data: {0}")
        return string.format(self.getDictionary)


class ERAIsf(ERAIgeneric):
    """
    Returns an object for ERA-Interim data that has methods for querying the
    ECMWF server for surface forecast variables (prec, swin, lwin).

    Args:
        target:    File name of the netcdf file to be created.

        variables:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].

        ERAgen:    ERAIgeneric() object with generic information on the area and
                   time span to be

    Example:
        from datetime import datetime
        date  = {'beg' : datetime(1994, 1, 1),
                 'end' : datetime(1994, 1, 1)}
        area  = {'north' :  40.0,
                 'south' :  40.0,
                 'west'  :  60.0,
                 'east'  :  60.0}
        variables  = ['prec','swin','lwin']
        directory = '/Users/stgruber/Desktop'
        ERAIsf = ERAIsf(date, area, variables, directory)
        ERAIsf.download()
    """
    def __init__(self, date, area, variables, directory):
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile = 'erai_sf' + self.getDstring() + '.nc'
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable names.
        # [m] total precipitation
        # [J m-2] short-wave downward
        # [J m-2] long-wave downward
        dpar = {'precipitation_amount'              : '228.128',
                'downwelling_shortwave_flux_in_air' : '169.128',
                'downwelling_longwave_flux_in_air'  : '175.128'}

        # translate variables into those present in ERA pl data
        self.TranslateCF2ERA(variables, dpar)
        print(self.param)


    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "00/12",
           'step'     : "3/6/9/12",
           'type'     : "fc",
           'param'    : self.param,
           'target'   : self.file_ncdf
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary)

class ERAIto(ERAIgeneric):
    """
    Returns an object for downloading and handling ERA-Interim
    topography (invariant).

    Args:

    Example:
        area  = {'north' :  40.0,
                 'south' :  45.0,
                 'west'  :  60.0,
                 'east'  :  65.0}
        directory = '/Users/stgruber/Desktop'
        ERAIto = ERAIto(area, directory)
        ERAIto.download()
    """
    def __init__(self, area, directory):
        self.area       = area
        self.date       = {'beg' : datetime(1979, 1, 1),
                           'end' : datetime(1979, 1, 1)}
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory,'erai_to.nc')

    def getDictionary(self):
        self.dictionary = {
           'levtype'  : "sfc",
           'time'     : "12",
           'step'     : "0",
           'type'     : "an",
           'param'    : "129.128/172.128", # geopotential and land-sea mask
           'target'   : self.file_ncdf
           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA-Interim air tenperature data: {0}")
        return string.format(self.getDictionary)


class ERAIdownload(GenericDownload):
    """
    Class for ERA-Interim data that has methods for querying
    the ECMWF server, returning all variables usually needed.

    Args:
        pfile: Full path to a Globsim Download Parameter file.

    Example:
        ERAd = ERAdownload(pfile)
        ERAd.retrieve()
    """

    def __init__(self, pfile):
        super().__init__(pfile)
        par = self.par
        self._set_data_directory("erai")

        # time bounds
        self.date  = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                      'end': datetime.strptime(par['end'], '%Y/%m/%d')}

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size']

    def retrieve(self):
        """
        Retrieve all required ERA-Interim data from MARS server.
        """
        # prepare time loop
        date_i = {}
        slices = floor(float((self.date['end'] - self.date['beg']).days)/
                       self.chunk_size)+1
        # topography
        top = ERAIto(self.area, self.directory)
        top.download()

        for ind in range (0, int(slices)):
            #prepare time slices
            date_i['beg'] = self.date['beg'] + timedelta(days =
                            self.chunk_size * ind)
            date_i['end'] = self.date['beg'] + timedelta(days =
                            self.chunk_size * (ind+1) - 1)
            if ind == (slices-1):
                date_i['end'] = self.date['end']

            #actual functions
            pl = ERAIpl(date_i, self.area, self.elevation,
                       self.variables, self.directory)
            sa = ERAIsa(date_i, self.area, self.variables, self.directory)
            sf = ERAIsf(date_i, self.area, self.variables, self.directory)

            #download from ECMWF server convert to netCDF
            ERAli = [pl, sa, sf]
            for era in ERAli:
                era.download()


        # report inventory
        self.inventory()


    def inventory(self):
        """
        Report on data avaialbe in directory: time slice, variables, area
        """
        print("\n\n\n")
        print("=== INVENTORY FOR GLOBSIM ERA-INTERIM DATA === \n")
        print("Download parameter file: \n" + self.pfile + "\n")
        # loop over filetypes, read, report
        file_type = ['erai_pl_*.nc', 'erai_sa_*.nc', 'erai_sf_*.nc', 'erai_t*.nc']
        for ft in file_type:
            infile = path.join(self.directory, ft)
            nf = len(filter(listdir(self.directory), ft))
            print(str(nf) + " FILE(S): " + infile)

            if nf > 0:
                # open dataset
                ncf = nc.MFDataset(infile, 'r')

                # list variables
                keylist = [str_encode(x) for x in ncf.variables.keys()]

                print("    VARIABLES:")
                print("        " + str(len(keylist)) +
                      " variables, inclusing dimensions")
                for key in keylist:
                    print("        " + ncf.variables[key].long_name)

                # time slice
                time = ncf.variables['time']
                tmin = nc.num2date(min(time[:]), time.units,
                                   calendar=time.calendar).strftime('%Y/%m/%d')
                tmax = nc.num2date(max(time[:]), time.units,
                                   calendar=time.calendar).strftime('%Y/%m/%d')
                print("    TIME SLICE")
                print("        " + str(len(time[:])) + " time steps")
                print("        " + tmin + " to " + tmax)

                # area
                lon = ncf.variables['longitude']
                lat = ncf.variables['latitude']
                nlat = str(len(lat))
                nlon = str(len(lon))
                ncel = str(len(lat) * len(lon))
                print("    BOUNDING BOX / AREA")
                print("        " + ncel + " cells, " + nlon +
                      " W-E and " + nlat + " S-N")
                print("        N: " + str(max(lat)))
                print("        S: " + str(min(lat)))
                print("        W: " + str(min(lon)))
                print("        E: " + str(max(lon)))

                ncf.close()

    def __str__(self):
        return "Object for ERA-Interim data download and conversion"
