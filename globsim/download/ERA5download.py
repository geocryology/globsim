
from __future__ import print_function

import cdsapi
import numpy as np
import netCDF4 as nc
import urllib3

from datetime import datetime, timedelta
from fnmatch import filter
from math import floor
from os import path, listdir

from globsim.common_utils import str_encode
from globsim.download.GenericDownload import GenericDownload
from globsim.meteorology import pressure_from_elevation

urllib3.disable_warnings()


class ERA5generic(object):
    """
    Parent class for other ERA5 classes.
    """
    CDS_DICT = {'130.128': 'temperature',
                '157.128': 'relative_humidity',
                '131.128': 'u_component_of_wind',
                '132.128': 'v_component_of_wind',
                '129.128': 'geopotential',
                '167.128': '2m_temperature',
                '168.128': '2m_dewpoint_temperature',
                '206.128': 'total_column_ozone',         # also consider surface_solar_radiation_downward_clear_sky
                '137.128': 'total_column_water_vapour',  # also consider surface_solar_radiation_downward_clear_sky
                '165.128': '10m_u_component_of_wind',
                '166.128': '10m_v_component_of_wind',
                '228.128': 'total_precipitation',
                '169.128': 'surface_solar_radiation_downwards',
                '175.128': 'surface_thermal_radiation_downwards',
                '172.128': 'land_sea_mask'}

    def areaString(self, area):
        """Converts numerical coordinates into string: North/West/South/East"""
        res = str(round(area['north'],2)) + "/"
        res += str(round(area['west'], 2)) + "/"
        res += str(round(area['south'],2)) + "/"
        res += str(round(area['east'], 2))
        return(res)

    def timeString(self, era5type):

        if era5type == 'ensemble_members':
            times = np.arange(0, 24, 3)
        elif era5type == 'reanalysis':
            times = np.arange(0, 24)

        times = [str(t).zfill(2) for t in times]
        times = [t + ':00' for t in times]

        return times

    def dateString(self, date):
        """Converts datetime objects into string"""
        res = (date['beg'].strftime("%Y-%m-%d") + "/" +
               date['end'].strftime("%Y-%m-%d"))
        return(res)

    def typeString(self, era5type):

        if era5type == 'ensemble_members':
            return 'era5_ens'
        elif era5type == 'reanalysis':
            return 'era5_rea'

    def getPressureLevels(self, elevation):
        """Restrict list of ERA5 pressure levels to be downloaded"""
        Pmax = pressure_from_elevation(elevation['min']) + 55
        Pmin = pressure_from_elevation(elevation['max']) - 55
        levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775,
                         800, 825, 850, 875, 900, 925, 950, 975, 1000])
        mask = (levs >= Pmin) * (levs <= Pmax)  # select
        levs = [str(levi) for levi in levs[mask]]
        return levs

    def getDictionaryGen(self, area, date):
        """
        Makes dictionary of generic variables for a server call
        """
        dictionary_gen = {
            'area'    : self.areaString(area),
            'date'    : self.dateString(date),
            'dataset' : "era5",
            'stream'  : "oper",
            'class'   : "ea",
            'format'  : "netcdf"}
        return dictionary_gen

    def getDstring(self):
        return ('_' + self.date['beg'].strftime("%Y%m%d")
                + "_to_"
                + self.date['end'].strftime("%Y%m%d"))

    def getDownFile(self, levStr):
        '''make outfile name with data level string'''
        typeStr = self.typeString(self.era5type)
        outfile = typeStr + '_' + levStr + self.getDstring() + '.nc'
        return outfile

    def download(self, storage='cds'):
        ''' 'patched' download function to allow download using cdsapi from CDS servers'''
        server = cdsapi.Client()

        query = self.getDictionary()
        target = query.pop('target')
        levtype = query.pop('levtype')

        # select dataset
        dataset = {'pl': 'reanalysis-era5-pressure-levels',
                   'sfc': 'reanalysis-era5-single-levels'
                   }[levtype]
        query = self.ECM2CDS(query)

        # launch download request
        print(server.info('=== ERA5 ({}API): START ACCESS ON {} ===='.format("CDS", storage.upper())))
        if path.isfile(target):
            print("WARNING: File '{}' already exists and was skipped".format(target))
        else:
            server.retrieve(dataset, query, target)
        print(server.info('=== ERA5 ({}API): END ACCESS ON {} ===='.format("CDS", storage.upper())))

    def ECM2CDS(self, query):
        ''' convert ECMWF query to CDS format '''

        # remove unnecessary keys
        for key in ['class', 'dataset', 'stream', 'type']:
            query.pop(key)

        # replace key('param'):val(str)  string with key('variables'):val(list)
        query['variable'] = [self.CDS_DICT[L] for L in query.pop('param').split('/')]

        return(query)

    def TranslateCF2ERA(self, variables, dpar):
        """
        Translate CF Standard Names into ERA5 code numbers.
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
                  "ERA5 data: {0}")
        return string.format(self.getDictionary)


class ERA5pl(ERA5generic):
    """Returns an object for ERA5 data that has methods for querying the
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
        ERA5pl = ERA5pl(date, area, elevation, variables, directory)
        ERA5pl.download()
    """
    def __init__(self, era5type, date, area, elevation, variables, directory):
        self.era5type = era5type
        self.date       = date
        self.area       = area
        self.elevation  = elevation
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory, self.getDownFile('pl'))

        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable names.
        dpar = {'air_temperature'   : '130.128',           # [K]
                'relative_humidity' : '157.128',           # [%]
                'wind_speed'        : '131.128/132.128'}   # [m s-1]

        # translate variables into those present in ERA pl data
        self.TranslateCF2ERA(variables, dpar)
        self.param += '/129.128'  # geopotential always needed [m2 s-2]

    def getDictionary(self):
        self.dictionary = {'product_type': self.era5type,
                           'levtype'  : "pl",
                           'levellist': self.getPressureLevels(self.elevation),
                           'time'     : self.timeString(self.era5type),
                           'step'     : "0",
                           'type'     : "an",
                           'param'    : self.param,
                           'target'   : self.file_ncdf
                           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA5 data: {0}")
        return string.format(self.getDictionary)


class ERA5sa(ERA5generic):
    """
    Returns an object for ERA5 data that has methods for querying the
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
        ERA5sa = ERA5sa(date, area, variables, directory)
        ERA5sa.download()
    """
    def __init__(self, era5type, date, area, variables, directory):
        self.era5type = era5type
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile         = self.getDownFile('sa')
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable names.
        dpar = {'air_temperature': '167.128',  # [K] 2m values
                'relative_humidity': '168.128',  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky':
                    '206.128/137.128',  # [kg m-2] Total column ozone
                                        # [kg m-2] Total column W vapor
                'wind_speed': '165.128/166.128'}   # [m s-1] 10m values

        # translate variables into those present in ERA pl data
        self.TranslateCF2ERA(variables, dpar)

    def getTime(self):

        times = np.arange(0, 24)
        times = [str(t).zfill(2) for t in times]
        times = [t + ':00:00' for t in times]
        times = '/'.join(times)

        return times

    def getDictionary(self):
        self.dictionary = {'product_type': self.era5type,
                           'levtype'  : "sfc",
                           'time'     : self.timeString(self.era5type),
                           'step'     : "0",
                           'type'     : "an",
                           'param'    : self.param,
                           'target'   : self.file_ncdf
                           }
        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("Class for ERA5 surface analysis data: {0}")
        return string.format(self.getDictionary)


class ERA5sf(ERA5generic):
    """
    Returns an object for ERA5 data that has methods for querying the
    ECMWF server for surface forecast variables (prec, swin, lwin).

    Args:
        target:    File name of the netcdf file to be created.

        variables:  List of variable(s) to download that can include one, several
                   , or all of these: ['airt', 'rh', 'geop', 'wind'].

        ERAgen:    ERA5generic() object with generic information on the area and
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
        ERA5sf = ERA5sf(date, area, variables, directory)
        ERA5sf.download()
    """
    def __init__(self, era5type, date, area, variables, directory):
        self.era5type = era5type
        self.date       = date
        self.area       = area
        self.directory  = directory
        outfile         = self.getDownFile('sf')
        self.file_ncdf  = path.join(self.directory, outfile)

        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable names.
        # [m] total precipitation
        # [J m-2] short-wave downward
        # [J m-2] long-wave downward
        dpar = {'precipitation_amount'              : '228.128',  # [kg m**-2 s**-1]
                'downwelling_shortwave_flux_in_air' : '169.128',  # [J m**-2]
                'downwelling_longwave_flux_in_air'  : '175.128'}  # [J m**-2]

        # translate variables into those present in ERA pl data
        self.TranslateCF2ERA(variables, dpar)

    def getStep(self):
        steps = np.arange(1, 13)
        steps = [str(s) for s in steps]
        steps = '/'.join(steps)

        return steps

    def getDictionary(self):

        self.dictionary = {'product_type': self.era5type,
                           'levtype'  : "sfc",
                           'time'     : self.timeString(self.era5type),
                           'step'     : self.getStep(),
                           'type'     : "fc",
                           'param'    : self.param,
                           'target'   : self.file_ncdf}

        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA5 air tenperature data: {0}")
        return string.format(self.getDictionary)


class ERA5to(ERA5generic):
    """
    Returns an object for downloading and handling ERA5
    topography (invariant).

    Args:

    Example:
        area  = {'north' :  40.0,
                 'south' :  45.0,
                 'west'  :  60.0,
                 'east'  :  65.0}
        directory = '/Users/stgruber/Desktop'
        ERA5to = ERA5to(area, directory)
        ERA5to.download()
    """
    def __init__(self, era5type, area, directory):
        self.era5type   = era5type
        self.area       = area
        self.date       = {'beg' : datetime(2017, 1, 1),
                           'end' : datetime(2017, 1, 1)}
        self.directory  = directory
        self.file_ncdf  = path.join(self.directory, self.typeString(self.era5type) + '_to.nc')

    def getDictionary(self):
        self.dictionary = {
           'product_type' : self.era5type,
           'levtype'  : "sfc",
           'time'     : "00:00:00",
           'step'     : "0",
           'type'     : "an",
           'param'    : "129.128/172.128",  # geopotential and land-sea mask
           'target'   : self.file_ncdf}

        self.dictionary.update(self.getDictionaryGen(self.area, self.date))
        return self.dictionary

    def __str__(self):
        string = ("List of variables to query ECMWF server for "
                  "ERA5 air tenperature data: {0}")
        return string.format(self.getDictionary)


class ERA5download(GenericDownload, ERA5generic):
    """
    Class for ERA5 data that has methods for querying
    the ECMWF server, returning all variables usually needed.

    Args:
        pfile: Full path to a Globsim Download Parameter file.
        api  : Which API to use. Either 'cds' (default) or 'ecmwf' (deprecated'
        storage : Which server to access data from. Either 'cds'
                    (default) or 'ecmwf' (deprecated). Note
                 that you can use the cds api to access the ecmwf storage (MARS)
    Example:
        ERAd = ERA5download(pfile)
        ERAd.retrieve()
    """

    def __init__(self, pfile, era5type='reanalysis', api='cds', storage='cds'):
        super().__init__(pfile)
        par = self.par
        self.era5type = era5type

        if self.era5type == 'reanalysis':
            self._set_input_directory("era5")
            self.topo_file = 'era5_rea_to.nc'
            
        elif self.era5type == 'ensemble_members':
            self._set_input_directory("era5ens")
            self.topo_file = 'era5_ens_to.nc'
            
        # time bounds
        self.date = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                     'end': datetime.strptime(par['end'], '%Y/%m/%d')}

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size']

        # set api and download server/storage
        self.api = api
        self.storage = storage

    def retrieve(self):
        """
        Retrieve all required ERA5 data from MARS/CDS server.
        """
        # prepare time loop
        date_i = {}
        slices = floor(float((self.date['end'] - self.date['beg']).days) /
                       self.chunk_size) + 1

        # topography
        if path.isfile(path.join(self.directory, self.topo_file)):
            print(f"WARNING: File {self.topo_file} already exists. Skipping.")
        else:
            top = ERA5to(self.era5type, self.area, self.directory)
            top.download(self.storage)

        for ind in range(0, int(slices)):
            # prepare time slices
            date_i['beg'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * ind))
            date_i['end'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * (ind + 1) - 1))
            if ind == (slices - 1):
                date_i['end'] = self.date['end']

            # actual functions
            pl = ERA5pl(self.era5type, date_i, self.area, self.elevation,
                        self.variables, self.directory)
            sa = ERA5sa(self.era5type, date_i, self.area,
                        self.variables, self.directory)
            sf = ERA5sf(self.era5type, date_i, self.area,
                        self.variables, self.directory)

            ERAli = [pl, sa, sf]
            for era in ERAli:
                era.download(self.storage)

        # report inventory
        self.inventory()

    def inventory(self):
        """
        Report on data avaialbe in directory: time slice, variables, area
        """
        print("\n\n\n")
        print("=== INVENTORY FOR GLOBSIM ERA5 DATA ===\n")
        print("Download parameter file: \n" + self.pfile + "\n")
        # loop over filetypes, read, report
        file_type = [self.typeString(self.era5type) + '_pl_*.nc',
                     self.typeString(self.era5type) + '_sa_*.nc',
                     self.typeString(self.era5type) + '_sf_*.nc',
                     self.typeString(self.era5type) + '_t*.nc']
        for ft in file_type:
            infile = path.join(self.directory, ft)
            nf = len(filter(listdir(self.directory), ft))
            print(str(nf) + " FILE(S): " + infile)

            if nf > 0:
                # open dataset
                ncf = nc.MFDataset(infile, 'r', aggdim='time')

                # list variables
                keylist = [str_encode(x) for x in ncf.variables.keys()]
                print("    VARIABLES:")
                print("        " + str(len(keylist)) +
                      " variables, including dimensions")
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
        return "Object for ERA5 data download and conversion"