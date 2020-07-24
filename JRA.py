#!/usr/bin/env python
# -*- coding: utf-8 -*- 
from __future__      import print_function

import urllib.request
import urllib.error
import http.cookiejar
import json
import glob
import tarfile
import time
import sys
import re

import netCDF4       as nc
import numpy         as np

from datetime            import datetime, timedelta
from globsim.generic     import ParameterIO, StationListRead, ScaledFileOpen, str_encode, series_interpolate, variables_skip, create_empty_netcdf, GenericDownload
from globsim.meteorology import LW_downward, pressure_from_elevation
from os                  import path, listdir, remove, makedirs
from math                import exp, floor, atan2, pi
from fnmatch             import filter
from scipy.interpolate   import interp1d

try:
    import ESMF
    
    # Check ESMF version.  7.0.1 behaves differently than 7.1.0r 
    ESMFv = int(re.sub("[^0-9]", "", ESMF.__version__))
    ESMFnew = ESMFv > 701  
except ImportError:
    print("*** ESMF not imported, interpolation not possible. ***")
    pass 



class RDA(object):
    
    def __init__(self, username, password):
        '''Return an object for RDA data sets, to submit subset requests 
        on select gridded data sets, and to check on the processing status of 
        any subset requests. Details could be found:
            https://www2.cisl.ucar.edu/data-portals/research-data-archive/command-line-subset-requests-%E2%80%93-rdams'''
        
        self.base = 'https://rda.ucar.edu/apps/'
        self.loginurl = 'https://rda.ucar.edu/cgi-bin/login'
        self.cookie_file = 'auth.rda_ucar_edu'
        self.username = username
        self.password = password
        
    def makeOpener(self, theurl):
        '''make the opener based on username and password'''
        
        passman = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        passman.add_password(None, theurl, self.username, self.password)
        authhandler = urllib.request.HTTPBasicAuthHandler(passman)
        opener = urllib.request.build_opener(authhandler)
        urllib.request.install_opener(opener)
        
        return opener
    
    def tarExtract(self, tarf, directory):
        '''extract the tar file'''
        
        filename, file_extension = path.splitext(tarf)
        
        if file_extension == '.tar':
            tar = tarfile.open(tarf)
            tar.extractall(path = directory)
            tar.close()
            remove(tarf)
    
    def urlOpen(self, theurl):
        '''open the url'''
        
        opener = self.makeOpener(theurl)
        request = urllib.request.Request(theurl)
        
        try:
            url = opener.open(request)
        except urllib.error.HTTPError as e:
            if e.code == 401:
                print('RDA username and password invalid. Please try again\n')
                opener = self.makeOpener(theurl)
                try:
                    url = opener.open(request)
                except urllib.error.HTTPError as e:
                    if e.code == 401:
                        print(
                            'RDA username and password invalid, or you are'\
                            'not authorized to access this dataset.\n')
                        print('Please verify your login information at'\
                              'http://rda.ucar.edu\n.')
                        sys.exit()
        
        return url
    
    def add_http_cookie(self, url, authstring):
        '''add_http_cookie(url,authstring): Get and add authentication cookie
        to http file download handler'''     
        
        cj = http.cookiejar.MozillaCookieJar(self.cookie_file)
        openrf = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(cj))
        frequest = urllib.request.Request(url, authstring)
        cj.add_cookie_header(frequest)
        response = openrf.open(frequest)
        openerf = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(cj))
        urllib.request.install_opener(openerf)
    
    def getHelp(self):
        '''get the help information'''
        
        theurl = self.base + 'help'
        url = self.urlOpen(theurl)
        
        print(url.read().decode())
        
    
    def getSummary(self, dsID = None):
        '''get the summary of given dataset'''
        
        print('\nGetting summary information. Please wait as this may take'\
              'awhile.\n')
        
        theurl = theurl = self.base + 'summary/' + dsID
        url = self.urlOpen(theurl)
        
        print(url.read().decode())
        
    def getMeta(self, dsID):
        '''get the metadata of gieven dataset'''
        
        theurl = self.base + 'metadata/' + dsID
        url = self.urlOpen(theurl)
        
        print(url.read().decode())
        
    def getParaSummary(self, dsID):
        '''submit a subset request control file.
           Subset request control files are built from the parameters dumped 
           out by the '-get_metadata <dsnnn.n>' option.'''
        
        print('\nGetting parameter summary. Please wait as this may take'\
              'awhile.\n')
        
        theurl = self.base + 'paramsummary/' + dsID
        url = self.urlOpen(theurl)
        
        print(url.read().decode())
        
    def getStatus(self):
        '''get the status of requests'''
        
        theurl = self.base + 'request'
        url = self.urlOpen(theurl)
        
        print(url.read().decode())
        
        
    def submit(self, controlparms):      
        '''submit JRA55 dataset request based on predescribed control
        parameters'''
        
        theurl = self.base + 'request'
        
        jsondata = '{'
        for k in list(controlparms.keys()):
            jsondata += '"' + k + '"' + ":" + '"' + controlparms[k] + '",'
        jsondata = jsondata[:-1]
        jsondata += '}'
        print('\n ====== Submitting Request ====\n')
        
        if len(jsondata) > 1:
            request = urllib.request.Request(
                theurl, jsondata.encode(), {'Content-type':'application/json'})
        else:
            request = urllib.request.Request(theurl)
        
        opener = self.makeOpener(theurl)
        try:
            url = opener.open(request)
        except urllib.error.HTTPError as e:
            if e.code == 401:
                print('RDA username and password invalid.  Please try again\n')
                (username, password) = get_userinfo()
                opener = self.makeOpener(theurl)
                try:
                    url = opener.open(request)
                except urllib.error.HTTPError as e:
                    if e.code == 401:
                        print(
                            'RDA username and password invalid, or you are'\
                            'not authorized to access this dataset.\n')
                        print('Please verify your login information at'\
                              'http://rda.ucar.edu\n.')
                        sys.exit()
        
        print(url.read().decode())
    
    def downloadSinglefile(self, remfile, outfile):
        '''download_file(remfile,outfile) : download a file from a remote
         server (remfile) to a local location (outfile)'''
        
        frequest = urllib.request.Request(remfile)
        fresponse = urllib.request.urlopen(remfile)
        with open(outfile, 'wb') as handle:
            handle.write(fresponse.read())
    
    def downloadFile(self, filelist, directory):
        '''download_files(filelist,directory): Download multiple files from the 
        rda server and save them to a local directory'''
        
        backslash = '/'
        filecount = 0
        percentcomplete = 0
        localsize = ''
        length = 0
        length = len(filelist)
        if not path.exists(directory):
            makedirs(directory)
        for key, value in filelist.items():
            downloadpath, localfile = key.rsplit("/", 1)
            outpath = directory + backslash + localfile
            percentcomplete = (float(filecount) / float(length))
            self.update_progress(percentcomplete, directory)
            if path.isfile(outpath):
                localsize = path.getsize(outpath)
                if(str(localsize) != value):
                    self.downloadSinglefile(key, outpath)
            elif(not path.isfile(outpath)):
                self.downloadSinglefile(key, outpath)
            self.tarExtract(outpath, directory)
                
    
            filecount = filecount + 1
            percentcomplete = (float(filecount) / float(length))
        self.update_progress(percentcomplete, directory)
        
    
    def update_progress(self, progress, outdir):
        
        barLength = 20  # Modify this to change the length of the progress bar
        status = ""
        if isinstance(progress, int):
            progress = float(progress)
        if not isinstance(progress, float):
            progress = 0
            status = "error: progress var must be float\r\n\n"
        if progress < 0:
            progress = 0
            status = "Halt...\r\n\n"
        if progress >= 1:
            progress = 1
            status = "Done...\r\n\n"
        block = int(round(barLength * progress))
        text = "\r ====== Downloading Request ======\n "
        sys.stdout.write(text)
        sys.stdout.flush()
        
    def getDSindex(self):
        '''get the index of submitted dataset index'''
        
        theurl = self.base + 'request/'
        url = self.urlOpen(theurl)
        
        authdata = 'email='+self.username+'&password='+self.password+'&action=login'
        authdata = authdata.encode()
        
        responses = url.read().decode()
        responses = responses.split('\n')
        responses = [item for item in responses if item.startswith("RequestIndex")]
        
        status = [item.split('- ')[1] for item in responses]
        status = np.where(np.asarray(status) == 'Online')[0]
        dsIndex = [item.split(':  ')[1] for item in responses]
        dsIndex = [item.split(', ')[0] for item in dsIndex]
        dsIndex = np.asarray(dsIndex)[status]
        
        return dsIndex

    
    def download(self, directory, ds):
        '''download all the data completed from the NCAR server'''
        
        if not isinstance(ds, str):
            ds = str(ds)
            
        theurl = self.base + 'request/' + ds + '/filelist'
        url = self.urlOpen(theurl)
        
        authdata = 'email='+self.username+'&password='+self.password+'&action=login'
        authdata = authdata.encode()
        
        jsonfilelist = url.read().decode()

        if not jsonfilelist[0] != "{":
            
            filelist = json.loads(jsonfilelist)

            # get cookie required to download data files
            self.add_http_cookie(self.loginurl, authdata)
        
            print("\n\nStarting Download.\n\n")
        
            self.downloadFile(filelist, directory)
                
    def purge(self, dsIndex):
        '''delete dataset from NCAR server based on given dsIndex'''
        
        for ds in dsIndex:
            if not isinstance(ds, str):
                ds = str(ds)
            theurl = self.base + 'request/' + ds
            
            opener = self.makeOpener(theurl)
            request = urllib.request.Request(theurl)
            request.get_method = lambda: 'DELETE'
            
            try:
                url = opener.open(request)
            except urllib.error.HTTPError as e:
                if e.code == 401:
                    print('RDA username and password invalid.  Please try again\n')
                    opener = self.makeOpener(theurl)
                    try:
                        url = opener.open(request)
                    except urllib.error.HTTPError as e:
                        if e.code == 401:
                            print(
                                'RDA username and password invalid, or you are not authorized to access this dataset.\n')
                            print('Please verify your login information at http://rda.ucar.edu\n.')
                            sys.exit()
            
            print(url.read().decode())   

class JRApl(object):
    
    def __init__(self, date, area, elevation, variables, rda):
        '''Returns an object for JRA55 data that has methods for querying the
        NCAR server for pressure level variables (prec, swin, lwin). '''
        
        self.date = date
        self.area = area
        self.elevation = elevation
        #self.directory = directory
        
        dpar = {'air_temperature'   : ['Temperature'],
                'relative_humidity' : ['Relative humidity'],
                'wind_speed'        : ['u-component of wind', 
                                       'v-component of wind']}
        
        self.param = self.getParam(dpar, variables)
        
    def getParam(self, dpar, variables):
        
        varlist = [] 
        for var in variables:
            varlist.append(dpar.get(var))

        varlist = [item for item in varlist if item is not None]
        varlist = [item for sublist in varlist for item in sublist]
        
        if len(varlist) > 0:
            varlist += ['Geopotential height', 'level', 'pressure']
        
        return varlist
    
    def makeDate(self):
        '''convert data format to NCAR RDA request'''
        
        beg = self.date['beg'].strftime('%Y%m%d%H%M')
        end = self.date['end'].strftime('%Y%m%d%H%M')
        dateRange = beg + '/to/' + end
        
        return dateRange

    def getPressureLevels(self):
        total_ele37 = [100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 
                       600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 
                       925, 950, 975, 1000]
        
        # flip max and min because 1000 is the bottom and 0 is the top
        elevationMax = pressure_from_elevation(self.elevation['min'])
        elevationMin = pressure_from_elevation(self.elevation['max'])
        
        minNum = min(total_ele37, key=lambda x:abs(x-elevationMin))
        maxNum = min(total_ele37, key=lambda x:abs(x-elevationMax))
        
        if (minNum > elevationMin and total_ele37.index(minNum) > 0 ):
            elevationMinRange = total_ele37.index(minNum) - 1
        else:
            elevationMinRange = total_ele37.index(minNum)
        
        if (maxNum < elevationMin and total_ele37.index(maxNum) < 36 ):
            elevationMaxRange = total_ele37.index(maxNum) - 1
        else:
            elevationMaxRange = total_ele37.index(maxNum)
            
        elevation = []
        for e in range(elevationMinRange, elevationMaxRange + 1):
            elevation.append(total_ele37[e])
        
        elevation = [str(ele) for ele in elevation]
        elevation = '/'.join(elevation)
        
        return elevation
        
    
    def getDictionary(self):
        self.dictionary = {
                'dataset': 'ds628.0',
                'date': self.makeDate(),
                'param': '/'.join(self.param),
                'level': 'Isobaric surface:'+self.getPressureLevels(),
                'oformat': 'netCDF',
                'nlat': str(self.area['north']),
                'slat': str(self.area['south']),
                'wlon': str(self.area['west']),
                'elon': str(self.area['east']),
                'product': 'Analysis',
                'compression': 'NN',
                'gridproj': 'latLon',
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'
                }
        
        return self.dictionary
               

class JRAsa(object):
    
    def __init__(self, date, area, variables, rda):
        '''Returns an object for JRA55 data that has methods for querying the
        NCAR server for surface forecast variables (prec, swin, lwin). '''
        
        self.date = date
        self.area = area
        #self.directory = directory

        dpar = {'air_temperature'  : ['Temperature'],
                'relative_humidity': ['Relative humidity'],
                'wind_speed'       : ['u-component of wind', 
                                      'v-component of wind']}
        
        self.param = self.getParam(dpar, variables)
        
    def getParam(self, dpar, variables):
        
        varlist = [] 
        for var in variables:
            varlist.append(dpar.get(var))

        varlist = [item for item in varlist if item is not None]
        varlist = [item for sublist in varlist for item in sublist]
        
        return varlist
    
    def makeDate(self):
        '''convert data format to NCAR RDA request'''
        
        beg = self.date['beg'].strftime('%Y%m%d%H%M')
        end = self.date['end'].strftime('%Y%m%d%H%M')
        dateRange = beg + '/to/' + end
        
        return dateRange
        
    
    def getDictionary(self):
        self.dictionary = {
                'dataset': 'ds628.0',
                'date': self.makeDate(),
                'param': '/'.join(self.param),
                'level': 'Specified height above ground:2/10',
                'oformat': 'netCDF',
                'nlat': str(self.area['north']),
                'slat': str(self.area['south']),
                'wlon': str(self.area['west']),
                'elon': str(self.area['east']),
                'product': 'Analysis',
                'compression': 'NN',
                'gridproj': 'latLon',
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'
                }
        
        return self.dictionary
    

class JRAsf(object):
    
    def __init__(self, date, area, variables, rda):
        '''Returns an object for JRA55 data that has methods for querying the
        NCAR server for surface forecast variables (prec, swin, lwin). '''
        
        self.date = date
        self.area = area
        
        dpar = {'precipitation_amount'              :
                    ['Total precipitation'],
                'downwelling_shortwave_flux_in_air' : 
                    ['Downward solar radiation flux'],
                'downwelling_longwave_flux_in_air'  : 
                    ['Downward longwave radiation flux'],
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': 
                    ['Clear sky downward solar radiation flux'],
                'downwelling_longwave_flux_in_air_assuming_clear_sky':
                    ['Clear sky downward longwave radiation flux']}
        
        self.param = self.getParam(dpar, variables)
        
    def getParam(self, dpar, variables):
        
        varlist = [] 
        for var in variables:
            varlist.append(dpar.get(var))

        varlist = [item for item in varlist if item is not None]
        varlist = [item for sublist in varlist for item in sublist]
        varlist.append('Pressure')
        
        return varlist
    
    def makeDate(self):
        '''convert data format to NCAR RDA request'''
        
        beg = self.date['beg'].strftime('%Y%m%d%H%M')
        end = self.date['end'].strftime('%Y%m%d%H%M')
        dateRange = beg + '/to/' + end
        
        return dateRange
        
    
    def getDictionary(self):
        self.dictionary = {
                'dataset': 'ds628.0',
                'date': self.makeDate(),
                'param': '/'.join(self.param),
                'level': 'Ground or water surface:0',
                'oformat': 'netCDF',
                'nlat': str(self.area['north']),
                'slat': str(self.area['south']),
                'wlon': str(self.area['west']),
                'elon': str(self.area['east']),
                'product': '3-hour Average (initial+0 to initial+3)',
                'compression': 'NN',
                'gridproj': 'latLon',
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'
                }
        
        return self.dictionary


class JRAdownload(GenericDownload):
    '''Return an objet to download JRA55 dataset client based on RDA'''
    
    #TODO add credential
    def __init__(self, pfile):
        super().__init__(pfile)
        par = self.par
        self._set_data_directory("jra55")
        
        self.dsID = 'ds628.0'

        self.__varCheck(par)
                 
        # time bounds
        self.date  = self.getDate(par)
        
        self.credential = path.join(par.credentials_directory, ".jrarc")
        #print(self.credential)
        self.account = open(self.credential, "r")
        self.inf = self.account.readlines()
        self.username = ''.join(self.inf[0].split())
        self.password = ''.join(self.inf[1].split()) 
            
        # chunk size for downloading and storing data [days]        
        self.chunk_size = par.chunk_size*2000
        
        self.ncfVar  = {
                'initial_time0_hours':   'time', 
                'initial_time0':         'time',
                'initial_time0_encoded': 'time',
                'lv_ISBL1':              'level', 
                'g0_lat_1':              'latitude', 
                'g0_lon_2':              'longitude',
                'g0_lat_2':              'latitude', 
                'g0_lon_3':              'longitude',
                'HGT_GDS0_ISBL':         'Geopotential height',
                'RH_GDS0_ISBL':          'Relative humidity',
                'VGRD_GDS0_ISBL':        'v-component of wind',
                'UGRD_GDS0_ISBL':        'u-component of wind',
                'TMP_GDS0_ISBL':         'Temperature',
                'TMP_GDS0_HTGL':         'Temperature',
                'VGRD_GDS0_HTGL':        'v-component of wind',
                'UGRD_GDS0_HTGL':        'u-component of wind',
                'RH_GDS0_HTGL':          'Relative humidity',
                'PRES_GDS0_SFC_ave3h':   'Pressure',
                'TPRAT_GDS0_SFC_ave3h':  'Total precipitation',
                'CSDSF_GDS0_SFC_ave3h':  'Clear sky downward solar radiation flux',
                'CSDLF_GDS0_SFC_ave3h':  'Clear sky downward longwave radiation flux',
                'DSWRF_GDS0_SFC_ave3h':  'Downward solar radiation flux',
                'DLWRF_GDS0_SFC_ave3h':  'Downward longwave radiation flux'}

    
    def __varCheck(self, par):
        '''convert one variable to a list'''
        
        if not isinstance(par.variables, (list,)):
            par.variables = [par.variables]
            
    def getDate(self, par):
        '''get download daterange'''
        
        dateRange = {'beg' : par.beg, 'end' : par.end}
        dateRange['end'] = dateRange['end'] + timedelta(hours=23)
        
        return dateRange
        
    def getDataLev(self, dsi):
        '''get data level of the download data set'''
        
        flist = glob.glob(path.join(self.directory, '*'+dsi+'*'))
        lev = path.basename(flist[0]).split('_')[1].split('.')[0]
        if lev == 'p125':
            dataLev = 'pl'
        elif lev == 'surf125':
            dataLev = 'sa'
        elif lev == 'phy2m125':
            dataLev = 'sf'
        
        return dataLev
    
    def getVars(self, dsi):
        '''get all download variable names'''
        
        flist = glob.glob(path.join(self.directory, '*'+dsi+'*'))
        varlist = []
        for f in flist:
            fname = path.basename(f)
            var = fname.split('_')[2]
            varlist.append(var.split('.')[0])
        
        variables = np.unique(varlist)
        
        return variables
        
    def getOutFile(self, ncf, dataLev):
        
        times = nc.num2date(ncf[self.timeName][:], 
                            units = ncf[self.timeName].units, 
                            calendar='standard')
        begStr = np.min(times).strftime('%Y%m%d')
        endStr = np.max(times).strftime('%Y%m%d')
        
        fileName = ['jra55', dataLev, begStr, 'to', endStr]
        fileName = '_'.join(fileName) + '.nc'
        fileName = path.join(self.directory, fileName)
        
        return fileName
    
    def getDimName(self, dataLev):
        '''get the dimension [time, level, latitude, longitude] name in the 
        original JRA55 ncf'''
        
        self.timeName = 'initial_time0_hours'
        self.levName = 'lv_ISBL1'
        
        if dataLev == 'pl':
            self.lonName = 'g0_lon_3'
            self.latName = 'g0_lat_2'
        elif dataLev in ['sa', 'sf']:
            self.lonName = 'g0_lon_2'
            self.latName = 'g0_lat_1'
        
    
    def makeNCF(self, dsi):
        
        variables = self.getVars(dsi)
        dataLev = self.getDataLev(dsi)
        self.getDimName(dataLev)
        
        varf = np.sort(glob.glob(path.join(self.directory, 
                                           '*' + variables[0]+'*')))
        ncf = nc.MFDataset(varf.tolist(), aggdim ='initial_time0_hours')
        
        if dataLev == 'pl':
            Levs = ncf['lv_ISBL1'][:].data
            
        Times = ncf[self.timeName][:]
        Lats  = ncf[self.latName][:].data
        Lons  = ncf[self.lonName][:].data
        
        file_new = self.getOutFile(ncf, dataLev)
        
        #initialize new data file and create group
        ncn = nc.Dataset(file_new, 'w', format='NETCDF4_CLASSIC')
        
        #make dimensions
        if dataLev == 'pl':
            Levs = ncf[self.levName][:].data
            ncn.createDimension('level', len(Levs))
            levels     = ncn.createVariable('level',  'i4',('level',))
            levels.long_name  = 'pressure level'
            levels.units      = 'mbar' 
            levels[:] = Levs
        ncn.createDimension('time',  len(Times))
        ncn.createDimension('latitude',  len(Lats))
        ncn.createDimension('longitude', len(Lons))
        
        #make dimension variables
        times      = ncn.createVariable('time',    'd',('time',))
        latitudes  = ncn.createVariable('latitude',   'f8', ('latitude',))
        longitudes = ncn.createVariable('longitude',  'f8', ('longitude',))
        
        times.standard_name = 'time'
        times.units     = ncf[self.timeName].units
        times.calendar  = 'standard'
        latitudes.standard_name  = ncf[self.latName].long_name
        latitudes.units      = ncf[self.latName].units 
        longitudes.standard_name = ncf[self.lonName].long_name
        longitudes.units     = ncf[self.lonName].units
        
        ncf.close()
        
        #assign dimensions
        times[:] = Times
        longitudes[:] = Lons
        latitudes[:] = Lats
        
        for vari in variables:
            flist = np.sort(glob.glob(path.join(self.directory, '*'+vari+'*')))
            ncf = nc.MFDataset(flist.tolist(), aggdim = self.timeName)
            for n, var in enumerate(ncf.variables.keys()):
                if variables_skip(self.ncfVar[var]):
                    continue                 
                print("VAR: ", var)
                if dataLev == 'pl':
                    vari = ncn.createVariable(self.ncfVar[var],'f4',
                                              ('time','level',
                                               'latitude','longitude',))
                    vari[:,:,:,:] = ncf[var][:,:,:,:]
                else:
                    vari = ncn.createVariable(self.ncfVar[var],'f4',
                                              ('time',
                                               'latitude','longitude'))
                    vari[:,:,:] = ncf[var][:,:,:]
                    
                vari.long_name = ncf[var].long_name
                vari.units     = ncf[var].units
            
            
            # clear
            ncf.close()
            for f in flist:
                remove(f)
        
        ncn.close()
        
    def requestClear(self, rda):
        '''clear online datasets before downloading'''
        
        print('\n======== Clear Online Datasets From NCAR Server ========\n')
        dsIndex = rda.getDSindex()
        if len(dsIndex) > 1:
            rda.purge(dsIndex)
        print('\n======== Online Dateset Cleared ========\n')
        
    def requestSubmit(self, rda):
        
        print('\n======== Submit Request ========\n')
        
        slices = floor(float((self.date['end'] - self.date['beg']).days)/
                       self.chunk_size) + 1
        
        dsN = 0
        for ind in range (0, int(slices)): 
            date_i = {}
            #prepare time slices   
            date_i['beg'] = self.date['beg'] + timedelta(days = 
                            self.chunk_size * ind)
            date_i['end'] = self.date['beg'] + timedelta(days = 
                            self.chunk_size * (ind+1) - 1)
            if ind == (slices-1):
                date_i['end'] = self.date['end']
            
            pl = JRApl(date_i, self.area, self.elevation, 
                       self.variables, rda)
            sa = JRAsa(date_i, self.area, self.variables, rda)
            sf = JRAsf(date_i, self.area, self.variables, rda)
            
            # get download data level
            JRAli = []
            for jrai in [pl, sa, sf]:
                if len(jrai.param) > 0:
                    JRAli.append(jrai)
                    
            for jrai in JRAli:
                rda.submit(jrai.getDictionary())
                dsN += 1
        print('\n======== Submit Completed ========\n')  
        
        return dsN
        
    def requestDownload(self, rda, dsN):
        '''download all the request'''
        
        doneI = []
        while len(doneI) < dsN:
            
            print('\n======== Geting Available Dataset ========\n')
            dsIndex = rda.getDSindex()
            dsIndex = [item for item in dsIndex if item not in doneI]
            if len(dsIndex) > 0:
                for ds in dsIndex:
                    rda.download(self.directory, ds)
                    self.makeNCF(ds)
                doneI += dsIndex
            time.sleep(60*10)# check available data every 10 mins
        
        print('''\n======== Download Completed ========\n''')
        
    def retrieve(self):
        '''submit and download all the dataset'''
        
        print('''\n======== JRA55: STRAT ========\n''')
        
        rda = RDA(self.username, self.password) # initialize RDA server
        self.requestClear(rda) # clear online request
        dsN = self.requestSubmit(rda) # submit request
        self.requestDownload(rda, dsN) # download dataset
        
        print('''\n======== JRA55: STOP ========\n''')
        
         


class JRAinterpolate(object):
    """
    Algorithms to interpolate JRA55 netCDF files to station coordinates. 
    All variables retains their original units and time-steps. 
    
    Referenced from era_interim.py (Dr.Stephan Gruber): Class ERAinterpolate()     
    
    Args:
        ifile: Full path to a Globsim Interpolate Paramter file
        JRAinterpolate(ifile)


    Example:
        ifile = '/home/xquan/src/globsim/examples/par/examples.globsim_interpolate'
        JRAinterpolate(ifile)
      
    """

    def __init__(self, ifile):
        #read parameter file
        self.ifile = ifile
        par = ParameterIO(self.ifile)
        self.dir_inp = path.join(par.project_directory,'jra55') 
        self.dir_out = self.makeOutDir(par)
        self.variables = par.variables     
        self.list_name = par.station_list.split(path.extsep)[0]
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        
        #read station points 
        self.stations = StationListRead(self.stations_csv)  
        
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end + timedelta(days=1)}

        # chunk size: how many time steps to interpolate at the same time?
        # A small chunk size keeps memory usage down but is slow.
        self.cs  = int(par.chunk_size) * 200
    
    
    def makeOutDir(self, par):
        '''make directory to hold outputs'''
        
        dirIntp = path.join(par.project_directory, 'interpolated')
        
        if not (path.isdir(dirIntp)):
            makedirs(dirIntp)
            
        return dirIntp
    
    def JRAinterp2D(self, ncfile_in, ncf_in, points, tmask_chunk,
                    variables=None, date=None):    
        """
        Biliner interpolation from fields on regular grid (latitude, longitude) 
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all jra55 files).
          
        Args:
            ncfile_in: Full path to an JRA-55 derived netCDF file. This can
                        contain wildcards to point to multiple files if temporal
                        chunking was used.
              
            ncf_in: A netCDF4.MFDataset derived from reading in JRA-55 
                    multiple files (def JRA2station())
            
              
            points: A dictionary of locations. See method StationListRead in
                    generic.py for more details.
        
            variables:  List of variable(s) to interpolate such as 
                        [air_temperature, easteard_wind, northward_wind, relative_humidy, surface_temperature, 
                        downwelling_shortwave_flux_in_air, downwelling_longwave_flux_in_air,
                        downwelling_shortwave_flux_in_air_assuming_clear_sky, 
                        downwelling_longwave_flux_in_air_assuming_clear_sky].
                        Defaults to using all variables available.
        
            date: Directory to specify begin and end time for the derived time 
                  series. Defaluts to using all times available in ncfile_in.
              
        Example:
            from datetime import datetime
            date  = {'beg' : datetime(2008, 1, 1),
                      'end' : datetime(2008,12,31)}
            variables  = [air_temperature, easteard_wind, northward_wind]       
            stations = StationListRead("points.csv")      
            JRA2station('jra_surf_*.nc', 'jra_sa_inter.nc', stations, 
                        variables=variables, date=date)        
        """   
        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # get spatial dimensions
        if pl: # only for pressure level files
            lev  = ncf_in.variables['level'][:]
            nlev = len(lev)
              
        # test if time steps to interpolate remain
        nt = sum(tmask_chunk)
        if nt == 0:
            raise ValueError('No time steps from netCDF file selected.')
    
        # get variables
        varlist = [str_encode(x) for x in ncf_in.variables.keys()]
        varlist.remove('time')
        varlist.remove('latitude')
        varlist.remove('longitude')
        if pl: #only for pressure level files
            varlist.remove('level')
    
        #list variables that should be interpolated
        if variables is None:
            variables = varlist
        #test is variables given are available in file
        if (set(variables) < set(varlist) == 0):
            raise ValueError('One or more variables not in netCDF file.')           

        # Create source grid from a SCRIP formatted file. As ESMF needs one
        # file rather than an MFDataset, give first file in directory.
        ncsingle = filter(listdir(self.dir_inp), path.basename(ncfile_in))[0]
        ncsingle = path.join(self.dir_inp, ncsingle)
        sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)

        # create source field on source grid
        if pl: #only for pressure level files
            sfield = ESMF.Field(sgrid, name='sgrid',
                                staggerloc=ESMF.StaggerLoc.CENTER,
                                ndbounds=[len(variables), nt, nlev])
        else: # 2D files
            sfield = ESMF.Field(sgrid, name='sgrid',
                                staggerloc=ESMF.StaggerLoc.CENTER,
                                ndbounds=[len(variables), nt])

        # assign data from ncdf: (variable, time, latitude, longitude) 
        for n, var in enumerate(variables):
            if pl: # only for pressure level files
                if ESMFnew:
                    sfield.data[:,:,n,:,:] = ncf_in.variables[var][
                            tmask_chunk,:,:,:].transpose((3,2,0,1)) 
                else:
                    sfield.data[n,:,:,:,:] = ncf_in.variables[var][
                            tmask_chunk,:,:,:].transpose((0,1,3,2)) 
            else:
                if ESMFnew:
                    sfield.data[:,:,n,:] = ncf_in.variables[var][
                            tmask_chunk,:,:].transpose((2,1,0))
                else:
                    sfield.data[n,:,:,:] = ncf_in.variables[var][
                            tmask_chunk,:,:].transpose((0,2,1))

        # create locstream, CANNOT have third dimension!!!
        locstream = ESMF.LocStream(len(self.stations), 
                                   coord_sys=ESMF.CoordSys.SPH_DEG)
        locstream["ESMF:Lon"] = list(self.stations['longitude_dd'])
        locstream["ESMF:Lat"] = list(self.stations['latitude_dd'])

        # create destination field
        if pl: # only for pressure level files
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt, nlev])
        else:
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt])    

        # regridding function, consider ESMF.UnmappedAction.ERROR
        regrid2D = ESMF.Regrid(sfield, dfield,
                                regrid_method=ESMF.RegridMethod.BILINEAR,
                                unmapped_action=ESMF.UnmappedAction.IGNORE,
                                dst_mask_values=None)
                  
        # regrid operation, create destination field (variables, times, points)
        dfield = regrid2D(sfield, dfield)        
        sfield.destroy() #free memory                  

        return dfield, variables

    def JRA2station(self, ncfile_in, ncfile_out, points,
                    variables = None, date = None):
        
        """
        Biliner interpolation from fields on regular grid (latitude, longitude) 
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all JRA55 files). The type 
        of variable and file structure are determined from the input.
        
        This function creates an empty of netCDF file to hold the interpolated 
        results, by calling self.netCDF_empty(). Then, data is 
        interpolated in temporal chunks and appended. The temporal chunking can 
        be set in the interpolation parameter file.
        
        Args:
        ncfile_in: Full path to an JRA-55 derived netCDF file. This can
                   contain wildcards to point to multiple files if temporal
                  chunking was used.
            
        ncfile_out: Full path to the output netCDF file to write.     
        
        points: A dictionary of locations. See method StationListRead in
                generic.py for more details.
    
        variables:  List of variable(s) to interpolate such as 
                    [air_temperature, easteard_wind, northward_wind, 
                    relative_humidy, surface_temperature, 
                    downwelling_shortwave_flux_in_air, 
                    downwelling_longwave_flux_in_air,
                    downwelling_shortwave_flux_in_air_assuming_clear_sky, 
                    downwelling_longwave_flux_in_air_assuming_clear_sky].
                    Defaults to using all variables available.
    
        date: Directory to specify begin and end time for the derived time 
                series. Defaluts to using all times available in ncfile_in.
        
        cs: chunk size, i.e. how many time steps to interpolate at once. This 
            helps to manage overall memory usage (small cs is slower but less
            memory intense).          
        """
        
        # read in one type of mutiple netcdf files       
        ncf_in = nc.MFDataset(ncfile_in, 'r', aggdim ='time')
        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # build the output of empty netCDF file
        rootgrp = create_empty_netcdf(ncfile_out, self.stations, ncf_in, 
                                      time_units='hours since 1800-01-01 00:00:0.0') 
        rootgrp.source = 'JRA55, interpolated bilinearly to stations'
        rootgrp.close()  
        
        # open the output netCDF file, set it to be appendable ('a')
        ncf_out = nc.Dataset(ncfile_out, 'a')

        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]
        #"hours since 1900-01-01 00:00:0.0"
        t_unit = ncf_in.variables['time'].units 
        try :
            t_cal = ncf_in.variables['time'].calendar
        except AttributeError :  # attribute doesn't exist
            t_cal = u"gregorian" # standard
        time = nc.num2date(nctime, units = t_unit, calendar = t_cal)
        
        # detect invariant files (topography etc.)
        if len(time) ==1:
            invariant=True
        else:
            invariant=False   
        
        # restrict to date/time range if given
        if date is None:
            tmask = time < datetime(3000, 1, 1)
        else:
            tmask = (time < date['end']) * (time >= date['beg'])
                              
        # get time indices
        time_in = nctime[tmask]     

        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter  = len(time_in) // self.cs
        niter += ((len(time_in) % self.cs) > 0)

        # loop over chunks
        for n in range(niter):
            # indices
            beg = n * self.cs
            # restrict last chunk to lenght of tmask plus one (to get last time)
            end = min(n*self.cs + self.cs, len(time_in)) -1
            
            # time to make tmask for chunk 
            beg_time = nc.num2date(time_in[beg], units=t_unit, calendar=t_cal)
            if invariant:
                # allow topography to work in same code, len(nctime) = 1
                end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                #end = 1
            else:
                end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)
                
            #'<= end_time', would damage appending
            tmask_chunk = (time <= end_time) * (time >= beg_time)
            if invariant:
                # allow topography to work in same code
                tmask_chunk = [True]
                 
            # get the interpolated variables
            dfield, variables = self.JRAinterp2D(ncfile_in, ncf_in, 
                                                     self.stations, tmask_chunk,
                                                     variables=None, date=None) 

            # append time
            ncf_out.variables['time'][:] = np.append(ncf_out.variables['time'][:], 
                                                     time_in[beg:end+1])
                                  
            #append variables
            for i, var in enumerate(variables):
                if variables_skip(var):
                    continue
                                                              
                if pl:
                    lev = ncf_in.variables['level'][:]
                    
                    if ESMFnew:
                        ncf_out.variables[var][beg:end+1,:,:] = dfield.data[:,i,:,:].transpose((1,2,0))
                    else:
                        # dimension: time, level, latitude, longitude
                        ncf_out.variables[var][beg:end+1,:,:] = dfield.data[i,:,:,:]      
                else:
                    if ESMFnew:
                        ncf_out.variables[var][beg:end+1,:] = dfield.data[:,i,:].transpose((1,0))
                    else:
                        # time, latitude, longitude
                        ncf_out.variables[var][beg:end+1,:] = dfield.data[i,:,:]
                    
                    
        #close the file
        ncf_in.close()
        ncf_out.close()

   
    def levels2elevation(self, ncfile_in, ncfile_out):    
        """
        Linear 1D interpolation of pressure level data available for individual
        stations to station elevation. Where and when stations are below the 
        lowest pressure level, they are assigned the value of the lowest 
        pressure level.
        
        """
        # open file 
        ncf = nc.MFDataset(ncfile_in, 'r', aggdim='time')
        height = ncf.variables['height'][:]
        nt = len(ncf.variables['time'][:])
        nl = len(ncf.variables['level'][:])
        
        # list variables
        varlist = [str_encode(x) for x in ncf.variables.keys()]
        for V in ['time', 'station', 'latitude', 'longitude', 'level', 'height']:
            varlist.remove(V)

        # === open and prepare output netCDF file =============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'JRA55, interpolated (bi)linearly to stations'
        rootgrp.featureType = "timeSeries"

        # dimensions
        station = rootgrp.createDimension('station', len(height))
        time    = rootgrp.createDimension('time', nt)

        # base variables
        time           = rootgrp.createVariable('time',     'i4',('time'))
        time.long_name = 'time'
        time.units     = 'hours since 1800-01-01 00:00:0.0'
        time.calendar  = 'gregorian'
        station             = rootgrp.createVariable('station',  'i4',
                                                     ('station'))
        station.long_name   = 'station for time series data'
        station.units       = '1'
        latitude            = rootgrp.createVariable('latitude', 'f4',
                                                     ('station'))
        latitude.long_name  = 'latitude'
        latitude.units      = 'degrees_north'    
        longitude           = rootgrp.createVariable('longitude','f4',
                                                     ('station'))
        longitude.long_name = 'longitude'
        longitude.units     = 'degrees_east'  
        height           = rootgrp.createVariable('height','f4',
                                                  ('station'))
        height.long_name = 'height_above_reference_ellipsoid'
        height.units     = 'm'  
       
        # assign base variables
        time[:] = ncf.variables['time'][:]
        station[:]   = ncf.variables['station'][:]
        latitude[:]  = ncf.variables['latitude'][:]
        longitude[:] = ncf.variables['longitude'][:]
        height[:]    = ncf.variables['height'][:]
        
        # create and assign variables from input file
        for var in varlist:
            vname = str_encode(ncf.variables[var].long_name)
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))    
            tmp.long_name = str_encode(ncf.variables[var].long_name)
            tmp.units     = str_encode(ncf.variables[var].units)
            
        # add air pressure as new variable
        var = 'air_pressure'
        varlist.append(var)
        tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))    
        tmp.long_name = str_encode(var)
        tmp.units     = str_encode('hPa') 
        # end file preparation ================================================
    
                                                                                                
        # loop over stations
        for n, h in enumerate(height): 
            # convert geopotential [millibar] to height [m]
            # shape: (time, level)
            ele = ncf.variables['Geopotential height'][:,:,n]
            
            # difference in elevation. 
            # level directly above will be >= 0
            dele = ele - h
            # vector of level indices that fall directly above station. 
            # Apply after ravel() of data.
            va = np.argmin(dele + (dele < 0) * 100000, axis=1) 
            # mask for situations where station is below lowest level
            mask = va < (nl-1)
            va += np.arange(ele.shape[0]) * ele.shape[1]
            
            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask # +1 when OK, +0 when below lowest level
            
            # weights
            wa = np.absolute(dele.ravel()[vb]) 
            wb = np.absolute(dele.ravel()[va])
            
            wt = wa + wb
            
            wa /= wt # Apply after ravel() of data.
            wb /= wt # Apply after ravel() of data.
            
            #loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                if var == 'air_pressure':
                    # pressure [Pa] variable from levels, shape: (time, level)
                    data = np.repeat([ncf.variables['level'][:]],
                                      len(time),axis=0).ravel()                 
                else:
                    #read data from netCDF
                    data = ncf.variables[var][:,:,n].ravel()
            
                multvawa = np.multiply(data[va], wa)
                multvbwb = np.multiply(data[vb], wb)
                ipol = multvawa + multvbwb
                rootgrp.variables[var][:,n] = ipol # assign to file   
    
        rootgrp.close()
        ncf.close()
        # closed file =========================================================


    def TranslateCF2short(self, dpar):
        """
        Map CF Standard Names into short codes used in JRA-55 netCDF files.
        """
        varlist = [] 
        for var in self.variables:
            varlist.append(dpar.get(var))
        # drop none
        varlist = [item for item in varlist if item is not None]      
        # flatten
        varlist = [item for sublist in varlist for item in sublist]         
        return(varlist) 

    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to 
        the more generically JRA-like interpolation functions.
        """                       

        # === 2D Interpolation for Surface  Data ===    
        # dictionary to translate CF Standard Names into JRA55
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['surface_temperature'], # [K] 2m values
                'relative_humidity' : ['relative_humidity'], # [%]                                                       
                'wind_speed' : ['eastward_wind', 'northward_wind']} # [m s-1] 2m & 10m values   
        varlist = self.TranslateCF2short(dpar)
        self.JRA2station(path.join(self.dir_inp,'jra55_sa_*.nc'),
                            path.join(self.dir_out,'jra_sa_' + 
                                      self.list_name + '.nc'), self.stations,
                                      varlist, date = self.date)
        
        # 2D Interpolation for Radiation Data 
        # dictionary to translate CF Standard Names into JRA55
        # pressure level variable keys.       
        dpar = {'downwelling_shortwave_flux_in_air': ['downwelling_shortwave_flux_in_air'], # [W/m2] short-wave downward
                'downwelling_longwave_flux_in_air' : ['downwelling_longwave_flux_in_air'], # [W/m2] long-wave downward
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['downwelling_shortwave_flux_in_air_assuming_clear_sky'], # [W/m2] short-wave downward assuming clear sky
                'downwelling_longwave_flux_in_air_assuming_clear_sky': ['downwelling_longwave_flux_in_air_assuming_clear_sky'],
                'precipitation_amount' : ['total_precipitation']} # [W/m2] long-wave downward assuming clear sky
        varlist = self.TranslateCF2short(dpar)                           
        self.JRA2station(path.join(self.dir_inp,'jra55_sf_*.nc'), 
                          path.join(self.dir_out,'jra_sf_' + 
                                    self.list_name + '.nc'), self.stations,
                                    varlist, date = self.date)          

        # 2D Interpolation for Pressure-Level, Analyzed Meteorological DATA
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['air_temperature'],  # [K]
                'relative_humidity' : ['relative_humidity'], # [%]
                'wind_speed'        : ['eastward_wind', 'northward_wind']}  # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('geopotential_height')
        self.JRA2station(path.join(self.dir_inp,'jra55_pl_*.nc'), 
                          path.join(self.dir_out,'jra_pl_' + 
                                    self.list_name + '.nc'), self.stations,
                                    varlist, date = self.date)  
                                                                                          
        # 1D Interpolation for Pressure Level Data 
        self.levels2elevation(path.join(self.dir_out,'jra_pl_' + 
                                        self.list_name + '.nc'), 
                              path.join(self.dir_out,'jra_pl_' + 
                                        self.list_name + '_surface.nc'))


   
class JRAscale(object):
    """
    Class for JRA-55 data that has methods for scaling station data to
    better resemble near-surface fluxes.
    
    Processing kernels have names in UPPER CASE.
       
    Args:
        sfile: Full path to a Globsim Scaling Parameter file. 
              
    Example:          
        JRAd = JRAscale(sfile) 
        JRAd.process()
    """
        
    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        par = ParameterIO(self.sfile)
        self.intpdir = path.join(par.project_directory, 'interpolated')
        self.scdir = self.makeOutDir(par)
        self.list_name = par.station_list.split(path.extsep)[0]
        
        # read kernels
        self.kernels = par.kernels
        if not isinstance(self.kernels, list):
            self.kernels = [self.kernels]
            
        # input file names
        self.nc_pl = nc.Dataset(path.join(self.intpdir,'jra_pl_' + 
                                self.list_name + '_surface.nc'), 'r')
        self.nc_sa = nc.Dataset(path.join(self.intpdir,'jra_sa_' + 
                                self.list_name + '.nc'), 'r')
        self.nc_sf = nc.Dataset(path.join(self.intpdir,'jra_sf_' + 
                                self.list_name + '.nc'), 'r')
                               
        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(par, 'jra55')
            
        # time vector for output data
        # get time and convert to datetime object
        nctime = self.nc_pl.variables['time'][:]
        self.t_unit = self.nc_pl.variables['time'].units
        self.t_cal  = self.nc_pl.variables['time'].calendar
        time = nc.num2date(nctime, units = self.t_unit, calendar = self.t_cal)
        
        #number of time steps
        self.nt = floor((max(time)-min(time)).total_seconds()/3600/par.time_step)+1
        self.time_step = par.time_step
                              
        # vector of output time steps as datetime object
        self.times_out    = [min(time) + timedelta(hours=x*par.time_step) 
                                for x in range(0, self.nt)]
        # vector of output time steps as written in ncdf file
        self.times_out_nc = nc.date2num(self.times_out, units = self.t_unit, 
                                        calendar = self.t_cal)
        # get the station file
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        #read station points 
        self.stations = StationListRead(self.stations_csv)
        
        
    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one 
        variable and adds it to the netCDF file.
        """    
        self.rg = ScaledFileOpen(self.output_file, self.nc_pl, 
                                 self.times_out_nc, self.nc_pl['time'].units)
        
        # add station names to netcdf
        # first convert to character array
        names_out = nc.stringtochar(np.array(self.stations['station_name'], 'S32'))
        
        # create space in the netcdf
        nchar        = self.rg.createDimension('name_strlen', 32) 
        st           = self.rg.createVariable('station_name', "S1", 
                                              ('station', 'name_strlen'))
        st.standard_name = 'platform_name'
        st.units     = ''
        
        # add data
        st[:] = names_out
        
        # iterate through kernels and start process
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                print(kernel_name)
                getattr(self, kernel_name)()
            
        # self.conv_geotop()
            
        # close netCDF files   
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
        
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


    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """        
        # add variable to ncdf file
        vn = 'PRESS_JRA55_Pa_pl' # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))    
        var.long_name = 'air_pressure JRA-55 pressure levels only'
        var.units     = 'Pa'
        var.standard_name = 'surface_air_pressure'
        
        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:].astype(np.int64)  
        values  = self.nc_pl.variables['air_pressure'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            #scale from hPa to Pa 
            self.rg.variables[vn][:, n] = series_interpolate(self.times_out_nc, 
                                        time_in, values[:, n]) * 100

    def AIRT_C_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """        
        vn = 'AIRT_JRA55_C_pl' # variable name
        var           = self.rg.createVariable(vn,'f4',('time','station'))    
        var.long_name = 'air_temperature JRA55 pressure levels only'
        var.units     = 'degrees_C'
        var.standard_name = 'air_temperature'
        
        # interpolate station by station
        time_in = self.nc_pl.variables['time'][:]
        values  = self.nc_pl.variables['Temperature'][:]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n])-273.15            
                
 
    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """   
        
        # add variable to ncdf file
        vn = 'AIRT_JRA55_C_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '2_metre_temperature JRA55 surface only'
        var.units     = 'degrees_C'
        var.standard_name = 'air_temperature'
        
        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['Temperature'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n])-273.15  
    def RH_per_sur(self):
        """
        Relative Humidity derived from surface data, exclusively.
        """        
        # add variable to ncdf file
        vn = 'RH_JRA55_per_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'relative humidity JRA55 surface only'
        var.units     = 'percent'
        var.standard_name = 'relative_humidity'
        
        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['Relative humidity'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n])  

    def WIND_sur(self):
        """
        Wind at 10 metre derived from surface data, exclusively.
        """   
        
        # add variable to ncdf file
        vn = '10 metre U wind component' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '10 metre U wind component'
        var.units     = str_encode(self.nc_sa.variables['u-component of wind'].units)  
        
        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['u-component of wind'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n]) 

        # add variable to ncdf file
        vn = '10 metre V wind component' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '10 metre V wind component'
        var.units     = str_encode(self.nc_sa.variables['v-component of wind'].units)  
        
        # interpolate station by station
        time_in = self.nc_sa.variables['time'][:]
        values  = self.nc_sa.variables['v-component of wind'][:]                   
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n]) 

        # add variable to ncdf file
        vn = 'WSPD_JRA55_ms_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '10 metre wind speed JRA55 surface only'
        var.units     = 'm s-1'  
        var.standard_name = 'wind_speed'
        
        # add variable to ncdf file
        vn = 'WDIR_JRA55_deg_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = '10 metre wind direction JRA55 surface only'
        var.units     = 'degree'  
        var.standard_name = 'wind_from_direction'
                                
        # convert
        # u is the ZONAL VELOCITY, i.e. horizontal wind TOWARDS EAST.
        # v is the MERIDIONAL VELOCITY, i.e. horizontal wind TOWARDS NORTH.
        V = self.rg.variables['10 metre V wind component'][:]
        U = self.rg.variables['10 metre U wind component'][:] 

        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            WS = np.sqrt(np.power(V,2) + np.power(U,2))
            WD = [atan2(V[i, n], U[i, n])*(180/pi) + 
                  180 for i in np.arange(V.shape[0])]
            self.rg.variables['WSPD_JRA55_ms_sur'][:, n] = WS
            self.rg.variables['WDIR_JRA55_deg_sur'][:,n] = WD
    

    def SW_Wm2_sur(self):
        """
        solar radiation downwards derived from surface data, exclusively.
        """   
        
        # add variable to ncdf file
        vn = 'SW_JRA55_Wm2_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Surface solar radiation downwards JRA-55 surface only'
        var.units     = 'W m-2'
        var.standard_name = 'surface_downwelling_shortwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:]
        values  = self.nc_sf.variables['Downward solar radiation flux'][:]                                
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n]) 

    def LW_Wm2_sur(self):
        """
        Long-wave radiation downwards derived from surface data, exclusively.
        """   
        
        # add variable to ncdf file
        vn = 'LW_JRA55_Wm2_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Surface thermal radiation downwards JRA-55 surface only'
        var.units     = 'W m-2' 
        var.standard_name = 'surface_downwelling_longwave_flux'

        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:]
        values  = self.nc_sf.variables['Downward longwave radiation flux'][:]                                
        for n, s in enumerate(self.rg.variables['station'][:].tolist()):  
            self.rg.variables[vn][:, n] = np.interp(self.times_out_nc, 
                                                    time_in, values[:, n]) 

    def PREC_mm_sur(self):
        """
        Precipitation derived from surface data, exclusively.
        Convert unit: mm/day to mm/time_step (hours)
        """   
        
        # add variable to ncdf file
        vn = 'PREC_JRA55_mm_sur' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Total precipitation JRA55 surface only'
        var.units     = 'kg m-2 s-1'
        var.standard_name = 'precipitation_amount'
        
        # interpolate station by station
        time_in = self.nc_sf.variables['time'][:]
        values  = self.nc_sf.variables['Total precipitation'][:]/24/3600 #[mm/h]
        for n, s in enumerate(self.rg.variables['station'][:].tolist()): 
            f = interp1d(time_in, values[:, n], kind = 'linear')
            self.rg.variables[vn][:, n] = f(self.times_out_nc) * 3600           

    def LW_Wm2_topo(self):
        """
        Long-wave radiation downwards [W/m2]
        https://www.geosci-model-dev.net/7/387/2014/gmd-7-387-2014.pdf
        """             
        # get sky view, and interpolate in time
        N = np.asarray(list(self.stations['sky_view'][:]))

        # add variable to ncdf file
        vn = 'LW_JRA55_Wm2_topo' # variable name
        var           = self.rg.createVariable(vn,'f4',('time', 'station'))    
        var.long_name = 'Incoming long-wave radiation JRA-55 surface only'
        var.units     = str_encode('W m-2') 
        var.standard_name = 'surface_downwelling_longwave_flux'

        # compute                            
        for i in range(0, len(self.rg.variables['RH_JRA55_per_sur'][:])):
            for n, s in enumerate(self.rg.variables['station'][:].tolist()):
                LW = LW_downward(self.rg.variables['RH_JRA55_per_sur'][i, n],
                     self.rg.variables['AIRT_JRA55_C_sur'][i, n]+273.15, N[n])
                self.rg.variables[vn][i, n] = LW

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg]
        https://crudata.uea.ac.uk/cru/pubs/thesis/2007-willett/2INTRO.pdf
        '''
        print("Warning: SH_JRA_kgkg_sur is not defined. Specific humidity data are not currently available")
