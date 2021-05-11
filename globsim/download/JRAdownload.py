#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import glob
import http.cookiejar
import json
import netCDF4 as nc
import numpy as np
import sys
import tarfile
import time
import urllib.error
import urllib.request

from datetime import datetime, timedelta
from math import floor
from os import path, remove, makedirs

from globsim.common_utils import variables_skip
from globsim.download.GenericDownload import GenericDownload
from globsim.meteorology import pressure_from_elevation


def get_userinfo():
    return None, None


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
            tar.extractall(path=directory)
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
                            'RDA username and password invalid, or you are'
                            'not authorized to access this dataset.\n')
                        print('Please verify your login information at'
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
        _ = openrf.open(frequest)
        openerf = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(cj))
        urllib.request.install_opener(openerf)

    def getHelp(self):
        '''get the help information'''

        theurl = self.base + 'help'
        url = self.urlOpen(theurl)

        print(url.read().decode())

    def getSummary(self, dsID=None):
        '''get the summary of given dataset'''

        print('\nGetting summary information. Please wait as this may take'
              'a while.\n')

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

        print('\nGetting parameter summary. Please wait as this may take'
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
                            'RDA username and password invalid, or you are'
                            'not authorized to access this dataset.\n')
                        print('Please verify your login information at'
                              'http://rda.ucar.edu\n.')
                        sys.exit()

        print(url.read().decode())

    def downloadSinglefile(self, remfile, outfile):
        '''download_file(remfile,outfile) : download a file from a remote
         server (remfile) to a local location (outfile)'''

        _ = urllib.request.Request(remfile)
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

        authdata = 'email=' + self.username + '&password=' + self.password + '&action=login'
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

        authdata = 'email=' + self.username + '&password=' + self.password + '&action=login'
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
                    print('RDA username and password invalid. Please try again\n')
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

        minNum = min(total_ele37, key=lambda x:abs(x - elevationMin))
        maxNum = min(total_ele37, key=lambda x:abs(x - elevationMax))

        if (minNum > elevationMin and total_ele37.index(minNum) > 0):
            elevationMinRange = total_ele37.index(minNum) - 1
        else:
            elevationMinRange = total_ele37.index(minNum)

        if (maxNum < elevationMin and total_ele37.index(maxNum) < 36):
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
                'level': 'Isobaric surface:' + self.getPressureLevels(),
                'oformat': 'netCDF',
                'nlat': str(self.area['north']),
                'slat': str(self.area['south']),
                'wlon': str(self.area['west']),
                'elon': str(self.area['east']),
                'product': 'Analysis',
                'compression': 'NN',
                'gridproj': 'latLon',
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'}

        return self.dictionary


class JRAsa(object):

    def __init__(self, date, area, variables, rda):
        '''Returns an object for JRA55 data that has methods for querying the
        NCAR server for surface forecast variables (prec, swin, lwin). '''

        self.date = date
        self.area = area

        dpar = {'air_temperature'  : ['Temperature'],
                'relative_humidity': ['Relative humidity'],
                'specific_humidity': ['Specific humidity'],
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
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'}

        return self.dictionary


class JRAsf(object):

    def __init__(self, date, area, variables, rda):
        '''Returns an object for JRA55 data that has methods for querying the
        NCAR server for surface forecast variables (prec, swin, lwin). '''

        self.date = date
        self.area = area

        dpar = {'precipitation_amount':
                    ['Total precipitation'],
                'downwelling_shortwave_flux_in_air':
                    ['Downward solar radiation flux'],
                'downwelling_longwave_flux_in_air':
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
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'}

        return self.dictionary


class JRAdownload(GenericDownload):
    '''Return an objet to download JRA55 dataset client based on RDA'''

    # TODO: add credential
    def __init__(self, pfile):
        super().__init__(pfile)
        par = self.par
        self.__set_input_directory("jra55")

        self.dsID = 'ds628.0'

        self.__varCheck(par)

        # time bounds
        self.date  = self.getDate(par)

        self.credential = path.join(par['credentials_directory'], ".jrarc")
        self.account = open(self.credential, "r")
        self.inf = self.account.readlines()
        self.username = ''.join(self.inf[0].split())
        self.password = ''.join(self.inf[1].split())

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size'] * 2000

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
                'SPFH_GDS0_HTGL':        'Specific humidity',
                'PRES_GDS0_SFC_ave3h':   'Pressure',
                'TPRAT_GDS0_SFC_ave3h':  'Total precipitation',
                'CSDSF_GDS0_SFC_ave3h':  'Clear sky downward solar radiation flux',
                'CSDLF_GDS0_SFC_ave3h':  'Clear sky downward longwave radiation flux',
                'DSWRF_GDS0_SFC_ave3h':  'Downward solar radiation flux',
                'DLWRF_GDS0_SFC_ave3h':  'Downward longwave radiation flux'}

    def __varCheck(self, par):
        '''convert one variable to a list'''

        if not isinstance(par['variables'], (list,)):
            par['variables'] = [par['variables']]

    def getDate(self, par):
        '''get download daterange'''

        dateRange = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                     'end': datetime.strptime(par['end'], '%Y/%m/%d')}
        dateRange['end'] = dateRange['end'] + timedelta(hours=23)

        return dateRange

    def getDataLev(self, dsi):
        '''get data level of the download data set'''

        flist = glob.glob(path.join(self.directory, '*' + dsi + '*'))
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

        flist = glob.glob(path.join(self.directory, '*' + dsi + '*'))
        varlist = []
        for f in flist:
            fname = path.basename(f)
            var = fname.split('_')[2]
            varlist.append(var.split('.')[0])

        variables = np.unique(varlist)

        return variables

    def getOutFile(self, ncf, dataLev):

        times = nc.num2date(ncf[self.timeName][:],
                            units=ncf[self.timeName].units,
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
                                           f'*{variables[0]}*')))
        ncf = nc.MFDataset(varf.tolist(), aggdim='initial_time0_hours')

        if dataLev == 'pl':
            Levs = ncf['lv_ISBL1'][:].data

        Times = ncf[self.timeName][:]
        Lats  = ncf[self.latName][:].data
        Lons  = ncf[self.lonName][:].data

        file_new = self.getOutFile(ncf, dataLev)

        # initialize new data file and create group
        ncn = nc.Dataset(file_new, 'w', format='NETCDF4_CLASSIC')

        # make dimensions
        if dataLev == 'pl':
            Levs = ncf[self.levName][:].data
            ncn.createDimension('level', len(Levs))
            levels     = ncn.createVariable('level', 'i4',('level',))
            levels.long_name  = 'pressure level'
            levels.units      = 'mbar'
            levels[:] = Levs
        ncn.createDimension('time', len(Times))
        ncn.createDimension('latitude', len(Lats))
        ncn.createDimension('longitude', len(Lons))

        # make dimension variables
        times      = ncn.createVariable('time', 'd',('time',))
        latitudes  = ncn.createVariable('latitude', 'f8', ('latitude',))
        longitudes = ncn.createVariable('longitude', 'f8', ('longitude',))

        times.standard_name = 'time'
        times.units     = ncf[self.timeName].units
        times.calendar  = 'standard'
        latitudes.standard_name  = ncf[self.latName].long_name
        latitudes.units      = ncf[self.latName].units
        longitudes.standard_name = ncf[self.lonName].long_name
        longitudes.units     = ncf[self.lonName].units

        ncf.close()

        # assign dimensions
        times[:] = Times
        longitudes[:] = Lons
        latitudes[:] = Lats

        for vari in variables:
            flist = np.sort(glob.glob(path.join(self.directory, f'*{vari}*')))
            ncf = nc.MFDataset(flist.tolist(), aggdim=self.timeName)
            for n, var in enumerate(ncf.variables.keys()):
                if variables_skip(self.ncfVar[var]):
                    continue
                print("VAR: ", var)
                if dataLev == 'pl':
                    vari = ncn.createVariable(self.ncfVar[var], 'f4',
                                              ('time', 'level',
                                               'latitude', 'longitude',))
                    vari[:,:,:,:] = ncf[var][:,:,:,:]
                else:
                    vari = ncn.createVariable(self.ncfVar[var],'f4',
                                              ('time',
                                               'latitude', 'longitude'))
                    vari[:,:,:] = ncf[var][:,:,:]

                vari.long_name = ncf[var].long_name
                vari.units     = ncf[var].units

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

        slices = floor(float((self.date['end'] - self.date['beg']).days) /
                       self.chunk_size) + 1

        dsN = 0
        for ind in range(0, int(slices)):
            date_i = {}
            # prepare time slices
            date_i['beg'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * ind))
            date_i['end'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * (ind + 1) - 1))
            if ind == (slices - 1):
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
            time.sleep(60 * 10)  # check available data every 10 mins

        print('''\n======== Download Completed ========\n''')

    def retrieve(self):
        '''submit and download all the dataset'''

        print('''\n======== JRA55: STRAT ========\n''')

        rda = RDA(self.username, self.password)  # initialize RDA server
        self.requestClear(rda)  # clear online request
        dsN = self.requestSubmit(rda)  # submit request
        self.requestDownload(rda, dsN)  # download dataset

        print('''\n======== JRA55: STOP ========\n''')

