
import sys
import urllib.request
import urllib.error
import urllib.parse
import os
import getpass
import http.cookiejar
import json
import tarfile
import time

import numpy as np

from math import exp, floor
from datetime     import datetime, timedelta
from generic import ParameterIO


class RDA(object):
    
    def __init__(self, username, password):
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
        
        filename, file_extension = os.path.splitext(tarf)
        
        if file_extension == '.tar':
            tar = tarfile.open(tarf)
            tar.extractall(path = directory)
            tar.close()
            os.remove(tarf)
    
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
        
        theurl = self.base + 'request'
        
        jsondata = '{'
        for k in list(controlparms.keys()):
            jsondata += '"' + k + '"' + ":" + '"' + controlparms[k] + '",'
        jsondata = jsondata[:-1]
        jsondata += '}'
        print('\nSubmitting request. Please wait as this may take awhile.\n')
        
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
        if not os.path.exists(directory):
            os.makedirs(directory)
        for key, value in filelist.items():
            downloadpath, localfile = key.rsplit("/", 1)
            outpath = directory + backslash + localfile
            percentcomplete = (float(filecount) / float(length))
            self.update_progress(percentcomplete, directory)
            if os.path.isfile(outpath):
                localsize = os.path.getsize(outpath)
                if(str(localsize) != value):
                    self.downloadSinglefile(key, outpath)
            elif(not os.path.isfile(outpath)):
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
        text = "\rDownloading Request to './{0}' directory."\
        "Download Progress: [{1}] {2}% {3}".format(
            outdir, "="*block + " " * (barLength-block), progress*100, status)
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

    
    def download(self, directory, dsIndex):
        
        for ds in dsIndex:
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
        NCAR server for surface forecast variables (prec, swin, lwin). '''
        
        self.date = date
        self.area = area
        self.elevation = elevation
        self.variables = variables
        #self.directory = directory
        
        dpar = {'air_temperature'    : 'Temperature',
                'relative_humidity'  : 'Relative humidity',
                'eastward_wind'      : 'u-component of wind',
                'northward_wind'     : 'v-component of wind',
                'geopotential_height': 'Geopotential height',
                'pressure_level'     : 'level'}
        
        self.param = self.getParam(dpar)
        
    def getParam(self, dpar):
        
        varlist = [] 
        for var in self.variables:
            varlist.append(dpar.get(var))

        varlist = [item for item in varlist if item is not None]
        
        return varlist
    
    def makeDate(self):
        '''convert data format to NCAR RDA request'''
        
        beg = self.date['beg'].strftime('%Y%m%d%H%M')
        end = self.date['end'].strftime('%Y%m%d%H%M')
        dateRange = beg + '/to/' + end
        
        return dateRange
    
    def getPressure(self, elevation):
        g  = 9.80665   #Gravitational acceleration [m/s2]
        R  = 8.31432   #Universal gas constant for air [N.m /(mol.K)]    
        M  = 0.0289644 #Molar mass of Earth's air [kg/mol]
        P0 = 101325    #Pressure at sea level [Pa]
        T0 = 288.15    #Temperature at sea level [K]
        #http://en.wikipedia.org/wiki/Barometric_formula
        return P0 * exp((-g * M * elevation) / (R * T0)) / 100 #[hPa] or [bar]
    
    def getPressureLevels(self):
        total_elevations = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 
                            175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 
                            600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 
                            925, 950, 975, 1000]
        
        # flip max and min because 1000 is the bottom and 0 is the top
        elevationMax = self.getPressure(self.elevation['min'])
        elevationMin = self.getPressure(self.elevation['max'])
        
        minNum = min(total_elevations, key=lambda x:abs(x-elevationMin))
        maxNum = min(total_elevations, key=lambda x:abs(x-elevationMax))
        
        if (minNum > elevationMin and total_elevations.index(minNum) > 0 ):
            elevationMinRange = total_elevations.index(minNum) - 1
        else:
            elevationMinRange = total_elevations.index(minNum)
        
        if (maxNum < elevationMin and total_elevations.index(maxNum) < 36 ):
            elevationMaxRange = total_elevations.index(maxNum) - 1
        else:
            elevationMaxRange = total_elevations.index(maxNum)
            
        elevation = []
        for e in range(elevationMinRange, elevationMaxRange + 1):
            elevation.append(total_elevations[e])
        
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
        self.variables = variables
        #self.directory = directory

        dpar = {'2-metre_air_temperature'  : 'Temperature',
                '2-metre_relative_humidity': 'Relative humidity',
                '10-metre_eastward_wind'   : 'u-component of wind',
                '10-metre_northward_wind'  : 'v-component of wind'}
        
        self.param = self.getParam(dpar)
        
    def getParam(self, dpar):
        
        varlist = [] 
        for var in self.variables:
            varlist.append(dpar.get(var))

        varlist = [item for item in varlist if item is not None]
        
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
        self.variables = variables
        
        dpar = {'precipitation_amount'              :
                    'Total precipitation',
                'downwelling_shortwave_flux_in_air' : 
                    'Downward solar radiation flux',
                'downwelling_longwave_flux_in_air'  : 
                    'Downward longwave radiation flux',
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': 
                    'Clear sky downward solar radiation flux',
                'downwelling_longwave_flux_in_air_assuming_clear_sky':
                    'Clear sky downward longwave radiation flux'}
        
        self.param = self.getParam(dpar)
        
    def getParam(self, dpar):
        
        varlist = [] 
        for var in self.variables:
            varlist.append(dpar.get(var))

        varlist = [item for item in varlist if item is not None]
        
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


class JRADownload(object):
    
    #TODO add credential
    def __init__(self, pfile, username, password):
        self.dsID = 'ds628.0'
        self.pfile = pfile
        #self.credential = path.join(par.credentials_directory, ".jrarc")
        self.username = username
        self.password = password
        par = ParameterIO(self.pfile)
        
        # assign bounding box
        self.area  = {'north':  par.bbN,
                      'south':  par.bbS,
                      'west' :  par.bbW,
                      'east' :  par.bbE}
                 
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}

        # elevation
        self.elevation = {'min' : par.ele_min, 
                          'max' : par.ele_max}
        
        # data directory for JRA-55  
        self.directory = par.project_directory + '/jra55'
     
        # variables
        self.variables = par.variables
            
        # chunk size for downloading and storing data [days]        
        self.chunk_size = par.chunk_size     
        
    
    
    def retrieve(self):
        
        date_i = {}
        slices = floor(float((self.date['end'] - self.date['beg']).days)/
                       self.chunk_size)+1
        
        rda = RDA(self.username, self.password)
        
        # submit request
        dsN = 0
        for ind in range (0, int(slices)): 
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
            
            JRAli = [pl, sa, sf]
            for jrai in JRAli:
                rda.submit(jrai.getDictionary())
                dsN += 1
        
        # download dataset
        
        doneI = []
        while len(doneI) < dsN:
            print('\nGeting available dataset... Please wait as this may take'\
                  ' awhile\n')
            dsIndex = rda.getDSindex()
            dsIndex = [item for item in dsIndex if item not in doneI]
            if len(dsIndex) > 1:
                rda.download(self.directory, dsIndex)
                doneI += dsIndex
            time.sleep(10)

        #dsIndex = rda.getDSindex()
        #rda.download(self.directory, dsIndex)
        
###############################################################################
###############################################################################
# -----------------------------------------------------------------------------
# --------------------------------- TEST --------------------------------------
# -----------------------------------------------------------------------------
        
# ==== SETTING-UP =============================================================
#pfile = 'C:/OneDrive/GitHub/globsim/examples/par/examples.globsim_download'
pfile = '/Users/bincao/OneDrive/GitHub/globsim/examples/par/examples.globsim_download'
username = 'caobin198912@outlook.com'
password = 'caobin1989.12.24'

# ==== RUN ====================================================================

jra = JRADownload(pfile, username, password)
#jra.retrieve()

###############################################################################
###############################################################################

date_i = {}
slices = floor(float((jra.date['end'] - jra.date['beg']).days)/
               jra.chunk_size)+1

rda = RDA(jra.username, jra.password)

# submit request
dsN = 0
for ind in range (0, int(slices)): 
    #prepare time slices   
    date_i['beg'] = jra.date['beg'] + timedelta(days = 
                    jra.chunk_size * ind)
    date_i['end'] = jra.date['beg'] + timedelta(days = 
                    jra.chunk_size * (ind+1) - 1)
    if ind == (slices-1):
        date_i['end'] = jra.date['end']
    
    pl = JRApl(date_i, jra.area, jra.elevation, 
               jra.variables, rda)
    sa = JRAsa(date_i, jra.area, jra.variables, rda)
    sf = JRAsf(date_i, jra.area, jra.variables, rda)
    
    JRAli = [pl, sa, sf]
    for jrai in JRAli:
        rda.submit(jrai.getDictionary())
        dsN += 1

# download dataset

doneI = []
while len(doneI) < dsN:
    print(len(doneI))
    #print('\nGeting available dataset... Please wait as this may take'\
          #' awhile\n')
    dsIndex = rda.getDSindex()
    dsIndex = [item for item in dsIndex if item not in doneI]
    if len(dsIndex) > 1:
        rda.download(jra.directory, dsIndex)
        doneI += dsIndex
    time.sleep(10)