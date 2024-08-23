import urllib.error
import urllib.request
from os import remove, makedirs, path
import sys
import tarfile
import http.cookiejar
import json
import urllib.error
import urllib.request
import logging
import numpy as np

logger = logging.getLogger(__name__)
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
                logger.error('RDA username and password invalid. Please try again\n')
                opener = self.makeOpener(theurl)
                try:
                    url = opener.open(request)
                except urllib.error.HTTPError as e:
                    if e.code == 401:
                        logger.error(
                            'RDA username and password invalid, or you are'
                            'not authorized to access this dataset. Please verify your login information at http://rda.ucar.edu')
                        sys.exit()
            else:
                logger.error(f"Received response: '{e}' when opening url {theurl}")
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

        logger.debug(url.read().decode())

    def getSummary(self, dsID=None):
        '''get the summary of given dataset'''

        logger.info('Getting summary information. Please wait as this may take a while.')

        theurl = theurl = self.base + 'summary/' + dsID
        url = self.urlOpen(theurl)

        logger.info(url.read().decode())

    def getMeta(self, dsID):
        '''get the metadata of gieven dataset'''

        theurl = self.base + 'metadata/' + dsID
        url = self.urlOpen(theurl)

        logger.info(url.read().decode())

    def getParaSummary(self, dsID):
        '''submit a subset request control file.
           Subset request control files are built from the parameters dumped
           out by the '-get_metadata <dsnnn.n>' option.'''

        logger.info('Getting parameter summary. Please wait as this may take a while.')

        theurl = self.base + 'paramsummary/' + dsID
        url = self.urlOpen(theurl)

        logger.info(url.read().decode())

    def getStatus(self):
        '''get the status of requests'''

        theurl = self.base + 'request'
        url = self.urlOpen(theurl)

        logger.info(url.read().decode())

    def submit(self, controlparms):
        '''submit JRA55 dataset request based on predescribed control
        parameters'''
        theurl = self.base + 'request'

        jsondata = '{'
        for k in list(controlparms.keys()):
            jsondata += '"' + k + '"' + ":" + '"' + controlparms[k] + '",'
        jsondata = jsondata[:-1]
        jsondata += '}'
        logger.info('Submitting request')

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
                logger.error('RDA username and password invalid.  Please try again\n')
                (username, password) = get_userinfo()
                opener = self.makeOpener(theurl)
                try:
                    url = opener.open(request)
                except urllib.error.HTTPError as e:
                    if e.code == 401:
                        logger.error('RDA username and password invalid, or you are not authorized to access this dataset.')
                        logger.error('Please verify your login information at http://rda.ucar.edu.')
                        sys.exit()
            else:
                logger.error(f"Received response: '{e}' when opening url {theurl}")
                sys.exit()
        logger.info(url.read().decode())

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

            logger.info("Starting Download.")

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
                    logger.error('RDA username and password invalid. Please try again')
                    opener = self.makeOpener(theurl)
                    try:
                        url = opener.open(request)
                    except urllib.error.HTTPError as e:
                        if e.code == 401:
                            logger.error('RDA username and password invalid, or you are not authorized to access this dataset.')
                            logger.error('Please verify your login information at http://rda.ucar.edu.')
                            sys.exit()

            logger.info(url.read().decode())