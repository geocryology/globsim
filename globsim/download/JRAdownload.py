#!/usr/bin/env python
# -*- coding: utf-8 -*-
import glob

import netCDF4 as nc
import numpy as np
import re
import tarfile
import time
import logging



from datetime import datetime, timedelta
from math import floor
from os import path, remove
from pathlib import Path
from . import rdams_client as rc

from globsim.common_utils import variables_skip
from globsim.download.GenericDownload import GenericDownload
from globsim.meteorology import pressure_from_elevation
from globsim.download.RDA import Rdams 
from globsim.download.jra_dict_formatters import J55DictFormatter

logger = logging.getLogger('globsim.download')

def get_userinfo():
    return None, None



class JRAdownload(GenericDownload):
    '''Return an objet to download JRA55 dataset client based on RDA'''
    JRA_VERSION = 'jra55'
    API = Rdams
    DICT_FORMATTER = J55DictFormatter
    dsID = 'ds628.0'
    timeName = 'initial_time0_hours'
    levName = 'lv_ISBL1'
    

    def __init__(self, pfile):
        super().__init__(pfile)
        self.retry_delay_min = 1
        par = self.par
        self._set_input_directory(self.JRA_VERSION)
        
        self.api = self.get_accessor()

        self.__varCheck(par)

        # time bounds
        self.date  = getDate(par)

        self.credential = path.join(par['credentials_directory'], ".jrarc")
        self.account = open(self.credential, "r")
        self.inf = self.account.readlines()
        self.username = ''.join(self.inf[0].split())
        self.password = ''.join(self.inf[1].split())

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size'] * 2000

    @property
    def ncVar(self):
        return  {
                'initial_time0_hours':   'time',
                'initial_time0':         'time',
                'initial_time0_encoded': 'time',
                'lv_ISBL1':              'level',
                'g0_lat_1':              'latitude',
                'g0_lon_2':              'longitude',
                'g0_lat_2':              'latitude',
                'g0_lon_3':              'longitude',
                'GP_GDS0_SFC':           'Geopotential',
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
                'DLWRF_GDS0_SFC_ave3h':  'Downward longwave radiation flux',
                'DSWRF_GDS0_NTAT_ave3h':  'skip',}
        
    def get_accessor(self):
        cred = path.join(self.par['credentials_directory'], "rdams_token.txt")
        return self.API(cred)
    
    def __varCheck(self, par):
        '''convert one variable to a list'''

        if not isinstance(par['variables'], (list,)):
            par['variables'] = [par['variables']]

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
        elif lev == 'gp':
            dataLev = 'to'
        else:
            raise KeyError(f"File name unknown for data {dsi}")

        return dataLev

    def getVars(self, dsi):
        '''get all download variable names'''

        flist = glob.glob(path.join(self.directory, '*' + dsi + '*'))
        varlist = []
        for f in flist:
            fname = path.basename(f)
            try:
                var = fname.split('_')[2]
                varlist.append(var.split('.')[0])
            except IndexError:  # one variable (i.e. 'gp' geopotential)
                var = re.search(r"\.\d{3}_([a-zA-Z]*)\.", fname).group(1)
                varlist.append(var)

        variables = np.unique(varlist)

        return variables

    def getOutFile(self, ncf, dataLev):

        times = nc.num2date(ncf[self.timeName][:],
                            units=ncf[self.timeName].units,
                            calendar='standard')
        begStr = np.min(times).strftime('%Y%m%d')
        endStr = np.max(times).strftime('%Y%m%d')

        fileName = [self.JRA_VERSION, dataLev, begStr, 'to', endStr]
        fileName = '_'.join(fileName) + '.nc'
        fileName = path.join(self.directory, fileName)

        return fileName

    def getDimName(self, dataLev):
        '''get the dimension [time, level, latitude, longitude] name in the
        original JRA55 ncf'''

        if dataLev == 'pl':
            self.lonName = 'g0_lon_3'
            self.latName = 'g0_lat_2'
        elif dataLev in ['sa', 'sf', 'to']:
            self.lonName = 'g0_lon_2'
            self.latName = 'g0_lat_1'

    def extract_downloaded_tar_files(self, dsi):
        '''find downloaded tar files'''

        tarf = np.sort(list(Path(self.directory).glob(f"*{dsi}*.tar")))
        for f in tarf:
            if Path(f).suffix == '.tar':
                tar = tarfile.open(f)
                tar.extractall(path=self.directory, filter='fully_trusted')
                tar.close()
                remove(f)
        return tarf
    
    def get_nc_files(self, variable, dsi):
        return np.sort(list(Path(self.directory).glob(f"*{variable}*{dsi}*.nc")))
    
    def makeNCF(self, dsi):

        variables = self.getVars(dsi)
        dataLev = self.getDataLev(dsi)
        self.getDimName(dataLev)
        
        self.extract_downloaded_tar_files(dsi)
        import pdb;pdb.set_trace()
        
        nc_template_files = self.get_nc_files(variables[0], dsi)
        
        ncf = nc.MFDataset(nc_template_files.tolist(), aggdim='initial_time0_hours')

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
            
            flist = np.sort(list(Path(self.directory).glob(f"*_p125*{vari}*{dsi}*.nc")))
            ncf = nc.MFDataset(flist.tolist(), aggdim=self.timeName)
            for n, var in enumerate(ncf.variables.keys()):
                if variables_skip(self.ncfVar[var]):
                    continue
                logger.info(f"Creating variable: {var}")
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

    def requestClear(self):
        '''clear online datasets before downloading'''

        logger.info('Clear Online Datasets From NCAR Server')
        dsIndex = self.api.get_status()['data']
        if len(dsIndex) > 1:
            for req in dsIndex:
                i = req['request_index']
                self.api.purge_request(str(i))
                logger.info(f'Cleared request {i}')

    def requestSubmit(self):

        logger.info('Submit Request')

        slices = floor(float((self.date['end'] - self.date['beg']).days) /
                       self.chunk_size) + 1
        request_list = []
        for ind in range(0, int(slices)):
            date_i = {}
            # prepare time slices
            date_i['beg'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * ind))
            date_i['end'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * (ind + 1) - 1))
            if ind == (slices - 1):
                date_i['end'] = self.date['end']
            
            Formatter = self.DICT_FORMATTER
            self.formatter = Formatter(date=date_i,
                                       area=self.area,
                                       elevation=self.elevation,
                                       variables=self.variables)
            

            pl = self.formatter.get_pl_dict(self.variables) 
            sa = self.formatter.get_sa_dict(self.variables) 
            sf = self.formatter.get_sf_dict(self.variables) 
            to = self.formatter.get_to_dict()
            
            for request_dict in [to, sa, sf, pl]:
                if len(request_dict['param']) == 0:
                    continue
                
                summary = self.api.get_param_summary(self.dsID)['data']['data']
                request_dict['param'] = self.formatter.param_descriptions_to_name(request_dict['param'], summary)
                logger.debug("Requesting data for: " + str(request_dict))
                req = self.api.submit_json(request_dict)
        
                if req['http_response'] != 200:
                    logger.error(f"Error in request {request_dict}")
                else:
                    logger.info(f"Request {req['data']['request_id']} submitted")
                
                # append request ID
                try:
                    request_list.append(req['data']['request_id'])
                except KeyError:
                    import pdb;pdb.set_trace()
                    print("WHATAS THIS ERRORRR")
                    print(req)
                    print(req['data'])
                    print(req['data']['request_id'])
        
        logger.info('Request submitted')

        return set(request_list)

    def update_status(self, request_list, failed=None, done=None):
        done = done if done is not None else []
        failed = failed if failed is not None else []
        to_get = []
        logger.info('Geting Active Requets from NCAR Server')
        active_requests = self.api.get_status()['data']
        
        if request_list is None:
            request_list = [str(r['request_index']) for r in active_requests]
                            
        for rid in request_list:
            if (rid in done) or (rid in failed):
                continue
            if rid not in [str(r['request_index']) for r in active_requests]:
                logger.error(f"Request {rid} not found in active requests")
                failed.append(rid)

        for r in active_requests:
            ix = str(r['request_index'])
            if (ix in done) or (ix in failed):
                continue
            elif str(r['request_index']) in request_list:
                to_get.append(r)
        
        return to_get, failed, done

    def requestDownload(self, request_list=None):
        '''download all the requests'''
        failed=None
        done=None
        to_get = [None]

        while len(to_get) > 0:
            to_get, failed, done = self.update_status(request_list, failed=failed, done=done)

            print(f"remaining requests :{[r['request_id'] for r in to_get]}")

            for ds in to_get:
                if (ds['request_index'] in done) or(ds['request_index'] in failed):
                    continue

                if ds['status'] == 'Error':
                    logger.error(f"Error in request {ds['request_index']}")
                    done.append(ds['request_index'])
                    
                elif ds['status'] == 'Queued for Processing':
                    logger.info(f"Request {ds['request_index']} Queued for Processing")
                    continue
                
                elif ds['status'] == 'Processing':
                    logger.info(f"Request {ds['request_index']} Processing")
                    continue

                elif ds['status'] == 'Set for Purge':
                    logger.info(f"Request {ds['request_index']} Set for Purge")
                    done.append(ds['request_index'])

                elif ds['status'] == 'Completed':
                    logger.info(f"Request {ds['request_index']} Complete")
                    res = self.api.download(ds['request_index'], path.join(self.directory, self.JRA_VERSION))
                    
                    # untar / extract
                    files = [file['wfile'] for file in res['data']['web_files']]
                    self.makeNCF(str(ds['request_index']))
                    done.append(ds['request_index'])

                else:
                    logger.error(f"Unknown status {ds['status']} for request {ds['request_index']}")
                    done.append(ds['request_index'])
                    continue
            
            if len(done) + len(failed) == len(to_get):
                break
            #if self.retry_delay_min < 10:
            #    self.retry_delay_min *= 1.6
            self.retry_delay_min = (1/6)
            logger.info(f"Waiting {self.retry_delay_min} minutes before checking again")

            time.sleep(60 * self.retry_delay_min)  # check available data every 10 mins

        logger.info('''Download Completed''')

    def retrieve(self):
        '''submit and download all the dataset'''

        logger.info(f'''Starting {self.JRA_VERSION} download''')

        #self.requestClear()  # clear online request
        request_list = self.requestSubmit()  # submit request
        request_list = None
        self.requestDownload(request_list)  # download dataset

        logger.info(f'''{self.JRA_VERSION} Complete''')


def getDate(par):
    '''get download daterange'''

    dateRange = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                 'end': datetime.strptime(par['end'], '%Y/%m/%d')}
    dateRange['end'] = dateRange['end'] + timedelta(hours=23)

    return dateRange



