#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import logging

from datetime import datetime, timedelta
from math import floor
from os import path, remove
from typing import Optional
from pathlib import Path
from requests import exceptions as req_exceptions

from globsim.download.GenericDownload import GenericDownload
from globsim.download.RDA import Rdams 
from globsim.download.jra_dict_formatters import J55DictFormatter
from globsim.download.JraDownloadHandler import J55DownloadHandler
logger = logging.getLogger('globsim.download')

def get_userinfo():
    return None, None


class JRAdownload(GenericDownload):
    '''Return an objet to download JRA55 dataset client based on RDA'''
    JRA_VERSION = 'jra55'
    API = Rdams
    DICT_FORMATTER = J55DictFormatter
    FILE_HANDLER = J55DownloadHandler
    dsID = 'ds628.0'
    timeName = 'initial_time0_hours'
    levName = 'lv_ISBL1'
    
    def __init__(self, pfile):
        super().__init__(pfile)
        self.retry_delay_min = 0.2
        self.retry_delay_current = 0.2
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
        self.chunk_size = par['chunk_size'] * 365
        
    def get_accessor(self):
        cred = path.join(self.par['credentials_directory'], "rdams_token.txt")
        return self.API(cred)
    
    def __varCheck(self, par):
        '''convert one variable to a list'''

        if not isinstance(par['variables'], (list,)):
            par['variables'] = [par['variables']]

    def requestClear(self,
                     only_with_status:Optional[list[str]]=None,
                     skip_with_status:Optional[list[str]]=None):
        '''clear online datasets 
        Parameters
        ----------
        only_with_status : list
            Only clear requests with status matching item of list

        Known status codes as of 2024-10:
        'Completed', 'Queued for Processing'
        '''
        logger.info('Clear Online Datasets From NCAR Server')
        dsIndex = self.api.get_status()['data']
        if len(dsIndex) > 1:
            for req in dsIndex:
                status = req['status']
                i = req['request_index']
                skip = False

                if skip_with_status and (status in skip_with_status):
                    skip = True
                if only_with_status and (status not in only_with_status):
                    skip = True

                if skip:
                    logger.debug(f"Skip purge_request for request {i} (Status = {status})")
                    continue
                else:
                    self.api.purge_request(str(i))
                    logger.info(f'Cleared request {i}')

    def _prepare_requests(self) -> list:
        slices = floor(float((self.date['end'] - self.date['beg']).days) /
                       self.chunk_size) + 1
        
        request_chunks = []

        for chunk_index in range(0, int(slices)):
            date_i = {}
            # prepare time slices
            date_i['beg'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * chunk_index))
            date_i['end'] = (self.date['beg']
                             + timedelta(days=self.chunk_size * (chunk_index + 1) - 1))
            if chunk_index == (slices - 1):
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

            if len(request_chunks) == 0:
                r = zip([to, sa, sf, pl], ['to','sa','sf','pl'])
            else:
                r = zip([sa, sf, pl], ['sa','sf','pl'])
            
            # check if file already exists
            request_chunks.append([i for i in r if not self._chunk_file_exists(*i)])
            
        return request_chunks

    def _chunk_file_exists(self, request_dict, filetype) -> bool:
        output_file = Path(self.directory,
                           f"{self.JRA_VERSION}_{filetype}_{request_dict['date'][:8]}_to_{request_dict['date'][-12:-4]}.nc")
        if output_file.exists():
            logger.info(f"File {output_file.name} already exists. Skipping")
            return True
        else:
            return False

    def submit_chunk(self, chunk, retries=3) -> dict:
        """submit request to the server
        
        Returns
        -------
        dict
            dictionary of submitted requests and their file type (e.g. 'to', 'sa', 'sf', 'pl')
            for example, {'1234':'to', '1235':'sa', '1236':'sf', '1237':'pl'}
        """
        logger.info('Submitting Requests')

        submitted_requests = dict()

        for (request_dict, filetype) in chunk:
            if len(request_dict['param']) == 0:
                continue

            # check if file already exists
            output_file = Path(self.directory,
                                f"{self.JRA_VERSION}_{filetype}_{request_dict['date'][:8]}_to_{request_dict['date'][-12:-4]}.nc")
            if output_file.exists():
                logger.info(f"File {output_file.name} already exists. Skipping")
                continue

            summary = self.api.get_param_summary(self.dsID)['data']['data']
            request_dict['param'] = self.formatter.param_descriptions_to_name(request_dict['param'], summary)
            logger.debug("Requesting data for: " + str(request_dict))
            
            successful = False
            ntries = 0
            while not successful and ntries < retries:
                ntries += 1
                try:
                    req = self.requestSubmit(request_dict)
                except req_exceptions.JSONDecodeError as e:
                    logger.error(f"JSONDecodeError in request {request_dict} (tries = {ntries})")
                    logger.error(e)
                    continue
                
                if req['http_response'] != 200:
                    logger.error(f"Error in request {request_dict} (tries = {ntries})\n(http_response={req['http_response']})")
                    for e in req.get('error_messages', []):
                        logger.error(e)
                    time.sleep(1)
                    continue

                elif not 'request_id' in req['data'].keys():
                    logger.error(f"Missing request_id in response for {request_dict} (tries = {ntries})\n(data={req['data']})")
                    time.sleep(1)
                    continue
                else:
                    successful=True
                    logger.info(f"Request {req['data']['request_id']} submitted successfully")
            
            # append request ID
            rix = str(req['data']['request_id'])
            try:
                submitted_requests[rix] = filetype
            except KeyError:
                logger.error(f"Error in request {request_dict}")
                logger.error(f"Error messages: {req.get('error_messages')}")
                continue

        return submitted_requests
    
    def requestSubmit(self, request_dict:dict) -> dict:
        """submit request to the server
        
        Returns
        -------
        dict
            dictionary of submitted requests and their file type (e.g. 'to', 'sa', 'sf', 'pl')
            for example, {'1234':'to', '1235':'sa', '1236':'sf', '1237':'pl'}
        """
        return self._submit_request(request_dict)  # append _test for server-free testing

    def _submit_request(self, request_dict:dict) -> dict:
        """ Actual server request """
        req = self.api.submit_json(request_dict)
        return req
    
    def _submit_request_test(self, request_dict:dict) -> dict:
        """ Dummy request for testing """
        req = {'http_response':200, 'data':{'request_id': str(time.monotonic())}}  
        print("temporary dummy request", request_dict)
        return req
    
    def update_status(self, submitted_requests:Optional[dict], failed=None, done=None):
        ''' 
        
        '''
        done = done if done is not None else []
        failed = failed if failed is not None else []
        to_get = []
        logger.info('Geting Active Requets from NCAR Server')
        active_requests = self.api.get_status()['data']  # type: dict
        
        if submitted_requests is None:
            submitted_requests = {str(r['request_index']):None for r in active_requests}

        for rix in submitted_requests.keys():  # Omit requests that 
            if (str(rix) in done) or (str(rix) in failed):
                continue
            if str(rix) not in [str(r['request_index']) for r in active_requests]:
                import pdb;pdb.set_trace()
                logger.error(f"Request {rix} not found in active requests")
                failed.append(str(rix))

        for r in active_requests:  # Omit requests that are already done or failed
            rix = r['request_index']
            if (str(rix) in done) or (str(rix) in failed):
                continue
            elif str(rix) in submitted_requests.keys():
                to_get.append(r)
        
        return to_get, failed, done

    def _request_download(self, submitted_requests:Optional[dict]=None):
        '''download all the requests
        request_id = 'USER14123'
        request_index = 14123
        '''
        failed=None
        done=None
        to_get = [None]

        while len(to_get) > 0:
            to_get, failed, done = self.update_status(submitted_requests, failed=failed, done=done)

            remaining = {r['request_id']: r['status'] for r in to_get}
            logger.info(f"remaining requests :{remaining}")
            logger.debug(f"done: {done} ")
            logger.debug(f"failed: {failed} ")
            
            for ds in to_get:
                rix = str(ds['request_index'])
                if (rix in done) or(rix in failed):
                    continue

                if ds['status'] == 'Error':
                    logger.error(f"Error in request {rix}")
                    failed.append(rix)
                    continue
                    
                elif ds['status'] == 'Queued for Processing':
                    logger.info(f"Request {rix} Queued for Processing")
                    continue
                
                elif ds['status'] == 'Processing':
                    logger.info(f"Request {rix} Processing")
                    continue

                elif ds['status'] == 'Set for Purge':
                    logger.info(f"Request {rix} Set for Purge")
                    done.append(rix)
                    continue

                elif ds['status'] == 'Completed':
                    logger.info(f"Request {rix} Complete")
                    res = self.api.download(ds['request_index'], path.join(self.directory, self.JRA_VERSION))
                    
                    # untar / extract
                    Handler = self.FILE_HANDLER
                    handler = Handler()
                    if submitted_requests is not None:
                        filetype = submitted_requests[rix]
                    else:
                        filetype = None
                    handler.make_globsim_dataset(self.directory, rix, filetype=filetype)
                    done.append(rix)
                    self.api.purge_request(rix)

                else:
                    logger.error(f"Unknown status {ds['status']} for request {ds['request_index']}")
                    failed.append(rix)
                    self.api.purge_request(str(ds['request_index']))
                    continue
                
            if (len(to_get) == 1):
                last_ri = str(to_get[0]['request_index'])
                if (last_ri in done) or (last_ri in failed):
                    break
                
            if self.retry_delay_current < 10:
                self.retry_delay_current *= 1.2
            
            logger.debug(f"Waiting {self.retry_delay_current:0.2f} minutes before next check")
            time.sleep(60 * self.retry_delay_current)  # check available data every 10 mins

        logger.info("All requests processed")
        self.retry_delay_current = self.retry_delay_min * 1

        if done is not None:
            for rix in done:
                logger.info(f"{rix}: Done")
        if failed is not None:
            for rix in failed:
                logger.error(f"{rix}: Failed ()")
            
        logger.info('''Download Completed''')

    def _request_download_test(self, submitted_requests:Optional[dict]=None):
        """ Dummy request function for testing """
        if not submitted_requests:
            return
        print("temporary request download")
        for r in submitted_requests.keys():
            print(f"downloaded {r}")
            time.sleep(1)
        print("downloaded all requests")

    def requestDownload(self, submitted_requests:Optional[dict]=None):
        """ """
        self._request_download(submitted_requests) # append _test for testing

    def ensure_requests_fewer_than(self, N=7, timeout=600):
        """ Make sure there are fewer than N requests on the server """
        active_requests = self.api.get_status()['data']

        # first, just try clearing those that haven't finished yet
        if len(active_requests) > N:
            logger.error(f"Found {len(active_requests)} active requests (>7). Will clear them before proceeding")
            self.requestClear(only_with_status=['Completed', 'Error'])  
        
        time.sleep(10)  # allow time for purge to register
        
        # if that's not enough, purge them all (including those that are queued for processing)
        active_requests = [r for r in self.api.get_status()['data'] if r['status'] not in ['Purged', 'Set for Purge']]
        if len(active_requests) > N:
            self.requestClear()  
        
        self.wait_for_purge(timeout=timeout)
        
    def wait_for_purge(self, timeout=3600):
        ''' wait until there are no requests on the server that await purging'''
        n_to_purge = 1
        s = 60
        while (n_to_purge > 0):
            logger.info(f"{n_to_purge} requests remain to purge on server. Waiting {int(s/60)} minutes")
            if timeout < 0:
                raise RuntimeError("Waiting for purge on server timed out. Try again later. ")

            active_requests = self.api.get_status()['data']
            n_to_purge = len([r for r in active_requests if r['status'] in ['Purged', 'Set for Purge']])
            s *= 1.2

            time.sleep(s)
            timeout -= s

    def retrieve(self):
        '''submit and download all the dataset'''

        logger.info(f'''Starting {self.JRA_VERSION} download''')
        chunked_requests = self._prepare_requests()
        
        concurrent_chunks = 1

        while chunked_requests:
            submitted_requests = {}
            self.ensure_requests_fewer_than(7)

            while chunked_requests:
                for i in range(concurrent_chunks):
                    chunk = chunked_requests.pop(0)
                    submitted_requests.update(self.submit_chunk(chunk))
            
                self.requestDownload(submitted_requests)

        logger.info(f'''{self.JRA_VERSION} Complete''')


def getDate(par):
    '''get download daterange'''

    dateRange = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                 'end': datetime.strptime(par['end'], '%Y/%m/%d')}
    dateRange['end'] = dateRange['end'] + timedelta(hours=23)

    return dateRange



