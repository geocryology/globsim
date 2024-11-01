#!/usr/bin/env python
"""List dataset metadata, subset data subset requests,
check on request status.

Usage:
```
rdams-client.py -get_summary <dsnnn.n>
rdams-client.py -get_metadata <dsnnn.n> <-f>
rdams-client.py -get_param_summary <dsnnn.n> <-f>
rdams-client.py -submit [control_file_name]
rdams-client.py -get_status <RequestIndex> <-proc_status>
rdams-client.py -download [RequestIndex]
rdams-client.py -globus_download [RequestIndex]
rdams-client.py -get_control_file_template <dsnnn.n>
rdams-client.py -help
```
"""
__version__ = '3.0.0'
__author__ = 'Doug Schuster (schuster@ucar.edu), Riley Conroy (rpconroy@ucar.edu)'

import sys
import os
import requests
import json
import argparse
import logging

from pathlib import Path

logger = logging.getLogger(__name__)


class Rdams(object):
    BASE_URL = 'https://rda.ucar.edu/api/'
    ENVAUTHFILE = 'GLOBSIM_RDA_AUTH_FILE'

    def __init__(self, auth_file):
        self.token = None
        if auth_file is None:
            auth_file = self.look_for_auth_file()
        
        logger.info(f"Using credential file {auth_file}")
        self.DEFAULT_AUTH_FILE = auth_file
        
    def look_for_auth_file(self):
        env_auth = os.environ.get(self.ENVAUTHFILE, None)
        if env_auth is not None:
            return env_auth

        default_auth = Path("~", 'rdams_token.txt').expanduser()
        if default_auth.exists():
            return str(default_auth)
        
        else:
            print(f"Could not find an authentication for RDAMS file. Please set the environment variable {self.ENVAUTHFILE} or create a file at {default_auth}")
            print("See 'https://rda.ucar.edu/accounts/profile/' to get a token.")
            sys.exit(1)

    def query(self, args=None):
        """Perform a query based on command line like arguments.

        Args:
            args (list): argument list of querying commands.

        Returns:
            (dict): Output of json decoded API query.

        Example:
            ```
            >>> query(['-get_status', '123456'])

            >>> query(['-get_metadata', 'ds083.2'])
            ```
        """
        parser = self.get_parser()
        if args is None or len(args) == 0:
            parser.parse_args(['-h'])
        args = parser.parse_args(args)
        args_dict = args.__dict__
        func,params = self.get_selected_function(args_dict)
        if args_dict['outdir'] and func==self.download:
            out_dir = args_dict['outdir']
            return func(params, out_dir)
        result = func(params)
        if not args.noprint:
            print(json.dumps(result, indent=3))
        return result

    def add_ds_str(self, ds_num):
        """Adds 'ds' to ds_num if needed.
        Throws error if ds number isn't valid.
        """
        ds_num = ds_num.strip()
        if ds_num[0:2] != 'ds':
            ds_num = 'ds' + ds_num
        if len(ds_num) != 7:
            print("'" + ds_num + "' is not valid.")
            sys.exit()
        return ds_num

    def get_userinfo(self):
        """Get token from command line."""
        print('Please visit https://rda.ucar.edu/accounts/profile/ to access token.')
        token = input("Paste that token here: ")
        self.write_token_file(token)
        return token

    def write_token_file(self, token, token_file):
        """Write token to a file."""
        with open(token_file, "w") as fo:
            fo.write(token)

    def read_token_file(self, token_file):
        """Read user information from token file.

        Args:
            token_file (str): location of token file.

        Returns:
            (str): token
        """
        with open(token_file, 'r') as f:
            token = f.read()
        return token.strip()

    def read_control_file(self, control_file):
        """Reads control file, and return python dict.

        Args:
            control_file (str): Location of control file to parse.
                    Or control file string.

        Returns:
            (dict) python dict representing control file.
        """
        control_params = {}
        if os.path.exists(control_file):
            myfile = open(control_file, 'r')
        else:
            myfile = control_file.split('\n')

        for line in myfile:
            line = line.strip()
            if line.startswith('#') or line == "":
                continue
            li = line.rstrip()
            (key, value) = li.split('=', 2)
            control_params[key] = value

        # Handle empty params
        if 'param' in control_params and control_params['param'].strip() == '':
            all_params = self.get_all_params(control_params['dataset'])
            control_params['param'] = '/'.join(all_params)

        try:
            myfile.close()
        except:
            pass
        return control_params

    def get_parser(self, ):
        """Creates and returns parser object.

        Returns:
            (argparse.ArgumentParser): Parser object from which to parse arguments.
        """
        description = "Queries NCAR RDA REST API."
        parser = argparse.ArgumentParser(prog='rdams', description=description)
        parser.add_argument('-noprint', '-np',
                action='store_true',
                required=False,
                help="Do not print result of queries.")
        parser.add_argument('-outdir', '-od',
                nargs='?',
                required=False,
                help="Change the output directory of downloaded files")
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-get_summary', '-gsum',
                type=str,
                metavar='<dsid>',
                required=False,
                help="Get a summary of the given dataset.")
        group.add_argument('-get_metadata', '-gm',
                type=str,
                metavar='<dsid>',
                required=False,
                help="Get metadata for a given dataset.")
        group.add_argument('-get_param_summary', '-gpm',
                type=str,
                metavar='<dsid>',
                required=False,
                help="Get only parameters for a given dataset.")
        group.add_argument('-submit', '-s',
                type=str,
                metavar='<control file>',
                required=False,
                help="Submit a request using a control file.")
        group.add_argument('-get_status', '-gs',
                type=str,
                nargs='?',
                const='ALL',
                metavar='<Request Index>',
                required=False,
                help="Get a summary of the given dataset.")
        group.add_argument('-download', '-d',
                type=str,
                required=False,
                metavar='<Request Index>',
                help="Download data given a request id.")
        group.add_argument('-get_filelist', '-gf',
                type=str,
                required=False,
                metavar='<Request Index>',
                help="Query the filelist for a completed request.")
        group.add_argument('-globus_download', '-gd',
                type=str,
                required=False,
                metavar='<Request Index>',
                help="Start a globus transfer for a give request index.")
        group.add_argument('-get_control_file_template', '-gt',
                type=str,
                metavar='<dsid>',
                required=False,
                help="Get a template control file used for subsetting.")
        group.add_argument('-purge', # Sorry no -p
                type=str,
                metavar='<Request Index>',
                required=False,
                help="Purge a request.")
        return parser

    def check_status(self, ret):
        """Checks that status of return object.

        Exits if a 401 status code.

        Args:
            ret (response.Response): Response of a request.
            token_file (str) : password file. Will remove if auth incorrect

        Returns:
            None
        """
        if ret.status_code == 401: # Not Authorized
            print(ret.content)
            exit(1)

    def check_file_status(self, filepath, filesize):
        """Prints file download status as percent of file complete.

        Args:
            filepath (str): File being downloaded.
            filesize (int): Expected total size of file in bytes.

        Returns:
            None
        """
        sys.stdout.write('\r')
        sys.stdout.flush()
        size = int(os.stat(filepath).st_size)
        percent_complete = (size/filesize)*100
        sys.stdout.write('%.3f %s' % (percent_complete, '% Completed'))
        sys.stdout.flush()

    def download_files(self, filelist, out_dir='./', retries=3, cookie_file=None):
        """Download files in a list.

        Args:
            filelist (list): List of web files to download.
            out_dir (str): directory to put downloaded files

        Returns:
            None
        """
        for _file in filelist:
            tries = 0
            while tries < retries:
                tries += 1
                try:
                    self._download_file(_file, out_dir)
                except Exception as e:
                    logger.error("Problem downloading file (attempt {tries}): {e}")
            
    def _download_file(self, _file, out_dir):
        file_base = os.path.basename(_file)
        out_file = out_dir + file_base
        print('Downloading',file_base)
        header = requests.head(_file, allow_redirects=True, stream=True)
        filesize = int(header.headers['Content-Length'])
        req = requests.get(_file, allow_redirects=True, stream=True)
        with open(out_file, 'wb') as outfile:
            chunk_size=1048576
            for chunk in req.iter_content(chunk_size=chunk_size):
                outfile.write(chunk)
                if chunk_size < filesize:
                    self.check_file_status(out_file, filesize)
        self.check_file_status(out_file, filesize)

    def encode_url(self, url, token):
        return url + '?token=' + token

    def get_authentication(self):
        """Attempts to get authentication.

        Args:
            token_file (str): location of password file.

        Returns:
            (tuple): token
        """
        token_file=self.DEFAULT_AUTH_FILE

        if self.token is not None:
            return self.token
        
        elif os.path.isfile(token_file) and os.path.getsize(token_file) > 0:
            token = self.read_token_file(token_file)
            self.token = token
        
        else:
            token = self.get_userinfo()
            self.token = token
        
        return token


    def get_summary(self, ds):
        """Returns summary of dataset.

        Args:
            ds (str): Datset id. e.g. 'ds083.2'

        Returns:
            dict: JSON decoded result of the query.
        """
        url = self.BASE_URL + 'summary/'
        url += ds

        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()

    def get_metadata(self, ds):
        """Return metadata of dataset.

        Args:
            ds (str): Datset id. e.g. 'ds083.2'

        Returns:
            dict: JSON decoded result of the query.
        """
        url = self.BASE_URL + 'metadata/'
        url += ds

        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()

    def get_all_params(self, ds):
        """Return set of parameters for a dataset.

        Args:
            ds (str): Datset id. e.g. 'ds083.2'

        Returns:
            set: All unique params in dataset.
        """
        res = self.get_param_summary(ds)
        res_data = res['data']['data']
        param_names = set()
        for param in res_data:
            param_names.add(param['param'])
        return param_names


    def get_param_summary(self, ds):
        """Return summary of parameters for a dataset.

        Args:
            ds (str): Datset id. e.g. 'ds083.2'

        Returns:
            dict: JSON decoded result of the query.
        """
        url = self.BASE_URL + 'paramsummary/'
        url += ds

        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()


    def submit_json(self, json_file):
        """Submit a RDA subset or format conversion request using json file or dict.

        Args:
            json_file (str): json control file to submit.
                    OR
                    Python dict to submit.

        Returns:
            dict: JSON decoded result of the query.
        """
        if type(json_file) is str:
            assert os.path.isfile(json_file)
            with open(json_file) as fh:
                control_dict = json.load(fh)
        else:
            assert type(json_file) is dict
            control_dict = json_file

        url = self.BASE_URL + 'submit/'

        token = self.get_authentication()
        ret = requests.post(self.encode_url(url,token), json=control_dict)

        self.check_status(ret)
        
        return ret.json()

    def submit(self, control_file_name):
        """Submit a RDA subset or format conversion request.
        Calls submit json after reading control_file

        Args:
            control_file_name (str): control file to submit.

        Returns:
            dict: JSON decoded result of the query.
        """
        _dict = self.read_control_file(control_file_name)
        return self.submit_json(_dict)


    def get_status(self, request_idx=None):
        """Get status of request.
        If request_ix not provided, get all open requests

        Args:
            request_idx (str, Optional): Request Index, typcally a 6-digit integer.

        Returns:
            dict: JSON decoded result of the query.
        """
        if request_idx is None:
            request_idx = 'ALL'
        url = self.BASE_URL + 'status/'
        url += str(request_idx)


        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()

    def get_filelist(self, request_idx):
        """Gets filelist for request

        Args:
            request_idx (str): Request Index, typically a 6-digit integer.

        Returns:
            dict: JSON decoded result of the query.
        """
        url = self.BASE_URL + 'get_req_files/'
        url += str(request_idx)

        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)

        return ret.json()


    def download(self, request_idx, out_dir='./'):
        """Download files given request Index

        Args:
            request_idx (str): Request Index, typically a 6-digit integer.

        Returns:
            None
        """
        ret = self.get_filelist(request_idx)
        if len(ret['data']) == 0:
            return ret

        filelist = ret['data']['web_files']

        token = self.get_authentication()

        web_files = list(map(lambda x: x['web_path'], filelist))

        # Only download unique files.
        self.download_files(set(web_files), out_dir)
        return ret

    def globus_download(self, request_idx):
        """Begin a globus transfer.

        Args:
            request_ix (str): Request Index, typically a 6-digit integer.

        Returns:
            dict: JSON decoded result of the query.
        """
        url = self.BASE_URL + 'request/'
        url += request_idx
        url += '-globus_download'

        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()

    def get_control_file_template(self, ds):
        """Write a control file for use in subset requests.

        Args:
            ds (str): datset id. e.g. 'ds083.2'

        Returns:
            dict: JSON decoded result of the query.
        """
        url = self.BASE_URL + 'control_file_template/'
        url += ds

        token = self.get_authentication()
        ret = requests.get(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()

    def write_control_file_template(self, ds, write_location='./'):
        """Write a control file for use in subset requests.

        Args:
            ds (str): datset id. e.g. 'ds083.2'
            write_location (str, Optional): Directory in which to write.
                    Defaults to working directory

        Returns:
            dict: JSON decoded result of the query.
        """
        _json = self.get_control_file_template(ds)
        control_str = _json['data']['template']

        template_filename = write_location + self.add_ds_str(ds) + '_control.ctl'
        if os.path.exists(template_filename):
            print(template_filename + " already exists.\nExiting")
            exit(1)
        with open(template_filename, 'w') as fh:
            fh.write(control_str)

        return _json

    def purge_request(self, request_idx):
        """Write a control file for use in subset requests.

        Args:
            ds (str): datset id. e.g. 'ds083.2'
            write_location (str, Optional): Directory in which to write.
                    Defaults to working directory

        Returns:
            None
        """
        url = self.BASE_URL + 'purge/'
        url += request_idx

        token = self.get_authentication()
        ret = requests.delete(self.encode_url(url,token))

        self.check_status(ret)
        return ret.json()

    def get_selected_function(self, args_dict):
        """Returns correct function based on options.
        Args:
            options (dict) : Command with options.

        Returns:
            (function): function that the options specified
        """
        # Maps an argument to function call
        action_map = {
                'get_summary' : self.get_summary,
                'get_metadata' : self.get_metadata,
                'get_param_summary' : self.get_param_summary,
                'submit' : self.submit,
                'get_status' : self.get_status,
                'download' : self.download,
                'get_filelist' : self.get_filelist,
                'globus_download' : self.globus_download,
                'get_control_file_template' : self.write_control_file_template,
                'purge' : self.purge_request
                }
        for opt,value in args_dict.items():
            if opt in action_map and value is not None:
                return (action_map[opt], value)

class RdamsWrapper(Rdams):
    
    def __init__(self, auth_file):
        super().__init__(auth_file)

    # wrappers for the functions
    def purge(self, ds):
        self.purge_request(ds)

    def download(self, directory, ds):
        super().download(ds, directory)
