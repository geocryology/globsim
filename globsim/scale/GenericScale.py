from __future__ import print_function

import tomlkit
import warnings
import logging

from os import path, makedirs
from pathlib import Path

from globsim.common_utils import StationListRead

# handle python 3 string types
try:
    basestring
except NameError:
    basestring = str

logger = logging.getLogger('globsim.scale')


class GenericScale:

    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        with open(self.sfile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config['scale']
        self.intpdir = path.join(par['project_directory'], 'interpolated')
        self.output_dir = self.make_output_directory(par)
        self.list_name  = path.basename(path.normpath(par['station_list'])).split(path.extsep)[0]

        # get the station file
        self.stations_csv = path.join(par['project_directory'],
                                      'par', par['station_list'])
        # read station points
        self.stations = StationListRead(self.stations_csv)

        # read kernels
        self.kernels = par['kernels']
        if not isinstance(self.kernels, list):
            self.kernels = [self.kernels]

    def getOutNCF(self, par, data_source_name):
        """make out file name"""

        timestep = str(par['time_step']) + 'h'
        snowCor  = 'scf' + str(par['scf'])
        src = '_'.join(['scaled', data_source_name, timestep, snowCor])

        src = src + '.nc'
        fname = path.join(self.output_dir, src)

        return fname

    def make_output_directory(self, par):
        """make directory to hold outputs"""
        output_dir = None
        
        logger.debug(f"Attempting to create directory: {output_dir}")
        
        if par.get('output_directory'):
            try:
                test_path = Path(par.get('output_directory'))
            except TypeError:
                msg = "You provided an output_directory for scaled files that does not exist. Saving files to project directory"
                logger.error(msg)
                warnings.warn(msg)
            
            if test_path.is_dir():
                output_dir = Path(test_path, "scaled")
            else:
                logger.warning("You provided an output_directory for scaled files that was not understood. Saving files to project directory.")
                
        if not output_dir:
            output_dir = path.join(par['project_directory'], 'scaled')

        if not Path(output_dir).is_dir():
            makedirs(output_dir)

        return output_dir
