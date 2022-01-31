from __future__ import print_function

import tomlkit
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
            self.par = config['scale']
        self.set_parameters(self.par)

    def set_parameters(self, par):
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
        
        # should file be overwritten - default to false
        try:
            self._overwrite_output = par['overwrite']
        except KeyError as e:
            self._overwrite_output = False
        finally:
            logger.debug("Scale configured to overwrite output files")

        # read snow correction info
        try:
            self.scf = par['scf']
        except KeyError as e:
            logger.warning("Missing snow correction factor in control file. Reverting to default (1).")
            self.scf = 1
        finally:
            logger.info(f"Snow correction factor for scaling set to {self.scf}")

    def getOutNCF(self, par, data_source_name):
        """make out file name"""

        timestep = str(par['time_step']) + 'h'
        snowCor  = 'scf' + str(self.scf)
        src = '_'.join(['scaled', data_source_name, timestep, snowCor])

        src = src + '.nc'
        output_file = Path(self.output_dir, src)

        if output_file.is_file():
            if self._overwrite_output:
                output_file.unlink()
                logger.info(f"Removed existing output file {output_file}")
            else:
                logger.error("Scaled file output already exists")
                raise FileExistsError(f"Output file {output_file} exists. Remove file or set 'overwite=true' in control file.")

        return output_file.resolve()

    def make_output_directory(self, par):
        """make directory to hold outputs"""
        output_dir = None
        
        if par.get('output_directory'):
            try:
                test_path = Path(par.get('output_directory'))
            except TypeError:
                msg = "You provided an output_directory for scaled files that does not exist. Saving files to project directory"
                logger.warning(msg)
            
            if test_path.is_dir():
                output_dir = Path(test_path, "scaled")
            else:
                logger.warning("You provided an output_directory for scaled files that was not understood. Saving files to project directory.")
                
        if not output_dir:
            logger.debug(f"Attempting to create directory: {output_dir}")
            output_dir = path.join(par['project_directory'], 'scaled')

        if not Path(output_dir).is_dir():
            makedirs(output_dir)

        return output_dir
