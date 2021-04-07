from __future__ import print_function

from os import path, makedirs

import tomlkit
import re

from globsim.common_utils import StationListRead
# handle python 3 string types
try:
    basestring
except NameError:
    basestring = str



class GenericScale:

    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        with open(self.sfile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config['scale']
        self.intpdir = path.join(par['project_directory'], 'interpolated')
        self.scdir = self.makeOutDir(par)
        self.list_name = par['station_list'].split(path.extsep)[0]

        # get the station file
        self.stations_csv = path.join(par['project_directory'],
                                      'par', par['station_list'])
        # read station points
        self.stations = StationListRead(self.stations_csv)

        # read kernels
        self.kernels = par['kernels']
        if not isinstance(self.kernels, list):
            self.kernels = [self.kernels]

    def getOutNCF(self, par, src, scaleDir='scale'):
        """make out file name"""

        timestep = str(par['time_step']) + 'h'
        src = '_'.join(['scaled', src, timestep])

        src = src + '.nc'
        fname = path.join(self.scdir, src)

        return fname

    def makeOutDir(self, par):
        """make directory to hold outputs"""

        dirSC = path.join(par['project_directory'], 'scaled')

        if not (path.isdir(dirSC)):
            makedirs(dirSC)

        return dirSC