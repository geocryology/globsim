import tomlkit
import numpy as np

from os import path, makedirs


class GenericDownload(object):
    """
    Generic functionality for download classes
    """

    def __init__(self, pfile):
        # read parameter file
        self.pfile = pfile
        with open(self.pfile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config['download']

        self._set_elevation(par)
        self._set_area(par)
        self._check_area(par)

        self.variables = par['variables']

    def _check_area(self, par):
        if (par['bbN'] < par['bbS']) or (par['bbE'] < par['bbW']):
            raise ValueError("Bounding box is invalid: {}".format(self.area))

        if (np.abs(par['bbN'] - par['bbS']) < 1.5) or (np.abs(par['bbE'] - par['bbW']) < 1.5):
            raise ValueError("Download area is too small to conduct interpolation.")

    def _set_area(self, par):
        self.area = {'north': par['bbN'],
                     'south': par['bbS'],
                     'west': par['bbW'],
                     'east': par['bbE']}

    def _set_elevation(self, par):
        self.elevation = {'min': par['ele_min'],
                          'max': par['ele_max']}

    def _set_data_directory(self, name):
        self.directory = path.join(self.par['project_directory'], name)
        if not path.isdir(self.directory):
            makedirs(self.directory)
