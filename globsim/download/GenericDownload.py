import numpy as np
import tomlkit
import logging

from os import path, makedirs
from abc import ABC, abstractmethod

from globsim.meteorology import pressure_from_elevation

logger = logging.getLogger("globsim.download")


class GenericDownload(ABC):
    """
    Generic functionality for download classes
    """
    retry_delay_min = 0.5
    
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

    @abstractmethod
    def retrieve(self):
        raise NotImplementedError("This method must be implemented by a subclass.")

    def _check_area(self, par):
        if (par['bbN'] < par['bbS']) or (par['bbE'] < par['bbW']):
            raise ValueError("Bounding box is invalid: {}".format(self.area))

        if (np.abs(par['bbN'] - par['bbS']) < 1.5) or (np.abs(par['bbE'] - par['bbW']) < 1.5):
            raise ValueError("Download area is too small to conduct interpolation.")

    def _set_area(self, par):
        logger.info(f"Bounding box: {par['bbS']} (S) to {par['bbN']} (N) and {par['bbW']} (W) to {par['bbE']} (E)")
        self.area = {'north': par['bbN'],
                     'south': par['bbS'],
                     'west': par['bbW'],
                     'east': par['bbE']}

    def _set_elevation(self, par):
        self.elevation = {'min': par['ele_min'],
                          'max': par['ele_max']}

    def _set_input_directory(self, name):
        self.directory = path.join(self.par['project_directory'], name)
        if not path.isdir(self.directory):
            makedirs(self.directory)
        
        if self.par.get("output_directory") is not None:
            self.output_directory = self.par.get("output_directory")
            if not path.isdir(self.output_directory):
                makedirs(self.output_directory)
        else:
            self.output_directory = self.directory

    @staticmethod
    def getPressureLevels(elevations: list, min, max):
        # flip max and min because 1000 is the bottom and 0 is the top
        elevationMax = pressure_from_elevation(min)
        elevationMin = pressure_from_elevation(max)

        minNum = min(elevations, key=lambda x:abs(x - elevationMin))
        maxNum = min(elevations, key=lambda x:abs(x - elevationMax))

        if (minNum > elevationMin and elevations.index(minNum) > 0):
            elevationMinRange = elevations.index(minNum) - 1
        else:
            elevationMinRange = elevations.index(minNum)

        if (maxNum < elevationMin and elevations.index(maxNum) < 36):
            elevationMaxRange = elevations.index(maxNum) - 1
        else:
            elevationMaxRange = elevations.index(maxNum)

        elevation = []
        for e in range(elevationMinRange, elevationMaxRange + 1):
            elevation.append(elevations[e])

        elevation = [str(ele) for ele in elevation]
        elevation = '/'.join(elevation)

        return elevation