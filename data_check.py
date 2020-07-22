#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Copyright Xiaojing Quan & Stephan Gruber
# ==============================================================================
# A script for checking the continuity of Time with the given chunk for downloaded netCDF
# files from ERA_Interim and MERRA-2,
#
# -- Step 1: read in downloaded netCDF files with a wildcard expression;
# -- Step 2: check the time coverage of variables to be continuous;
# -- Step 3: if there was any missing one, add to list the missing series of time
#
# Referenced from data_checker.py (Dr.Stephan Gruber): Class MERRA2DataCheck()
# ===============================================================================
from __future__        import print_function
from datetime          import datetime, timedelta, date
from os                import path
#from netCDF4           import Dataset, MFDataset
from generic           import ParameterIO, StationListRead


import numpy as np
import csv
import netCDF4 as nc
import math
import time as tc
import sys


class DataCheck(object):
    """
        To check  the continuity of the time coverage from all the downloaded
        netCDF files with a wildward expression

        Args:
        ifile: Full path to interpolation parameter file.

    """

    def __init__(self, ifile, varF):
        self.ifile = ifile                      # read parameter file for interpolaiton
        par = ParameterIO(self.ifile)           # Reads generic par files and makes values available as dictionary.

        # read file, access first reanalysis variable, access data_directory, access N of bounding box
        self.dir_inp = path.join(par.project_directory, varF)
        self.dir_out = path.join(par.project_directory, 'station')
        self.variables = par.variables

        # self.list_name = par.list_name
        self.stations_csv = path.join(par.project_directory, 'par', par.station_list)

        # read station points
        self.stations = StationListRead(self.stations_csv)

        # time bounds
        self.date = {'beg': par.beg, 'end': par.end}

    def findStep(self, file_in):
        '''
        Function for determining time step to be used to ensure all time data is present.
        '''
        f = file_in.split('/')
        f = f[-1]
        stepDic = {'merra_pl': 6, 'merra_sa': 1, 'merra_sf': 1, 'erai_sf_': 3,
                   'erai_pl_': 6, 'erai_sa_': 6}
        if f[:4] == 'era5':
            step = 1
        elif f[:5] == 'jra55':
            step = 6
        else:
            step = stepDic.get(f[:8])
        return step

    def checkVar(self, ):
        '''
        Takes list of missing times
        want to check individual variables for time gaps - long periods of data drought
        '''

    def DataReadin(self, file_in):
        """
        To read in all the downloaded netCDF files with a wildward expression

        Args:
            ncfile_in: Full path to an reanalysis derived netCDF file. This can
                        contain wildcards to point to multiple files if temporal
                        chunking was used.

            ncfile_out: Full path to the output netCDF file to write.
          """
        # open netcdf file handle, can be one file of several with wildcards
        # noinspection PyBroadException
        try:
            ncf = nc.MFDataset(file_in, 'r', aggdim='time')
        except IndexError:
            print('%s does not exist' % file_in)
            return 0
        # get variable 'time'
        nctime = ncf.variables['time'][:]
        t_unit = ncf.variables['time'].units
        try:
            t_cal = ncf.variables['time'].calendar
        except AttributeError:  # Attribute doesn't exist
            t_cal = u"gregorian"  # or standard
        STEP = self.findStep(file_in)
        # STEP DEFINED BY DATA TIME FREQUENCY (6HR, 3HR, 1HR)
        try:
            missingIndex = list(set(range(nctime[0], nctime[-1] + 1, STEP)) - set(nctime))
        except TypeError:
            missingIndex = list(set(range(int(nctime[0]), int(nctime[-1]) + 1, STEP)) - set(nctime))

        if missingIndex == []:
            print("NO MISSING TIME FOR", file_in)
        else:
            missingIndex.sort()
            missingTime = nc.num2date(missingIndex, units=t_unit, calendar=t_cal)

            print('FOR', file_in)

            print("Number of Missing Indices of Time : ", len(missingIndex))

            print("Unit of Time: ", t_unit)

            print('Missing Time List:', list(missingTime))

    def process(self,varF):
        """
        combine the specific type of single or mutiple netCDF files
        """
        # === Surface  Data, Pressure-level Data, Radiation Data ===
        if varF == 'merra2':
            varF = 'merra'
        for d in ['sa', 'pl', 'sf']:
            self.DataReadin(path.join(self.dir_inp, '%s_%s_*.nc' % (varF, d)))