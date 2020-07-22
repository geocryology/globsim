##!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic functions will be placed here for testing new parameterizations.
#
# Thrown together by Hannah Macdonell (June, 2020)
#
# Functions:
#   netInfo(ncFile)
#   checkData()
#   newSfile(PATH,PROJ)     -f file path to project directory
# ===============================================================================

import argparse
import netCDF4
from netCDF4 import Dataset, MFDataset
import h5py
import numpy as np
import os.path
from globsim.data_check import DataCheck, GlobsimScale
import unittest
from osgeo import gdal, osr, ogr
import matplotlib.pyplot as plt
import shapefile

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="What are we testing?", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', default=None, type=str, help="Project name.")
    parser.add_argument('-d', default=None, type=str, nargs="*", help="ERAI, ERA5, MERRA, or JRA")
    parser.add_argument('-m', default=None, type=str, help="Run testScale().")
    parser.add_argument('-c', action="store_true", default=False, help="Run checkDownloads().")

    args = parser.parse_args()
    PROJ = args.f
    PATH = os.path.dirname(os.path.abspath('globsim'))
    if PROJ: PATH = PATH+'/'+PROJ
    MET = args.m
    if args.d: DATA = args.d[0]

class testingFunctions():
    '''
    self.
        - file path
        - project name

    '''

    def __init__(self):
        # read parameter file
        self.pfile = pfile
        par = ParameterIO(self.pfile)

        # assign bounding box
        self.area = {'north': par.bbN, 'south': par.bbS, 'west': par.bbW, 'east': par.bbE}

        # sanity check to make sure area is good
        if (par.bbN < par.bbS) or (par.bbE < par.bbW):
            raise Exception("Bounding box is invalid: {}".format(self.area))
        parser = argparse.ArgumentParser(description="What are we testing?",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-f', default=None, type=str, help="File path to project directory.")

        args = parser.parse_args()

        self.prjpath = args.f #File path to project
        self.prj = self.prjpath.split('/')[-1] #Name of project

    def checkData(self):
        """
        A file that calls data_check.py and checks for all reanalysis data for missing time.
        """
        ifile = '%s/%s.globsim_interpolate' % self.prjpath, self.prj
        for i in ['merra2', 'jra55', 'erai', 'era5']:
            print('------------Checking %s data files------------' % i)
            varF = i
            dCheck = DataCheck(ifile, varF)
            dCheck.process(varF)

    def newSfile(self):
        '''
        If -m argument was given, testScale is called to write a .globsim_scale function with that
        new kernal scaling function inserted.
        '''
        MET = raw_input("Enter new met kernal name : ")
        #TODO: check dictionary for met variable
        fin = open(("%s/par/%s.globsim_scale" % (self.prjpath, self.prj)), "rt")
        fout = open(("%s/par/%s-%s.globsim_scale" % (self.prjpath, self.prj, MET)), "wt")
        for line in fin:
            fout.write(line.replace("kernels =", "kernels = "+MET+', '))
        fin.close()
        fout.close()
        return [fin.name, fout.name]

    def testScale(self):
        '''
        Runs scaling with and without tested scaling function.
        '''
        ls = newSfile(self.prjpath, self.prj)  # newSfile returns list of two par file paths
        GlobsimScale(ls[0], ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA)
        GlobsimScale(ls[1], ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA)

# function to create the mask of your shapefile
    def makeMask(self, lon, lat, res):
        source_ds = ogr.Open(shapefile)
        source_layer = source_ds.GetLayer()

        # Create high res raster in memory
        mem_ds = gdal.GetDriverByName('MEM').Create('', lon.size, lat.size, gdal.GDT_Byte)
        mem_ds.SetGeoTransform((lon.min(), res, 0, lat.max(), 0, -res))
        band = mem_ds.GetRasterBand(1)

        # Rasterize shapefile to grid
        gdal.RasterizeLayer(mem_ds, [1], source_layer, burn_values=[1])

        # Get rasterized shapefile as numpy array
        array = band.ReadAsArray()

        # Flush memory file
        mem_ds = None
        band = None
        return array

    def runMask(self):
        # set the data directories
        datadir = self.prjpath
        shapefile = raw_input('Input path to shp file: ')
        ncs = raw_input('Input path to nc file: ')
        ifile = ncs.split('/')[-1]
        # read the netcdf data file
        nc = netCDF4.Dataset(ncs, 'r')

        # get the precipitation
        pr = nc.variables['pr'][:]

        # show the precipitation
        plt.imshow(pr)
        plt.show()
        # get the longitude information
        lons = nc.variables['lon'][:]
        # get the latitude information
        lats = nc.variables['lat'][:]
        # calculate the cellsize
        cellsize = lons[:][1] - lons[:][0]

        # create the mask
        mask = makeMask(lons, lats, cellsize)

        # show the mask
        plt.imshow(mask)
        plt.show()

        # mask the precipitation data
        precip = np.ma.masked_where(mask == 0, pr)

        plt.imshow(precip)
        plt.show()

        # print some stats
        print('she functions')
        # print np.min(precip), np.mean(precip), np.max(precip)
