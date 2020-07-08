##!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic functions will be placed here for testing new parameterizations.
#
# Thrown together by Hannah Macdonell (June, 2020)
#   Lines that work:
#       python test.py -f JUNE11_BC -c
#       python test.py -f JUNE11_BC -m TOPO
#   -f PROJECT NAME - used for newSfile, testData and checkDownloads
#   -c EMPTY - calls checkDownloads() function. You can specify project with -f
#   -m METHOD - used for testScale, specifying which kernal we're testing
#   -d DATA - which reanalysis data are we using?
# ===============================================================================



import argparse
import netCDF4
import numpy as np
import sys
import os
import os.path
#from globsim.globsim_main import GlobsimScale
from globsim.data_check import DataCheck
import unittest

# =====download the wanted variable from multiple re-analysis date sources======

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

def netInfo(ncFile):
    '''
    f = netCDF object that holds metadata.
        #print(temp.dimensions)     ('time', 'level', 'lats', 'lons')
        #print(temp.shape)          (8, 13, 9, 7)
        #print(f.variables.keys())  (['RH', 'H', 'T', 'U', 'V', 'time', 'level', 'latitude', 'longitude'])
    '''
    f = netCDF4.Dataset(ncFile)
    temp = f.variables['T']
    temp = temp[:]
    depth = f.variables['level']
    dpth = depth[:]

    x = f.variables['latitude']
    y = f.variables['longitude']

    xx, yy = x[:], y[:]
    print('shape of temp variable: %s' % repr(temp.shape))
    tempslice = temp[0, dpth > 400, yy > yy.max() / 2, xx > xx.max() / 2]
    print('shape of temp slice: %s' % repr(tempslice.shape))

def checkData():
    """
    A file that calls data_check.py and checks for all reanalysis data for missing time.
    """
    ifile = 'examples/Example1/par/examples.globsim_interpolate'
    for i in ['merra2', 'jra55', 'erai', 'era5']:
        print('------------Checking %s data files------------' % i)
        varF = i
        DC = DataCheck(ifile, varF)
        DC.process(varF)

def newSfile():
    '''
    If -m argument was given, testScale is called to write a .globsim_scale function with that
    new kernal scaling function inserted.
    '''
    fin = open(("%s/par/%s.globsim_scale" % (PATH, PROJ)), "rt")
    fin = open(("%s/par/%s-%s.globsim_scale" % (PATH, PROJ,MET)), "wt")
    for line in fin:
        fout.write(line.replace("kernels =", "kernels = "+MET+', '))
    fin.close()
    fout.close()
    return [fin.name, fout.name]

def testScale():
    '''
    Runs scaling with and without tested scaling function.
    '''
    ls = newSfile() #newSfile returns list of two par file paths
    GlobsimScale(ls[0], ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA)
    GlobsimScale(ls[1], ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA)


'''
other test.py stuff:
print(os.path.exists('/home/hma000/globsim/examples/Example1/merra2/merra_pl_20170701_to_20170702.nc'))
print(globsimScaled2Pandas('/home/hma000/globsim/examples/Example1/merra2/merra_pl_20170701_to_20170702.nc', 1))
if args.m: testScale()
if args.c: checkDownloads()
netInfo(args.f)
if file_in == path.join(self.dir_inp, '%s_pl_*.nc' % varD):


print(os.path.join(PATH,'par/%s.globsim_download' % PROJ))

varF = 'merra2'
varD = 'MERRA'
for dataType in ['sa', 'pl', 'sf']:
    print(os.path.join(os.path.abspath('globsim'),'%s_%s_*.nc' % (varF, dataType)))

for i in [['merra2','MERRA'],['jra55','JRA'],['erai','ERAI'],['era5','ERA5']]:
    varF, varD = i[0],i[1]
    print(varF+' - '+varD)
'''

