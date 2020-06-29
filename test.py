##!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic functions will be placed here for testing new parameterizations.
#
# Thrown together by Hannah Macdonell (June, 2020)
#   Lines that work:
#       python test.py -f JUNE11_BC -c
#       python test.py -f JUNE11_BC -m TOPO
#   -f PROJECT NAME - used for newSfile, testData and checkData
#   -c EMPTY - calls checkData() function. You can specify project with -f
#   -m METHOD - used for testScale, specifying which kernal we're testing
#   -d DATA - which reanalysis data are we using?
# ===============================================================================
import argparse
import netCDF4
import numpy as np
import os
import os.path

from globsim.globsim_main import GlobsimScale

# =====download the wanted variable from multiple re-analysis date sources======

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="What are we testing?", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', default=None, type=str, help="Project name.")
    parser.add_argument('-d', default=None, type=str, nargs="*", help="ERAI, ERA5, MERRA, or JRA")
    parser.add_argument('-m', default=None, type=str, help="Run testScale().")
    parser.add_argument('-c', action="store_true", default=False, help="Run checkData().")

    args = parser.parse_args()
    PROJ = args.f
    PATH = os.path.dirname(os.path.abspath('globsim'))
    if PROJ:
        PATH = PATH+'/'+PROJ
    MET = args.m
    DATA = args.d[0]

def netInfo():
    '''
    f = netCDF object that holds metadata.
    :return:
    '''
    f = netCDF4.Dataset()
    print(f.variables.keys())
    print(f)

def checkData():
    """
    Eventually:
        - Check all data in all downloaded files
        - If files match your query, copy them over into cur project instead of downloading new
        - Run download for what's missing
    Right now:
        - Searches specified project and creates list of all data files
    TODO: - Possible fun time transition for easy reading:
        dt_stamp.strftime('%Y/%m/%d') # should return: '2019/01/09' using import datetime
    """
    ls =[]
    for root, dirs, files in os.walk(PATH):
        for file in files:
            if file.endswith(".nc"):
                #TODO - run netInfo(file) to retrieve data stephan ideally wants
                f = file.strip('.nc').split('_')
                try:
                    ls.append([f[0], f[1], f[2], f[4]])
                except:
                    IndexError
    for i in ls: print(i)

def newSfile():
    '''
    If -m argument was given, testScale is called to write a .globsim_scale function with that
    new kernal scaling function inserted.
    '''
    print(PATH)
    fin = open(PATH+'/par/'+PROJ+".globsim_scale", "rt")
    fout = open(PATH+'/par/'+PROJ+'-'+MET+".globsim_scale", "wt")
    for line in fin:
        fout.write(line.replace("kernels =", "kernels = "+MET+', '))
    fin.close()
    fout.close()


def testScale():
    '''
    Runs selected scale with and without tested scaling function.
    :return:
    '''
    newSfile()
    path = PATH+'/par/'+PROJ
    os.system('python globsim_scale.py -d '+DATA+' -f '+path+'.globsim_scale')
    os.system('python globsim_scale.py -d '+DATA+' -f '+path+'-'+MET+'.globsim_scale')

if args.f and args.c: checkData()
if args.c and not args.f: checkData(False)
if args.m: testScale()

