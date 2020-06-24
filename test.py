##!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic functions will be placed here for testing new parameterizations.
#
# Thrown together by Hannah Macdonell (June, 2020)
#
# ===============================================================================
import argparse
import datetime
import os
import os.path
from os import path

from globsim.globsim_main import GlobsimScale

# =====download the wanted varialble from multiple re-analysis date sources======


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="What are we testing?", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', default=None, type=str, help="project file path i.e.globsim/examples/JUNE11_BC ")
    parser.add_argument('-d', default=None, type=str, nargs="*", help="What data sources should run? ERAI, ERA5, MERRA, JRA")
    parser.add_argument('-m', default=None, type=str, help="Which method are we testing?")
    parser.add_argument('-c', action="store_true", default=False, help="Run data checker")

    args = parser.parse_args()
    PROJ = args.f
    MET = args.m
    PATH = "/home/hma000/examples/" + PROJ +'/'
    data = ['MERRA', 'ERA', 'ERA5', 'JRA55']

    ERAI = True if args.d is None or "ERAI" in args.d else False
    ERA5 = True if args.d is None or "ERA5" in args.d else False
    JRA = True if args.d is None or "JRA" in args.d else False
    MERRA = True if args.d is None or "MERRA" in args.d else False

    #temp = True if args.v is None or "temp" in args.v else False
    # TODO: add lines as more variable scaling functions are implemented


def checkData():
    """
    Data checks to insure there is interpolated data present in given project directory before running scale.
-> checks all proj folders unless proj is specified
    specify project name
    checks folders and returns what you have
    i.e. MERRA
    """

    for i in ['merra2', 'erai', 'era5', 'jra55']:
        # Returns list of file names from data directory
        #print(str(path.exists("/home/hma000/examples/JUNE11_BC/erai")))
        #print(path.exists("/home/hma000/examples/JUNE11_BC"))
        print(PATH)
        print(path.exists(PATH))
        print(os.listdir(PATH))
        if path.exists(PATH+i):
            print('yup')
            ls = [os.listdir(PATH+i + "*.txt")]
            print(ls)
    #dt_stamp.strftime('%Y/%m/%d') # should return: '2019/01/09'

def testScale():
    '''
    First bit: writes new par file containing tested scale method.
    '''
    fin = open(PATH+PROJ+".globsim_scale", "rt")
    fout = open(PATH+MET+".globsim_scale", "wt")
    for line in fin:
        print(line)
        fout.write(line.replace("kernels =", "kernels = TOPOSCALE,"))
    fin.close()
    fout.close()

if args.c == True: checkData()
'''
    if not args.c:checkData()
    else: checkData()
    
'''
