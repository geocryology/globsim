##!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic functions will be placed here for testing new parameterizations.
#
# Thrown together by Hannah Macdonell (June, 2020)
#
# ===============================================================================
import argparse
import os
import subprocess
from globsim_main import GlobsimDownload

# =====download the wanted varialbe from multiple re-analysis date sources======
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get MET data", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', default=None, type=str, help="file path to download parameter file")
    parser.add_argument('-d', default=None, type=str, nargs="*",
                        help="What data sources should run? ERAI, ERA5, MERRA, JRA")
    parser.add_argument('-v', default=None, type=str, nargs="*", help="Which variable are we testing?")
    # TODO: add argument specifying method (REDCAPP vs. TopoScale)

    args = parser.parse_args()
    pfile = args.f

    ERAI = True if args.d is None or "ERAI" in args.d else False
    ERA5 = True if args.d is None or "ERA5" in args.d else False
    JRA = True if args.d is None or "JRA" in args.d else False
    MERRA = True if args.d is None or "MERRA" in args.d else False

    temp = True if args.v is None or "temp" in args.v else False  # TODO: add lines as more variable scaling functions are implemented


def runDownload(parFile):
    '''
    Function downloads data twice; once with and once without implementing scaling functions.

    :return: rn - returns parfile updated, but also runs download and will put in place two downloaded files


    subprocess.call("python globsim_download.py -d MERRRA, ERA5, JRA -f " + pfile, shell=True)
    parFile = parFile + '/fn'
    subprocess.call("python globsim_download.py -d MERRRA, ERA5, JRA -f ./examples/JUNE11_BC/par/JUNE11_BC.globsim_download", shell=True)
    return(parFile)
    '''
    print("okokokok")
    return (parFile)


runDownload(pfile)
'''
Download reanalysis data for area specified in par twice.
runrunrunrun
   python globsim_download.py -d MERRA -f ./examples/JUNE11_BC/par/JUNE11_BC.globsim_download
   Multiple: python3 opt/globsim/globsim_download.py -d ERA5, MERRA -f ./examples/JUNE11_BC/par/JUNE11_BC.globsim_download

    python test.py -d MERRA -f ./example/file

   If it could: 
   Down load data you requested w/o testing functions
   -> WITH testing functions
   -> Then compare the two to observations. 

   i.e. test.temp.topoScale()

   --> create class for each variable? within that, each function 
   Seperate classes for each variable 
   - bc they will be implemented and tested systematically: one at a time 
   '''
