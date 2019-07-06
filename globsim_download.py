#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic classes, methods, functions used for more than one reanalysis.
#
#
# (C) Copyright Stephan Gruber (2017)
#         
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
from globsim.globsim_main import GlobsimDownload

import argparse
import configparser
import time

#=====download the wanted varialbe from multiple re-analysis date sources======
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get MET data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f',             default=None, type=str, help="file path to download parameter file")
    parser.add_argument('-d',             default=None, type=str, nargs="*", help="What data sources should run? ERAI, ERA5, MERRA, JRA")
    parser.add_argument('-r', '--retry',  default=1,    type=int, help="Number of times to re-launch download if it crashes.")
    parser.add_argument('-m', '--multi', action='store_true',    help="Download all data sources simultaneously")
    
    args = parser.parse_args()
    
    pfile = args.f
    
    ERAI =  True if args.d is None or "ERAI"  in args.d else False
    ERA5 =  True if args.d is None or "ERA5"  in args.d else False
    JRA =   True if args.d is None or "JRA"   in args.d else False
    MERRA = True if args.d is None or "MERRA" in args.d else False
    
    r_max = args.retry
    i = 0
    
    if r_max <= 1:
        GlobsimDownload(pfile, ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA, multithread=args.multi)
    else:
        while i < r_max:
            try:
                GlobsimDownload(pfile, ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA, multithread=args.multi)
            except Exception as e:
                print(e)
            time.sleep(360)
            i += 1