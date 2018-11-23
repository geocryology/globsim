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
import argparse
import configparser
from globsim import GlobsimScale


#===scale the variables from mutiple re-analysis data at stations===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scale reanalysis data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--f',    default=None, type=str, help="file path to download parameter file")
    parser.add_argument('--d',    default=None, nargs="*", type=str, help="What data sources should run? ERAI, ERA5, MERRA, JRA")
    
    
    args = parser.parse_args()

    sfile = args.f
    
    ERAI =  True if args.d is None or "ERAI"  in args.d else False
    ERA5 =  True if args.d is None or "ERA5"  in args.d else False
    JRA =   True if args.d is None or "JRA"   in args.d else False
    MERRA = True if args.d is None or "MERRA" in args.d else False
    
    GlobsimScale(sfile, ERAI=ERAI, ERA5=ERA5, JRA=JRA, MERRA=MERRA)
    





