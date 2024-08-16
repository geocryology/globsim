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

from globsim.globsim_main import GlobsimScale


def main(args):
        sfile = args.f

        ERA5    = True if args.d is None or "ERA5"     in args.d else False
        ERA5ENS = True if args.d is None or "ERA5ENS"  in args.d else False
        JRA     = True if args.d is None or "JRA"      in args.d else False
        JRA3Q   = True if args.d is None or "JRA3Q"    in args.d else False
        MERRA   = True if args.d is None or "MERRA"    in args.d else False
        
        if sum([ ERA5, ERA5ENS, JRA, MERRA, JRA3Q]) > 0:
            GlobsimScale(sfile, ERA5=ERA5, ERA5ENS=ERA5ENS, 
                        JRA=JRA, MERRA=MERRA, JRA3Q=JRA3Q)

        
        else: print("Failed! Reanalysis source should be ERA5, MERRA, JRA, JRA3Q please check")

#===scale the variables from mutiple re-analysis data at stations===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scale reanalysis data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', default=None, type=str, 
                        help="file path to download parameter file")
    parser.add_argument('-d', default=None, nargs="*", type=str, 
                        help="What data sources should run?  ERA5, ERA5ENS, MERRA, JRA, JRA3Q")
    
    
    args = parser.parse_args()

    main(args)
