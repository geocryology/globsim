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

from globsim.globsim_main import GlobsimInterpolateStation


def main(args):
    pfile = args.f

    ERA5     =  bool(args.d is None or "ERA5"  in args.d)
    ERA5ENS =  bool(args.d is None or "ERA5ENS"  in args.d)
    JRA     =  bool(args.d is None or "JRA"   in args.d)
    JRA3Q   = bool(args.d is None or "JRA3Q"   in args.d)
    JRA3QG   = bool(args.d is None or "JRA3QG"   in args.d)
    MERRA     = bool(args.d is None or "MERRA" in args.d)

    if sum([ERA5, ERA5ENS, JRA, MERRA, JRA3Q, JRA3QG]) > 0:

        GlobsimInterpolateStation(pfile,
                                  ERA5=ERA5,
                                  ERA5ENS=ERA5ENS,
                                  JRA=JRA,
                                  MERRA=MERRA,
                                  JRA3Q=JRA3Q,
                                  JRA3QG=JRA3QG,
                                  **vars(args))

    else: print("Failed! Reanalysis source should be  ERA5, MERRA, JRA, JRA3Q, please check")


# ===interpolate the variables from multiple re-analysis data to individual stations===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interpolate reanalysis data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f',    default=None, type=str, 
                        help="file path to download parameter file")
    parser.add_argument('-d',    default=None, nargs="*", type=str, 
                        help="What data sources should run?  ERA5, MERRA, JRA, JRA3Q, JRA3QG")

    args = parser.parse_args()

    main(args)
