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
from globsim.registry import parse_enabled
from globsim.registry import all_keys

def main(args):
    enabled = parse_enabled(args.d)
    if not any(enabled.values()):
        print("No reanalysis source selected. Use -d to specify one.")
        return
    GlobsimInterpolateStation(args.f, **enabled, **vars(args))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interpolate reanalysis data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f',    default=None, type=str, 
                        help="file path to download parameter file")
    parser.add_argument('-d',    default=None, nargs="*", type=str.upper, choices=all_keys(),
                        help="What data sources should run?")

    args = parser.parse_args()

    main(args)
