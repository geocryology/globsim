"""globsim.globsim: provides entry point main()."""


__version__ = "1.0.1"


import sys
import argparse

"""
from .era_interim import ERAIdownload, ERAIinterpolate, ERAIscale
from .era5 import ERA5download, ERA5interpolate, ERA5scale
from .merra_2 import MERRAdownload, MERRAinterpolate, MERRAscale
from .JRA import JRAdownload, JRAinterpolate, JRAscale
"""

from .globsim_download import main as globsim_download

def main():
    parser = argparse.ArgumentParser(description="GlobSim: meteorological reanalysis for point-scale simulation. Find out more at https://globsim.readthedocs.io/en/latest",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('action', default=None, type=str, choices=["download", "interpolate", "scale"],
                        help="Which GlobSim action should be run.")

    parser.add_argument("-f", "--config-file", default=None, type=str,
                        dest='file', help="Path to GlobSim configuration file.")

    parser.add_argument("-d", default=None, nargs="*", type=str, choices=["ERAI", "ERA5", "ERA5ENS", "MERRA", "JRA"],
                        dest='data', help="What data sources should run?")

    parser.add_argument('-r', '--retry',  default=1,    type=int, help="Number of times to re-launch download if it crashes (globsim_download only) ")
    parser.add_argument('-m', '--multi', action='store_true',    help="Download all data sources simultaneously (globsim_download only) ")

    parser.add_argument("-v", "--version", action='version', version=f"GlobSim version {__version__}")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    print(f"Running globsim_{args.action} for {', '.join(args.data)}")
