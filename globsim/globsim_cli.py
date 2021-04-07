"""globsim.globsim: provides entry point main()."""


import sys
import argparse

from globsim import globsim_download, globsim_scale, globsim_interpolate
from globsim._version import __version__

action_dict = {'download': globsim_download,
               'interpolate': globsim_interpolate,
               'scale': globsim_scale}


def main():
    parser = argparse.ArgumentParser(description="GlobSim: meteorological reanalysis for point-scale simulation. Find out more at https://globsim.readthedocs.io/en/latest",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v", "--version", action='version', version=f"GlobSim version {__version__}")

    parser.add_argument('action', default=None, type=str, choices=["download", "interpolate", "scale"],
                        help="Which GlobSim action should be run.")

    parser.add_argument("-f", "--config-file", default=None, type=str,
                        dest='f', help="Path to GlobSim configuration file.")

    parser.add_argument("-d", default=None, nargs="*", type=str, choices=["ERAI", "ERA5", "ERA5ENS", "MERRA", "JRA"],
                        dest='d', help="What data sources should run?")

    parser.add_argument('-r', '--retry',  default=1,    type=int, help="Number of times to re-launch download if it crashes (globsim_download only) ")
    parser.add_argument('-m', '--multi', action='store_true',    help="Download all data sources simultaneously (globsim_download only) ")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    else:
        args = parser.parse_args()
        print(f"Running globsim_{args.action}")
        task = action_dict[args.action]
        task.main(args)

