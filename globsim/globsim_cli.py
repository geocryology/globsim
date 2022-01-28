"""globsim.globsim: provides entry point main()."""

import logging
import sys
import argparse

from pathlib import Path
from globsim import globsim_download, globsim_scale, globsim_interpolate
from globsim._version import __version__

def configure_logging(args):
    #logging.basicConfig(format='%(asctime)s  %(asctime)s ')
    logger.setLevel(args.level)

    console_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', datefmt="%H:%M:%S")
    # console_formatter = logging.Formatter('%(message)s', datefmt="%H:%M:%S")

    ch = logging.StreamHandler()
    ch.setLevel(args.level)
    ch.setFormatter(console_formatter)
    logger.addHandler(ch)

    if args.logfile:
        file_formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
        logfile = args.logfile  # TODO: write logfile to project directory if missing
        fh = logging.FileHandler(logfile)
        fh.setLevel(args.level)
        fh.setFormatter(file_formatter)
        logger.addHandler(fh)

action_dict = {'download': globsim_download,
               'interpolate': globsim_interpolate,
               'scale': globsim_scale}


logger = logging.getLogger("globsim")

def main():
    parser = argparse.ArgumentParser(description="GlobSim: meteorological reanalysis for point-scale simulation. Find out more at https://globsim.readthedocs.io/en/latest",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--version", action='version', version=f"GlobSim version {__version__}")

    parser.add_argument('action', default=None, type=str, choices=["download", "interpolate", "scale"],
                        help="Which GlobSim action should be run.")

    parser.add_argument("-f", "--config-file", default=None, type=str,
                        dest='f', help="Path to GlobSim configuration file.")

    parser.add_argument("-d", default=None, nargs="*", type=str, choices=["ERAI", "ERA5", "ERA5ENS", "MERRA", "JRA"],
                        dest='d', help="What data sources should run?")

    parser.add_argument('-r', '--retry',  default=1,    type=int, help="Number of times to re-launch download if it crashes (globsim_download only) ")
    parser.add_argument('-m', '--multi', action='store_true',    help="Download all data sources simultaneously (globsim_download only) ")
    parser.add_argument("-v", "--verbose", nargs="?", default=logging.INFO, const=logging.DEBUG, dest='level',
                        help="Show detailed output")
    parser.add_argument("-L", "--logfile", default=False, nargs="?", const=True, dest='logfile',
                        help="Whether to save output to a logfile. If a path is provided, it will be used instead of the default path")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    else:
        args = parser.parse_args()
        configure_logging(args)
        logger.info(f"Running globsim {args.action}")
        task = action_dict[args.action]
        task.main(args)

