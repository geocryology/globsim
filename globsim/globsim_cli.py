"""globsim.globsim: provides entry point main()."""

import logging
import sys
import argparse

from globsim import globsim_convert, globsim_download, globsim_scale, globsim_interpolate
from globsim.LazyLoader import LazyLoader
from globsim.globsim_convert import export_styles
from globsim._version import __version__

gsview = LazyLoader("globsim.view.view_main")

def configure_logging(args: argparse.Namespace):
    # logging.basicConfig(format='%(asctime)s  %(asctime)s ')
    try:
        level = args.level
    except AttributeError:
        level = logging.INFO

    logger.setLevel(level)
    console_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', datefmt="%H:%M:%S")
    # console_formatter = logging.Formatter('%(message)s', datefmt="%H:%M:%S")

    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(console_formatter)
    logger.addHandler(ch)

    if getattr(args, "logfile", None):
        file_formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
        logfile = args.logfile  # TODO: write logfile to project directory if missing
        fh = logging.FileHandler(logfile)
        fh.setLevel(level)
        fh.setFormatter(file_formatter)
        logger.addHandler(fh)


action_dict = {'download': globsim_download,
               'interpolate': globsim_interpolate,
               'scale': globsim_scale,
               'convert': globsim_convert}


logger = logging.getLogger("globsim")


def main():
    mainparser = argparse.ArgumentParser(description="GlobSim: meteorological reanalysis for point-scale simulation." +
                                                     "Find out more at https://globsim.readthedocs.io/en/latest",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = mainparser.add_subparsers()

    download = subparsers.add_parser("download", description="Download reanalysis data from a provider")
    interpolate = subparsers.add_parser("interpolate")
    scale = subparsers.add_parser("scale")
    convert = subparsers.add_parser("convert")
    view = subparsers.add_parser("view")

    mainparser.add_argument("--version", action='version', version=f"GlobSim version {__version__}")

    for parser in [download, interpolate, scale]:

        parser.add_argument("-f", "--config-file",
                            default=None, type=str, required=True, dest='f',
                            help="Path to GlobSim configuration file.")

        parser.add_argument("-d",
                            default=None, nargs="*",required=True, type=str.upper,
                            choices=["ERA5", "ERA5ENS", "MERRA", "JRA", "JRA3Q", "JRA3QG"],
                            dest='d', help="What data sources should run?")

    for parser in [download, interpolate, scale, convert]:
        parser.add_argument("-v", "--verbose", nargs="?", default=logging.INFO, const=logging.DEBUG, dest='level',
                            help="Show detailed output")
        parser.add_argument("-L", "--logfile", default=False, nargs="?", const=True, dest='logfile',
                            help="Whether to save output to a logfile. If a path is provided, it will be used instead of the default path")

    download.add_argument('-r', '--retry',  default=1, type=int, help="Number of times to re-launch download if it crashes")
    download.add_argument('-m', '--multi', action='store_true',    help="Download all data sources simultaneously ")

    interpolate.add_argument('--skip-sa', dest='skip_sa', action='store_true', help='skip creation of "sa" file')
    interpolate.add_argument('--skip-sf', dest='skip_sf', action='store_true', help='skip creation of "sf" file')
    interpolate.add_argument('--skip-pl', dest='skip_pl', action='store_true', help='skip creation of "pl" file')
    interpolate.add_argument('--resume', dest='resume', action='store_true', help='pick up interpolation where it left off')
    interpolate.add_argument('--skip-checks', dest='skip_checks', action='store_true', help='skip data integrity checks')
    
    download.set_defaults(func=globsim_download.main)
    interpolate.set_defaults(func=globsim_interpolate.main)
    scale.set_defaults(func=globsim_scale.main)
    convert.set_defaults(func=globsim_convert.main)

    convert.add_argument('-f', "--file", dest='file', default=None, required=True, type=str,
                         help="Path to scaled Globsim *.nc file")
    convert.add_argument('-F', "--format", dest='format', default=None, required=True, type=str, choices=export_styles.keys(),
                         help=f"What kind of output to generate. Chosen from {export_styles.keys()}")
    convert.add_argument('-o', "--output", dest='output', default=None, required=True, type=str,
                         help="Output directory to write new files")
    convert.add_argument('-s', "--site", dest='site', default=None, nargs="*", type=str,
                         help="(optional) The name of the site you want to export. If not provided, all sites will be exported")
    convert.add_argument('-p', "--profile", dest='profile', default=None, type=str, help="Path to an 'export profile' TOML file (geotop only) ")

    view.set_defaults(func=gsview.main_args)
    view.add_argument("file", nargs="?", type=str, help="file to plot")
    view.add_argument("--file", dest='file', type=str, help="file to plot")
    view.add_argument("reanalysis",  type=str, choices=('era5','jra3qg','merra'), nargs='?', default=None,
                        help="if file is a TOML file, specify the reanalysis to plot.")
    view.add_argument("ftype",  type=str, choices=('sa','pl','sf'), nargs='?', default=None,
                        help="if file is a TOML file, specify the type of file to plot.")
    view.add_argument("-v", "--var", type=str, dest='variable', help="variable to plot")
    view.add_argument("-a", "--agg", choices=["1h", "6h", "D", "ME", "YE"], dest='aggregate', default="ME", help="aggregate data")
    view.add_argument("-o", "--output", type=str, dest='output', help="output directory")
   
    if len(sys.argv) == 1:
        mainparser.print_help(sys.stderr)
        sys.exit(1)

    else:
        args = mainparser.parse_args()
        configure_logging(args)
        args.func(args)

