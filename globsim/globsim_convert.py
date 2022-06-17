from globsim.exporttools import globsim_to_geotop, globsim_to_classic_met

import argparse

export_styles = {"geotop": globsim_to_geotop,
                 "classicmet": globsim_to_classic_met}


def main(args):
    source_file = args.file
    dest_directory = args.output
    export_type = args.format
    site = args.site
    export_profile = args.profile

    if export_type == "geotop":
        globsim_to_geotop(ncd=source_file,
                          out_dir=dest_directory,
                          export_profile=export_profile,
                          site=site)

    elif export_type == "classicmet":
        globsim_to_classic_met(ncd=source_file,
                               out_dir=dest_directory,
                               site=site)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export globsim data to desired format",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', "--file", dest='file', default=None, required=True, type=str, help="Path to scaled Globsim file")
    parser.add_argument('-F', "--format", dest='format', default=None, required=True, type=str, choices=export_styles.keys(), help=f"What kind of output to generate. Chosen from {export_styles.keys()}")
    parser.add_argument('-o', "--output", dest='output', default=None, type=str, help="Output directory to write new files")
    parser.add_argument('-p', "--profile", dest='profile', default=None, type=str, help="Path to an 'export profile' TOML file (geotop only) ")
    parser.add_argument('-s', "--site", dest='site', default=None, nargs="*", type=str, help="(optional)")

    args = parser.parse_args()

    main(args)
