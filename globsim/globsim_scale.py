import argparse
from globsim.globsim_main import GlobsimScale
from globsim.registry import parse_enabled
from globsim.registry import all_keys

def main(args):
    enabled = parse_enabled(args.d)
    if not any(enabled.values()):
        print("No reanalysis source selected. Use -d to specify one.")
        return
    GlobsimScale(args.f, **enabled)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scale reanalysis data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', default=None, type=str, 
                        help="file path to download parameter file")
    parser.add_argument('-d', default=None, nargs="*", type=str.upper, choices=all_keys(),
                        help="What data sources should run?")
    
    
    args = parser.parse_args()

    main(args)
