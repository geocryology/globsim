from globsim.globsim_main import GlobsimDownload
from globsim.registry import parse_enabled
from globsim.registry import all_keys

import argparse
import time


def main(args):
    enabled = parse_enabled(args.d)
    r_max = args.retry

    for attempt in range(max(r_max, 1)):
        try:
            GlobsimDownload(args.f, multithread=args.multi, **enabled)
            break  # success — stop retrying
        except Exception as e:
            print(e)
            if attempt < r_max - 1:
                time.sleep(360)


# =====download the wanted variable from multiple re-analysis date sources======
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get MET data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f',             default=None, type=str, help="file path to download parameter file")
    parser.add_argument('-d',             default=None, type=str.upper, nargs="*", choices=all_keys(), help="What data sources should run?")
    parser.add_argument('-r', '--retry',  default=1,    type=int, help="Number of times to re-launch download if it crashes.")
    parser.add_argument('-m', '--multi', action='store_true',    help="Download all data sources simultaneously")

    args = parser.parse_args()

    main(args)
