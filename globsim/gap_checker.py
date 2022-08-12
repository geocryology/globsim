import argparse
import netCDF4 as nc
import numpy as np
import re
import scipy.stats as st

from typing import Optional


def check_file_time_integrity(pattern: "str", interval: "Optional[float]" = None) -> "tuple[np.ndarray,np.ndarray,np.ndarray]":
    ncdf = nc.MFDataset(pattern, aggdim='time')
    check_mf_time_integrity(ncdf, interval)


def check_mf_time_integrity(mfdataset: "nc.MFDataset", interval: "Optional[float]" = None) -> "tuple[np.ndarray,np.ndarray,np.ndarray]":
    t = find_time(mfdataset)

    time = mfdataset[t]
    
    result = check_time_integrity(time, interval)


def check_time_integrity(time: "nc.Variable", interval: "Optional[float]" = None) -> "tuple[np.ndarray,np.ndarray,np.ndarray]":
    if not interval:
        interval = estimate_interval(time[:])
    
    gaps = find_gaps(time, interval)
    
    gap_starts = time[:][gaps]
    gap_ends = time[:][gaps + 1]

    if not hasattr(time, 'units'):
        raise AttributeError(f"Time variable {time} missing units")

    calendar = time.calendar if hasattr(time, 'calendar') else 'standard'

    pytime = nc.num2date(time[:], time.units, calendar)

    return pytime, gap_starts, gap_ends


def find_gaps(var: "nc.Variable", interval: float):
    """ Find gaps in a variable with regular interval """
    diff = np.diff(var[:])

    gaps = np.where(diff != interval)[0]
    
    return gaps

#def report_gaps(pytime, gap_starts, gap_ends, interpolate_start, interpolate_end):


def find_time(ncdf: "nc.Dataset") -> str:
    # Search by unit name
    unit_pattern = re.compile(r"^\w* since ")
    candidates = [n for n,v in ncdf.variables.items() if unit_pattern.match(v.units)]

    if len(candidates) > 1:
        raise ValueError("Multiple time variables found")

    elif len(candidates) == 1:
        return candidates[0]

    # Search by axis attribute
    if ncdf.get_variables_by_attributes(axis='T'):
        return ncdf.get_variables_by_attributes(axis='T')

    # Search by name
    if any([re.search('time', y, re.IGNORECASE) for y in ncdf.variables]):
        return 'time'

    raise ValueError("No suitable time variable found")


def estimate_interval(i: np.ndarray) -> float:
    mode = st.mode(np.diff(i))
    
    if not len(mode.mode) == 1:
        raise ValueError("Could not guess interval (multiple modes found)")
    
    return mode.mode


def main(args):
    gap = find_gaps(args.pattern, args.interval)
    gap = find_gaps("/fs/yedoma/usr-storage/nbr512/esj-2022-merra/merra2/MERRA2_400.inst3_3d_asm_Np.*")
    gap = find_gaps("/fs/yedoma/data/globsim/n60/era5/era5_sa_*")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for gaps",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-I', '--interval', dest='interval', type=float,
                        help="Spacing between successive times")
    parser.add_argument('-D', '--dimension')
    parser.add_argument('-p', type=str, dest='pattern',
                        help="Pattern to match")

    args = parser.parse_args()

    main(args)
