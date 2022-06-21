"""
Export functions to convert globsim output to other file types
"""
from os import path, makedirs
from pathlib import Path

import logging
import netCDF4 as nc
import numpy as np
import pandas as pd
import pkg_resources
import tomlkit
import shutil

from globsim.common_utils import variables_skip

logger = logging.getLogger("globsim.convert")


def globsimScaled2Pandas(ncdf_in, station_nr):
    """
    Read a scaled (or interpolated) globsim netCDF file and return all values
    for one station as a Pandas data frame.

    ncdf_in: full path to a globsim netCDF (by station)

    station_nr: station_number, as given in the stations .csv file to identify
                the station.

    """
    # open file
    ncf = nc.Dataset(ncdf_in, 'r')

    # station mask
    sm = ncf.variables['station'][:] == int(station_nr)
    # list variables
    varlist = [x.encode('UTF8') for x in ncf.variables.keys()]

    # get and convert time
    time = ncf.variables['time'][:]
    t_unit = ncf.variables['time'].units
    t_cal = ncf.variables['time'].calendar
    time = nc.num2date(time, units=t_unit, calendar=t_cal)

    # make data frame with time
    df = pd.DataFrame(data=time, columns=['time'])
    # add variables
    for var in varlist:
        if variables_skip(var):
            continue
        data = ncf.variables[var][:,sm]
        df = pd.concat([df, pd.DataFrame(data=data, columns=[var])], axis=1)

    return df


def globsim_to_classic_met(ncd, out_dir, site=None):
    """
    @args
    ncd : str
        path to output globsim dataset
    out_dir : str
        where files should be saved
    site : int or str
        index of site (zero-indexed) or site name if only one site is desired.
        For all sites, leave as None

    Returns
    list : paths of output files
    """
    # open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)

    # find number of stations
    nstn = len(n['station'][:])

    # get date / time columns
    time = nc.num2date(n['time'][:],
                       units=n['time'].units,
                       calendar=n['time'].calendar)
    time_step = (time[1] - time[0]).seconds  # in second

    HH = [x.timetuple().tm_hour for x in time]
    MM = [x.timetuple().tm_min for x in time]
    DDD = [x.timetuple().tm_yday for x in time]
    YYYY = [x.timetuple().tm_year for x in time]

    TIME = np.stack((HH, MM, DDD, YYYY))

    # get shortwave
    SW = "SW_sur"
    SW = n[SW][:]

    # get longwave
    LW = "LW_sur"
    LW = n[LW][:]

    # get precip
    PREC = "PREC_sur"
    PREC = n[PREC][:]  # Defaults to mm/s (CLASSIC-compatible)

    # get temp
    AIRT = "AIRT_sur"
    AIRT = n[AIRT][:]

    # get specific humidity
    SH = "SH_sur"
    SH = n[SH][:]

    # get wind speed
    WSPD = "WSPD_sur"
    WSPD = n[WSPD][:]

    # get pressure
    PRESS = "PRESS_pl"
    PRESS = n[PRESS][:]

    # get site names
    NAMES = nc.chartostring(n['station_name'][:])

    data = np.stack((SW, LW, PREC, AIRT, SH, WSPD, PRESS))

    # write output files
    files = []
    for i in range(nstn):
        if (site is None) or (site == i) or (site == NAMES[i]):
            # massage data into the right shape
            out_array = data[:, :, i]
            out_array = np.concatenate((TIME, out_array), 0)
            out_array = np.transpose(out_array)

            # get station name
            st_name = NAMES[i]

            filename = "{}_{}.MET".format(i, st_name)
            savepath = path.join(out_dir, filename)
            files.append(savepath)
            # create file
            np.savetxt(savepath, out_array,
                       newline='\n',  # windows line endings breaks the converter
                       fmt=[" %02u", "%02u", " %03u",
                            " %2u", "%8.2f", "%8.2f",
                            "%13.4E", "%8.2f", "%11.3E",
                            "%7.2f", "%11.2f"])
    return files


def globsim_to_geotop(ncd, out_dir, export_profile=None, site=None, start=None, end=None) -> "list[str]":
    """
    Export a scaled globsim file to a GEOtop-style text file

    Parameters
    ----------
    ncd : str or Dataset
        netcdf dataset or path to dataset
    site : str or int
        site name or index
    export_profile : str, optional
        path to TOML file that provides configuration information. If not provided, a default configuration is created in your home
        directory in a .globsim folder (i.e. '~/.globsim/geotop_profile.toml'). For example configuration files see the globsim/data/ folder
        in the package.

    Returns
    -------
    list : list of file paths to created files
    """
    # Open geotop profile
    if export_profile is None:

        export_profile = Path("~/.globsim/geotop_profile.toml").expanduser()

        if not Path(export_profile).is_file():
            default = pkg_resources.resource_filename("globsim", "data/geotop_profile_default.toml")
            
            if not export_profile.parent.is_dir():
                makedirs(export_profile.parent)
            
            shutil.copy(default, export_profile)
            logger.warning(f"Created default geotop export profile: {export_profile}")

    with open(export_profile) as p:
        profile = tomlkit.loads(p.read())
        logger.info(f"Loaded geotop export profile from {export_profile}")

    # open netcdf if string provided
    if type(ncd) is str:
        ncd = nc.Dataset(ncd)

    logger.debug(f"Read file {ncd.filepath}")

    # find number of stations
    nstn = len(ncd['station'][:])

    # get date / time column
    time = nc.num2date(ncd['time'][:],
                       units=ncd['time'].units,
                       calendar=ncd['time'].calendar)

    time = [x.strftime('%d/%m/%Y %H:%M') for x in time]
    time = pd.DataFrame(time)

    output_dict = {}

    for out_var, cfg in profile.items():
        var_name = cfg.get("input")
        scale_factor = cfg.get("scale_factor", 1)
        offset = cfg.get("offset", 0)
        
        try:
            arr = ncd[var_name][:] * scale_factor + offset
            output_dict[out_var] = arr
        
        except IndexError:
            logger.error(f"Scaled netCDF file has no variable '{var_name}'")

    # get site names
    NAMES = nc.chartostring(ncd['station_name'][:])

    # combine data variables into array
    data = np.stack([arr for arr in output_dict.values()])

    # write output files
    files = []
    for i in range(nstn):
        if (site is None) or (site == i) or (site == NAMES[i]):
            # massage data into the right shape
            out_df = pd.DataFrame(np.transpose(data[:, :, i]))
            out_df = pd.concat([time, out_df], axis=1)
            out_df.columns = ["Date"] + [name for name in output_dict.keys()]

            # get station name
            st_name = NAMES[i]

            # prepare paths
            filename = "{}-{}_Forcing_0001.txt".format(i, st_name)
            savepath = path.join(out_dir, filename)
            files.append(savepath)

            # create file
            logger.info(f"writing file '{savepath}'")
            out_df.to_csv(savepath, index=False, float_format="%10.5f")

    return files
