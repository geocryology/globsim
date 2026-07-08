"""
Export functions to convert globsim output to other file types
"""
from os import path, makedirs
from pathlib import Path

import datetime
import logging
from typing import Union
import netCDF4 as nc
import numpy as np
import pandas as pd
import pkg_resources
import tomlkit
import shutil
import bisect 
from cfunits import Units
import bisect 

from globsim.common_utils import variables_skip, get_scaled_site_names

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


def globsim_to_classic_met(ncd, out_dir, site: int|str|list[str|int]|None=None, export_profile=None, start=None, end=None,
                           file_vars=["sw", "lw", "precip", "airt", "sh", "wspd", "press"]):
    """
    @args
    ncd : str
        path to output globsim dataset
    out_dir : str
        where files should be saved
    site : int or str
        index of site (zero-indexed) or site name if only one site is desired.
        For all sites, leave as None
    export_profile : str or None
        path to a TOML export-profile file.  When *None* a default profile
        is created / used at ``~/.globsim/classic_met_profile.toml``.

    Returns
    list : paths of output files
    """
    if isinstance(site, (int, str)):
        site = [site]  
    # open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)
    else:
        n = ncd

    # find number of stations
    nstn = len(n['station'][:])

    # get date / time columns
    time = nc.num2date(n['time'][:],
                       units=n['time'].units,
                       calendar=n['time'].calendar)
    if isinstance(time, datetime.datetime):
        time = [time]
        
    time_slice = time_slice_index(time, start=start, end=end)
    sliced_time = time[time_slice]

    HH = [x.timetuple().tm_hour for x in sliced_time]
    MM = [x.timetuple().tm_min for x in sliced_time]
    DDD = [x.timetuple().tm_yday for x in sliced_time]
    YYYY = [x.timetuple().tm_year for x in sliced_time]

    TIME = np.stack((HH, MM, DDD, YYYY))
    
    # Open classic profile
    profile = get_export_profile("classic_met", export_profile=export_profile)
    
    data = []

    for var in file_vars:
        # get variable
        var_profile = profile.get(var, {})
        
        if not var_profile:
            logger.critical(f"No profile information found for variable '{var}' in export profile. Skipping.")
            continue
        
        input_var_name = var_profile.get("input")

        if input_var_name not in n.variables:
            logger.critical(f"Variable '{input_var_name}' not found in netCDF file. Skipping {var}.")
            continue

        nc_var = n[input_var_name]
        output_units = var_profile.get("output_units")
        var_data = nc_var[time_slice]
        
        if output_units:
            var_data = Units.conform(var_data, Units(nc_var.units), Units(output_units))
        
        else:
            scale_factor = var_profile.get("scale_factor", 1)
            offset = var_profile.get("offset", 0)
            var_data = var_data * scale_factor + offset  # apply scale factor and offset if output_units not provided

        data.append(var_data)

    # get site names
    names = get_scaled_site_names(n)

    STN_ID  = n['station'][:].astype('str')  # integer values of station index
    NAMES = [f"{i}_{n}" for n, i in zip(names, STN_ID)]  

    data = np.stack(data)

    # write output files
    files = []
    for i in range(nstn):
        if (site is None) or (names[i] in site) or (str(STN_ID[i]) in site):
            print(f"Processing site {i} ({STN_ID[i]} - {names[i]})")
            # massage data into the right shape
            out_array = data[:, :, i]
            out_array = np.concatenate((TIME, out_array), 0)
            out_array = np.transpose(out_array)

            # get station name
            st_name = NAMES[i]

            filename = f"{st_name}.MET"
            savepath = path.join(out_dir, filename)
            files.append(savepath)
            # create file
            np.savetxt(savepath, out_array,
                       newline='\n',  # windows line endings breaks the converter
                       fmt=[" %02u", "%02u", " %03u",
                            " %2u", "%8.2f", "%8.2f",
                            "%13.4E", "%8.2f", "%11.3E",
                            "%7.2f", "%11.2f"])
            logger.info(f"Created classic MET file '{savepath}' for site '{st_name}'")
    return files


def globsim_to_svs2(ncd, out_dir, site=None, export_profile=None, start=None, end=None):
    """
    Export a scaled globsim file to SVS2-style text files

    Parameters
    ----------
    ncd : str or Dataset
        netcdf dataset or path to dataset
    site : str or int
        site name or index
    export_profile : str or None
        path to a TOML export-profile file.  When *None* a default profile
        is created / used at ``~/.globsim/svs2_profile.toml``.

    Returns
    -------
    list : list of file paths to created files
    """
    profile_file = get_export_profile_file("svs2", export_profile=export_profile)
    return globsim_to_classic_met(ncd=ncd, out_dir=out_dir, site=site, export_profile=profile_file,
                                  start=start, end=end,
                                  file_vars=["FSIN", "FLIN", "PRE", "TA", "QA", "UV", "PRES", "PRERN", "PRESNO"])


def globsim_to_geotop(ncd, out_dir, site=None, export_profile=None, start=None, end=None) -> "list[str]":
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
    profile = get_export_profile("geotop", export_profile=export_profile)

    # open netcdf if string provided
    if type(ncd) is str:
        ncd = nc.Dataset(ncd)

    logger.debug(f"Read file {ncd.filepath}")

    # find number of stations
    nstn = len(ncd['station_name'][:])

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
            arr = ((ncd["SW_topo_diffuse"][:] + ncd["SW_topo_direct"][:]) if var_name == "SW_topo" else ncd[var_name][:]) * scale_factor + offset
            output_dict[out_var] = arr
        
        except IndexError:
            logger.error(f"Scaled netCDF file has no variable '{var_name}'")

    # get site names
    names = get_scaled_site_names(ncd)

    # combine data variables into array
    data = np.stack([arr for arr in output_dict.values()])

    # write output files
    files = []
    for i in range(nstn):
        if (site is None) or (site == i) or (site == names[i]):
            # massage data into the right shape
            out_df = pd.DataFrame(np.transpose(data[:, :, i]))
            out_df = pd.concat([time, out_df], axis=1)
            out_df.columns = ["Date"] + [name for name in output_dict.keys()]

            # get station name
            st_name = names[i]

            # prepare paths
            filename = "{}-{}_Forcing_0001.txt".format(i, st_name)
            savepath = path.join(out_dir, filename)
            files.append(savepath)

            # create file
            logger.info(f"writing file '{savepath}'")
            out_df.to_csv(savepath, index=False, float_format="%10.5f")

    return files


def globsim_to_freethaw(ncd, out_dir, site=None, export_profile=None, start=None, end=None) -> "list[str]":
    """
    Export a scaled globsim file to a Freethaw-style (OMS3) text file

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
    # open netcdf if string provided
    if type(ncd) is str:
        ncd = nc.Dataset(ncd)

    logger.debug(f"Read file '{ncd.filepath()}'")

    # find number of stations
    nstn = len(ncd['station'][:])

    # get date / time column
    time = nc.num2date(ncd['time'][:],
                       units=ncd['time'].units,
                       calendar=ncd['time'].calendar)

    time = [x.strftime('%Y-%m-%d %H:%M') for x in time]  # TODO: this could be faster
    time = pd.DataFrame(time)

    # write headers lines
    headers = [f"@Table,table,\n",
               f"Created,{datetime.datetime.now().isoformat()[:19]},\n",
               "Author,Globsim,\n",
               "ID,,0\n",
               "Type,Date,double\n",
               "@Header,timestamp,temperature\n",
               "Format,yyyy-MM-dd HH:mm,,\n"]

    AIRT = "AIRT_sur"
    AIRT = ncd[AIRT][:]

    # get site names
    NAMES = nc.chartostring(ncd['station_name'][:])

    # combine data variables into array
    data = AIRT

    # write output files
    files = []
    for i in range(nstn):
        if (site is None) or (site == i) or (site == NAMES[i]):
            # massage data into the right shape
            out_df = pd.DataFrame(data[:,i])
            out_df = pd.concat([time, out_df], axis=1)
            out_df.insert(loc=0, column="Format",value="")

            # get station name
            st_name = NAMES[i]

            # prepare paths
            filename = "{}-{}_Temperature.csv".format(i, st_name)
            savepath = path.join(out_dir, filename)
            files.append(savepath)

            # create file
            logger.info(f"writing file '{savepath}'")
            with open(savepath, 'w') as file:
                file.writelines(headers)
                file.write(out_df.to_csv(header=False, index=False))

    return files


def get_export_profile_file(model, export_profile:Union[None, Path, str]=None) -> Union[None, Path, str]:
    if export_profile is None:

        export_profile = Path(f"~/.globsim/{model}_profile.toml").expanduser()

        if not Path(export_profile).is_file():
            default = pkg_resources.resource_filename("globsim", f"data/{model}_profile_default.toml")
            
            if not Path(default).is_file():
                logger.critical(f"Default export profile for model '{model}' not found in package resources. This is likely a bug. Expected path: {default}")
                raise FileNotFoundError(f"Default export profile for model '{model}' not found in package resources. This is likely a bug. Expected path: {default}")
            
            if not export_profile.parent.is_dir():
                makedirs(export_profile.parent)
            
            shutil.copy(default, export_profile)
            logger.warning(f"Created default {model} export profile: {export_profile}")
    
    return export_profile


def get_export_profile(model, export_profile=None) -> "tomlkit.TOMLDocument":
    
    export_profile_file = get_export_profile_file(model, export_profile=export_profile)
    
    with open(export_profile_file) as p:
        profile = tomlkit.loads(p.read())
        logger.info(f"Loaded {model} export profile from {export_profile_file}")

    return profile



def time_slice_index(times_as_dates, start=None, end=None) -> slice:
    """
    Return a slice object that can be used to index a time series between two dates

    Parameters
    ----------
    times_as_dates : list of datetime
        list of datetime objects, MUST BE SORTED
    start : str, optional
        start date in YYYY-MM-DD format. If not provided, the first date in the list is used.
    end : str, optional
        end date in YYYY-MM-DD format. If not provided, the last date in the list is used.

    Returns
    -------
    slice : slice object that can be used to index the time series
    """

    if start is None:
        start_index = 0
    else:
        start = pd.to_datetime(start)
        # bisect_left finds the first index where t >= start
        start_index = bisect.bisect_left(times_as_dates, start)
    
    if end is None:
        end_index = len(times_as_dates)
    else:
        end = pd.to_datetime(end)
        # bisect_right finds the first index where t > end
        end_index = bisect.bisect_right(times_as_dates, end)
    
    return slice(start_index, end_index)

