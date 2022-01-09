#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic classes, methods, functions used for more than one reanalysis.
#
#
# (C) Copyright Stephan Gruber (2017)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ===============================================================================
from __future__ import print_function

from datetime import datetime, timedelta
from os import mkdir, path

import pandas as pd
import netCDF4 as nc
import numpy as np

import glob

# handle python 3 string types
try:
    basestring
except NameError:
    basestring = str


def variables_skip(variable_name):
    """
    Which variable names to use? Drop the ones that are dimensions.
    """
    skip = 0
    dims = ('time', 'number', 'level',
            'latitude', 'longitude', 'station', 'height')
    if variable_name in dims:
        skip = 1
    return skip


def StationListRead(sfile):
    """
    Reads ASCII station list and returns a pandas dataframe.

    # read station list
    stations = StationListRead('examples/par/examples_list1.globsim_interpolate')
    print(stations['station_number'])
    """
    # read file; allow for comma or semicolon as delimiter
    raw = pd.read_csv(sfile)
    raw = raw.rename(columns=lambda x: x.strip())
    print('Keys: '+ raw.keys())
    if len(raw.keys()) < 2:
        raw = pd.read_csv(sfile, sep=';')
        raw = raw.rename(columns=lambda x: x.strip())
    print('Keys: '+ raw.keys())
    return(raw)


def convert_cummulative(data):
    """
    Convert values that are serially cummulative, such as precipitation or
    radiation, into a cummulative series from start to finish that can be
    interpolated on for sacling.
    data: 1-dimensional time series
    """
    # get increment per time step
    diff = np.diff(data)
    diff = np.concatenate(([data[0]], diff))

    # where new forecast starts, the increment will be smaller than 0
    # and the actual value is used
    mask = diff < 0
    diff[mask] = data[mask]

    # get full cummulative sum
    return np.cumsum(diff, dtype=np.float64)


def cummulative2total(data, time):
    """
    Convert values that are serially cummulative, such as precipitation or
    radiation, into a cummulative series from start to finish that can be
    interpolated on for sacling.
    data: 1-dimensional time series
    """
    # get increment per time step
    diff = np.diff(data)
    diff = np.concatenate(([data[0]], diff))

    # where new forecast starts, the increment will be smaller than 0
    # and the actual value is used

    mask = [timei.hour in [3, 15] for timei in time]
    diff[mask] = data[mask]

    mask = diff < 0
    diff[mask] = 0

    return diff


def get_begin_date(par, data_folder, match_strings):
    """ Get the date to begin downloading when some files already exist

    Parameters
    ----------
    par : dict
        download section of configuration file as read in by tomlkit
    data_folder : str
        name of subdirectory containing data files. Examples: merra2, era5
    match_strings : list
        list of glob-style strings to check. Examples ["merra_pl*", "merra_sa*","merra_sf*"]
    Returns
    -------
    datetime
        datetime object corresponding to the desired begin date (replaces par['beg'])
    This makes an inventory of all the files that have been downloaded so far and
    returns the next date to begin downloading.  If all match_strings are downloaded up to the same
    day, then the following day is returned. Otherwise, the
    """
    directory = par['project_directory']
    print("Searching for existing files in directory")

    if not all([len(glob.glob(path.join(directory, data_folder, s))) > 0 for s in match_strings]):
        print("No existing files found. Starting download from {}".format(par['beg'].strftime("%Y-%m-%d")))
        begin_date = datetime.strptime(par['beg'], '%Y/%m/%d')
    else:
        datasets = [nc.MFDataset(path.join(directory, data_folder, s)) for s in match_strings]
        dates = [nc.num2date(x['time'][:], x['time'].units, x['time'].calendar) for x in datasets]

        latest = [max(d) for d in dates]
        latest = [dt.replace(hour=0, minute=0, second=0, microsecond=0) for dt in latest]
        latest_complete = min(latest)
        begin_date = latest_complete + timedelta(days=1)

        print("Found some files in directory. Beginning download on {}".format(begin_date.strftime("%Y-%m-%d")))

    return(begin_date)


def series_interpolate(time_out, time_in, value_in, cum=False):
    """
    Interpolate single time series. Convenience function for usage in scaling
    kernels.
    time_out: Array of times [s] for which output is desired. Integer.
    time_in:  Array of times [s] for which value_in is given. Integer.
    value_in: Value time series. Must have same length as time_in.
    cum:      Is valiable serially cummulative like LWin? Default: False.
    """
    time_step_sec = time_out[1] - time_out[0]

    # convert to continuous cummulative, if values are serially cummulative
    if cum:
        value_in = convert_cummulative(value_in)

    # interpolate
    vi = np.interp(time_out, time_in, value_in)

    # convert from cummulative to normal time series if needed
    if cum:
        vi = np.diff(vi) / time_step_sec
        vi = np.float32(np.concatenate(([vi[0]], vi)))

    return vi


def str_encode(value, encoding="UTF8"):
    """
    handles encoding to allow compatibility between python 2 and 3
    specifically with regards to netCDF variables. Python 2 imports
    variable names as unicode, whereas python 3 imports them as str.
    """
    if type(value) == str:
        return(value)
    else:
        return(value.encode(encoding))


def create_globsim_directory(target_dir, name):
    """
    creates globsim directory
    """
    # create top-level
    TL = path.join(target_dir, name)
    mkdir(TL)

    # create subdirectories
    mkdir(path.join(TL, "eraint"))
    mkdir(path.join(TL, "Grib"))
    mkdir(path.join(TL, "jra55"))
    mkdir(path.join(TL, "merra2"))
    mkdir(path.join(TL, "par"))
    mkdir(path.join(TL, "scale"))
    mkdir(path.join(TL, "station"))
    mkdir(path.join(TL, "era5"))

    return(True)


def get_begin_date(par, data_folder, match_strings):
    """ Get the date to begin downloading when some files already exist

    Parameters
    ----------
    par : dict
        download section of configuration file as read in by tomlkit
    data_folder : str
        name of subdirectory containing data files. Examples: merra2, era5
    match_strings : list
        list of glob-style strings to check. Examples ["merra_pl*", "merra_sa*","merra_sf*"]

    Returns
    -------
    datetime
        datetime object corresponding to the desired begin date (replaces par['beg'])

    This makes an inventory of all the files that have been downloaded so far and
    returns the next date to begin downloading.  If all match_strings are downloaded up to the same
    day, then the following day is returned. Otherwise, the
    """
    directory = par['project_directory']
    print("Searching for existing files in directory")

    if not all([len(glob.glob(path.join(directory, data_folder, s))) > 0 for s in match_strings]):
        print("No existing files found. Starting download from {}".format(par['beg']))
        return par['beg']

    datasets = [nc.MFDataset(path.join(directory, data_folder, s)) for s in match_strings]
    dates = [nc.num2date(x['time'][:], x['time'].units, x['time'].calendar) for x in datasets]

    latest = [max(d) for d in dates]
    latest = [dt.replace(hour=0, minute=0, second=0, microsecond=0) for dt in latest]
    latest_complete = min(latest)

    begin_date = latest_complete + timedelta(days=1)

    print("Found some files in directory. Beginning download on {}".format(begin_date.strftime("%Y-%m-%d")))
    return(begin_date)
