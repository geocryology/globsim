"""
Export functions to convert globsim output to other file types
"""
from os import path

import pandas as pd
import netCDF4 as nc
import numpy as np
import logging

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


def globsim2CLASS(ncdf_globsim, met_class, station_nr):
    """
    Convert globsim scaled netCDF to CLASS-CTEM .met file.

    ncdf_globsim: full path to a globsim scaled netCDF (by station)

    met_class: full path to the CLASS-CTEM met file to write.

    station_nr: station_number, as given in the stations .csv file to identify
                the station.

    The columns in CLASS-CTEM MET files are:
    1)  Hour
    2)  Minute
    3)  Day of year
    4)  Year YYYY
    5)  Shortwave Radiation (W/m2)
    6)  Longwave Radiation (W/m2)
    7)  Precip (mm/s)
    8)  Temp.(°C)
    9)  Specific Humidity (Kg/Kg)
    10) Wind Speed (m/s)
    11) Pressure (Pa)

    """

    # columns to export
    columns = ['time', 'SW_sur', 'LW_sur', 'PREC_sur',
               'AIRT_sur', 'SH_sur', 'WSPD_sur',
               'AIRT_pl']

    # output ASCII formatting
    formatters = {"time": "  {:%H %M  %j  %Y}".format,
                  "SW_ERA_Wm2_sur": "{:8.2f}".format,
                  "LW_ERA_Wm2_sur": "{:8.2f}".format,
                  "PREC_ERA_mmsec_sur": "{:13.4E}".format,
                  "AIRT_ERA_C_sur": "{:8.2f}".format,
                  "SH_ERA_kgkg_sur": "{:11.3E}".format,
                  "WSPD_ERA_ms_sur": "{:7.2f}".format,
                  "AIRT_PRESS_Pa_pl": "{:11.2f}".format}

    # get data
    df = globsimScaled2Pandas(ncdf_globsim, station_nr)

    # convert precipitation
    df["PREC_ERA_mmsec_sur"] = df["PREC_ERA_mm_sur"] / 1800.0

    # write FORTRAN formatted ASCII file
    with open(met_class, 'w') as f:
        f.write(' ')
        f.write(df.to_string(columns=columns,
                formatters=formatters,
                header=False, index=False))
    f.close()


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


def globsim_to_geotop(ncd, out_dir, site=None, start=None, end=None):
    """
    @args
    ncd: netcdf dataset
    site: site name or index
    """
    # open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)
    
    logger.debug(f"Reading file {n.filepath}")
    
    # find number of stations
    nstn = len(n['station'][:])

    # get date / time column
    time = nc.num2date(n['time'][:],
                       units=n['time'].units,
                       calendar=n['time'].calendar)

    time = [x.strftime('%d/%m/%Y %H:%M') for x in time]
    time = pd.DataFrame(time)

    # get precip
    PREC = "PREC_sur"
    PREC = n[PREC][:]  * 3600  # Convert from mm/s to mm/hr

    # get wind velocity
    WSPD = "WSPD_sur"
    WSPD = n[WSPD][:]

    # get wind direction
    WDIR = "WDIR_sur"
    WDIR = n[WDIR][:]

    # get windx and windy

    # get RH
    RH = "RH_sur"
    RH = n[RH][:]
    RH *= 0.01  # Convert RH to dimensionless fraction for geotop

    # get air temp
    AIRT = "AIRT_sur"
    AIRT = n[AIRT][:]

    # get dew temp (missing?)

    # get air pressure
    PRESS = "PRESS_pl"
    PRESS = n[PRESS][:]
    PRESS *= 1e-5      # convert to bar for geotop

    # get shortwave solar global (direct / diffuse missing?)
    SW = "SW_sur"
    SW = n[SW][:]

    # get longwave incoming
    LW = "LW_sur"
    LW = n[LW][:]

    # get site names
    NAMES = nc.chartostring(n['station_name'][:])

    # combine data variables into array
    data = np.stack((PREC, WSPD, WDIR, RH, AIRT, PRESS, SW, LW))

    # write output files
    files = []
    for i in range(nstn):
        if (site is None) or (site == i) or (site == NAMES[i]):
            # massage data into the right shape
            out_df = pd.DataFrame(np.transpose(data[:, :, i]))
            out_df = pd.concat([time, out_df], axis=1)
            out_df.columns = ["Date", "IPrec", "WindVelocity", "WindDirection", "RH",
                              "AirTemp", "AirPress", "SWglobal", "LWin"]

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
