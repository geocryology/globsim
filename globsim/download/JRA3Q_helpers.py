
import base64
import numpy as np

import xarray as xr
import pygrib

from pathlib import Path

from globsim.download.ERA3Q_dl import GetAccessor
# workflow ->

# make an empty netcdf file with the correct dimensions, unlimited time dimension
# netcdf file is chunked according to 'chunks' argument, or maybe just yearly files

# download a grib file, open it with pygrib, get the data, close the file
# write the data to the netcdf file
# repeat for all grib files for a particular time step

# 


def build_dir_structure(project_dir):
    pass


def download_constant(access, dest_dir):
    url, path, name = url_LL125_surf()
    path = str(Path(dest_dir, path))
    access.dl(url, path, name)


def download_daily_gribs(access, dest_dir, year, month, day):
    for hour in range(0, 24, 6):
        # surface
        url, path, name = url_anl_surf125(year, month, day, hour)
        path = str(Path(dest_dir, path[1:]))
        print(f"Downloading {name} {year}-{month}-{day} {hour}")
        access.dl(url, path, name, absolute_path=True)

        # pressure-level
        for var in ['tmp', 'spfh', 'ugrd', 'vgrd']:
            url, path, name = url_anl_p125(year, month, day, hour, var)
            path = str(Path(dest_dir, path[1:]))
            print(f"Downloading {name} ({var}) {year}-{month}-{day} {hour}")
            access.dl(url, path, name, absolute_path=True)

    # forecast
    for hour in range(0, 24, 1):
        url, path, name = url_fcst_phy2m125(year, month, day, hour)
        path = str(Path(dest_dir, path[1:]))
        print(f"Downloading {name} {year}-{month}-{day} {hour}")
        access.dl(url, path, name, absolute_path=True)


def url_LL125_surf():
    url = r"https://data.diasjp.net/dl/storages/downloadCmd/L0pSQTNRL0NvbnN0L0xMMTI1X3N1cmYuZ3JpYjI="
    p1 = "/JRA3Q/Const/"
    p2 = "LL125_surf.grib2"
    
    return url, p1, p2


def url_anl_surf125(year, month, day, hour):
    """
    1 Total column vertically-integrated water vapour
    2 Potential temperature
    3 Water equivalent of accumulated snow depth (deprecated)
    4 Surface pressure
    5 unknown
    6 unknown
    7 Pressure reduced to MSL
    8 2 metre temperature
    9 unknown
    10 2 metre specific humidity
    11 2 metre relative humidity
    12 10 metre U wind component
    13 10 metre V wind component
    """
    YYYY = "{:04}".format(year)
    MM = "{:02}".format(month)
    DD = "{:02}".format(day)
    HH = "{:02}".format(hour)

    p1 = f"/JRA3Q/Hist/Daily/anl_surf125/{YYYY}{MM}/"
    p2 = f"anl_surf125.{YYYY}{MM}{DD}{HH}"
    
    coded = base64.b64encode((p1 + p2).encode("utf-8")).decode("UTF-8")
    
    url = f"https://data.diasjp.net/dl/storages/downloadCmd/{coded}"

    return url, p1, p2
    

def url_anl_p125(year, month, day, hour, var):
    """
    Parameters
    ----------
    var : str
        tmp: temperature, spfh: specific humidity, ugrd: u wind, vgrd: v wind,
    """
    YYYY = "{:04}".format(year)
    MM = "{:02}".format(month)
    DD = "{:02}".format(day)
    HH = "{:02}".format(hour)

    p1 = f"/JRA3Q/Hist/Daily/anl_p125/{YYYY}{MM}/"
    p2 = f"anl_p125_{var}.{YYYY}{MM}{DD}{HH}"

    coded = base64.b64encode((p1 + p2).encode("utf-8")).decode("UTF-8")
    
    # url = f"https://data.diasjp.net/dl/storages/file/{coded}"
    url = f"https://data.diasjp.net/dl/storages/downloadCmd/{coded}"

    return url, p1, p2


def url_fcst_phy2m125(year, month, day, hour):
    """
    1 unknown
    2 unknown
    3 Latent heat net flux
    4 Instantaneous surface sensible heat flux
    5 Convective precipitation rate
    6 Total precipitation rate
    7 Total snowfall rate water equivalent
    8 Large scale precipitation rate
    9 Instantaneous moisture flux
    10 Momentum flux, u component
    11 Momentum flux, v component
    12 Surface pressure
    13 Downward short-wave radiation flux
    14 Upward short-wave radiation flux
    15 Downward short-wave radiation flux, clear sky
    16 Upward short-wave radiation flux, clear sky
    17 Downward long-wave radiation flux
    18 Upward long-wave radiation flux
    19 Downward long-wave radiation flux, clear sky
    20 unknown
    21 unknown
    22 unknown
    23 unknown
    24 Downward short-wave radiation flux
    25 Upward short-wave radiation flux
    26 Upward short-wave radiation flux, clear sky
    27 Upward long-wave radiation flux
    28 Net long-wave radiation flux, clear sky
    """
    YYYY = "{:04}".format(year)
    MM = "{:02}".format(month)
    DD = "{:02}".format(day)
    HH = "{:02}".format(hour)

    p1 = f"/JRA3Q/Hist/Daily/fcst_phy2m125/{YYYY}{MM}/"
    p2 = f"fcst_phy2m125.{YYYY}{MM}{DD}{HH}"

    coded = base64.b64encode((p1 + p2).encode("utf-8")).decode("UTF-8")

    url = f"https://data.diasjp.net/dl/storages/downloadCmd/{coded}"

    return url, p1, p2

