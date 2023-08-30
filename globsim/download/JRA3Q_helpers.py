
import base64
import logging
import numpy as np
import pygrib

from pathlib import Path
from typing import Optional

logger = logging.getLogger("globsim.download")

# workflow ->

# make an empty netcdf file with the correct dimensions, unlimited time dimension
# netcdf file is chunked according to 'chunks' argument, or maybe just yearly files

# download a grib file, open it with pygrib, get the data, close the file
# write the data to the netcdf file
# repeat for all grib files for a particular time step

# 


def download_constant(access, dest_dir) -> Path:
    url, path, name = url_LL125_surf()
    dirpath = Path(dest_dir, path[1:])
    fname = Path(dirpath, name)

    if Path(dirpath, name).exists():
        logger.warning(f"Already exists: {fname}")
    else:
        logger.info(f"Downloading {name}")
        access.dl(url, str(dirpath), name, absolute_path=True)

    return fname


def download_daily_gribs(access, dest_dir, year, month, day) -> "tuple[list[Path], list[Path], list[Path]]":
    sa, pl, sf = [], [], []

    for hour in range(0, 24, 6):
        # surface
        url, path, name = url_anl_surf125(year, month, day, hour)
        dirpath = Path(dest_dir, path[1:])
        fname = Path(dirpath, name)
        sa.append(fname)
        if Path(dirpath, name).exists():
            logger.debug(f"Already exists: {fname}")
        else:
            logger.info(f"Downloading {name} {year}-{month}-{day} {hour}")
            access.dl(url, str(dirpath), name, absolute_path=True)

        # pressure-level
        for var in ['tmp', 'spfh', 'rh', 'ugrd', 'vgrd', 'hgt']:
            url, path, name = url_anl_p125(year, month, day, hour, var)
            dirpath = Path(dest_dir, path[1:])
            fname = Path(dirpath, name)
            pl.append(fname)
            if Path(dirpath, name).exists():
                logger.debug(f"Already exists: {fname}")
            else:
                logger.info(f"Downloading {name} ({var}) {year}-{month}-{day} {hour}")
                access.dl(url, str(dirpath), name, absolute_path=True)

    # forecast
    for hour in range(0, 24, 1):
        url, path, name = url_fcst_phy2m125(year, month, day, hour)
        dirpath = Path(dest_dir, path[1:])
        fname = Path(dirpath, name)
        sf.append(fname)
        if Path(dirpath, name).exists():
            logger.debug(f"Already exists: {fname}")
        else:
            logger.info(f"Downloading {name} {year}-{month}-{day} {hour}")
            access.dl(url, str(dirpath), name, absolute_path=True)

    return sa, sf, pl


def download_monthly_gribs(access, dest_dir, year, month):
    if month in [1, 3, 5, 7, 8, 10, 12]:
        days = 31
    elif month in [4, 6, 9, 11]:
        days = 30
    elif month == 2:
        if year % 4 == 0:
            days = 29
        else:
            days = 28
    else:
        raise ValueError(f"Invalid month: {month}")
    
    for day in range(1, days + 1):
        try:
            sa, sf, pl = download_daily_gribs(access, dest_dir, year, month, day)
        except Exception as e:
            print(f"Error downloading {year}-{month}-{day}: {e}")


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
        tmp: temperature, spfh: specific humidity,  rh: relative humidity, ugrd: u wind, vgrd: v wind, hgt: geopotential
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


class GribSubsetter:
    DEFAULT_LEV_HPA = [0.1, 0.3, 1.0, 3.0, 7.0,
                       1, 2, 3, 5, 7, 10, 20, 30,
                       40, 50, 60, 70, 85, 100,
                       125, 150, 175, 200, 225,
                       250, 300, 350, 400, 450,
                       500, 550, 600, 650, 700,
                       750, 775, 800, 825, 850,
                       875, 900, 925, 950, 975,
                       1000]

    def __init__(self, lon_min, lon_max, lat_min, lat_max,
                 levels:Optional[list] = None):
        self.lon_min = np.floor((lon_min % 360) / 1.25) * 1.25
        self.lon_max = np.ceil((lon_max % 360) / 1.25) * 1.25
        self.lat_min = np.floor((lat_min) / 1.25) * 1.25
        self.lat_max = np.ceil((lat_max) / 1.25) * 1.25
        self.levels = self.DEFAULT_LEV_HPA if levels is None else levels

    @property
    def lev_min(self):
        return min(self.levels)

    @property
    def lev_max(self):
        return max(self.levels)

    @property
    def lats(self):
        return np.arange(self.lat_min, self.lat_max + 1.25, 1.25)

    @property
    def lons(self):
        return np.arange(self.lon_min, self.lon_max + 1.25, 1.25)

    def subset(self, record) -> tuple:
        vals, lats, lons = record.data(lat1=self.lat_min, lat2=self.lat_max, lon1=self.lon_min, lon2=self.lon_max)
        return vals, lats, lons

    def xsubset(self, var) -> tuple:
        if 'isobaricInhPa' not in var.dims:
            ss = var.sel(latitude=slice(self.lat_max, self.lat_min),  # indexed from N to S
                         longitude=slice(self.lon_min % 360, self.lon_max % 360))
        else:
            ss = var.sel(latitude=slice(self.lat_max, self.lat_min),  # indexed from N to S
                         longitude=slice(self.lon_min % 360, self.lon_max % 360),
                         isobaricInhPa=slice(self.lev_max, self.lev_min))  # levels positive down

        vals, lats, lons = ss.values, ss.latitude.values, ss.longitude.values
        return vals, lats, lons


