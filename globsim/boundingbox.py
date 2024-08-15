import numpy as np
import netCDF4 as nc
from typing import Union


class BoundingBox:

    def __str__(self):
        return f"BoundingBox({self.xmin},{self.xmax},{self.ymin},{self.ymax})"

    def __repr__(self):
        return f"BoundingBox({self.xmin},{self.xmax},{self.ymin},{self.ymax})"

    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def to_numpy(self) -> "np.ndarray":
        return np.array([self.xmin, self.xmax, self.ymin, self.ymax])

    def widen(self, amount:float, about: str = 'centre') -> None:
        if about.lower() == 'centre':
            self.xmin -= amount / 2
            self.xmax += amount / 2

        elif about.lower() == 'left':
            self.xmax += amount

        elif about.lower() == 'right':
            self.xmin -= amount

    def heighten(self, amount:float, about:str = 'centre') -> None:
        if about.lower() == 'centre':
            self.ymin -= amount / 2
            self.ymax += amount / 2

        elif about.lower() == 'left':
            self.ymax += amount

        elif about.lower() == 'right':
            self.ymin -= amount

    def contains_bbox(self, bbox: "BoundingBox") -> bool:
        contains = ((self.xmin % 360 < bbox.xmin % 360) and
                    (self.xmax % 360 > bbox.xmax % 360) and
                    (self.ymin % 360 < bbox.ymin % 360) and
                    (self.ymax % 360 > bbox.ymax % 360))
        return contains

    def within_bbox(self, bbox: "BoundingBox") -> bool:
        return bbox.contains_bbox(self)


def stations_bbox(stations) -> BoundingBox:
    bbox = BoundingBox(xmin=stations['longitude_dd'].describe()["min"],
                       xmax=stations['longitude_dd'].describe()["max"],
                       ymin=stations['latitude_dd'].describe()["min"],
                       ymax=stations['latitude_dd'].describe()["max"])
    return bbox


def netcdf_bbox(ncf: "Union[str, nc.Dataset]") -> BoundingBox:

    if isinstance(ncf, str):
        ncf = nc.Dataset(ncf)
    if isinstance(ncf, nc.Dataset):
        x,y = _nc_get_lat_lon(ncf)
        xval = x[:]
        yval = y[:]
    
    else:  # XARRAY
        x,y = _xr_get_lat_lon(ncf)
        xval = x.values
        yval = y.values

    bbx = BoundingBox(np.min(xval), np.max(xval), np.min(yval), np.max(yval))

    return bbx


def _xr_get_lat_lon(xrf):
    # TODO: implement this
    return xrf['longitude'], xrf['latitude']

    x = xrf.filter_by_attrs(axis='X')
    if not x:
        x = xrf.filter_by_attrs(standard_name='longitude')
    if not x:
        x = xrf.filter_by_attrs(long_name='longitude')
    if not x:
        x = xrf['longitude'] if 'longitude' in xrf.coords else None
    if not x:
        raise KeyError("Could not find x-coordinate")

    y = xrf.filter_by_attrs(axis='Y')
    if not y:
        y = xrf.filter_by_attrs(standard_name='latitude')
    if not y:
        y = xrf.filter_by_attrs(long_name='latitude')
    if not y:
        y = xrf['latitude'] if 'latitude' in xrf.coords else None
    if not y:
        raise KeyError("Could not find y-coordinate")

    return x, y

def _nc_get_lat_lon(ncf):
    x = ncf.get_variables_by_attributes(axis='X')
    if not x:
        x = ncf.get_variables_by_attributes(standard_name='longitude')
    if not x:
        x = ncf.get_variables_by_attributes(long_name='longitude')
    if not x:
        x = ncf['longitude'] if 'longitude' in ncf.variables else None
    if not x:
        raise KeyError("Could not find x-coordinate")

    y = ncf.get_variables_by_attributes(axis='Y')
    if not y:
        y = ncf.get_variables_by_attributes(standard_name='latitude')
    if not y:
        y = ncf.get_variables_by_attributes(long_name='latitude')
    if not y:
        y = ncf['latitude'] if 'latitude' in ncf.variables else None
    if not y:
        raise KeyError("Could not find y-coordinate")
    
    return x,y
"""

from collections import namedtuple
import ESMF
import numpy as np


BoundingBox = namedtuple('BoundingBox', ['xmin', 'xmax', 'ymin', 'ymax'])


def create_stations_bbox(stations) -> BoundingBox:
    # get max/min of points lat/lon from self.stations

    stations_bbox = BoundingBox(xmin=stations['longitude_dd'].describe()["min"],
                                xmax=stations['longitude_dd'].describe()["max"],
                                ymin=stations['latitude_dd'].describe()["min"],
                                ymax=stations['latitude_dd'].describe()["max"])
    # add generous buffer
    buffer = 2.0  # assume degrees here
    return BoundingBox(stations_bbox.xmin - buffer,
                       stations_bbox.xmax + buffer,
                       stations_bbox.ymin - buffer,
                       stations_bbox.ymax + buffer)
"""

"""
def grid_from_bbox(latitudes: np.ndarray,
                   longitudes: np.ndarray,
                   bbox: BoundingBox) -> ESMF.Grid:

    valid_lat = latitudes[np.where((latitudes >= bbox.ymin) & (latitudes <= bbox.ymax))[0]]
    valid_lon = longitudes[np.where((longitudes >= bbox.xmin) & (longitudes <= bbox.xmax))[0]]

    grid = ESMF.Grid(max_index=np.array([len(valid_lon), len(valid_lat)]))
    grid.coords[0][0] = np.repeat(valid_lon[np.newaxis, :], len(valid_lat), axis=1).ravel().reshape(len(valid_lon),len(valid_lat))
    grid.coords[0][1] = np.repeat(valid_lat[np.newaxis, :], len(valid_lon), axis=0)
    
    return grid
"""
