import numpy as np
import logging
from typing import Union
from numbers import Number

logger = logging.getLogger(__name__)


def ele_interpolate(elevation: np.ndarray, h: float, nl: int):
        """ 
        Parameters
        ----------
        elevation : np.ndarray
            elevation of data at (time, level) in meters. Must be monotonically increasing/decreasing along axis 1
        h : float
            elevation of station
        nl : int
            number of levels in the dataset

        Returns
        -------
        elev_diff : np.ndarray
            difference in elevation between station and pressure level, level directly above will be >= 0
        va : np.ndarray
            array of indices corresponding to nearest level above station
        vb : np.ndarray 
            array of indices corresponding to nearest level below station
        """
        # difference in elevation, level directly station above will be >= 0
        elev_diff = elevation - h
        i_ravel = np.arange(elevation.shape[0]) * elevation.shape[1]  # indices for first level in each time 
        if np.all(np.diff(elevation, axis=1) > 0):  # elevations monotonically increasing
            inverted = True
        elif np.all(np.diff(elevation, axis=1) < 0):  # elevations monotonically decreasing 
            inverted = False
        else:
            raise ValueError("Elevation must be monotonically decreasing or increasing")
        if inverted:
            va = np.argmin(elev_diff + (elev_diff < 0) * 1e6, axis=1) # level indices that are directly above station.
            mask = (va != 0)  # station is above lowest level
            va += i_ravel  # will be added to raveled data
            vb = va - mask  # next-lowest station index when va != 0
        else:
            # vector of level indices that fall directly above station.
            # Apply after ravel() of data.
            va = np.argmin(elev_diff + (elev_diff < 0) * 100000, axis=1)
            # mask for situations where station is below lowest level
            mask = va < (nl - 1)
            va += i_ravel

            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask  # +1 when OK, +0 when below lowest level
        
        return elev_diff, va, vb


def extrapolate_below_grid(elevation: np.ndarray, data:np.ndarray, h: float):
        """ 
        Parameters
        ----------
        elevation : np.ndarray
            elevation of data at (time, level) in meters. Must be monotonically increasing/decreasing along axis 1
        data : np.ndarray
            data to interpolate.
        h : float
            elevation of station

        Returns
        -------
        epol : np.ndarray
            extrapolated data for station below lowest level (masked where station is above lowest level)
        """
        nl = elevation.shape[1]
        i_ravel = np.arange(elevation.shape[0]) * nl  # indices for first level in each time 

        if np.all(np.diff(elevation, axis=1) > 0): 
            inverted = True   # elevations monotonically increasing
        elif np.all(np.diff(elevation, axis=1) < 0):  
            inverted = False  # elevations monotonically decreasing 
        else:
            raise ValueError("Elevation must be monotonically decreasing or increasing")

        if inverted:
            delta_s = elevation[:, 0] - h
            delta_L = elevation[:, 1] - elevation[:, 0]
            i_lowest_diff = i_ravel
            i_lowest_data = i_ravel
            R = -(delta_s / delta_L)
        
        else:
            delta_s = elevation[:, -1] - h
            delta_L = elevation[:, -2] - elevation[:, -1]
            i_lowest_diff = i_ravel + nl - 2
            i_lowest_data = i_ravel + nl - 1
            R = delta_s / delta_L

        below_lowest = np.where(np.min(elevation, axis=1) > h)[0]  # indices of times where station is below lowest level
        delta_V = np.diff(data)[i_lowest_diff]  # difference between levels at each time [nt,]
        epol = np.ma.MaskedArray(data=data[i_lowest_data] + R * delta_V, mask=True)  # extrapolated values
        epol.mask[below_lowest] = False  

        return epol


def calculate_weights(elev_diff, va, vb) -> tuple:
    """
    Parameters
    ----------
    elev_diff : np.ndarray
        difference in elevation between station and , level directly above will be >= 0
    va : np.ndarray 
        array of indices corresponding to nearest level above station
    vb : np.ndarray
        array of indices corresponding to nearest level below station
    
    Returns
    -------
    wa : np.ndarray
        weight for level above station
    wb : np.ndarray
        weight for level below station
    """
    wa = np.absolute(elev_diff.ravel()[vb])
    wb = np.absolute(elev_diff.ravel()[va])

    # guard if wa + wb == 0
    both_zero = wa + wb == 0
    wa[both_zero] = 1
    wb[both_zero] = 1

    wt = wa + wb
    wa /= wt  # Apply after ravel() of data.
    wb /= wt  # Apply after ravel() of data.

    return wa, wb


def interpolate_level_data(data: np.ndarray, elevation: np.ndarray, h: Union[float, int, np.ndarray]):
    """ Linear interpolation of data to station elevation
    
    Parameters
    ----------
    data : np.ndarray
        data to interpolate. Must be (time, level) or (time, level, station)
    elevation : np.ndarray
        elevation of data in meters. Same dimensions as data. Must be monotonically decreasing along axis 1
    h : float
        elevation of station

    Returns
    ------- 
    interp : np.ndarray
        interpolated data

    Examples
    --------
    
    """
    nl = data.shape[1]

    if (data.ndim == 3) and isinstance(h, np.ndarray):
        interp = np.zeros(shape=data.shape[slice(0,None,2)])  # type: np.ndarray
        for i, h_i in enumerate(h):
            interp[:, i] = interpolate_level_data(data[:, :, i], elevation[:, :, i], h_i)
    
    elif (data.ndim == 2) and isinstance(h, Number):
        elev_diff, va, vb = ele_interpolate(elevation, h, nl)
        wa, wb = calculate_weights(elev_diff, va, vb)
        interp = data.ravel()[va] * wa + data.ravel()[vb] * wb
    
    else:
        raise ValueError(f"Data (dim: {data.shape}) and elevation (dim: {data.shape}) must have same dimensions, and h ({type(h)}) must be same type as  a float or np.ndarray")
    
    return interp
