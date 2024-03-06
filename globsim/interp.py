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
            elevation of data at (time, level) in meters. Must be monotonically decreasing along axis 1
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

        # difference in elevation, level directly above will be >= 0
        if np.all(np.diff(elevation, axis=1) > 0):  # monotonically increasing
            elev_diff = -(elevation - h)
        elif np.all(np.diff(elevation, axis=1) < 0):  # monotonically decreasing 
            elev_diff = elevation - h
        else:
            raise ValueError("Elevation must be monotonically decreasing or increasing")
        
        # vector of level indices that fall directly above station.
        # Apply after ravel() of data.
        va = np.argmin(elev_diff + (elev_diff < 0) * 100000, axis=1)
        # mask for situations where station is below lowest level
        mask = va < (nl - 1)
        va += np.arange(elevation.shape[0]) * elevation.shape[1]

        # Vector level indices that fall directly below station.
        # Apply after ravel() of data.
        vb = va + mask  # +1 when OK, +0 when below lowest level

        return elev_diff, va, vb


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
