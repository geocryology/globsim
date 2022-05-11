from globsim.meteorology import emissivity_clear_sky, boltzmann
from typing import Union
import numpy as np


def emissivity_all_sky(lw: "Union[float, np.ndarray]",
                       t: "Union[float, np.ndarray]") -> "Union[float, np.ndarray]":
    """ Eq (2) from Fiddes & Gruber 
    lw : float
        longwave radiation [W m-2] 
    t : float
        temperature [K] 
    """
    e_all_sky = lw / (boltzmann * np.power(t, 4))

    return e_all_sky


def emissivity_cloud(e_cl: "Union[float, np.ndarray]", 
                     e_as: "Union[float, np.ndarray]") -> "Union[float, np.ndarray]":
    """ Estimate of cloud-based component of emissivity
    e_cl : float
        clear-sky component of emissivity
    e_as : float
        all-sky component of emissivity
    """
    delta_e = e_as - e_cl

    return delta_e


def lw_down_toposcale(t_sub: "Union[float, np.ndarray]", 
                      rh_sub: "Union[float, np.ndarray]",
                      t_sur: "Union[float, np.ndarray]",
                      rh_sur: "Union[float, np.ndarray]",
                      lw_sur: "Union[float, np.ndarray]") -> "Union[float, np.ndarray]":
    """ Eq (3) From Fiddes & Gruber

    Parameters
    ----------
    t_sub : float
        temperature [K] at subgrid level
    rh_sub : float
        relative humidity [1] (unitless) at subgrid level
    t_sur : float
        surface temperature [K] (GRID temperature in Fiddes & Gruber)
    rh_sur : float
        relative humidity [1] (unitless) at surface level (GRID level)
    lw_sur : float
        total downwelling LW radiation [W m-2] at surface level (GRID level)
    """
    clear_sky_emissivity_sub = emissivity_clear_sky(rh_sub, t_sub)
    clear_sky_emissivity_sur = emissivity_clear_sky(rh_sur, t_sur)
    all_sky_emissivity_sur = emissivity_all_sky(lw_sur, t_sur)

    delta_e = emissivity_cloud(clear_sky_emissivity_sur, all_sky_emissivity_sur)

    lw_d_sub = (clear_sky_emissivity_sub + delta_e) * boltzmann * np.power(t_sub, 4)

    return lw_d_sub
