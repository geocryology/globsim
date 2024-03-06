import numpy as np
import netCDF4 as nc
from globsim.interp import interpolate_level_data


def delta_T_c(T_sa: np.ndarray, airT_pl: np.ndarray, elevation: np.ndarray, h_sur: np.ndarray):
    """
    Calculate land-surface influences on surface air temperature at elevation of coarse scale
    
    Parameters
    ----------
    T_sa : np.ndarray
        Surface air temperature [K] with dimensions (time, station)
    airT_pl : np.ndarray
        Air temperature at pressure levels [K] with dimensions (time, level, station)
    elevation : np.ndarray
        Elevation at pressure levels [m] with dimensions (time, level, station)
    h_sur : np.ndarray
        Station or fine-scale topography elevation [m] with dimensions (station,)

    Returns
    ------- 
    Delta_T_c : np.ndarray
        Land-surface influences on surface air temperature at elevation of coarse scale (K) with dimensions (time, station)
    """
    T_pl_c = interpolate_level_data(airT_pl, elevation, h_sur)

    Delta_T_c = T_sa - T_pl_c  # eq. (4)  

    return Delta_T_c

def add_var_delta_T(group: nc.Dataset):
    """ Add REDCAPP delta T to the netCDF file"""
    vn  = 'AIRT_redcapp_DeltaT'  # variable name
    var = group.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Land-surface influences on surface air temperature at elevation of coarse scale'
    var.comment = 'Cao et al. (2017) 10.5194/gmd-10-2905-2017'
    var.units     = 'degree_C'
    
    return group[vn]
