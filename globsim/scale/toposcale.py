from globsim.meteorology import emissivity_clear_sky, boltzmann


def emissivity_all_sky(lw_grid, t_grid):
    """ Eq (2) from Fiddes & Gruber """
    e_all_sky = lw_grid / (boltzmann * t_grid ** 4)

    return e_all_sky


def emissivity_cloud(e_cl_grid, e_as_grid):
    """ Estimate of cloud-based component of emissivity """
    delta_e = e_as_grid - e_cl_grid

    return delta_e


def lw_down_toposcale(t_sub, rh_sub, t_sur, rh_sur, lw_sur):
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
        total downwelling LW radiation at surface level
    """
    clear_sky_emissivity_sub = emissivity_clear_sky(rh_sub, t_sub)
    clear_sky_emissivity_sur = emissivity_clear_sky(rh_sur, t_sur)
    all_sky_emissivity_sur = emissivity_all_sky(lw_sur, t_sur)

    delta_e = emissivity_cloud(clear_sky_emissivity_sur, all_sky_emissivity_sur)

    lw_d_sub = (clear_sky_emissivity_sub + delta_e) * boltzmann * t_sub ** 4

    return lw_d_sub
