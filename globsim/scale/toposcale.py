from globsim.meteorology import emissivity_clear_sky, boltzmann
from typing import Union
from pysolar import solar, radiation
import numpy as np
import warnings


def emissivity_all_sky(lw: "Union[float, np.ndarray]",
                       t: "Union[float, np.ndarray]") -> "Union[float, np.ndarray]":
    """ Eq (2) from Fiddes & Gruber
    lw : float
        longwave radiation [W m-2]
    t : float
        temperature [K]
    """
    e_all_sky = lw / (boltzmann * np.power(t, 4, dtype='float64'))

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
                      lw_sur: "Union[float, np.ndarray]",
                      v_d: "Union[float, np.ndarray]"=1.0) -> "Union[float, np.ndarray]":
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
    v_d : float
        Sky view factor [1] (unitless) between 0.0 and 1.0
    """
    if not (0.0 <= np.atleast_1d(v_d)).all() and (np.atleast_1d(v_d) <= 1.0).all():
        raise ValueError(f"Sky view factor must be between 0 and 1 (not {v_d})")

    clear_sky_emissivity_sub = emissivity_clear_sky(rh_sub, t_sub)
    clear_sky_emissivity_sur = emissivity_clear_sky(rh_sur, t_sur)

    all_sky_emissivity_sur = emissivity_all_sky(lw_sur, t_sur)

    delta_e = emissivity_cloud(clear_sky_emissivity_sur, all_sky_emissivity_sur)

    lw_d_sub = (clear_sky_emissivity_sub + delta_e) * boltzmann * np.power(t_sub, 4, dtype='float64')

    return lw_d_sub * v_d


def sw_down_diffuse_toposcale(sw_down_sur, v_d):
    """ Eq (10) from Fiddes & Gruber 
    sw_down_sur : float
        diffuse SW radiation at surface (GRID) level
    v_d : float
        sky view factor at site (SUB)
    """
    return sw_down_sur * v_d


def solar_zenith(lat, lon, time, elevation=None):
    """ 
    lat : float
        latitude
    lon : float
        longitude
    time : datetime
        datetime with time zone
    
    Example
    -------
    solar_zenith(52, -117, datetime.datetime(year=2020, month=5, day=10, hour=18, tzinfo=datetime.timezone.utc))
    """
    solar_elev = solar.get_altitude(lat, lon, time, elevation)
    zenith = 90 - solar_elev
    return zenith


def delta_m(delta_z, zenith_angle):
    """ Eq (4) from Fiddes and Gruber 
    zenith_angle : float
        angle in degrees between vertical and incoming solar radiation
    """
    d_m = delta_z * np.cos(np.radians(zenith_angle))
    return d_m


def broadband_attenuation(sw_d, sw_toa, zenith_angle):
    """ Based on eq(5) in Fiddes and Gruber """
    if zenith_angle > 70:
        warnings.warn(f"Calculated zenith angle of {zenith_angle} is large. Beet-Lambert approximation [m = 1 / cos(zenith)] may not be valid.")
    k = - np.log(sw_d / sw_toa) / (1 / np.cos(np.radians(zenith_angle)))
    return k


def eq6(sw_dir, k, delta_z, zenith_angle, clearness_index=None):
    """ Eq (6) from Fiddes & Gruber """
    # TODO: constrain to reasonable values? github.com/joelfiddes/toposcale/blob/593447cd9de5b7f96a2d9d503d32816fd7bc5285/src/solar.r#L71
    ratio = 1 - np.exp(-k * delta_z * np.cos(np.radians(zenith_angle)))
    ratio = 1 - np.exp(-k * delta_z * np.cos(np.radians(zenith_angle)))
    print(ratio)
    delta_sw = ratio * sw_dir
    
    if clearness_index is None:
        warnings.warn("air-mass corrected clearness index is unknown. Shortwave radiation may be incorrect")
    elif clearness_index < 0.65:
        warnings.warn(f"air-mass corrected clearness index is {clearness_index}. Correction validity is questionable.")
    else:
        pass
    
    return delta_sw


def eq6mod(sw_dir, k, delta_z, zenith_angle, clearness_index=None):
    """ Eq (6) from Fiddes & Gruber """
    # TODO: constrain to reasonable values? github.com/joelfiddes/toposcale/blob/593447cd9de5b7f96a2d9d503d32816fd7bc5285/src/solar.r#L71
    delta_m = 
    ratio = 1 - np.exp(-k * / np.cos(np.radians(zenith_angle)))
    ratio = 1 - np.exp(-k * delta_z * np.cos(np.radians(zenith_angle)))
    print(ratio)
    delta_sw = ratio * sw_dir
    
    if clearness_index is None:
        warnings.warn("air-mass corrected clearness index is unknown. Shortwave radiation may be incorrect")
    elif clearness_index < 0.65:
        warnings.warn(f"air-mass corrected clearness index is {clearness_index}. Correction validity is questionable.")
    else:
        pass
    
    return delta_sw


def clearness_index(sw, sw_toa):
    kt = sw / sw_toa  # clearness index
    mask = np.logical_or(np.isinf(kt), np.less(kt, 0))
    kt[mask] = 0
    kt[kt > 1] = 0.8  # upper limit from topooscale R package
    return kt


def illumination_angle(zenith_angle, solar_azimuth, slope_angle, aspect):
    """ Eq(7) """
    sin = np.sin
    cos = np.cos
    θz = np.radians(zenith_angle)
    φ0 = np.radians(solar_azimuth)
    S = np.radians(slope_angle)
    A = np.radians(aspect)
    i_sub = np.arccos(cos(θz) * cos(S) + sin(θz) * sin(S) * cos(φ0 - A))
    
    return np.degrees(i_sub)


def sw_partition(sw, sw_toa):
    """
    return diffuse, direct"""
    kt = clearness_index(sw=sw, sw_toa=sw_toa)

    kd = 0.952 - 1.041 * np.exp(-np.exp(2.3 - 4.702 * kt))
    
    return kd, 1 - kd


def sw_toa(latitude, longitude, elevation, date):
    """ top-of-atmosphere solar shortwave

    Returns
    -------
    _type_
        _description_
    
    Examples
    --------
    import datetime
    time = [datetime.datetime(year=2020, month=5, day=10, hour=h, tzinfo=datetime.timezone.utc) for h in range(0,24)]
    toa = sw_toa((52,)*24, (-117,)*24, (800,)*24, time)
    toa = sw_toa(45, -110, 400, time[0])
    https://www.sciencedirect.com/topics/earth-and-planetary-sciences/insolation
    https://en.wikipedia.org/wiki/Solar_irradiance
    https://www.itacanet.org/the-sun-as-a-source-of-energy/part-2-solar-energy-reaching-the-earths-surface/#2.1.-The-Solar-Constant

    Seems wrong... TSI = 1360, 'apparent_ET_flux' only goes up to 1225 or so.
    Consider switching to Q = S_0 (1 + 0.034 cost(2 pi n / 365.25)) where n = julian day S0 = TSI
    """
    
    date = np.atleast_1d(date)
    args = list(zip(np.atleast_1d(latitude),
               np.atleast_1d(longitude),
               date,
               np.atleast_1d(elevation)))
            
    # TSI = np.array([radiation.get_apparent_extraterrestrial_flux(t.timetuple().tm_yday) for t in date])
    # warnings.warn("top-of-atmosphere irradiance values are supect. Seems to low")
    TSI = 1366  # better fit with ERA5 tisr
    alt = np.array([solar.get_altitude(x[0], x[1], x[2], x[3]) for x in args])

    zenith = 90 - alt
    zenith[zenith > 90] = 90

    toa = TSI * np.cos(np.radians(zenith))
    toa[toa < 1e-2] = 0

    return toa
    

def sw_direct_sub(sw_direct, incidence_sub, incidence_grid, shadow_mask=None):
    """ Topographically-corrected SW direct Eq (9) Fiddes & Gruber """
    
    if shadow_mask is None:
        shadow_mask = np.ones_like(sw_direct, dtype=bool)

    sw_dir_sub = sw_direct * (incidence_sub / incidence_grid)
    sw_dir_sub[shadow_mask] = 0

    return sw_dir_sub
    

if __name__ == "__main__":
    # worked example - longwave
    t = 273 + 0  # from reanalysis
    t0 = 273 + 15  # from reanalysis
    rh = 0.2  # from reanalysis
    rh0 = 0.6  # from reanalysis
    shortwave = 500  # from reanalysis
    svf = 1  # from horizon photo / provided in siteslist

    lw_down_toposcale(t_sub=t, rh_sub=rh, t_sur=t0, rh_sur=rh0, lw_sur=shortwave, v_d=svf)

    # worked example - shorty  wave
    import datetime
    lat = 62  # from siteslist
    lon = -114  # from siteslist
    elevation = 1420  # from siteslist
    grid_elev = 100 #  derive from geopotential E [m] = G [m2 s-2] / 9.80665 [m s-2]
    
    time = datetime.datetime(year=2022, month=5, day=9, hour=12, tzinfo=datetime.timezone(datetime.timedelta(seconds=60*60*-6)))

    zenith = solar_zenith(lat=lat, lon=lon, time=time, elevation=elevation)
    dm = delta_m(delta_z=grid_elev, zenith_angle=zenith)
    toa = sw_toa(latitude=lat, longitude=lon, elevation=grid_elev, date=time)  # shouldn't have to provide elevation here.... 
    kt = clearness_index(shortwave, toa)
    diffuse, direct = sw_partition(shortwave, toa)
    
    k = broadband_attenuation(sw_d=shortwave * direct, sw_toa=toa, zenith_angle=zenith)
    
    delta_sw = eq6(sw_dir=shortwave * direct, k=k, delta_z=elevation-grid_elev, zenith_angle=zenith, clearness_index=kt)
    

""" 
lat = 51.3
lon = -117
when = [datetime.datetime(year=2022, month=5, day=9, hour=h, tzinfo=datetime.timezone(datetime.timedelta(seconds=60*60*-6)) ) for h in range(0,24)]
az = [solar.get_azimuth(lat, lon, w ) for w in when]
alt = [solar.get_altitude(lat, lon, w ) for w in when]
"""
