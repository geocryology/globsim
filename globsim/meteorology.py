import numpy as np


boltzmann = 5.67 * 10**(-8)  # J/s/m/K4 Stefan-Boltzmann constant


def satvapp_kPa_fT(T):
    '''
    Saturation water vapour pressure [kPa] following the Tetens formula, Eq 4.2
    in Stull, Practical Meteorology.

    T: Temperature [C]
    '''
    e0 = 0.6113  # [kPa]
    b = 17.2694  # fitting constant
    T1 = 273.15  # [K]
    T2 = 35.86   # [K]
    T += T1
    return e0 * np.exp((b * (T - T1)) / (T - T2))


def vapp_kPa_fTd(Td):
    '''
    Water vapour pressure [hPa] derived from dewpoint temperature [C]. Taken
    from www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    where it is attributed to (Bolton 1980)
    https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2

    Td: Dew point temperature [C]

    '''
    # (Bolton 1980)
    # https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
    return 6.112 * np.exp((17.67 * Td) / (Td + 243.5))

    # https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
    # return 6.112 * np.power(10,(7.5 * Td)/(237.3 + Td))


def spec_hum_kgkg(Td, Pr):
    '''
    Specific humidity [Kg/Kg]. Eq 4.7 in Stull, Practical Meteorology.
    Td: Dewpoint temperature [C]
    Pr:  Air pressure [Pa]
    '''
    E = 0.622  # density of vater vapour / density of dry air
    e = vapp_kPa_fTd(Td) / 10
    P = Pr / 1000.  # from Pa to kPa
    spec_hum = E * e / (P - e * (1 - E))
    return spec_hum


def water_vap_pressure(RH,T):
    '''
    water vapour pressure [unit:1], Eq C9,C10 in Fiddes and Gruber (2014)
    https://doi.org/10.5194/gmd-7-387-2014
    RH: relative humidity (%)
    Tair: air temperature (kelvin)
    '''
    es0 = 6.11  # reference saturation vapour pressure at 0ºC
    T0 = 273.15  # Kelvin
    lv = 2.5 * 1000000  # latent heat of vaporization of water
    Rv = 461.5  # gas constant for water vapour
    es = es0 * np.exp((lv) / Rv * (1 / T0 - 1 / T))
    pv = (RH * es) / 100
    return pv


def emissivity_clear_sky(RH,T):
    '''
    clear sky emissivity, Eq(1) in Fiddes and Gruber (2014)W
    https://doi.org/10.5194/gmd-7-387-2014
   
    Parameters
    ----------
    RH : float or array
        water vapour pressure [%]
    T : float or array
        air temperature (kelvin)
    
    Returns 
    -------
    float or array : clear sky emissivity [1] (unitless)
    '''
    pv = water_vap_pressure(RH, T)
    x1 = 0.43
    x2 = 5.7
    e_clear = 0.23 + x1 * (pv / T)**(1 / x2)
    return e_clear


def pressure_from_elevation(elevation):
    """Convert elevation into air pressure using barometric formula"""
    g = 9.80665    # Gravitational acceleration [m/s2]
    R = 8.31432    # Universal gas constant for air [N·m /(mol·K)]
    M = 0.0289644  # Molar mass of Earth's air [kg/mol]
    P0 = 101325    # Pressure at sea level [Pa]
    T0 = 288.15    # Temperature at sea level [K]
    # http://en.wikipedia.org/wiki/Barometric_formula
    return P0 * np.exp((-g * M * elevation) / (R * T0)) / 100  # [hPa] or [bar]


def LW_downward(RH,T,N):
    '''
    incoming longware radiation [W/m2], Eq(14) in Fiddes and Gruber (2014)
    https://doi.org/10.5194/gmd-7-387-2014
    e_clear: clear sky emissivity
    N: cloud cover
    T: air temperature
    '''
    e_clear = emissivity_clear_sky(RH,T)
    p1 = 6
    p2 = 4
    e_as = 0.979
    lw = e_clear * (1 - N**p1) + (e_as * (N**p2)) * boltzmann * T**4
    return lw


def relhu_approx_lawrence(t: np.ndarray, td: np.ndarray) -> np.ndarray:
    '''
    https://doi.org/10.1175/BAMS-86-2-225
    t : temperature [C]
    td : dewpoint temperature [C]
    returns : relative humidity [%]
    '''
    t = np.atleast_1d(t)
    td = np.atleast_1d(td)

    relhu = 100 - 5 * (t - td)
    relhu = np.where(relhu < 0 , 0, relhu)
    relhu = np.where(relhu > 100 , 100, relhu)
    
    return relhu


def rh_liston(t: np.ndarray, td: np.ndarray) -> np.ndarray:
    """
    t : temperature [C]
    td : dewpoint temperature [C]
    returns : relative humidity [%]
    """
    t = np.atleast_1d(t)
    td = np.atleast_1d(td)
    
    a = 611.21
    b = 17.502
    c = 240.97
    
    e = a * np.exp(b * td / (c + td))
    es = a * np.exp(b * t / (c + t))
    rh = 100 * e / es
    
    return rh
