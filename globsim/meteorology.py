import numpy as np
import pandas as pd


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


def rh_lawrence(t: np.ndarray, td: np.ndarray) -> np.ndarray:
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
    https://doi.org/10.1175/JHM486.1
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


def degree_days(df: pd.DataFrame) -> pd.DataFrame:
    """ 

    df : pd.DataFrame
        Dataframe with a time index and first column temperatures
    """
    daily = df.resample("D").mean()  # ensure daily means
    df['degree_days'] = daily['T']


def forteau_2021_conductivity(rho, T):
    """
    Forteau et al. (2021) thermal conductivity model for snow

    rho : float or array
        Density of snow [kg m-3]
    T : float or array
        Temperature of snow [K]
    """
    # thermal conductivity [W m-1 K-1]
    if np.max(T) < 50:
        print("Warning: temperatures should be in degrees K")

    fitted = {223.0: (+2.564, -0.059, +0.0205),
              248.0: (+2.172, +0.015, +0.0252),
              263.0: (+1.985, +0.073, +0.0336),
              268.0: (+1.883, +0.107, +0.0386),
              273.0: (+1.776, +0.147, +0.0455),}
    
    closest_value = min(list(fitted.keys()), key=lambda x: abs(T - x))
    
    c1, c2, c3 = fitted[closest_value]
    rho_i = 917  # density of ice [kg m-3]
    vfi: float = rho / rho_i
    
    k = c1 * (vfi)**2 + c2 * (vfi) + c3
    
    return k


def swe_accumulation(times, temps, total_precip) -> pd.Series:
    """ 
    Calculate snow water equivalent accumulation over a time period

    Parameters
    ----------
    times : array-like
        Times of observations
    temps : array-like
        Temperatures at each time [C]
    total_precip : array-like
        Precipitation at each time [mm]

    Returns
    -------
    swe : float
        Accumulated snow water equivalent over the time period
    """
    df = pd.DataFrame(data={"T":temps, "P":total_precip}, index=times)
    df['P'] = df['P'].clip(lower=0)  # remove negative precipitation
    df['P'] = df['P'].fillna(0)  # remove NaNs
    precip_fractions = precip_fraction(df.index, df['T'], df['P'])
    df['P_liq'] = df['P'] * precip_fractions['liquid']
    df['P_solid'] = df['P'] * precip_fractions['solid']
    
    return df["P_solid"]


def precip_fraction(times, temps, total_precip, method="static", **kwargs) -> pd.DataFrame:
    """
    Calculate the solid (snow) and liquid (rain) fraction of precipitation over a time period
    """
    df = pd.DataFrame(data={"T":temps, "P":total_precip}, index=times)
    # select and apply appropriate method here
    if method == "static":
        rain_frac = lambda t: np.heaviside(t - 2.2, 1)
    else:
        raise NotImplementedError
    
    df['liquid'] = rain_frac(df['T'])
    df['solid'] = 1 - df['liquid']

    return df[['liquid', 'solid']]


def snowmelt(times, temps) -> pd.DataFrame:
    """
    Calculate potential snow melt
    
    snow melt is calculated according to a constant degree-day factor of 
    3.0 mm oC-1 d-1 (https://core.ac.uk/reader/14923120)
    Returns
    --------
    melt : pd.DataFrame
        Data frame with 1 column of melt [mm]
    """
    df = pd.DataFrame(data = temps, index=times, columns = ["T"])
    df.resample("D").mean()  # ensure daily means
    melt = df['T'] * -3.0
    melt = melt.clip(upper=0)  # remove negative melt

    return pd.DataFrame(data={"melt_mm": melt.values}, index=df.index)


def snow_model(times, temps, total_precip):
    """
    Calculate snow water equivalent accumulation and melt over a time period

    Parameters
    ----------
    times : array-like
        Times of observations
    temps : array-like
        Temperatures at each time [C]
    precip : array-like
        Precipitation at each time [mm]

    Returns
    -------
    swe : float
        Accumulated snow water equivalent over the time period
    """
    df = pd.DataFrame(data={"T":temps, "P":total_precip}, index=times)
    # 

    accum = swe_accumulation(times, temps, total_precip)
    potential_melt = snowmelt(times, temps)
    
    df['accumulation'] = accum
    df['potential_melt'] = potential_melt

    df['snowpack_SWE'] = accumulate(accum.values, potential_melt.values)
    
    df['age_of_snowpack'] = df['snowpack_SWE'].groupby((df['snowpack_SWE'] == 0).cumsum()).cumcount()
    df['snowcover'] = (df['snowpack_SWE']>0).astype('int')

    return df


def thermal_transmissivity(swe, t, alpha: float, beta:float, nosnow=15):
    """ 
    Calculate heat transfer coefficient for snowpack
    Parameters
    ----------
    swe : float
        Snow water equivalent [mm]
    t : float
        Age of snowpack [days]
    alpha : float
        Constant [W m-3 C-1]
    beta : float   
        Constant [days]
    """
    if (max(t) >= beta):
        raise ValueError(f"beta ({beta}) must be greater than the longest snowpack duration ({max(t)} days)")
    
    swe = np.atleast_1d(swe)
    t = np.atleast_1d(t)
    Hs = np.empty_like(swe)
    Hs[swe == 0] = np.nan

    Hs[~np.isnan(Hs)] = 1 / (swe[~np.isnan(Hs)] * alpha * (1 - (t[~np.isnan(Hs)] / beta)))
    Hs[np.isnan(Hs)] = nosnow
    
    return Hs


def accumulate(accumulation, potential_melt, init=0.0) -> np.ndarray:
    out = np.zeros_like(accumulation, dtype='float64')
    for i, (acc, melt) in enumerate(zip(accumulation, potential_melt)):
        if i == 0:
            out[i] = max(0.0, (init + acc + melt))
        else:
            out[i] = max(0.0, (out[i - 1] + acc + melt))
    return out


if __name__ == "__main__":
    import pandas as pd
    w19 = pd.read_csv(r"C:\Users\Nick\Downloads\en_climate_daily_NT_2204101_2019_P1D.csv")
    w20 = pd.read_csv(r"C:\Users\Nick\Downloads\en_climate_daily_NT_2204101_2020_P1D.csv")
    w21 = pd.read_csv(r"C:\Users\Nick\Downloads\en_climate_daily_NT_2204101_2021_P1D.csv")
    w22 = pd.read_csv(r"C:\Users\Nick\Downloads\en_climate_daily_NT_2204101_2022_P1D.csv")
    EC = pd.concat([w19, w20, w21, w22])
    EC['Date/Time'] = pd.to_datetime(EC['Date/Time'])
    EC['Mean Temp (°C)'].interpolate(inplace=True)
    EC["Total Precip (mm)"].interpolate(inplace=True)
    EC.set_index("Date/Time", inplace=True)
    EC = EC[(EC.index.year >= 2020 )| (EC.index.month >= 8)]
    times = EC.index
    temps = EC["Mean Temp (°C)"].values
    total_precip = EC["Total Precip (mm)"].values
    snow = snow_model(times, temps, total_precip)
    htc = HTC_transfer(snow['snowpack_SWE'], snow['age_of_snowpack'], 0.14, 365)


    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('Date')
    ax1.set_ylabel('SWE', color="blue")

    ax1.plot(snow.index, snow['snowpack_SWE'], color="blue")
    ax2.set_ylabel('htc', color="blue")
    ax2.plot(snow.index, 1/htc, color="black")
    fig.tight_layout()
    plt.show()

    ## plotting
    fig, ax = plt.subplots()
    ax.plot(EC['Snow on Grnd (cm)'], label="Env. Can.", color="black")

    rho = 100
    ax.plot(snow['snowpack_SWE'] * 100 / (rho), label=f"model (rho={rho} kg m-3)", color="blue")
    ax.set_ylabel("Snow depth [cm]")
    ax.set_xlabel("Date")
    ax.set_title("Yellowknife")
    plt.legend()
    plt.show()