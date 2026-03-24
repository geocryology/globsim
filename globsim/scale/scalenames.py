import enum


class ScaleNames(enum.Enum):
    """Canonical variable names for scaling."""
    time = "time"
    pressure = "pressure"
    temperature = "temperature"
    dewpoint = "dewpoint"
    u_wind = "u_wind"
    v_wind = "v_wind"
    geopotential = "geopotential"
    elevation = "elevation"
    ozone = "ozone"
    water_vapour = "water_vapour"
    rh = "relative_humidity"
    specific_humidity = "specific_humidity"
    sw_down_flux = "shortwave_down_flux"
    sw_down_accumulated = "shortwave_down_accumulated"
    lw_down_flux = "longwave_down_flux"
    lw_down_accumulated = "longwave_down_accumulated"
    precipitation_rate = "precipitation_rate"  # instantaneous precipitation rate
    precipitation_total = "total_precipitation"  # cumulative precipitation from start of forecast

    
