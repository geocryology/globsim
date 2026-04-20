"""
Export functions to convert globsim scaled netCDF to CLASSIC meteorological
forcing netCDF files.

CLASSIC expects one netCDF file per meteorological variable (e.g. shortwave,
longwave, precipitation, temperature, humidity, wind, pressure).  An export
profile (TOML) maps each output variable to a globsim scaled input variable
and carries metadata (units, long_name, output filename, scale_factor, offset).
"""
from os import path, makedirs
from pathlib import Path
from globsim import __version__

import logging
import datetime
import netCDF4 as nc
import numpy as np
import pkg_resources
import shutil
import tomlkit

logger = logging.getLogger("globsim.convert")


# ---------------------------------------------------------------------------
# Helper – create a single CLASSIC-style netCDF file
# ---------------------------------------------------------------------------

def _fractional_day_to_cftime(value):
    ymd = int(value)
    frac = value - ymd
    
    year = ymd // 10000
    month = (ymd % 10000) // 100
    day = ymd % 100
    
    base_date = datetime.datetime(year, month, day)
    
    #Calculate total seconds and round to the nearest minut
    total_seconds = round(frac * 86400 / 60) * 60
    
    return base_date + datetime.timedelta(seconds=total_seconds)


def _format_cftime_fractional_day(t) ->float:
    """
    Format a cftime datetime as YYYYMMDD.FFFF
    where FFFF is the fraction of the day (7 decimals).
    """
    sec = (
        t.hour * 3600 +
        t.minute * 60 +
        t.second +
        t.microsecond / 1e6
    )
    
    frac = round(sec / 86400.0, 7)
    
    return float(f"{t.year:04d}{t.month:02d}{t.day:02d}") + frac


def _create_classic_nc(filepath, time_values,
                       lat, lon, var_name, var_data, var_units, var_long_name, title=None, extra_timestep=True):
    """Write a single CLASSIC meteorological forcing netCDF file.

    Parameters
    ----------
    filepath : str
        Full path for the output netCDF file.
    time_values : array-like
        Numeric time values (in *time_units* / *time_calendar*).
    lat, lon : float
        Station latitude / longitude.
    var_name : str
        NetCDF variable name (e.g. ``'Fss'``).
    var_data : np.ndarray
        1-D array of values with shape ``(ntime,)``.
    var_units : str
        CF units string for the variable.
    var_long_name : str
        Human-readable description of the variable.
    extra_timestep : bool
        Whether to add an extra timestep at the end of the file with the same value as the last real timestep. 
        This is a workaround for CLASSIC expecting one more timestep than the number of time values in the input file.
    """
    ds = nc.Dataset(filepath, "w", format="NETCDF3_CLASSIC")

    if extra_timestep:
        dt_seconds = round((time_values[-1] - time_values[-2]) * 24 * 3600, 0)
        dt_delta = datetime.timedelta(seconds=dt_seconds)
        
        current_last_date = _fractional_day_to_cftime(time_values[-1])
        
        target_date = current_last_date + datetime.timedelta(days=1)
        target_date = target_date.replace(hour=0, minute=0, second=0, microsecond=0)

        while current_last_date < target_date:
            current_last_date += dt_delta
            next_time = _format_cftime_fractional_day(current_last_date)
            
            time_values = np.append(time_values, next_time)
            var_data = np.append(var_data, var_data[-1])
    
    # dimensions
    nt = len(time_values)
    ds.createDimension("time", nt)  
    ds.createDimension("lat", 1)
    ds.createDimension("lon", 1)

    # coordinate variables
    t_var = ds.createVariable("time", "f8", ("time",))
    t_var.units = "day as YYYYMMDD.FFFF"
    t_var.calendar = "proleptic_gregorian"
    t_var.axis = "T"
    t_var.long_name = "time"
    t_var[:] = time_values

    lat_var = ds.createVariable("lat", "f8", ("lat",))
    lat_var.units = "degrees_north"
    lat_var.long_name = "Latitude"
    lat_var.standard_name = "latitude"
    lat_var.axis = "Y"
    lat_var[:] = [lat]

    lon_var = ds.createVariable("lon", "f8", ("lon",))
    lon_var.units = "degrees_east"
    lon_var.standard_name = "longitude"
    lon_var.long_name = "Longitude"
    lon_var.axis = "X"
    lon_var[:] = [lon]

    # data variable
    v = ds.createVariable(var_name, "f4", ("time", "lat", "lon"),
                          fill_value=1.0e38)
    v.units = var_units
    v.long_name = var_long_name
    v[:, 0, 0] = var_data

    # global attributes
    ds.title = title if title is not None else var_long_name
    ds.units = var_units
    ds.grid_type = "gaussian"
    ds._FillValue = v._FillValue
    ds.Conventions = "CF-1.6"
    ds.history = f"Created by globsim (v{__version__}) CLASSIC export"

    ds.close()
    logger.info(f"Wrote {filepath}")


# ---------------------------------------------------------------------------
# Main export routine
# ---------------------------------------------------------------------------

def globsim_to_classic(ncd, out_dir, site=None, export_profile=None):
    """Export a scaled globsim file to CLASSIC meteorological forcing netCDF
    files (one file per variable per site).

    Parameters
    ----------
    ncd : str or netCDF4.Dataset
        Path to scaled globsim netCDF file, or an already-opened Dataset.
    out_dir : str
        Directory where the output netCDF files will be written.
    site : str, int, list, or None
        Site name(s) or index to export.  May be a single value or a list
        of site names (as returned by CLI ``nargs="*"``).  ``None`` exports
        all sites.
    export_profile : str or None
        Path to a TOML export-profile file.  When *None* a default profile
        is created / used at ``~/.globsim/classic_profile.toml``.

    Returns
    -------
    list[str]
        Paths of all created netCDF files.
    """

    # --- load export profile ------------------------------------------------
    if export_profile is None:
        export_profile = Path("~/.globsim/classic_profile.toml").expanduser()

        if not Path(export_profile).is_file():
            default = pkg_resources.resource_filename(
                "globsim", "data/classic_profile_default.toml")

            if not export_profile.parent.is_dir():
                makedirs(export_profile.parent)

            shutil.copy(default, export_profile)
            logger.warning(
                f"Created default CLASSIC export profile: {export_profile}")

    with open(export_profile) as p:
        profile = tomlkit.loads(p.read())
        logger.info(f"Loaded CLASSIC export profile from {export_profile}")

    # --- open scaled netCDF -------------------------------------------------
    if isinstance(ncd, str):
        ncd = nc.Dataset(ncd)

    # time metadata
    time_values = ncd["time"][:]
    time_units = ncd["time"].units
    time_calendar = ncd["time"].calendar
    converted_time = nc.num2date(time_values, units=time_units, calendar=time_calendar)
    formatted_time = [_format_cftime_fractional_day(t) for t in converted_time]

    # station names
    try:
        raw = ncd["station_name"][:]
        try:
            names = nc.chartostring(raw)
        except (ValueError, TypeError):
            names = np.array(raw).astype("str")
    except KeyError:
        raw = ncd["station"][:]
        try:
            names = nc.chartostring(raw)
        except (ValueError, TypeError):
            names = np.array(raw).astype("str")

    nstn = len(names)

    # station coordinates
    try:
        lats = ncd["latitude"][:]
        lons = ncd["longitude"][:]
    except KeyError:
        lats = np.zeros(nstn)
        lons = np.zeros(nstn)

    # ensure output directory exists
    if not path.isdir(out_dir):
        makedirs(out_dir)

    # --- iterate sites ------------------------------------------------------
    files = []
    for i in range(nstn):
        if site is not None:
            if isinstance(site, list):
                if (i not in site) and (names[i] not in site):
                    continue
            elif (site != i) and (site != names[i]):
                continue

        st_name = names[i]
        site_dir = path.join(out_dir, str(st_name))
        if not path.isdir(site_dir):
            makedirs(site_dir)

        lat = float(lats[i])
        lon = float(lons[i])

        # --- iterate variables from profile ---------------------------------
        for out_var, cfg in profile.items():
            input_name = cfg.get("input")
            scale_factor = cfg.get("scale_factor", 1)
            offset = cfg.get("offset", 0)
            output_file = cfg.get("output_file", f"metVar_{out_var}.nc")
            units = cfg.get("units", "")
            long_name = cfg.get("long_name", out_var)
            title = cfg.get("title", None)

            try:
                data = ncd[input_name][:, i] * scale_factor + offset
            except (IndexError, KeyError):
                logger.warning(
                    f"Variable '{input_name}' not found in scaled file – "
                    f"skipping {out_var} for site {st_name}")
                continue

            filepath = path.join(site_dir, output_file)
            _create_classic_nc(
                filepath=filepath,
                time_values=formatted_time,
                lat=lat,
                lon=lon,
                var_name=out_var,
                var_data=data,
                var_units=units,
                var_long_name=long_name,
                title=title
            )
            files.append(filepath)

    return files
