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

import logging
import netCDF4 as nc
import numpy as np
import pkg_resources
import shutil
import tomlkit

logger = logging.getLogger("globsim.convert")


# ---------------------------------------------------------------------------
# Helper – create a single CLASSIC-style netCDF file
# ---------------------------------------------------------------------------

def _create_classic_nc(filepath, time_values, time_units, time_calendar,
                       lat, lon, var_name, var_data, var_units, var_long_name):
    """Write a single CLASSIC meteorological forcing netCDF file.

    Parameters
    ----------
    filepath : str
        Full path for the output netCDF file.
    time_values : array-like
        Numeric time values (in *time_units* / *time_calendar*).
    time_units : str
        CF-compliant time units string, e.g. ``'minutes since 1970-01-01 …'``.
    time_calendar : str
        CF calendar name, e.g. ``'standard'``.
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
    """
    ds = nc.Dataset(filepath, "w", format="NETCDF4_CLASSIC")

    # dimensions
    ds.createDimension("time", None)  # unlimited
    ds.createDimension("lat", 1)
    ds.createDimension("lon", 1)

    # coordinate variables
    t_var = ds.createVariable("time", "f8", ("time",))
    t_var.units = time_units
    t_var.calendar = time_calendar
    t_var.axis = "T"
    t_var.long_name = "time"
    t_var[:] = time_values

    lat_var = ds.createVariable("lat", "f8", ("lat",))
    lat_var.units = "degrees_north"
    lat_var.long_name = "latitude"
    lat_var.axis = "Y"
    lat_var[:] = [lat]

    lon_var = ds.createVariable("lon", "f8", ("lon",))
    lon_var.units = "degrees_east"
    lon_var.long_name = "longitude"
    lon_var.axis = "X"
    lon_var[:] = [lon]

    # data variable
    v = ds.createVariable(var_name, "f4", ("time", "lat", "lon"),
                          fill_value=1.0e20)
    v.units = var_units
    v.long_name = var_long_name
    v[:, 0, 0] = var_data

    # global attributes
    ds.Conventions = "CF-1.6"
    ds.history = "Created by globsim CLASSIC export"

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
    site : str, int, or None
        Site name or index to export.  ``None`` exports all sites.
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
        if (site is not None) and (site != i) and (site != names[i]):
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
                time_values=time_values,
                time_units=time_units,
                time_calendar=time_calendar,
                lat=lat,
                lon=lon,
                var_name=out_var,
                var_data=data,
                var_units=units,
                var_long_name=long_name,
            )
            files.append(filepath)

    return files
