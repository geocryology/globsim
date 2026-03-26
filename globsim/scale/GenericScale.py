import enum
import logging
import numpy as np
import netCDF4 as nc
import pandas as pd
import pytz
import tomlkit

from cfunits import Units
from datetime import datetime
from math import floor
from os import path, makedirs
from pathlib import Path
from pysolar.solar import get_azimuth_fast
from scipy.interpolate import interp1d
from typing import Callable

from globsim.constants import G
from globsim.common_utils import StationListRead, series_interpolate
from globsim.scale.scalenames import ScaleNames as SN
from globsim.scale.toposcale import lw_down_toposcale, elevation_corrected_sw, solar_zenith, shading_corrected_sw_direct, illumination_angle
from globsim.nc_elements import new_scaled_netcdf

import globsim.dreamit as dreamit
import globsim.meteorology as met
import globsim.redcapp as redcapp
import globsim.scale.kernel_templates as kt


logger = logging.getLogger('globsim.scale')


class GenericScale:
    REANALYSIS: str = ''

    SCALING = {"sf": {},
               "sa": {},
               "pl": {},
               "to": {},
               "pl_sur": {}}
    
    VARNAMES ={"sf": {},  # translates scaling canonical names to input variable names
               "sa": {},
               "pl": {},
               "to": {},
               "pl_sur": {}}
    
    CONVERTERS = {}  # translator methods for e.g. accumulations to rates

    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        with open(self.sfile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = config['scale']
        self.set_parameters(self.par)
        self.scaled_t_units = 'minutes since 1970-01-01 00:00:00+00:00'
        self.scaled_t_cal = 'standard'
        # check if output file exists and remove if overwrite parameter is set
        self.output_file = self.getOutNCF(self.par, f'{self.REANALYSIS}')

    def upscale(self, time_in, values:np.ndarray, times_out):
        """"""
        if time_in.shape == times_out.shape and np.allclose(time_in, times_out):
            return values
        
        axis = [i for i,v in enumerate(values.shape) if v == time_in.shape[0]][0]
        if not axis:
            axis=0
        fv = (values.take(indices=0, axis=axis), values.take(indices=-1, axis=axis))
        f = interp1d(time_in, values, kind='linear', bounds_error=False, fill_value=fv, axis=axis)

        if np.ma.isMaskedArray(times_out):
            times_out = times_out.data

        return f(times_out)
    
    def _convert(self, file:str, name:str|enum.Enum, data:np.ndarray, nc_var: nc.Dataset, _slice=None) -> tuple[np.ndarray, str]:
        """Apply physics conversion if registered, return (data, units_str)."""
        key = (file, name)
        if key in self.CONVERTERS:
            method = getattr(self, self.CONVERTERS[key])
            return method(data, nc_var, _slice)
        
        # No conversion needed — units are whatever the file says
        return data, nc_var.units
    
    def process(self):
        """Run all relevant processes and save data. Each kernel processes one variable and adds it to the netCDF file."""
        if not path.isdir(path.dirname(self.output_file)):
            makedirs(path.dirname(self.outfile))
        
        self.set_valid_stations()
        valid_indices = self.valid_stations['nc_index']
        self.rg = new_scaled_netcdf(ncfile_out=self.output_file, 
                                    nc_interpol=self.nc_pl_sur,
                                    times_out=self.times_out_nc, 
                                    t_unit=self.scaled_t_units,
                                    valid_indices=valid_indices,
                                    station_names=self.valid_stations['station_name'],)
        
        elev = np.squeeze(self.get_values("to", SN.elevation, units='m'))
        self.add_grid_elevation(self.rg, elev[valid_indices.values])
        
        self.run_kernels()
        logger.info(f"Created scaled output file {self.output_file}")

        self.close_files()
    
    def close_files(self):
        self.rg.close()
        self.nc_pl_sur.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()
        self.nc_pl.close()
    
    def get_name(self, file:str, name:str):
        """ Globsim translation of variable names. If no translation is found, returns the original name. """
        if name in self.VARNAMES.get(file, {}).keys():
            return self.VARNAMES[file][name]
        
        else:
            logger.debug(f"No name translation found for variable '{name}' in file '{file}'. Using name as provided.")

        return name
    
    def get_attr(self, file:str, name:str, attr:str):
        f = self.get_file(file)
        n = self.get_name(file, name)
        return f.variables[n].getncattr(attr)

    def get_values(self, file:str, name:str|enum.Enum, _slice=None, units:str=None):
        """Get values for a given variable and file. Handles scaling and attribute retrieval.
        If units is provided, bypasses internal scaling and converts based on netcdf attributes"""
        f = self.get_file(file)
        n = self.get_name(file, name)
        
        if n not in f.variables:
            raise KeyError(f"Variable '{n}' not found in '{file}' file. Available variables: {list(f.variables.keys())}")
        nc_var = f.variables[n]

        v = nc_var[_slice] if _slice else nc_var[:]  # raw data

        v, effective_units = self._convert(file, name, v, nc_var, _slice)  # physical conversion (e.g. rate to accumulation)

        if (units is not None) and (effective_units != units):  # Equivalent units conversion
            v = Units.conform(v, Units(effective_units), Units(units), inplace=True)

        elif name in self.SCALING.get(file).keys():  # (Legacy) Scaling factors available
            scale, offset = self.SCALING.get(file).get(name)
            v *= scale
            v += offset

        return v
    
    def get_station_values(self, file:str, name:str, station_ix:int, preserve_dims:bool=False, units=None) -> np.ndarray:
        """ Get station values for a given variable and station index. Handles slicing for 2D and 3D variables. """
        if preserve_dims:
            station_ix = slice(station_ix, station_ix + 1)

        if file in ['sa', 'sf', 'to', 'pl_sur']:  #  time, station
            _slice = (slice(None), station_ix)
        elif file == 'pl':  #  time, level, station
            _slice = (slice(None), slice(None), station_ix)
        
        return self.get_values(file, name, _slice=_slice, units=units)

    def get_station_values_at(self, file: str, name: str, station_ix: int,
                          target_file: str, preserve_dims: bool = False,
                          units: str = None) -> np.ndarray:
        """Get station values interpolated to the time grid of target_file.    
        If file and target_file have the same timestep, no interpolation occurs."""
        values = self.get_station_values(file, name, station_ix,
                                        preserve_dims=preserve_dims, units=units)
        
        source_time = self._file_time_numeric(file)
        target_time = self._file_time_numeric(target_file)
        
        # Fast path: same time grid, no work needed
        if source_time.shape == target_time.shape and np.allclose(source_time, target_time):
            return values
        
        return self.upscale(source_time, values, target_time)

    def _file_time_numeric(self, file: str) -> np.ndarray:
        """Get time array for a file in a common numeric representation (int64)."""
        f = self.get_file(file)
        time_var = f.variables['time']
        # Convert to a common epoch so times from different files are comparable
        time_dates = nc.num2date(time_var[:], time_var.units, time_var.calendar)
        return nc.date2num(time_dates, units=self.scaled_t_units, calendar=self.scaled_t_cal).astype(np.int64)

    def set_valid_stations(self):
        ipl_station_ix=self.nc_pl_sur['station'][:]
        ipl_station_lon=self.get_values('pl_sur', SN.longitude)
        ipl_station_lat=self.get_values('pl_sur', SN.latitude)
        ipl_station_elev=self.get_values('pl_sur', SN.elevation)

        try:
            ipl_station_names=nc.chartostring(self.nc_pl['station_name'][:])
        except IndexError:
            logger.warning("No station_name variable in interpolated netCDF.")
            ipl_station_names = None
        
        interpolated_stations = pd.DataFrame(data={'station_number':ipl_station_ix,
                                                   'longitude_dd': ipl_station_lon, 
                                                   'latitude_dd': ipl_station_lat, 
                                                   'station_name': ipl_station_names,
                                                   'elevation_m':ipl_station_elev})
        interpolated_stations['nc_index'] = interpolated_stations.index

        stations = self.stations.copy()
        stations['siteslist_index'] = stations.index
        
        if (~interpolated_stations['station_name'].isna().any()):
            station_df = stations.merge(interpolated_stations, on=['station_number', 'station_name'], 
                                        how='inner', suffixes=('_scale', '_interpolate'))
        else:
            logger.warning("One or more station names in interpolated netCDF are missing. Matching stations based on station number and coordinates only.")
            station_df = stations.merge(interpolated_stations, on='station_number', how='inner', suffixes=('_scale', '_interpolate'))

        station_df['lon_matches'] = np.isclose(station_df['longitude_dd_scale'] % 360, station_df['longitude_dd_interpolate'] % 360, atol=1e-6)
        station_df['lat_matches'] = np.isclose(station_df['latitude_dd_scale'] % 360, station_df['latitude_dd_interpolate'] % 360, atol=1e-6)
        station_df['elev_matches'] = np.isclose(station_df['elevation_m_scale'], station_df['elevation_m_interpolate'], atol=1e-4)
        station_df['coordinates_match'] = station_df['lon_matches'] & station_df['lat_matches'] & station_df['elev_matches']
        
        if ipl_station_names is None:
            station_df['station_name_scale'] = None
            station_df['name_matches'] = station_df['station_name_scale'] == station_df['station_name_interpolate']

            for _, row in station_df[~station_df['coordinates_match']].iterrows():
                logger.warning(f"Station {row['station_number']} ({row['station_name_scale']})" \
                            f"has mismatched coordinates between station list" \
                            f"({row['longitude_dd_scale']}, {row['latitude_dd_scale']}) and interpolated netCDF" \
                            f"({row['longitude_dd_interpolate']}, {row['latitude_dd_interpolate']}).")

            for _, row in station_df[~station_df['name_matches']].iterrows():
                if row['station_name_interpolate'] is not None:
                    logger.warning(f"Station {row['station_number']} ({row['station_name_scale']}) " \
                                f"has mismatched names between station list ({row['station_name_scale']})" \
                                f"and interpolated netCDF ({row['station_name_interpolate']}).")
        else:
            station_df['name_matches'] = True  
            logger.info(f"Found {len(station_df)} valid stations (out of {len(self.stations)} in siteslist and {len(interpolated_stations)} in interpolated netCDF) " \
                        f"with matching station numbers and coordinates between station list and interpolated netCDF. ") 

        self.valid_stations = station_df[station_df['coordinates_match']].copy()
        self.nstation = self.valid_stations.shape[0]


    def iterate_stations(self):
        """Iterate through stations, returning siteslist index and interpolated file index."""
        for row in self.valid_stations.itertuples():
            yield row.siteslist_index, row.nc_index

    def get_time_step(self, file:str) -> int:
        """ Returns time step of file in seconds """
        time_var = self.get_file(file).variables['time']
        time_in = time_var[:]
        t1 = nc.num2date(time_in[1], units=time_var.units, calendar=time_var.calendar)
        t0 = nc.num2date(time_in[0], units=time_var.units, calendar=time_var.calendar)
        return (t1 - t0).seconds 

    def get_file(self, file:str) -> "nc.Dataset":
        if file == "sa":
            f = self.nc_sa
        elif file == "sf":
            f = self.nc_sf
        elif file == "pl_sur":
            f = self.nc_pl_sur
        elif file == "pl":
            f = self.nc_pl
        elif file == "to":
            f = self.nc_to
        else:
            raise ValueError("sa, sf, to, or pl_sur")
        
        return f  

    def set_parameters(self, par):
        self.interp_dir = path.join(par['project_directory'], 'interpolated')
        self.output_dir = self.make_output_directory(par)
        self.list_name  = path.basename(path.normpath(par['station_list'])).split(path.extsep)[0]

        # get the station file
        self.stations_csv = path.join(par['project_directory'],
                                      'par', par['station_list'])
        # read station points
        self.stations = StationListRead(self.stations_csv)

        # read kernels
        self.kernels = par['kernels']
        if not isinstance(self.kernels, list):
            self.kernels = [self.kernels]

        # should file be overwritten - default to false
        try:
            self._overwrite_output = par['overwrite']
        except KeyError:
            logger.warning("Missing overwrite parameter in control file. Reverting to default (ovewrite = true).")
            self._overwrite_output = False
        finally:
            logger.debug(f"Overwriting of output files set to '{self._overwrite_output}'")

        # read snow correction info
        try:
            self.scf = par['scf']
        except KeyError:
            logger.warning("Missing snow correction factor parameter in control file. Reverting to default (scf = 1).")
            self.scf = 1
        finally:
            logger.debug(f"Snow correction factor for scaling set to {self.scf}")

        # read RH approximation
        try:
            rhf = par['rh_approximation']
        except KeyError:
            logger.warning("Missing relative humidity approximation choice in control file (rh_approximation). Reverting to default ('rh_liston').")
            rhf = 'rh_liston'
        finally:
            self._rh_function_name = rhf
            logger.debug(f"Using relative humidity approximation {rhf}")

    def getOutNCF(self, par, data_source_name):
        """make out file name"""

        timestep = str(par['time_step']) + 'h'
        snowCor  = 'scf' + str(self.scf)
        src = '_'.join(['scaled', data_source_name, timestep, snowCor])

        src = src + '.nc'
        output_file = Path(self.output_dir, src)

        if output_file.is_file():
            if self._overwrite_output:
                output_file.unlink()
                logger.info(f"Removed existing output file {output_file}")
            else:
                logger.error("Scaled file output already exists")
                raise FileExistsError(f"Output file {output_file} exists. Remove file or set 'overwite=true' in control file.")

        return output_file.resolve()

    def get_slope(self) -> "np.ndarray":
        if hasattr(self, "__slope"):
            return getattr(self, "__slope")

        if 'slope' in self.stations.columns:
            slope = self.stations['slope'].to_numpy(dtype='float32')
        else:
            logger.warning("No 'slope' column in siteslist. Assuming horizontal surface (00)")
            slope = np.zeros_like(self.stations['longitude_dd'].values)
        
        self.__slope = np.atleast_1d(slope)  # Cache for later use

        return self.__slope

    def get_aspect(self) -> "np.ndarray":
        if hasattr(self, "__aspect"):
            return getattr(self, "__aspect")

        if 'aspect' in self.stations.columns:
            aspect = self.stations['aspect'].to_numpy(dtype='float32')
        else:
            logger.warning("No 'aspect' column in siteslist. Assuming north aspect (000)")
            aspect = np.zeros_like(self.stations['longitude_dd'].values)
        
        self.__aspect = np.atleast_1d(aspect)  # Cache for later use

        return self.__aspect

    def get_sky_view(self) -> "np.ndarray":
        if hasattr(self, "__skyview"):
            return getattr(self, "__skyview")

        if 'sky_view' in self.stations.columns:
            svf = self.stations['sky_view'].to_numpy(dtype='float32')
        else:
            logger.warning("No 'sky_view' column in siteslist. Assuming SVF = 1.0")
            svf = np.ones_like(self.stations['longitude_dd'].values)
        
        self.__skyview = np.atleast_1d(svf)  # Cache for later use

        return self.__skyview
    
    def get_hypsometry(self) -> "np.ndarray":
        if hasattr(self, "__hypsometry"):
            return getattr(self, "__hypsometry")

        if 'hypsometry' in self.stations.columns:
            hyps = self.stations['hypsometry'].to_numpy(dtype='float32')
        else:
            logger.warning("No 'hypsometry' column in siteslist. Assuming h = 0 (top of ridge)")
            hyps = np.ones_like(self.stations['longitude_dd'].values)
        
        self.__hypsometry = np.atleast_1d(hyps)  # Cache for later use

        return self.__hypsometry

    def make_output_directory(self, par):
        """make directory to hold outputs"""
        output_dir = None

        if par.get('output_directory'):
            try:
                test_path = Path(par.get('output_directory'))
            except TypeError:
                msg = "You provided an output_directory for scaled files that does not exist. Saving files to project directory"
                logger.warning(msg)

            if test_path.is_dir():
                output_dir = Path(test_path, "scaled")
            else:
                logger.warning("You provided an output_directory for scaled files that was not understood. Saving files to project directory.")

        if not output_dir:
            logger.debug(f"Attempting to create directory: {output_dir}")
            output_dir = path.join(par['project_directory'], 'scaled')

        if not Path(output_dir).is_dir():
            makedirs(output_dir)

        return output_dir

    def set_time_scale(self, time_variable, time_step):
        nctime = time_variable[:]

        t_unit = time_variable.units
        t_cal  = time_variable.calendar
        self.scaled_tstep_h = time_step
        self.min_time = nc.num2date(min(nctime), units=t_unit, calendar=t_cal)

        max_time = nc.num2date(max(nctime), units=t_unit, calendar=t_cal)
        t1 = nc.num2date(nctime[1], units=t_unit, calendar=t_cal)
        t0 = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
        interval_in = (t1 - t0).seconds
        
        # number of time steps
        self.nt = floor((max_time - self.min_time).total_seconds() / (3600 * time_step)) + 1
        logger.debug(f"Output time array has {self.nt} elements between "
                     f"{self.min_time.strftime('%Y-%m-%d %H:%M:%S')} and "
                     f"{max_time.strftime('%Y-%m-%d %H:%M:%S')}"
                     f" (time step of {self.scaled_tstep_h} hours)")

    @staticmethod
    def build_datetime_array(start_time: datetime, timestep_in_hours: int, num_times: int, output_units:str, output_calendar:str):
        time_array = np.arange(num_times, dtype='float64') * timestep_in_hours * 3600

        datetime_array = nc.num2date(time_array,
                                     units=f"seconds since {start_time.strftime('%Y-%m-%d %H:%M:%S.%f')}",
                                     calendar="gregorian")

        result = nc.date2num(datetime_array,
                             units=output_units,
                             calendar=output_calendar)

        return result

    def run_kernels(self):
        for kernel_name in self.kernels:
            if hasattr(self, kernel_name):
                logger.info(f"running scaling kernel: '{kernel_name}'")
                kernel = getattr(self, kernel_name)
                _ = kernel()
            else:
                logger.error(f"Missing kernel {kernel_name}")

    def _rh(self) -> Callable:
        rh_function = getattr(met, self._rh_function_name)
        return rh_function
    
    def add_grid_elevation(self, ncf: "nc.Dataset", data: np.ndarray) -> None:
        """Add station elevation to the netCDF file"""
        elev_var = ncf.createVariable('grid_elevation', 'f4', ('station',))
        elev_var.units = 'm'
        elev_var.long_name = 'grid elevation'
        elev_var.comment = 'Elevation of the grid cell at the station location'
        elev_var[:] = data

        self.warn_station_elevation(data)

    def input_times_in_output_units(self, ncf):
        """
        Convert time in input data to time in Globsim.
        """
        raw = ncf['time'][:].astype(np.int64)
        time = nc.num2date(raw, units=ncf['time'].units, calendar=ncf['time'].calendar)
        converted = nc.date2num(time, units=self.scaled_t_units, calendar=self.scaled_t_cal)
        return converted.astype(np.int64)

    def warn_station_elevation(self, data: "np.ndarray") -> None:
        """Warn if there are stations with elevation below grid level"""
        stn_elev = self.stations['elevation_m']
        stn_name = self.stations['station_name']

        for i, (stn, grid) in enumerate(zip(stn_elev, data)):
            if stn < grid:
                logger.warning(f" {stn_name[i]} site elevation ({stn} m) is below reanalysis grid elevation ({grid} m). Results may be unreliable.")

    def PRESS_Pa_pl(self):
        """
        Surface air pressure from pressure levels.
        """
        vn = kt.PRESS_Pa_pl(self.rg, self.NAME)
        
        time_in = self.get_values("pl_sur", SN.time).astype(np.int64)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("pl_sur", SN.pressure, interp_ix, units="Pa") 
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc,
                                                             time_in,
                                                             values)
            
    def _geopotential_to_m(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        units = Units(nc_var.units) / Units("m s-2")
        return data / G, units.units
    
    def AIRT_C_pl(self):
        """
        Air temperature derived from pressure levels, exclusively.
        """
        vn = kt.AIRT_C_pl(self.rg, self.NAME)
        time_in = self.get_values("pl_sur","time")

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("pl_sur", SN.temperature, interp_ix, units="degree_C")
            self.rg.variables[vn][:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)

    def AIRT_C_sur(self):
        """
        Air temperature derived from surface data, exclusively.
        """
        vn = kt.AIRT_C_sur(self.rg, self.NAME)
        var = self.rg.variables[vn] 
        time_in = self.get_values("sa", SN.time)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sa", SN.temperature, interp_ix, units=var.units)
            var[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)

    def AIRT_redcapp(self):
        """
        Air temperature derived from surface data and pressure level data as
        shown by the method REDCAPP Cao et al. (2017) 10.5194/gmd-10-2905-2017
        """
        logger.warning(f"Globsim implementation of REDCAPP only provides Delta_T_c")

        var = redcapp.add_var_delta_T(self.rg)
        time_in = self.input_times_in_output_units(self.nc_sa)

        for siteslist_ix, interp_ix in self.iterate_stations():
            T_sa  = self.get_station_values_at("sa", SN.temperature, interp_ix, "sf", preserve_dims=True, units="degree_K")  
            airT_pl = self.get_station_values_at("pl", SN.temperature, interp_ix, "sf", preserve_dims=True, units="degree_K")
            elevation = self.get_station_values_at("pl", SN.elevation, interp_ix, "sf", preserve_dims=True, units='m')
            h_sur = self.get_station_values("to", SN.elevation, interp_ix, preserve_dims=False, units='m')

            Delta_T_c = redcapp.delta_T_c(T_sa=T_sa, 
                                          airT_pl=airT_pl, 
                                          elevation=elevation,
                                          h_sur=h_sur)  

            values  = Delta_T_c[:, 0]
            var[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)      

    def PREC_mm_sur(self):
        """
        Precipitation derived from surface data, exclusively.
        Convert unit: to mm/s (kg m-2 s-1)
        """
        vn  = kt.PREC_mm_sur(self.rg, self.NAME)
        var = self.rg.variables[vn]
        time_in = self.get_values("sf", SN.time)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", SN.precipitation_rate, interp_ix, units=var.units)
            var[:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values) * self.scf
    
    def RH_per_sur(self):
        """
        Relative Humidity derived from surface data, exclusively.
        """
        vn = kt.RH_per_sur(self.rg, self.NAME)
        var = self.rg.variables[vn]
        time_in = self.get_values("sa", SN.time)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sa", SN.rh, interp_ix, units=var.units)
            var[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values) 
    
    def RH_per_pl(self):
        """
        Relative Humidity derived from pressure-level data, exclusively.
        """
        vn = kt.RH_per_pl(self.rg, self.NAME)
        var = self.rg.variables[vn]
        time_in = self.get_values("pl_sur", SN.time)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("pl_sur", SN.rh, interp_ix, units=var.units)
            self.rg.variables[vn][:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)   

    def SH_kgkg_sur(self):
        '''
        Specific humidity [kg/kg] derived from surface data, exclusively
        '''
        vn = kt.SH_kgkg_sur(self.rg, self.NAME) 
        var = self.rg.variables[vn]
        time_in = self.get_values("sa", SN.time)

        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sa", SN.specific_humidity, interp_ix, units=var.units)
            var[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)

    def SW_Wm2_sur(self):
        """
        solar radiation downwards derived from surface data, exclusively.
        """
        vn = kt.SW_Wm2_sur(self.rg, self.NAME)
        var = self.rg.variables[vn]
        time_in = self.get_values("sf", SN.time)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", SN.sw_down_flux, interp_ix, units=var.units)
            var[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)

    def LW_Wm2_sur(self):
        """
        Long-wave radiation downwards derived from surface data, exclusively.
        """
        vn = kt.LW_Wm2_sur(self.rg, self.NAME) 
        var = self.rg.variables[vn]
        time_in = self.get_values("sf","time")
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values  = self.get_station_values("sf", SN.lw_down_flux, interp_ix, units=var.units)
            var[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values)

    def WIND_sur(self):
        """
        Wind at 10 metre derived from surface data, exclusively.
        """
        vn_u, vn_v, vn_spd, vn_dir = kt.WIND_sur(self.rg, self.NAME)
        var_u = self.rg.variables[vn_u]
        var_v = self.rg.variables[vn_v]
        var_wspd = self.rg.variables[vn_spd]
        var_wdir = self.rg.variables[vn_dir]

        time_in = self.get_values("sa", SN.time)
        
        for siteslist_ix, interp_ix in self.iterate_stations():
            values_u  = self.get_station_values("sa", SN.u_wind, interp_ix, units=var_u.units)
            var_u[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values_u)
            values_v  = self.get_station_values("sa", SN.v_wind, interp_ix, units=var_v.units)
            var_v[:, siteslist_ix] = np.interp(self.times_out_nc, time_in, values_v)

        # convert
        # u is the ZONAL VELOCITY, i.e. horizontal wind TOWARDS EAST.
        # v is the MERIDIONAL VELOCITY, i.e. horizontal wind TOWARDS NORTH.
        V = self.rg.variables['10 metre V wind component'][:]
        U = self.rg.variables['10 metre U wind component'][:]
        
        WS = np.sqrt(np.power(V, 2) + np.power(U, 2))
        WD = 90 - (np.arctan2(V, U) * (180 / np.pi)) + 180
        WD = np.mod(WD, 360)

        var_wspd[:] = WS
        var_wdir[:] = WD

    def SW_Wm2_topo(self):
        """
        Short-wave downwelling radiation corrected using a modified version of TOPOscale.
        Partitions into direct and diffuse
        """
        vn_dir, vn_diff, vn_glob = kt.SW_Wm2_topo(self.rg, self.NAME)

        nc_time = self.input_times_in_output_units(self.nc_sf)
        py_time = nc.num2date(nc_time[:], self.scaled_t_units, self.scaled_t_cal, only_use_cftime_datetimes=False)
        py_time = np.array([pytz.utc.localize(t) for t in py_time])
        
        lat = self.get_values('pl_sur', SN.latitude)
        lon = self.get_values('pl_sur', SN.longitude)

        svf = self.get_sky_view()
        slope = self.get_slope()
        aspect = self.get_aspect()

        grid_elev = np.squeeze(self.get_values('to', SN.elevation, units='m'))
        station_elev = self.get_values("pl_sur", SN.elevation, units='m')

        for siteslist_ix, interp_ix in self.iterate_stations():
            zenith = solar_zenith(lat=lat[interp_ix], lon=lon[interp_ix], time=py_time)
            sw = self.get_station_values('sf', SN.sw_down_flux, interp_ix, preserve_dims=False, units="W m-2")
            diffuse, corrected_direct = elevation_corrected_sw(zenith=zenith,
                                                               grid_sw=sw,
                                                               lat=np.ones_like(sw) * lat[interp_ix],
                                                               lon=np.ones_like(sw) * lon[interp_ix],
                                                               time=py_time,
                                                               grid_elevation=np.ones_like(sw) * grid_elev[interp_ix],
                                                               sub_elevation=np.ones_like(sw) * station_elev[interp_ix])

            diffuse = diffuse * svf[siteslist_ix]  # apply sky-view factor

            if slope[siteslist_ix] != 0:
                azimuth = get_azimuth_fast(lat[interp_ix], lon[interp_ix], py_time)
                cos_i_sub = illumination_angle(zenith, azimuth, slope[siteslist_ix], aspect[siteslist_ix])
                cos_i_grid = np.cos(np.radians(zenith))
                corrected_direct = shading_corrected_sw_direct(corrected_direct, cos_i_sub, cos_i_grid)
                
                sensible_values_mask = np.where(cos_i_grid < 0.001, 0, 1) * np.where(corrected_direct > 1366, 0, 1)
                corrected_direct *= sensible_values_mask

            global_sw = diffuse + corrected_direct

            self.rg.variables[vn_dir][:, siteslist_ix] = series_interpolate(self.times_out_nc, nc_time, corrected_direct)
            self.rg.variables[vn_diff][:, siteslist_ix] = series_interpolate(self.times_out_nc, nc_time, diffuse)
            self.rg.variables[vn_glob][:, siteslist_ix] = series_interpolate(self.times_out_nc, nc_time, global_sw)


    def LW_Wm2_topo(self):
        """ Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data"""
        vn = kt.LW_Wm2_topo(self.rg, self.NAME)

        time_in = self.input_times_in_output_units(self.nc_sf)
        svf = self.get_sky_view()

        for siteslist_ix, interp_ix in self.iterate_stations():
            t_sub = self.get_station_values_at("pl_sur", SN.temperature, interp_ix, "sf", units="degree_K")  # [K]
            rh_sub = self.get_station_values_at("pl_sur", SN.rh, interp_ix, "sf", units="percent")  # [%]
            t_grid = self.get_station_values_at("sa", SN.temperature, interp_ix, "sf", units="degree_K")  # [K]
            rh_grid = self.get_station_values_at("sa", SN.rh, interp_ix, "sf", units="percent")  # [%]
            lw_grid  = self.get_station_values_at("sf", SN.lw_down_flux, interp_ix, "sf")

            lw_sub = lw_down_toposcale(t_sub=t_sub, rh_sub=rh_sub, t_sur=t_grid, rh_sur=rh_grid, lw_sur=lw_grid)

            values = lw_sub * svf[siteslist_ix]
            self.rg.variables[vn][:, siteslist_ix] = series_interpolate(self.times_out_nc, time_in, values)

    def AIRT_DReaMIT(self):
        """
        Air temperature derived from surface data, pressure level data, and
        dynamically-computed inversion metrics as shown by the method DReaMIT
        """
        kt.AIRT_DReaMIT(self.rg)

        nc_time = self.nc_sa.variables['time']
        time_in = self.input_times_in_output_units(self.nc_sa)

        hypsometry = self.get_hypsometry()
        list_params = dreamit.get_model_params(self.REANALYSIS)
        time_frac_year = dreamit.time_frac_year(nc_time)
        pl_height = self.get_values("pl_sur", SN.elevation, units='m')

        for siteslist_ix, interp_ix in self.iterate_stations():
            T_pl_in = self.get_station_values("pl", SN.temperature, interp_ix, preserve_dims=True, units='degree_K')
            h_pl_in = self.get_station_values("pl", SN.elevation, interp_ix, preserve_dims=True, units='m')
            T_pl_surface_in = self.get_station_values("pl_sur", SN.temperature, interp_ix, preserve_dims=True, units='degree_K')
            h_pl_surface_in = pl_height[interp_ix: interp_ix+1] # preserve dimension
            grid_elev_in = self.get_station_values("to", SN.elevation, interp_ix, preserve_dims=True, units='m')

            T_sur_in = self.get_station_values("sa", SN.temperature, interp_ix, preserve_dims=True, units='degree_K')            

            mtrcs = dreamit.dreamit_metrics(reanalysis=self.REANALYSIS,
                                            T_pl_in=T_pl_in,
                                            h_pl_in=h_pl_in,
                                            T_pl_surface_in=T_pl_surface_in,
                                            h_pl_surface_in=h_pl_surface_in,
                                            grid_elev_in=grid_elev_in)
            
            z_top_inversion_m, T_lapse_grid_K, T_lapse_station_K, lapse_Cperm = mtrcs



            AIRT_DReaMIT_K, beta_t_C = dreamit.dreamit_air_T(T_lapse_grid=T_lapse_grid_K,
                                                            T_lapse_station=T_lapse_station_K,
                                                            T_sur=T_sur_in,
                                                            time_frac_year=time_frac_year,
                                                            hyps=np.atleast_1d(hypsometry),
                                                            params=list_params)
            
            var_ztop = self.rg.variables['z_top_inversion_m']
            var_ztop[:, siteslist_ix] = np.interp(self.times_out_nc,
                                                  time_in, 
                                                  np.squeeze(z_top_inversion_m))
            
            var_tl_grid = self.rg.variables['T_lapse_grid_C']
            var_tl_grid[:, siteslist_ix] = np.interp(self.times_out_nc,
                                                     time_in, 
                                                     np.squeeze(T_lapse_grid_K) - 273.15)
            
            var_tl_stn = self.rg.variables['T_lapse_station_C']
            var_tl_stn[:, siteslist_ix] = np.interp(self.times_out_nc,
                                                    time_in, 
                                                    np.squeeze(T_lapse_station_K) - 273.15)
            
            var_lapse = self.rg.variables['lapse_Cperm']
            var_lapse[:, siteslist_ix] = np.interp(self.times_out_nc,
                                                   time_in, 
                                                   np.squeeze(lapse_Cperm))
            
            var_airt_drm = self.rg.variables['AIRT_DReaMIT_C']
            var_airt_drm[:, siteslist_ix] = np.interp(self.times_out_nc,
                                                      time_in, 
                                                      np.squeeze(AIRT_DReaMIT_K) - 273.15)
        
        var_beta = self.rg.variables['beta_t_C']
        var_beta[:] = np.interp(self.times_out_nc, time_in, beta_t_C[:])   

def _check_timestep_length(nctime: "nc.Variable", source:str) -> None:
    """ Ensure that input data has a consistent timestep

    An inconsistent timestep usually suggests that there is a data gap.

    Parameters
    ----------
    nctime : netCDF4.Variable
        A time variable with 'units' attribute in the form of 'X since Y'
    source : str
        Where does the data come from?
    """
    # Check for a single time step
    steps = np.unique(np.diff(nctime[:]))

    units = nctime.units.split(" ")[0] if hasattr(nctime, 'units') else ""

    if len(steps) != 1:
        logger.critical(f"Input time step length is not consistent for {source}. Found time steps: {list(steps)} ({units}).")
