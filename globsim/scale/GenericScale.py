import tomlkit
import logging
import numpy as np
import netCDF4 as nc
import pandas as pd

from datetime import datetime
from os import path, makedirs
from pathlib import Path
from math import floor
from typing import Callable
from scipy.interpolate import interp1d

from globsim.common_utils import StationListRead
import globsim.meteorology as met

logger = logging.getLogger('globsim.scale')


class GenericScale:
    REANALYSIS: str = ''

    SCALING = {"sf": {},
               "sa": {},
               "pl": {},
               "to": {},
               "pl_sur": {}}

    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        with open(self.sfile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = config['scale']
        self.set_parameters(self.par)

    @staticmethod
    def upscale(time_in, values:np.ndarray, times_out):
        axis = [i for i,v in enumerate(values.shape) if v == time_in.shape[0]][0]
        if not axis:
            axis=0
        fv = (values.take(indices=0, axis=axis), values.take(indices=-1, axis=axis))
        f = interp1d(time_in, values, kind='linear', bounds_error=False, fill_value=fv, axis=axis)
        # handle masked array
        if np.ma.isMaskedArray(times_out):
            times_out = times_out.data
        return f(times_out)
    
    def get_name(self, file:str, name:str):
        return name
    
    def get_values(self, file:str, name:str, _slice=None, attr=None):
        f = self.get_file(file)
        n = self.get_name(file, name)
        
        if attr is not None:
            v = f.variables[n].getncattr(attr)
        
        else:
            if _slice is None:
                v = f.variables[n][:]
            else:
                v = f.variables[n][_slice]
            
            if name in self.SCALING.get(file).keys():
                scale, offset = self.SCALING.get(file).get(name)
                v *= scale
                v += offset

        return v
    
    def get_station_values(self, file:str, name:str, station_ix:int, preserve_dims:bool=False) -> np.ndarray:
        """ Get station values for a given variable and station index. Handles slicing for 2D and 3D variables. """
        if preserve_dims:
            station_ix = slice(station_ix, station_ix + 1)

        if file in ['sa', 'sf', 'to', 'pl_sur']:  #  time, station
            _slice = (slice(None), station_ix)
        elif file == 'pl':  #  time, level, station
            _slice = (slice(None), slice(None), station_ix)
        
        return self.get_values(file, name, _slice=_slice)

    def set_valid_stations(self, interpolated_ncf: nc.Dataset):
        ipl_station_ix=self.nc_pl_sur['station'][:]
        ipl_station_lon=self.get_values('pl_sur', 'longitude')
        ipl_station_lat=self.get_values('pl_sur', 'latitude')
        ipl_station_elev=self.get_values('pl_sur', 'elevation')

        try:
            ipl_station_name=self.nc_pl_sur['station_name'][:]
        except IndexError:
            logger.warning("No station_name variable in interpolated netCDF.")
            ipl_station_name = None
        
        interpolated_stations = pd.DataFrame(data={'station_number':ipl_station_ix,
                                                   'longitude_dd': ipl_station_lon, 
                                                   'latitude_dd': ipl_station_lat, 
                                                   'station_name': ipl_station_name,
                                                   'elevation_m':ipl_station_elev})
        interpolated_stations['nc_index'] = interpolated_stations.index

        stations = self.stations.copy()
        stations['siteslist_index'] = stations.index
        
        if (interpolated_stations['station_name'] != None).all():
            station_df = stations.merge(interpolated_stations, on=['station_number', 'station_name'], 
                                        how='inner', suffixes=('_scale', '_interpolate'))
        else:
            logger.warning("One or more station names in interpolated netCDF are missing. Matching stations based on station number and coordinates only.")
            station_df = stations.merge(interpolated_stations, on='station_number', how='inner', suffixes=('_scale', '_interpolate'))
        
        station_df['lon_matches'] = np.isclose(station_df['longitude_dd_scale'], station_df['longitude_dd_interpolate'], atol=1e-6)
        station_df['lat_matches'] = np.isclose(station_df['latitude_dd_scale'], station_df['latitude_dd_interpolate'], atol=1e-6)
        station_df['elev_matches'] = np.isclose(station_df['elevation_m_scale'], station_df['elevation_m_interpolate'], atol=1e-4)
        station_df['coordinates_match'] = station_df['lon_matches'] & station_df['lat_matches'] & station_df['elev_matches']
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
            
        logger.info(f"Found {len(station_df)} valid stations with matching station numbers and coordinates between station list and interpolated netCDF. " \
                    f"{(~station_df['name_matches']).sum()} have mismatched or missing names.")
        
        self.valid_stations = station_df[station_df['coordinates_match']].copy()
        self.nstation = self.valid_stations.shape[0]


    def iterate_stations(self):
        """Iterate through stations, returning siteslist index and interpolated file index."""
        for row in self.valid_stations.itertuples():
            yield row.siteslist_index, row.nc_index

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
        self.intpdir = path.join(par['project_directory'], 'interpolated')
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

        self.t_unit = time_variable.units
        self.t_cal  = time_variable.calendar
        self.time_step = time_step
        self.min_time = nc.num2date(min(nctime), units=self.t_unit, calendar=self.t_cal)

        self.max_time = nc.num2date(max(nctime), units=self.t_unit, calendar=self.t_cal)
        t1 = nc.num2date(nctime[1], units=self.t_unit, calendar=self.t_cal)
        t0 = nc.num2date(nctime[0], units=self.t_unit, calendar=self.t_cal)
        self.interval_in = (t1 - t0).seconds
        
        # number of time steps
        self.nt = floor((self.max_time - self.min_time).total_seconds() / (3600 * time_step)) + 1
        logger.debug(f"Output time array has {self.nt} elements between "
                     f"{self.min_time.strftime('%Y-%m-%d %H:%M:%S')} and "
                     f"{self.max_time.strftime('%Y-%m-%d %H:%M:%S')}"
                     f" (time step of {self.time_step} hours)")

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
                getattr(self, kernel_name)()
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

    def warn_station_elevation(self, data: "np.ndarray") -> None:
        """Warn if there are stations with elevation below grid level"""
        stn_elev = self.stations['elevation_m']
        stn_name = self.stations['station_name']

        for i, (stn, grid) in enumerate(zip(stn_elev, data)):
            if stn < grid:
                logger.warning(f" {stn_name[i]} site elevation ({stn} m) is below reanalysis grid elevation ({grid} m). Results may be unreliable.")


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
