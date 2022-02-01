import tomlkit
import logging
import numpy as np
import netCDF4 as nc

from datetime import datetime
from os import path, makedirs
from pathlib import Path
from math import floor

from globsim.common_utils import StationListRead

logger = logging.getLogger('globsim.scale')


class GenericScale:

    def __init__(self, sfile):
        # read parameter file
        self.sfile = sfile
        with open(self.sfile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = config['scale']
        self.set_parameters(self.par)

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
