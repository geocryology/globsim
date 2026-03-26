import logging
import netCDF4 as nc

from pathlib import Path

from globsim.scale.GenericScale import GenericScale, _check_timestep_length

logger = logging.getLogger('globsim.scale')


class JRAscale(GenericScale):
    """
    Class for JRA-55 data that has methods for scaling station data to
    better resemble near-surface fluxes.

    Processing kernels have names in UPPER CASE.

    Args:
        sfile: Full path to a Globsim Scaling Parameter file.

    Example:
        JRAd = JRAscale(sfile)
        JRAd.process()
    """
    NAME = "JRA-generic"
    REANALYSIS = "jra-generic"

    VARNAMES = {}
    CONVERTERS = {}
    

    def __init__(self, sfile):
        super().__init__(sfile)
        par = self.par

        # input file names
        self.nc_pl_sur = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_pl_{self.list_name}_surface.nc'), 'r')
        self.nc_pl = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_pl_{self.list_name}.nc'), 'r')
        self.nc_sa = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_sa_{self.list_name}.nc'), 'r')
        self.nc_sf = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_sf_{self.list_name}.nc'), 'r')
        self.nc_to = nc.Dataset(Path(self.interp_dir, f'{self.REANALYSIS}_to_{self.list_name}.nc'), 'r')
        
        # Check data integrity
        _check_timestep_length(self.nc_pl_sur.variables['time'], "pl_sur")
        _check_timestep_length(self.nc_pl.variables['time'], "pl")
        _check_timestep_length(self.nc_sa.variables['time'], "sa")
        _check_timestep_length(self.nc_sf.variables['time'], "sf")

        # time vector for output data
        # get time and convert to datetime object
        self.set_time_scale(self.nc_pl_sur.variables['time'], par['time_step'])
        
        self.times_out_nc = self.build_datetime_array(start_time=self.min_time,
                                                      timestep_in_hours=self.scaled_tstep_h,
                                                      num_times=self.nt,
                                                      output_units=self.scaled_t_units,
                                                      output_calendar=self.scaled_t_cal)
        

class JRA55(JRAscale):
      NAME = "JRA-55"
      REANALYSIS = "jra55"
      CONVERTERS = {}
      VARNAMES = {}