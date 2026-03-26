import urllib3
import logging

from globsim.nc_elements import new_scaled_netcdf
from globsim.scale.ERA5scale import ERA5scale

urllib3.disable_warnings()

logger = logging.getLogger('globsim.scale')


class ERA5Escale(ERA5scale):
    """
    Class for ERA5 data that has methods for scaling station data to
    better resemble near-surface fluxes.

    Processing kernels have names in UPPER CASE.

    Args:
        sfile: Full path to a Globsim Scaling Parameter file.

    Example:
        ERAd = ERA5scale(sfile)
        ERAd.process()
    """
    REANALYSIS = 'era5_ens'
    NAME = "ERA-5 Ensemble"

    def __init__(self, sfile):
        super().__init__(sfile)
        self._ni = None

    @property
    def current_member(self) -> int:
        if self._ni is None:
            raise ValueError("Ensemble number must first be set")
        return self._ni

    @current_member.setter
    def current_member(self, ni:int):
        if ni in self.nc_sa['number']:
            self._ni = ni
        else:
            raise ValueError(f"Ensemble member number {ni} not in {self.nc_sa['number']}")

    def get_values(self, *args, **kwargs):
        ni = self.current_member  # get current ensemble member number
        # ncf.variables[varStr][:, ni, :]  # additional subsetting for ensemble member
        raise RuntimeError("get_values is not yet implemented for ERA5-ensemble.")

    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """

        stations = self.stations['station_name']
        # iterate thorugh kernels and start process

        for ni in self.nc_sa['number']:
            self.current_member = ni
            src = '{}_{}'.format(self.REANALYSIS, ni)
            self.output_file = self.getOutNCF(self.par, src)
            self.rg = new_scaled_netcdf(self.output_file,
                                        self.nc_pl_sur, self.times_out_nc,
                                        t_unit=self.scaled_t_units,
                                        station_names=stations,
                                        valid_indices=self.valid_indices)
            self.run_kernels()
            logger.info(f"Created scaled output file {self.output_file}")

        # close netCDF files
        self.rg.close()
        self.nc_pl_sur.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()
