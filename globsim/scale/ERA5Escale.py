# -*- coding: utf-8 -*-
#
# Methods for downloading ERA5 data from the ECMWF server for limited
# areas and limited times.
#
#  OVERALL WORKFLOW
#
#  See: https://github.com/geocryology/globsim/wiki/Globsim
#       and the examples directory on github
#
# ECMWF and netCDF information:
# https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
#
# For variable codes and units, see:
#     http://www.ecmwf.int/publications/manuals/d/gribapi/param/
#
# Check ECMWF job status: http://apps.ecmwf.int/webmars/joblist/
#
#  CF-Convention: this is how you check netcdf file conformity:
#  http://pumatest.nerc.ac.uk/cgi-bin/cf-checker-dev.pl
#
# ===============================================================================
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
    src = 'era5_ens'

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

    def getValues(self, ncf, varStr):  # must redefine
        ni = self.current_member
        return ncf.variables[varStr][:, ni, :]

    def process(self):
        """
        Run all relevant processes and save data. Each kernel processes one
        variable and adds it to the netCDF file.
        """

        stations = self.stations['station_name']
        # iterate thorugh kernels and start process

        for ni in self.nc_sa['number']:
            self.current_member = ni
            src = '{}_{}'.format(self.src, ni)
            self.output_file = self.getOutNCF(self.par, src)
            self.rg = new_scaled_netcdf(self.output_file,
                                        self.nc_pl, self.times_out_nc,
                                        t_unit=self.scaled_t_units,
                                        station_names=stations)
            self.indProcess()
            logger.info(f"Created scaled output file {self.output_file}")

        # close netCDF files
        self.rg.close()
        self.nc_pl.close()
        self.nc_sf.close()
        self.nc_sa.close()
        self.nc_to.close()
