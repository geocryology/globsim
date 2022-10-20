import numpy as np

from globsim.download.ERA5download import ERA5download, ERA5pl, ERA5sa, ERA5sf, ERA5to


class ERA5Epl(ERA5pl):

    @property
    def eratype(self):
        return 'ensemble_members'


class ERA5Eto(ERA5to):

    @property
    def eratype(self):
        return 'ensemble_members'


class ERA5Esf(ERA5sf):

    @property
    def eratype(self):
        return 'ensemble_members'


class ERA5Esa(ERA5sa):

    @property
    def eratype(self):
        return 'ensemble_members'


class ERA5Edownload(ERA5download):
    """
    Class for ERA5 data that has methods for querying
    the ECMWF server, returning all variables usually needed.

    Args:
        pfile: Full path to a Globsim Download Parameter file.
        api  : Which API to use. Either 'cds' (default) or 'ecmwf' (deprecated'
        storage : Which server to access data from. Either 'cds'
                    (default) or 'ecmwf' (deprecated). Note
                 that you can use the cds api to access the ecmwf storage (MARS)
    Example:
        ERAd = ERA5download(pfile)
        ERAd.retrieve()
    """
    def __init__(self, pfile):
        super().__init__(pfile)

    def typeString(self):
        return 'era5_ens'

    def timeString(self):
        times = np.arange(0, 24, 3)
        times = [str(t).zfill(2) for t in times]
        times = [t + ':00' for t in times]

        return times

    @property
    def topo_file(self):
        return 'era5_ens_to.nc'

    @property
    def input_directory(self):
        return "era5ens"

    def list_downloaders(self, date_i):
        pl = ERA5Epl(date_i, self.area, self.elevation, self.variables, self.directory)
        sa = ERA5Esa(date_i, self.area, self.variables, self.directory)
        sf = ERA5Esf(date_i, self.area, self.variables, self.directory)

        ERAli = [pl, sa, sf]
        return ERAli
