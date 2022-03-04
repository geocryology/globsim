import numpy as np
import re
import tomlkit
import warnings
import logging

from datetime import datetime, timedelta
from os import path, makedirs, listdir
from pathlib import Path
from fnmatch import filter as fnmatch_filter

from globsim.common_utils import StationListRead, str_encode

logger = logging.getLogger('globsim.interpolate')

try:
    import ESMF

    # Check ESMF version.  7.0.1 behaves differently than 7.1.0r
    ESMFv = int(re.sub("[^0-9]", "", ESMF.__version__))
    ESMFnew = ESMFv > 701

except ImportError:
    print("*** ESMF not imported, interpolation not possible. ***")
    pass


class GenericInterpolate:

    def __init__(self, ifile):
        # read parameter file
        self.ifile = ifile
        with open(self.ifile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config['interpolate']
        self.output_dir = self.make_output_directory(par)
        self.variables = par['variables']
        self.list_name = path.basename(path.normpath(par['station_list'])).split(path.extsep)[0]
        
        # read station points
        self.stations_csv = self.find_stations_csv(par)
        self.stations = StationListRead(self.stations_csv)

        # time bounds, add one day to par['end'] to include entire last day
        self.date = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                     'end': datetime.strptime(par['end'], '%Y/%m/%d') + timedelta(days=1)}

        # chunk size: how many time steps to interpolate at the same time?
        # A small chunk size keeps memory usage down but is slow.
        self.cs = int(par['chunk_size'])

    def find_stations_csv(self, par):
        if Path(par['station_list']).is_file():
            return Path(par['station_list'])
        
        elif Path(par['project_directory'], 'par', par['station_list']).is_file():
            return Path(par['project_directory'], 'par', par['station_list'])
        
        else:
            raise FileNotFoundError(f"Siteslist file {par['station_list']} not found.")

    def _set_input_directory(self, name):
        self.input_dir = path.join(self.par['project_directory'], name)

    def TranslateCF2short(self, dpar):
        """
        Map CF Standard Names into short codes used in netCDF files.
        """
        varlist = []
        for var in self.variables:
            varlist.append(dpar.get(var))
        # drop none
        varlist = [item for item in varlist if item is not None]
        # flatten
        varlist = [item for sublist in varlist for item in sublist]
        return(varlist)

    def interp2D(self, ncfile_in, ncf_in, points, tmask_chunk,
                 variables=None, date=None):
        """
        Bilinear interpolation from fields on regular grid (latitude, longitude)
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files

        Args:
            ncfile_in: Full path to an Era-Interim derived netCDF file. This can
                       contain wildcards to point to multiple files if temporal
                       chunking was used.

            ncf_in: A netCDF4.MFDataset derived from reading in Era-Interim
                    multiple files (def ERA2station())

            points: A dictionary of locations. See method StationListRead in
                    common_utils.py for more details.

            tmask_chunk:

            variables:  List of variable(s) to interpolate such as
                        ['r', 't', 'u','v', 't2m', 'u10', 'v10', 'ssrd', 'strd', 'tp'].
                        Defaults to using all variables available.

            date: Directory to specify begin and end time for the derived time
                  series. Defaluts to using all times available in ncfile_in.

        Example:
            from datetime import datetime
            date  = {'beg' : datetime(2008, 1, 1),
                     'end' : datetime(2008,12,31)}
            variables  = ['t','u', 'v']
            stations = StationListRead("points.csv")
            ERA2station('era_sa.nc', 'era_sa_inter.nc', stations,
                        variables=variables, date=date)
        """
        logger.debug("Starting 2d interpolation")
        
        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()
        ens = 'number' in ncf_in.dimensions.keys()

        # get spatial dimensions
        if pl:  # only for pressure level files
            lev = ncf_in.variables['level'][:]
            nlev = len(lev)
        if ens:
            num = ncf_in.variables['number'][:]

        # test if time steps to interpolate remain
        nt = sum(tmask_chunk)
        if nt == 0:
            raise ValueError('No time steps from netCDF file selected.')

        # get variables
        varlist = [str_encode(x) for x in ncf_in.variables.keys()]
        self.remove_select_variables(varlist, pl, ens=False)

        # list variables that should be interpolated
        if variables is None:
            variables = varlist
        # test is variables given are available in file
        if (set(variables) < set(varlist) == 0):
            raise ValueError('One or more variables not in netCDF file.')

        sgrid = self.create_source_grid(ncfile_in)

        # create source field on source grid
        if ens:
            sfield = []
            for ni in num:
                if pl:  # only for pressure level files
                    sfield.append(
                        ESMF.Field(sgrid, name='sgrid',
                                   staggerloc=ESMF.StaggerLoc.CENTER,
                                   ndbounds=[len(variables), nt, nlev]))
                else:  # 2D files
                    sfield.append(
                        ESMF.Field(sgrid, name='sgrid',
                                   staggerloc=ESMF.StaggerLoc.CENTER,
                                   ndbounds=[len(variables), nt]))

        else:
            if pl:  # only for pressure level files
                sfield = ESMF.Field(sgrid, name='sgrid',
                                    staggerloc=ESMF.StaggerLoc.CENTER,
                                    ndbounds=[len(variables), nt, nlev])
            else:  # 2D files
                sfield = ESMF.Field(sgrid, name='sgrid',
                                    staggerloc=ESMF.StaggerLoc.CENTER,
                                    ndbounds=[len(variables), nt])

        self.nc_data_to_source_field(variables, sfield, ncf_in,
                                     tmask_chunk, pl, ens)

        locstream = self.create_loc_stream()

        # create destination field
        if ens:
            dfield = []
            for ni in num:
                if pl:  # only for pressure level files
                    di = ESMF.Field(locstream, name='dfield',
                                    ndbounds=[len(variables), nt, nlev])
                else:
                    di = ESMF.Field(locstream, name='dfield',
                                    ndbounds=[len(variables), nt])
                dfield.append(self.regrid(sfield[ni], di))
        else:
            if pl:  # only for pressure level files
                dfield = ESMF.Field(locstream, name='dfield',
                                    ndbounds=[len(variables), nt, nlev])
            else:
                dfield = ESMF.Field(locstream, name='dfield',
                                    ndbounds=[len(variables), nt])

            dfield = self.regrid(sfield, dfield)
        
        logger.debug("Created destination field")

        return dfield, variables

    def make_output_directory(self, par):
        """make directory to hold outputs"""
        output_root = None
        
        if par.get('output_directory'):
            try:
                test_path = Path(par.get('output_directory'))
            except TypeError:
                warnings.warn("You provided an output_directory for interpolation that was not understood. Saving files to project directory.")
            
            if test_path.is_dir():
                output_root = Path(test_path, "interpolated")
            else:
                warnings.warn("You provided an output_directory for interpolation that does not exist. Saving files to project directory.")
        
        if not output_root:
            output_root = path.join(par['project_directory'], 'interpolated')

        if not Path(output_root).is_dir():
            makedirs(output_root)

        return output_root

    def create_source_grid(self, ncfile_in):
        # Create source grid from a SCRIP formatted file. As ESMF needs one
        # file rather than an MFDataset, give first file in directory.
        flist = np.sort(fnmatch_filter(listdir(self.input_dir),
                                       path.basename(ncfile_in)))
        ncsingle = path.join(self.input_dir, flist[0])
        sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)
        return sgrid

    def create_loc_stream(self):
        # CANNOT have third dimension!!!
        locstream = ESMF.LocStream(len(self.stations),
                                   coord_sys=ESMF.CoordSys.SPH_DEG)
        locstream["ESMF:Lon"] = list(self.stations['longitude_dd'])
        locstream["ESMF:Lat"] = list(self.stations['latitude_dd'])

        return locstream

    @staticmethod
    def regrid(sfield, dfield):
        # regridding function, consider ESMF.UnmappedAction.ERROR
        
        regrid2D = ESMF.Regrid(sfield, dfield,
                               regrid_method=ESMF.RegridMethod.BILINEAR,
                               unmapped_action=ESMF.UnmappedAction.IGNORE,
                               dst_mask_values=None)
        
        # regrid operation, create destination field (variables, times, points)
        logger.debug("Attempting to regrid")
        dfield = regrid2D(sfield, dfield)
        logger.debug("Regridding complete")

        sfield.destroy()  # free memory

        return dfield

    @staticmethod
    def nc_data_to_source_field(variables, sfield, ncf_in,
                                tmask_chunk, pl, ens):
        # assign data from ncdf: (variable, time, latitude, longitude)

        for n, var in enumerate(variables):
            if ens:
                for ni in ncf_in['number'][:]:
                    if pl:
                        # sfield.data [lon, lat, var, time, number, level]
                        # vi [time, number, level, lat, lon]
                        vi = ncf_in[var][tmask_chunk,ni,:,:,:]
                        sfield[ni].data[:,:,n,:,:] = vi.transpose((3,2,0,1))
                        logger.debug(f"Wrote {var} data to source field for regridding")
                    else:
                        vi = ncf_in[var][tmask_chunk,ni,:,:]
                        sfield[ni].data[:,:,n,:] = vi.transpose((2,1,0))
                        logger.debug(f"Wrote {var} data to source field for regridding")
            else:
                if pl:
                    vi = ncf_in[var][tmask_chunk,:,:,:]
                    sfield.data[:,:,n,:,:] = vi.transpose((3,2,0,1))
                    logger.debug(f"Wrote {var} data to source field for regridding")
                else:
                    vi = ncf_in[var][tmask_chunk,:,:]
                    sfield.data[:,:,n,:] = vi.transpose((2,1,0))
                    logger.debug(f"Wrote {var} data to source field for regridding")

    @staticmethod
    def remove_select_variables(varlist, pl, ens=False):
        varlist.remove('time')
        try:
            varlist.remove('latitude')
        except ValueError as e:
            print(e)
            print('continue')
            pass
        try:
            varlist.remove('longitude')
        except ValueError as e:
            print(e)
            print('continue')
            pass
        if pl:  # only for pressure level files
            varlist.remove('level')
        if ens:
            varlist.remove('number')

    @staticmethod
    def ele_interpolate(elevation, h, nl):
        # difference in elevation, level directly above will be >= 0
        elev_diff = elevation - h
        # vector of level indices that fall directly above station.
        # Apply after ravel() of data.
        va = np.argmin(elev_diff + (elev_diff < 0) * 100000, axis=1)
        # mask for situations where station is below lowest level
        mask = va < (nl - 1)
        va += np.arange(elevation.shape[0]) * elevation.shape[1]

        # Vector level indices that fall directly below station.
        # Apply after ravel() of data.
        vb = va + mask  # +1 when OK, +0 when below lowest level

        return elev_diff, va, vb

    @staticmethod
    def calculate_weights(elev_diff, va, vb):
        wa = np.absolute(elev_diff.ravel()[vb])
        wb = np.absolute(elev_diff.ravel()[va])
        wt = wa + wb
        wa /= wt  # Apply after ravel() of data.
        wb /= wt  # Apply after ravel() of data.

        return wa, wb