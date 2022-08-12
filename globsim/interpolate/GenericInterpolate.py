import numpy as np
import netCDF4 as nc
import re
import tomlkit
import warnings
import logging

from datetime import datetime, timedelta
from os import path, makedirs, listdir
from pathlib import Path
from fnmatch import filter as fnmatch_filter
from typing import Union

from globsim.common_utils import StationListRead, str_encode
from globsim.boundingbox import stations_bbox, netcdf_bbox, BoundingBox
from globsim.gap_checker import check_time_integrity
from globsim.decorators import check

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

    def __init__(self, ifile: str):
        # read parameter file
        self.ifile = ifile
        with open(self.ifile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config.get('interpolate')
        self.output_dir = self.make_output_directory(par)
        self.variables = par.get('variables')
        self.skip_checks = bool(par.get("skip_checks", False))
        self.list_name = path.basename(path.normpath(par.get('station_list'))).split(path.extsep)[0]
        
        # read station points
        self.stations_csv = self.find_stations_csv(par)
        self.stations = StationListRead(self.stations_csv)
        self.stations_bbox = stations_bbox(self.stations)

        # time bounds, add one day to par['end'] to include entire last day
        self.date = {'beg': datetime.strptime(par.get('beg'), '%Y/%m/%d'),
                     'end': datetime.strptime(par.get('end'), '%Y/%m/%d') + timedelta(days=1)}

        # chunk size: how many time steps to interpolate at the same time?
        # A small chunk size keeps memory usage down but is slow.
        self.cs = int(par.get('chunk_size'))

    def find_stations_csv(self, par):
        if Path(par.get('station_list')).is_file():
            return Path(par.get('station_list'))

        elif Path(par.get('project_directory'), 'par', par.get('station_list')).is_file():
            return Path(par.get('project_directory'), 'par', par.get('station_list'))

        else:
            raise FileNotFoundError(f"Siteslist file {par.get('station_list')} not found.")

    @check
    def validate_stations_extent(self, ncdf):
        try:
            if not netcdf_bbox(ncdf).contains_bbox(self.stations_bbox):
                logger.error("Station coordinates exceed downloaded extent")
                raise ValueError("Station coordinates exceed downloaded extent")
            else:
                logger.info("Stations within bounding box of dataset")
        except Exception:
            logger.error("Could not verify whether stations are within downloaded netcdf")

    def _set_input_directory(self, name):
        self.input_dir = path.join(self.par.get('project_directory'), name)

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

    def interp2D(self,  ncf_in, points, tmask_chunk: "np.ndarray",
                 variables=None, date=None):
        """
        Bilinear interpolation from fields on regular grid (latitude, longitude)
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files

        Args:
            ncfile_in: nc Dataset  OR Full path to an Era-Interim derived netCDF file. This can
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

        logger.debug(f"Starting 2d interpolation for chunks {np.min(np.where(tmask_chunk == True))} to {np.max(np.where(tmask_chunk == True))} of {len(tmask_chunk)} ")

        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()
        ens = 'number' in ncf_in.dimensions.keys()

        # get spatial dimensions
        if pl:  # only for pressure level files
            nlev = len(ncf_in.variables['level'][:])
        else:
            nlev = 1
        
        if ens:
            num = ncf_in.variables['number'][:]
        else:
            num = []

        # test if time steps to interpolate remain
        nt = tmask_chunk.sum()  # TODO: could this just be length? what is being tested here?
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

        sgrid = self.create_source_grid(ncf_in)

        # create source field(s) on source grid
        if ens:
            sfield = []
            for ni in num:
                if pl:  # only for pressure level files
                    sfield.append(create_field(sgrid, variables, nt, nlev))
                else:  # 2D files
                    sfield.append(create_field(sgrid, variables, nt))

            self.nc_ensemble_data_to_source_field(variables, sfield, ncf_in, tmask_chunk, pl)

        else:
            if pl:  # only for pressure level files
                sfield = create_field(sgrid, variables, nt, nlev)
            else:  # 2D files
                sfield = create_field(sgrid, variables, nt)
            
            #self.nc_data_subset_to_source_field(variables, sfield, ncf_in, tmask_chunk, pl)
            self.nc_data_to_source_field(variables, sfield, ncf_in, tmask_chunk, pl)

        locstream = self.create_loc_stream(points)

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

    def make_output_directory(self, par) -> str:
        """make directory to hold outputs"""
        output_root = None

        if par.get('output_directory'):
            try:
                test_path = Path(par.get('output_directory'))

                if test_path.is_dir():
                    output_root = Path(test_path, "interpolated")
                else:
                    warnings.warn("You provided an output_directory for interpolation that does not exist. Saving files to project directory.")

            except TypeError:
                warnings.warn("You provided an output_directory for interpolation that was not understood. Saving files to project directory.")

        if not output_root:
            output_root = path.join(par.get('project_directory'), 'interpolated')

        if not Path(output_root).is_dir():
            makedirs(output_root)

        return str(output_root)

    def create_source_grid(self, ncf_in: "nc.MFDataset") -> "ESMF.Grid":
        # Create source grid from a SCRIP formatted file. As ESMF needs one
        # file rather than an MFDataset, give first file in directory.
        #flist = np.sort(fnmatch_filter(listdir(self.input_dir),
        #                               path.basename(ncfile_in)))
        #ncsingle = path.join(self.input_dir, flist[0])
        ncsingle = ncf_in._files[0]
        sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)
        return sgrid

    def create_loc_stream(self, points) -> "ESMF.LocStream":
        # CANNOT have third dimension!!!
        locstream = ESMF.LocStream(len(points),
                                   coord_sys=ESMF.CoordSys.SPH_DEG)
        locstream["ESMF:Lon"] = list(points['longitude_dd'])
        locstream["ESMF:Lat"] = list(points['latitude_dd'])

        return locstream

    @staticmethod
    def regrid(sfield: "ESMF.Field", dfield: "ESMF.Field") -> "ESMF.Field":
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
    def nc_ensemble_data_to_source_field(variables, sfield_list: "list[ESMF.Field]", ncf_in,
                                         tmask_chunk, pl: bool):
        for n, var in enumerate(variables):
            for ni in ncf_in['number'][:]:
                if pl:
                    # sfield.data [lon, lat, var, time, number, level]
                    # vi [time, number, level, lat, lon]
                    vi = ncf_in[var][tmask_chunk,ni,:,:,:]
                    sfield_list[ni].data[:,:,n,:,:] = vi.transpose((3,2,0,1))

                else:
                    vi = ncf_in[var][tmask_chunk,ni,:,:]
                    sfield_list[ni].data[:,:,n,:] = vi.transpose((2,1,0))

            logger.debug(f"Wrote {var} data to source field for regridding")

    @staticmethod
    def nc_data_to_source_field(variables, sfield: "ESMF.Field", ncf_in,
                                tmask_chunk, pl: bool):
        # assign data from ncdf: (time, [level], latitude, longitude)
        # to sfield (longitude, latitude, variable, time, [level])
        tmin = np.min(np.where(tmask_chunk))
        tmax = np.max(np.where(tmask_chunk))
        for n, var in enumerate(variables):
            if pl:
                vi = ncf_in[var][tmin:tmax + 1,:,:,:]
                sfield.data[:,:,n,:,:] = vi.transpose((3,2,0,1))
            else:
                vi = ncf_in[var][tmin:tmax + 1,:,:]
                sfield.data[:,:,n,:] = vi.transpose((2,1,0))

            logger.debug(f"Wrote {var} data to source field for regridding")

    def nc_data_subset_to_source_field(self, variables, sfield: "ESMF.Field", ncf_in,
                                       tmask_chunk, pl: bool):
        """ assign data from ncdf: (time, [level], latitude, longitude)
                         to sfield (longitude, latitude, variable, time, [level])
        """

        # get indices of lat/lon
        lat = ncf_in.get_variables_by_attributes(standard_name='latitude', axis='Y')[0][:]
        lon = ncf_in.get_variables_by_attributes(standard_name='longitude', axis='X')[0][:]

        lat_ix = np.where((lat >= self.stations_bbox.ymin) & (lat <= self.stations_bbox.ymax))[0]
        lon_ix = np.where((lon >= self.stations_bbox.xmin) & (lon <= self.stations_bbox.xmax))[0]

        lat_slice = slice(min(lat_ix), max(lat_ix) + 1)
        lon_slice = slice(min(lon_ix), max(lon_ix) + 1)

        tmin = np.min(np.where(tmask_chunk))
        tmax = np.max(np.where(tmask_chunk))
        time_slice = slice(tmin, tmax + 1)

        for n, var in enumerate(variables):
            if pl:
                vi = ncf_in[var][time_slice, :, lat_slice, lon_slice]
                sfield.data[lon_slice, lat_slice, n, :, :] = vi.transpose((3,2,0,1))
            else:
                vi = ncf_in[var][time_slice, lat_slice, lon_slice]
                sfield.data[lon_slice, lat_slice, n, :] = vi.transpose((2,1,0))

            logger.debug(f"Wrote {var} data to source field for regridding")

    @staticmethod
    def remove_select_variables(varlist: list, pl: bool, ens: bool = False):
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

    @check
    def ensure_datset_integrity(self, time: "nc.Variable", interval: float):
        """ Perform basic pre-flight checks on (downloaded) input datasets before running interpolate"""
        # Check coverage
        interpolate_start = self.date['beg']
        interpolate_end = self.date['end'] - timedelta(days=1)  # take off last day added in __init__
        data_start = nc.num2date(time[0], time.units, time.calendar)
        data_end = nc.num2date(time[-1], time.units, time.calendar)

        if not (data_start <= interpolate_start):
            raise ValueError(f"Requested interpolation start ({interpolate_start}) is before data bounds ({data_start})")

        elif not (data_end >= interpolate_end):
            raise ValueError(f"Requested interpolation end ({interpolate_end}) is after data bounds ({data_end})")

        # Check gaps
        critical = False
        for gap_i, gap_start, gap_end in zip(*check_time_integrity(time, interval)):
            message = f"Data gap found at index {gap_i}. Missing data from {gap_start} to {gap_end}."
            
            if interpolate_start <= gap_start <= gap_end:
                critical=True
                logger.critical(message)
            else:
                logger.warning(message)
        
        if critical:
            raise ValueError("Data gaps found within interpolation bounds.")


def create_stations_bbox(stations) -> BoundingBox:
    # get max/min of points lat/lon from self.stations

    stations_bbox = BoundingBox(xmin=stations['longitude_dd'].describe()["min"],
                                xmax=stations['longitude_dd'].describe()["max"],
                                ymin=stations['latitude_dd'].describe()["min"],
                                ymax=stations['latitude_dd'].describe()["max"])
    # add generous buffer
    buffer = 2.0  # assume degrees here
    return BoundingBox(stations_bbox.xmin - buffer,
                       stations_bbox.xmax + buffer,
                       stations_bbox.ymin - buffer,
                       stations_bbox.ymax + buffer)


def grid_from_bbox(latitudes: np.ndarray,
                   longitudes: np.ndarray,
                   bbox: BoundingBox) -> ESMF.Grid:

    valid_lat = latitudes[np.where((latitudes >= bbox.ymin) & (latitudes <= bbox.ymax))[0]]
    valid_lon = longitudes[np.where((longitudes >= bbox.xmin) & (longitudes <= bbox.xmax))[0]]

    grid = ESMF.Grid(max_index=np.array([len(valid_lon), len(valid_lat)]))
    grid.coords[0][0] = np.repeat(valid_lon[np.newaxis, :], len(valid_lat), axis=1).ravel().reshape(len(valid_lon),len(valid_lat))
    grid.coords[0][1] = np.repeat(valid_lat[np.newaxis, :], len(valid_lon), axis=0)
    
    return grid


def create_field(sgrid: "ESMF.Grid", variables: list, nt: int, nlev:int = 1) -> "ESMF.Field":
    nvar = len(variables)
    try:
        if nlev > 1:
            field = ESMF.Field(sgrid, name='sgrid',
                               staggerloc=ESMF.StaggerLoc.CENTER,
                               ndbounds=[nvar, nt, nlev])
        else:
            field = ESMF.Field(sgrid, name='sgrid',
                               staggerloc=ESMF.StaggerLoc.CENTER,
                               ndbounds=[nvar, nt])
    except TypeError as e:
        msg = "Tried to create a ESMF.Field that was too big. Try reducing the chunk_size in your configuration."
        logger.error(f"{msg} Currently there are {nt} time-steps (chunk size) {nvar} variables and {nlev} levels on a {sgrid.size[0][0]}-by-{sgrid.size[0][1]} grid (total size of {nt * nvar * nlev * sgrid.size[0][0] * sgrid.size[0][1]})")
        raise Exception(msg).with_traceback(e.__traceback__)

    return field
