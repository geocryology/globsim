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
from globsim.interpolate.create_grid_helper import clip_grid_to_indices, clipped_grid_indices, get_buffered_slices

logger = logging.getLogger('globsim.interpolate')

try:
    import ESMF

    # Check ESMF version.  7.0.1 behaves differently than 7.1.0r
    ESMFv = int(re.sub("[^0-9]", "", ESMF.__version__))
    ESMFnew = ESMFv > 701

except ModuleNotFoundError:
        print("*** ESMF not imported, trying esmpy. ***")
        try:
            import esmpy as ESMF
        except ImportError:
            print('Could not import ESMF or esmpy')
            pass

class GenericInterpolate:

    def __init__(self, ifile: str, **kwargs):
        # read parameter file
        self.ifile = ifile
        with open(self.ifile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config.get('interpolate')
        self.output_dir = self.make_output_directory(par)
        self.variables = par.get('variables')
        self.skip_checks = kwargs.get('skip_checks', bool(par.get("skip_checks", False)))
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

        self._array = np.array([])  # recycled numpy array
        self._plarray = np.array([])

        self._skip_sa = kwargs.get('skip_sa', False)
        self._skip_sf = kwargs.get('skip_sf', False)
        self._skip_pl = kwargs.get('skip_pl', False)

    @property
    def vn_time(self):
        return 'time'

    @property
    def vn_level(self):
        return 'level'
    
    def find_stations_csv(self, par):
        if Path(par.get('station_list')).is_file():
            return Path(par.get('station_list'))

        elif Path(par.get('project_directory'), 'par', par.get('station_list')).is_file():
            return Path(par.get('project_directory'), 'par', par.get('station_list'))

        else:
            raise FileNotFoundError(f"Siteslist file {par.get('station_list')} not found.")

    @check
    def validate_stations_extent(self, ncdf):
        data_bbox = netcdf_bbox(ncdf)
        msg = f"Station coordinates {self.stations_bbox} exceed downloaded extent {data_bbox}"
        try:
            if not data_bbox.contains_bbox(self.stations_bbox):
                logger.error(msg)
                raise ValueError(msg)
            else:
                logger.info("Stations within bounding box of dataset")
        
        except ValueError:
            raise ValueError(msg)
        
        except KeyError:
            logger.error("Could not verify whether stations are within downloaded netcdf")

    def process(self):
        self._preprocess()

        if not self._skip_sa:
            self._process_sa()
        else:
            logger.info("skipping interpolation of _sa file")

        if not self._skip_sf:
            self._process_sf()
        else:
            logger.info("skipping interpolation of _sf file")

        if not self._skip_pl:
            self._process_pl()
        else:
            logger.info("skipping interpolation of _pl file")

    def _preprocess(self):
        pass

    def _process_sa(self):
        pass

    def _process_sf(self):
        pass

    def _process_pl(self):
        pass

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
                 sgrid: "ESMF.Grid",
                 lon_subset:slice, lat_subset:slice,
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

            tmask_chunk: boolean array of which indices along the time axis to include

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
        tbeg = nc.num2date(ncf_in['time'][np.where(tmask_chunk)[0][0]], ncf_in['time'].units, ncf_in['time'].calendar)
        tend = nc.num2date(ncf_in['time'][np.where(tmask_chunk)[0][-1]],  ncf_in['time'].units, ncf_in['time'].calendar)
        logger.debug(f"2d interpolation for period {tbeg} to {tend}")

        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dims.keys()
        ens = 'number' in ncf_in.dims.keys()

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
            
            self.nc_data_subset_to_source_field(variables, sfield, ncf_in, tmask_chunk, pl, lon_subset, lat_subset)

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
                                    ndbounds=[len(variables), nt, nlev])  # TODO: we can just get this from sfield.ndbounds
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
        # flist = np.sort(fnmatch_filter(listdir(self.input_dir),
        #                               path.basename(ncfile_in)))
        # ncsingle = path.join(self.input_dir, flist[0])
        ncsingle = None
        for v in ncf_in.variables:
            if ncsingle is not None:
                break
            try:
                ncsingle = ncf_in[v].encoding["source"]
            except Exception: 
                pass
            #ncf_in._files[0]
       # import pdb;pdb.set_trace()
        #sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)
        
        # sg = ESMF.Grid(max_index = np.array([205,45]), num_peri_dims=1, periodic_dim=0,staggerloc=ESMF.StaggerLoc.CENTER)
        if True:  # single-level
            template = nc.Dataset(ncsingle)
            lat = template.variables['latitude'][:]
            lon = template.variables['longitude'][:]

        sgrid = grid_create_from_coordinates_periodic(lon, lat)
        return sgrid

    def create_subset_source_grid(self, sgrid: "ESMF.Grid", bbox: BoundingBox) -> "tuple[ESMF.Grid, slice, slice]":
        lon_subset, lat_subset = clipped_grid_indices(sgrid, bbox)
        lon_slice, lat_slice = get_buffered_slices(sgrid, lon_subset, lat_subset)
        subset_grid = clip_grid_to_indices(sgrid, lon_subset, lat_subset)
        
        return subset_grid, lon_slice, lat_slice
    
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
                                       tmask_chunk, pl: bool, lon_slice, lat_slice):
        """ assign data from ncdf: (time, [level], latitude, longitude)
                         to sfield (longitude, latitude, variable, time, [level])
        """
        tmin = np.min(np.where(tmask_chunk))
        tmax = np.max(np.where(tmask_chunk))
        time_slice = slice(tmin, tmax + 1)
        dlon, dlat, dtime = [x.stop - x.start for x in [lon_slice, lat_slice, time_slice]]

        logger.info("Reading source data from netcdf file")
        t0 = datetime.now()

        for n, v in enumerate(variables):
            var = ncf_in[v]

            if pl:
                logger.debug(f"Reading {v} data from source.")  # NB: writing is almost instantaneous
                if self._plarray.shape[0] == dtime:
                    self._plarray[:] = var[time_slice, :, lat_slice, lon_slice].values
                else:
                    self._plarray = var[time_slice, :, lat_slice, lon_slice].values
                # logger.debug(f"Chunksize {(self._plarray.size * self._plarray.itemsize * 1e-6)} Megabytes")
                sfield.data[:, :, n, :, :] = self._plarray.transpose((3,2,0,1))

            else:
                logger.debug(f"Reading {v} data from source")
                try:
                    if self._array.shape == (dtime, dlat, dlon):
                        self._array[:] = var[time_slice, lat_slice, lon_slice].values
                    else:
                        self._array = var[time_slice, lat_slice, lon_slice].values
                except IndexError as e:
                    import pdb;pdb.set_trace()
                # logger.debug(f"Chunksize {(self._plarray.size * self._plarray.itemsize * 1e-6)} Megabytes")
                
                sfield.data[:, :, n, :] = self._array.transpose((2,1,0))

        filltime = (datetime.now() - t0).total_seconds()
        logger.info(f"Finished reading source data for chunk ({filltime} seconds)")

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
        if 'number' in varlist:
            varlist.remove('number')
        if 'expver' in varlist:
            varlist.remove('expver')

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

import esmpy
def grid_create_from_coordinates_periodic(longitudes, latitudes, lon_corners=False, lat_corners=False, corners=False, domask=False):
    """
    Create a 2 dimensional periodic Grid using the 'longitudes' and 'latitudes'.
    :param longitudes: longitude coordinate values at cell centers
    :param latitudes: latitude coordinate values at cell centers
    :param lon_corners: longitude coordinate values at cell corners
    :param lat_corners: latitude coordinate values at cell corners
    :param corners: boolean to determine whether or not to add corner coordinates to this grid
    :param domask: boolean to determine whether to set an arbitrary mask or not
    :return: grid
    """
    [lon, lat] = [0, 1]

    # create a grid given the number of grid cells in each dimension the center stagger location is allocated
    max_index = np.array([len(longitudes), len(latitudes)])
    grid = esmpy.Grid(max_index, num_peri_dims=1, staggerloc=[esmpy.StaggerLoc.CENTER])

    # set the grid coordinates using numpy arrays, parallel case is handled using grid bounds
    gridXCenter = grid.get_coords(lon)
    lon_par = longitudes[grid.lower_bounds[esmpy.StaggerLoc.CENTER][lon]:grid.upper_bounds[esmpy.StaggerLoc.CENTER][lon]]
    gridXCenter[...] = lon_par.reshape((lon_par.size, 1))

    gridYCenter = grid.get_coords(lat)
    lat_par = latitudes[grid.lower_bounds[esmpy.StaggerLoc.CENTER][lat]:grid.upper_bounds[esmpy.StaggerLoc.CENTER][lat]]
    gridYCenter[...] = lat_par.reshape((1, lat_par.size))

    # create grid corners in a slightly different manner to account for the bounds format common in CF-like files
    if corners:
        grid.add_coords([esmpy.StaggerLoc.CORNER])
        lbx = grid.lower_bounds[esmpy.StaggerLoc.CORNER][lon]
        ubx = grid.upper_bounds[esmpy.StaggerLoc.CORNER][lon]
        lby = grid.lower_bounds[esmpy.StaggerLoc.CORNER][lat]
        uby = grid.upper_bounds[esmpy.StaggerLoc.CORNER][lat]

        gridXCorner = grid.get_coords(lon, staggerloc=esmpy.StaggerLoc.CORNER)
        for i0 in range(ubx - lbx - 1):
            gridXCorner[i0, :] = lon_corners[i0+lbx, 0]
        gridXCorner[i0 + 1, :] = lon_corners[i0+lbx, 1]

        gridYCorner = grid.get_coords(lat, staggerloc=esmpy.StaggerLoc.CORNER)
        for i1 in range(uby - lby - 1):
            gridYCorner[:, i1] = lat_corners[i1+lby, 0]
        gridYCorner[:, i1 + 1] = lat_corners[i1+lby, 1]

    # add an arbitrary mask
    if domask:
        mask = grid.add_item(esmpy.GridItem.MASK)
        mask[:] = 1
        mask[np.where((1.75 <= gridXCenter.any() < 2.25) &
                      (1.75 <= gridYCenter.any() < 2.25))] = 0

    return grid