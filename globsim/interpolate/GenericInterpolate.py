import esmpy
import gc
import numpy as np
import netCDF4 as nc
import re
import tomlkit
import warnings
import logging
import sys
import psutil
import xarray as xr

from datetime import datetime, timedelta
from os import path, makedirs
from pathlib import Path

from globsim.common_utils import StationListRead, variables_skip
from globsim.boundingbox import stations_bbox, netcdf_bbox, BoundingBox
from globsim.gap_checker import check_time_integrity
from globsim.decorators import check
from globsim.interpolate.create_grid_helper import clip_grid_to_indices, clipped_grid_indices, get_buffered_slices

logger = logging.getLogger('globsim.interpolate')

import esmpy as ESMF
if logger.level < 10:  # (DEBUG)
    ESMF.Manager(debug=True)

class GenericInterpolate:
    REANALYSIS = ''
    SAFE_MEM_LIMIT_PERCENT = 90

    def __init__(self, ifile: str, **kwargs):
        # read parameter file
        self.ifile = ifile
        with open(self.ifile) as FILE:
            config = tomlkit.parse(FILE.read())
            self.par = par = config.get('interpolate')
        self.output_dir = self.make_output_directory(par)
        self.variables = par.get('variables')
        self.skip_checks = kwargs.get('skip_checks', par.get("skip_checks", False))
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
        self.resume = bool(self.read_and_report(kwargs, 'resume', False))
        self._reordered = bool(self.read_and_report(kwargs, 'reordered', False))
        self.extrapolate_below_grid = bool(self.read_and_report(kwargs, 'extrapolate_below_grid', True))

        # how globsim handles dimension order internally
        self._working_order_3d = (self.vn_time, 'latitude', 'longitude')
        self._working_order_4d = (self.vn_time, self.vn_level, 'latitude', 'longitude')
        self._working_order_5d = (self.vn_time, 'number', self.vn_level, 'latitude', 'longitude')
        
    @property
    def _downloaded_order_3d(self):
        if self._reordered:
            return ('latitude', 'longitude',self.vn_time)
        else:
            return (self.vn_time, 'latitude', 'longitude')
    
    @property
    def _downloaded_order_4d(self):
        if self._reordered:
            return ('latitude', 'longitude', self.vn_level, self.vn_time)
        else:
            return (self.vn_time, self.vn_level, 'latitude', 'longitude')
    
    @property
    def _downloaded_order_5d(self):
        if self._reordered:
            return ('latitude', 'longitude', 'number', self.vn_level, self.vn_time)
        else:
            return (self.vn_time, 'number', self.vn_level, 'latitude', 'longitude')
        
    def r3d(self, arr, slices={}):
        ''' transform 3d array dimensions from on-disk order to working order 
            used to handle cases where reanalysis dimensions change order    
        '''
        working_order = self._working_order_3d
        downloaded_order = self._downloaded_order_3d
        return reorder_and_slice_array(arr, downloaded_order, working_order, slices)

    def r4d(self, arr, slices={}):
        ''' read downloaded file and reorder dimensions '''
        working_order =  self._working_order_4d
        downloaded_order = self._downloaded_order_4d
        return reorder_and_slice_array(arr, downloaded_order, working_order, slices)
    
    def r5d(self, arr, slices={}):
        ''' read downloaded file and reorder dimensions '''
        working_order =  self._working_order_5d
        downloaded_order = self._downloaded_order_5d
        return reorder_and_slice_array(arr, downloaded_order, working_order, slices)
    

    def read_and_report(self, kwargs, name=None, default=None):
        value = kwargs.get(name, "MISSING FROM KWARGS")

        if value == "MISSING FROM KWARGS":
            value = self.par.get(name, "MISSING FROM TOML")
            if value == "MISSING FROM TOML":
                value = default
                setfrom = "DEFAULT"
            else:
                setfrom = "TOML"
        else:
            setfrom = "CLI"
        logger.debug(f"CONFIG ({setfrom}) {name}: {value}")
        return value
        
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

    def getOutFile(self, kind):
        """
        Get the output file name for the given kind of interpolation.
        """
        return path.join(self.output_dir, f'{self.REANALYSIS}_{kind}_{self.list_name}.nc')
    
    def process(self):
        t_start = datetime.now()
        self._preprocess()

        if not self._skip_sa:
            self.require_safe_mem_usage(logging.INFO)
            self._process_sa()
        else:
            logger.info("skipping interpolation of _sa file")

        if not self._skip_sf:
            self.require_safe_mem_usage(logging.INFO)
            self._process_sf()
        else:
            logger.info("skipping interpolation of _sf file")

        if not self._skip_pl:
            self.require_safe_mem_usage(logging.INFO)
            self._process_pl()
        else:
            logger.info("skipping interpolation of _pl file")

        duration = human_readable_time(datetime.now() - t_start)
        text = f"Interpolation complete in {duration[0]} Days, {duration[1]} Hours, {duration[2]}:{duration[3]}"
        logger.info(text)

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
                 variables=None, date=None) -> "tuple[ESMF.Field, list[str]]":
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
        #import pdb;pdb.set_trace()
        #tbeg = nc.num2date(ncf_in[self.vn_time][np.where(tmask_chunk)[0][0]], ncf_in[self.vn_time].units, ncf_in[self.vn_time].calendar)
        #tend = nc.num2date(ncf_in[self.vn_time][np.where(tmask_chunk)[0][-1]],  ncf_in[self.vn_time].units, ncf_in[self.vn_time].calendar)
        #logger.info(f"2d interpolation for period {tbeg} to {tend}")

        # is it a file with pressure levels?
        pl = self.vn_level in ncf_in.sizes.keys()
        ens = 'number' in ncf_in.sizes.keys()

        # get spatial dimensions
        if pl:  # only for pressure level files
            nlev = len(ncf_in.variables[self.vn_level][:])
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
        varlist = [x for x in ncf_in.variables.keys()]
        self.remove_select_variables(varlist, pl, ens=False)

        # list variables that should be interpolated
        if variables is None:
            variables = varlist
        # test is variables given are available in file
        if (set(variables) < set(varlist) == 0):
            raise ValueError('One or more variables not in netCDF file.')

        # create source field(s) on source grid
        sfield = self.create_source_field(sgrid, variables, nt, ncf_in, pl)
        self.nc_data_subset_to_source_field(variables, sfield, ncf_in, 
                                            tmask_chunk, pl, lon_subset, lat_subset)

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
        
        # clean up and GC
        del sfield, locstream, ncf_in
        gc.collect()

        logger.debug("Created destination field")

        return dfield, variables

    def create_source_field(self,
                            sgrid: ESMF.Grid, 
                            variables, 
                            nt,
                            ncf_in,
                            pl:bool):

        if pl:  # only for pressure level files
            nlev = ncf_in.variables[self.vn_level].shape[0]
            sfield = create_field(sgrid, variables, nt, nlev)
        else:  # 2D files
            sfield = create_field(sgrid, variables, nt)

        return sfield

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

        #sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)
        
        # sg = ESMF.Grid(max_index = np.array([205,45]), num_peri_dims=1, periodic_dim=0,staggerloc=ESMF.StaggerLoc.CENTER)

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
    
    def write_dfield_to_file(self, 
                             dfield, 
                             variables:list,
                             ncf_out,
                             beg:int,
                             end:int,
                             pl:bool):

        for i, var in enumerate(variables):
            if variables_skip(var):
                continue
                        
            else:
                if pl:
                    vi = dfield.data[:,i,:,:].transpose((1,2,0))
                    ncf_out.variables[var][beg:end + 1,:,:] = vi
                else:
                    vi = dfield.data[:,i,:].transpose((1,0))
                    ncf_out.variables[var][beg:end + 1,:] = vi

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

    def nc_data_to_source_field(self, variables, sfield: "ESMF.Field", ncf_in,
                                tmask_chunk, pl: bool):
        # assign data from ncdf: (time, [level], latitude, longitude)
        # to sfield (longitude, latitude, variable, time, [level])
        tmin = np.min(np.where(tmask_chunk))
        tmax = np.max(np.where(tmask_chunk))
        for n, var in enumerate(variables):
            if pl:
                vi = self.r4d(ncf_in[var], {'time':slice(tmin,tmax+1)})
                sfield.data[:,:,n,:,:] = vi.transpose((3,2,0,1))
            else:
                vi = self.r3d(ncf_in[var], {'time':slice(tmin,tmax+1)})
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

        logger.debug("Reading source data from netcdf file")
        t0 = datetime.now()

        for n, v in enumerate(variables):
            var = ncf_in[v]

            if pl:
                logger.debug(f"Reading {v} data from source.")
                if self._plarray.shape[0] == dtime:
                    self._plarray[:] = self.r4d(var, {self.vn_time:time_slice, 'latitude':lat_slice, 'longitude':lon_slice}).values
                else:
                    self._plarray = self.r4d(var, {self.vn_time:time_slice, 'latitude':lat_slice, 'longitude':lon_slice}).values
                # logger.debug(f"Chunksize {(self._plarray.size * self._plarray.itemsize * 1e-6)} Megabytes")
                sfield.data[:, :, n, :, :] = self._plarray.transpose((3,2,0,1))

            else:
                logger.debug(f"Reading {v} data from source")
                if self._array.shape == (dtime, dlat, dlon):
                    self._array[:] = self.r3d(var, {self.vn_time:time_slice, 'latitude':lat_slice, 'longitude':lon_slice}).values
                else:
                    self._array = self.r3d(var, {self.vn_time:time_slice, 'latitude':lat_slice, 'longitude':lon_slice}).values
                
                # logger.debug(f"Chunksize {(self._plarray.size * self._plarray.itemsize * 1e-6)} Megabytes")
                
                sfield.data[:, :, n, :] = self._array.transpose((2,1,0))

        filltime = (datetime.now() - t0).total_seconds()
        logger.debug(f"Finished reading source data for chunk ({filltime} seconds)")

    def remove_select_variables(self, varlist: list, pl: bool, ens: bool = False):
        varlist.remove(self.vn_time)
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
            varlist.remove(self.vn_level)
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

    def prefilter_mf_paths(self, pattern):
        return prefilter_mf_paths(pattern, self.date['beg'], self.date['end'])
    
    def completed_successfully(self, file:str) -> bool:
        """ Check if interpolation was successful """
        if not Path(file).is_file():
            return False
        with nc.Dataset(file) as ncfile:
            if 'globsim_interpolate_success' in ncfile.ncattrs():
                return bool(int(ncfile.getncattr('globsim_interpolate_success')))
            else:
                return False
            
    def require_safe_mem_usage(self, level=logging.DEBUG):
        """ Check if memory usage is safe. Kill globsim if not """
        mem = psutil.virtual_memory()
        logger.log(level, f"Memory usage: {mem.used / 1024**3:.2f} GB ({mem.percent}%)")
        safe = mem.percent < self.SAFE_MEM_LIMIT_PERCENT
        if not safe:
            logger.critical(f"Memory use exceeds safe limit of {self.SAFE_MEM_LIMIT_PERCENT}%. Exiting safely")
            sys.exit(1)

    def file_can_be_resumed(self, file:str, fail_on_missing=False):
        ''' check if file can be resumed '''
        errmsgs = []
    
        if not (Path(file).is_file() or fail_on_missing):
            return True
        
        with nc.Dataset(file) as ncfile:
            if not 'globsim_last_chunk_written' in ncfile.ncattrs():
                errmsgs.append("No 'globsim_last_chunk_written' attribute found in file.")
            if not 'globsim_chunk_size' in ncfile.ncattrs():
                errmsgs.append("No 'globsim_chunk_size' attribute found in file.")
            else:
                if ncfile.getncattr('globsim_chunk_size') != self.cs:
                    errmsgs.append(f"Chunk size in file ({ncfile.getncattr('globsim_chunk_size')}) does not match configuration ({self.cs}).")
            if not 'globsim_interpolate_start' in ncfile.ncattrs():
                errmsgs.append("No 'globsim_interpolate_start' attribute found in file.")
            else:
                if ncfile.getncattr('globsim_interpolate_start') != self.par['beg']:
                    errmsgs.append(f"Interpolation start in file ({ncfile.getncattr('globsim_interpolate_start')}) does not match configuration ({self.date['beg']}).")
            if not 'globsim_interpolate_end' in ncfile.ncattrs():
                errmsgs.append("No globsim_interpolate_end attribute found in file.")
            else:
                if ncfile.getncattr('globsim_interpolate_end') != self.par['end']:
                    errmsgs.append(f"Interpolation end in file ({ncfile.getncattr('globsim_interpolate_end')}) does not match configuration ({self.date['end']}).")
        
        if errmsgs:
            logger.error(f"Cannot resume interpolation from file {file}:")
            for msg in errmsgs:
                logger.error(msg)
            return False
        else:
            return True
    
    def require_file_can_be_resumed(self, file:str, fail_on_missing:bool = False):

        if not self.file_can_be_resumed(file, fail_on_missing=fail_on_missing):
            logger.critical(f"Cannot resume interpolation from file {file}. Exiting safely.")
            sys.exit(1)
        

def prefilter_mf_paths(pattern:str, beg: "datetime", end: "datetime") -> list:
    """ Filter out files that do not fall within date range """ 
    pattern = Path(pattern)
    files = list(Path(pattern.parent).glob(pattern.name))
    dates = [re.search(r'(\d{8})_to_(\d{8})', str(f)).groups() for f in files]
    filtered_files = []
    for file, (b, e) in zip(files, dates):
        filebeg = datetime.strptime(b, '%Y%m%d')

        try:
            fileend = datetime.strptime(e, '%Y%m%d')
        except ValueError:  # end on 31st
            try: 
                fileend = datetime.strptime(e[:-2] + "30", '%Y%m%d')
            except ValueError:
                try:
                    fileend = datetime.strptime(e[:-2] + "29", '%Y%m%d')
                except ValueError:
                    fileend = datetime.strptime(e[:-2] + "28", '%Y%m%d')


        if fileend < beg or filebeg > end:
            pass
        else:
            filtered_files.append(file)

    return sorted(filtered_files)

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


def reorder_and_slice_array(arr, array_order: tuple, desired_order: tuple, slice_boundaries: dict) -> np.ndarray:
    """
    Reorders and slices a NumPy or xarray array based on the provided dimension order and slice boundaries.
    
    Parameters:
        arr (np.ndarray or xr.DataArray): Input array (can be a NumPy array or an xarray DataArray).
        array_order (tuple): Tuple indicating the current order of dimensions (e.g., ('time', 'latitude', 'longitude')).
        desired_order (tuple): Tuple indicating the desired order of dimensions.
        slice_boundaries (dict): Dictionary with dimension names as keys and corresponding slice boundaries as values.
            For example: {'time': slice(1, 3), 'latitude': slice(None), 'longitude': slice(0, 5)}.

    Returns:
        np.ndarray: Reordered and sliced array (either numpy ndarray or xarray DataArray).
    """
    
    # Create a mapping from dimension name to index based on the array's current order
    order_mapping = {dim: idx for idx, dim in enumerate(array_order)}
    
    slices = []
    for dim in array_order:
        boundary = slice_boundaries.get(dim, slice(None))
        
        # If the boundary is an integer, convert it to a slice
        if isinstance(boundary, int):
            slices.append(slice(boundary, boundary + 1))  # Treat integer as a slice of length 1
        else:
            slices.append(boundary)  # Otherwise, use the boundary as is (could be a slice or None)
    
    # Apply the slices to the array (slicing first)
    sliced_arr = arr[tuple(slices)]  # This applies the slices
    
    # Get the new indices based on the desired order
    new_axes = [order_mapping[dim] for dim in desired_order]
    
    # Reorder the dimensions of the sliced array
    if isinstance(arr, xr.DataArray):
        # For xarray, use .transpose to reorder dimensions
        reordered_arr = sliced_arr.transpose(*desired_order)
    else:
        # For NumPy, use np.transpose to reorder dimensions
        reordered_arr = np.transpose(sliced_arr, axes=new_axes)

    return reordered_arr

def human_readable_time(delta: timedelta) -> tuple:
    """
    Convert a timedelta object into a tuple of days, hours, minutes, and seconds.
    
    Parameters:
        delta (timedelta): A timedelta object representing a duration of time.
    
    Returns:
        tuple: A tuple of integers representing the number of days, hours, minutes, and seconds in the timedelta.
    """
    # Calculate the total number of seconds in the timedelta
    total_seconds = delta.total_seconds()
    
    # Calculate the number of days, hours, minutes, and seconds
    days = int(total_seconds // (24 * 3600))
    hours = int((total_seconds % (24 * 3600)) // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = int(total_seconds % 60)
    
    return days, hours, minutes, seconds
