import logging
import netCDF4 as nc
import numpy as np
import tarfile
import xarray as xr

from pathlib import Path
from typing import Optional

from globsim.common_utils import variables_skip
from globsim.download.jra_dict_formatters import lookup_param

logger = logging.getLogger(__name__)


class NcarDownloadHandler:
    REANALYSIS = ''

    def __init__(self):
        self._keep_raw_files = False

    def keep_raw_files(self):
        self._keep_raw_files = True

    def format_output_filename(self, beg, end, filetype):
        return f"{self.REANALYSIS}_{filetype}_{beg.strftime(r'%Y%m%d')}_to_{end.strftime(r'%Y%m%d')}.nc"


class J3QDownloadHandler(NcarDownloadHandler):
    REANALYSIS = 'jra3q'
    DIMS = {'time': 'time', 
            'latitude': 'lat',
            'longitude': 'lon',
            'level': 'pressure_level'}
    DATASET_ID = 'd640000'
    
    def __init__(self):
        super().__init__()
    
    def extract_dimensions(self, nc_template_files):
        '''extract dimensions from the template files'''

        dataset = xr.open_mfdataset(nc_template_files, decode_times=False)
        Times = dataset[self.DIMS['time']]
        Lats = dataset[self.DIMS['latitude']]
        Lons = dataset[self.DIMS['longitude']]

        if self.DIMS['level'] in dataset.dims: 
            Levs = dataset[self.DIMS['level']]
        else:
            Levs = None
        
        dataset.close()

        return Times, Lats, Lons, Levs

    def variable_to_globsim_name(self, param):
        varname = param.split('-')[0]
        d = {'rh2m': 'Relative humidity',
             'spfh2m': 'Specific humidity',
             'spfh': 'Specific humidity',
            'tmp2m': 'Temperature',
            'ugrd10m': 'u-component of wind',
            'vgrd10m': 'v-component of wind',
            'lat': 'latitude',
            'lon': 'longitude',
            'time': 'time',
            'gp' : 'Geopotential',
            'hgt': 'Geopotential height',
            'rh': 'Relative humidity',
            'tmp': 'Temperature',
            'ugrd': 'u-component of wind',
            'vgrd': 'v-component of wind',
            'dswrf1have': 'Downward solar radiation flux',
            'dswrfcs1have': 'Clear sky downward solar radiation flux',
            'dlwrf1have': 'Downward longwave radiation flux',
            'dlwrfcs1have': 'Clear sky downward longwave radiation flux',
            'tprate1have': 'Total precipitation',
            'pres': 'Pressure',
            }
        return d.get(varname)

    def get_dim_names(self, *args, **kwargs):
        return self.DIMS
    
    def get_nc_files(self, directory, variable, dsi) -> list[Path]:
        def filter_func(f):
            return dsi in f.name and variable in f.name
        
        valid_files = filter(filter_func, Path(directory).glob('*nc'))

        return list(valid_files)
    
    def make_globsim_dataset(self, directory:str, request_id:str, filetype:Optional[str]=None):
        '''make a globsim dataset from JRA-3Q dataset'''

        # extract the downloaded tar files
        # get all variables associated with the dataset
        # return the extracted tar files   
        logger.debug(f"Creating {filetype if filetype is not None else '(unknown type)'} dataset from requeset {request_id}")
        
        extract_downloaded_tar_files(directory, request_id, (not self._keep_raw_files))
        
        variables = get_downloaded_variable_names(directory, request_id)
        
        if len(variables) == 0:
            raise ValueError(f"No variables found in dataset {request_id}")
        
        if filetype is None:
            try:
                filetype = determine_output_file_type(directory, request_id)
            except KeyError:
                logger.error(f"Dataset type not provided and not able to guess.  Create data from {request_id} manually.")
                return
            
        dims = self.get_dim_names(directory, request_id)
        
        nc_template_files = self.get_nc_files(directory, variables[0], request_id)
        
        Times, Lats, Lons, Levs = self.extract_dimensions(nc_template_files)
        
        beg, end = nc.num2date(Times.values[[0,-1]], units=Times.units, calendar=Times.calendar)
        file_new = Path(directory, self.format_output_filename(beg, end, filetype))

        if file_new.is_file():
            logger.warning(f"Partially completed file {file_new} already exists. ")
            ncn = nc.Dataset(file_new, 'a')
        else:
            ncn = new_jra_download_file(str(file_new), Times, Lats, Lons, Levs)
        
        for varname in variables:
            if variables_skip(varname):
                logger.debug(f"Skipping variable: {varname}")

            flist = self.get_nc_files(directory, varname, request_id)
            ncf = xr.open_mfdataset(flist, decode_times=False)
            
            output_name = self.variable_to_globsim_name(varname)

            if output_name is None:
                output_name = ''.join(ch for ch in varname if ch.isalnum())
                logger.warning(f"Variable {varname} not found in lookup table.  Using {varname}")

            logger.debug(f"Creating variable: {varname}")
            
            if filetype == 'pl':
                try:
                    vari = ncn.createVariable(output_name, 'f4',
                                                ('time', 'level',
                                                'latitude', 'longitude',))
                    copy_variable_in_chunks(ncf[varname], vari, 4)
                except RuntimeError:
                    logger.error(f"Variable {varname} already exists.")
                    vari = ncn[output_name]
                
            else:
                vari = ncn.createVariable(output_name,'f4',
                                            ('time',
                                            'latitude', 'longitude'))
                copy_variable_in_chunks(ncf[varname], vari, 4)
            
            vari.long_name = ncf[varname].long_name
            vari.units     = ncf[varname].units
            vari.orig_name = varname

            ncf.close()
            if not self._keep_raw_files:
                for f in flist:
                    f.unlink()

        ncn.close()


class J55DownloadHandler(J3QDownloadHandler):
    REANALYSIS = 'jra55'
    DATASET_ID = 'ds628.0'

    def __init__(self):
        super().__init__()
    
    def variable_to_globsim_name(self, param):
        varname = param.split('-')[0]
        return  {
                'initial_time0_hours':   'time',
                'initial_time0':         'time',
                'initial_time0_encoded': 'time',
                'lv_ISBL1':              'level',
                'g0_lat_1':              'latitude',
                'g0_lon_2':              'longitude',
                'g0_lat_2':              'latitude',
                'g0_lon_3':              'longitude',
                'GP_GDS0_SFC':           'Geopotential',
                'HGT_GDS0_ISBL':         'Geopotential height',
                'RH_GDS0_ISBL':          'Relative humidity',
                'VGRD_GDS0_ISBL':        'v-component of wind',
                'UGRD_GDS0_ISBL':        'u-component of wind',
                'TMP_GDS0_ISBL':         'Temperature',
                'TMP_GDS0_HTGL':         'Temperature',
                'VGRD_GDS0_HTGL':        'v-component of wind',
                'UGRD_GDS0_HTGL':        'u-component of wind',
                'RH_GDS0_HTGL':          'Relative humidity',
                'SPFH_GDS0_HTGL':        'Specific humidity',
                'PRES_GDS0_SFC_ave3h':   'Pressure',
                'TPRAT_GDS0_SFC_ave3h':  'Total precipitation',
                'CSDSF_GDS0_SFC_ave3h':  'Clear sky downward solar radiation flux',
                'CSDLF_GDS0_SFC_ave3h':  'Clear sky downward longwave radiation flux',
                'DSWRF_GDS0_SFC_ave3h':  'Downward solar radiation flux',
                'DLWRF_GDS0_SFC_ave3h':  'Downward longwave radiation flux',
                'DSWRF_GDS0_NTAT_ave3h':  'skip',}[varname]
    
    def get_dim_names(self, filetype, *args, **kwargs) -> dict:
    
        if filetype == 'pl':
            lon = 'g0_lon_3'
            lat = 'g0_lat_2'
        elif filetype in ['sa', 'sf', 'to']:
            lon = 'g0_lon_2'
            lat = 'g0_lat_1'
        
        return {'time': 'initial_time0_hours', 
                'latitude': lat,
                'longitude': lon,
                'level': 'lv_ISBL'}


def determine_output_file_type(directory:str, request_id:str) -> str:
    '''determine the output file type for data request
    
    returns one of 'pl' 'sa' 'sf' or 'to'
    '''
   # TODO: [NB] make this more robust- its super janky
    flist = list(Path(directory).glob(f'*{request_id}*nc'))
    product = flist[0].name.split('_')[1]

    if 'gp' in product:
        file_type = 'to'
    elif 'p125' in product:
        file_type = 'pl'
    elif 'surf' in product:
        file_type = 'sa'
    elif 'phy2m125' in product:
        file_type = 'sf'
        
    else:
        raise KeyError(f"File type unknown for data request {request_id} ({product})")

    return file_type


class J3QgDownloadHandler(J3QDownloadHandler):
    REANALYSIS = 'jra3qg'


def new_jra_download_file(filename:str,
                          Times:xr.DataArray,
                          Lats:xr.DataArray,
                          Lons:xr.DataArray,
                          Levs:Optional[xr.DataArray]=None) -> nc.Dataset:
    output_file = nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC')

    logger.debug(f"Creating new empty output file: {Path(filename).name}")
    
    if Levs is not None:
        output_file.createDimension('level', len(Levs))
        levels     = output_file.createVariable('level', 'i4',('level',))
        levels.long_name  = 'pressure level'
        levels.units      = 'mbar'
        levels[:] = Levs.values

    dim_time = output_file.createDimension('time', len(Times))
    dim_lat = output_file.createDimension('latitude', len(Lats))
    dim_lon = output_file.createDimension('longitude', len(Lons))
    
    # make dimension variables
    var_times      = output_file.createVariable('time', 'd',('time',))
    var_latitude  = output_file.createVariable('latitude', 'f8', ('latitude',))
    var_longitude = output_file.createVariable('longitude', 'f8', ('longitude',))

    var_times[:] = Times.values
    var_latitude[:] = Lats.values
    var_longitude[:] = Lons.values

    var_times.standard_name = 'time'
    var_times.units     = Times.units
    var_times.calendar  = 'standard'
    var_latitude.standard_name  = Lats.long_name
    var_latitude.units      = Lats.units
    var_longitude.standard_name = Lons.long_name
    var_longitude.units     = Lons.units

    return output_file
   

def extract_downloaded_tar_files(directory:str, request_id:str, remove_when_completed:bool=True) -> list[Path]:
    '''find downloaded tar files for a given dataset id and extract them'''

    tarf = list(Path(directory).glob(f"*{request_id}*.tar"))
    
    for f in tarf:
        if Path(f).suffix == '.tar':
            tar = tarfile.open(f)
            tar.extractall(path=directory, filter='fully_trusted')
            tar.close()
            
            if remove_when_completed:
                f.unlink()
    
    return tarf


def get_downloaded_variable_names(directory:str, request_id:str) -> list[str]:
    '''get all variables associated with a given dataset id'''

    flist = list(Path(directory).glob(f'*{request_id}*'))
    varlist = []
    for f in flist:
        fname = f.name
        for part in fname.split("."):
            if len(part.split("-")) == 4:
                # var = part.split("-")[0]
                var = part
                varlist.append(var)
    unique_variables = []
    
    for v in varlist:
        if v not in unique_variables:
            unique_variables.append(v)

    return unique_variables

def copy_variable_in_chunks(src_var, dst_var, max_mem_gb = 4):
    
    size_in_gb = src_var.data.nbytes * 1e-9
    n_time = src_var.shape[0]
    n_chunks = int(np.ceil(size_in_gb / max_mem_gb))
    chunk_size = int(n_time / n_chunks)
    chunk_dims = list(src_var.shape)
    chunk_dims[0] = chunk_size
    
    # arr = np.empty(shape=chunk_dims, dtype=src_var.dtype.str)
    for i in range(n_chunks):
        start = i * chunk_size
        end = (i+1) * chunk_size
        print(f"Compiling data {start} to {end - 1}", end='\r')
        S = slice(start, end)
        if len(dst_var.shape) == 4:
            dst_var[S,:,:,:] = src_var[S,:,:,:]
        else:   
            dst_var[S,:,:] = src_var[S,:,:]

    if n_chunks * chunk_size < n_time:
        start = n_chunks * chunk_size
        print(f"Compiling data {start} to {n_time}")
        if len(dst_var.shape) == 4:
            dst_var[start:,:,:,:] = src_var[start:,:,:,:]
        else:
            dst_var[start:,:,:] = src_var[start:,:,:]

if __name__ == "__main__":
    import argparse
    import logging

    L = logging.getLogger()
    L.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(epilog="Example usage:\npython JraDownloadHandler.py /path/to/directory -I 194323 -t pl --gauss")

    parser.add_argument('dir')
    parser.add_argument('-I', '--request-id', dest='rid', type=str)
    parser.add_argument('-t', '--type', choices=["to", "pl", "sa", "sf"], default='pl', dest='type')
    parser.add_argument('-K', '--keep-files', action='store_true', dest='keep_files')
    parser.add_argument('--gauss', action='store_true')
    args, _ = parser.parse_known_args()

    if args.gauss:
        j = J3QgDownloadHandler()
    else:
        j = J3QDownloadHandler()

    if args.keep_files:
        j.keep_raw_files()

    print("Starting dataset conversion")
    j.make_globsim_dataset(directory=args.dir, request_id=args.rid, filetype=args.type )

    # python JraDownloadHandler.py /path/to/dir -I 524213 -t pl --gauss