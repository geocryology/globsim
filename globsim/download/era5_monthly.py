import logging
from os import scandir
import re
import subprocess
import time
import zipfile

from datetime import datetime
from functools import partial
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from os import rename
from requests.exceptions import ChunkedEncodingError

from globsim.download.GenericDownload import GenericDownload
from globsim.download.era_helpers import make_monthly_chunks, Era5Request, Era5RequestParameters, era5_pressure_levels, cf_to_cds_pressure, cf_to_cds_single


logger = logging.getLogger("globsim.download")
    

class ERA5MonthlyDownload(GenericDownload):

    def __init__(self, pfile, ensemble=False):
        super().__init__(pfile)
        par = self.par

        # time bounds
        self.date = {'beg': datetime.strptime(par.get('beg'), '%Y/%m/%d'),  # type: ignore
                     'end': datetime.strptime(par.get('end'), '%Y/%m/%d')}  # type: ignore

        logger.warning("Using optimal monthly chunks.  Chunk size in parameter file ignored.")
        self.chunks = make_monthly_chunks(self.date['beg'], self.date['end'])

        # pressure level bounds
        self.levels = era5_pressure_levels(self.elevation['min'], self.elevation['max'])

        if ensemble:
            raise NotImplementedError("Not implemented yet for ensemble download")
            # self._set_input_directory("era5ens")
            # self.topo_file = 'era5_ens_to.nc'
            # self.product_type = 'something else'
        else:
            self._set_input_directory("era5")
            self.product_type = 'reanalysis'
            self.topo_file = 'era5_to.nc'
        
    def list_requests(self):
        variables = self.par.get('variables')

        area = [self.area["north"],
                self.area["west"],
                self.area["south"],
                self.area["east"]]

        requests_pl = []
        requests_sl = []

        to_param = Era5RequestParameters(product_type=self.product_type,
                                         variable=["land_sea_mask", "geopotential"],
                                         year="2000",
                                         month="01",
                                         day="01",
                                         time="00:00",
                                         area=area,
                                         format='netcdf',
                                         data_format= 'netcdf',
                                         download_format = 'unarchived')

        request_to = Era5Request('reanalysis-era5-single-levels', self.directory, to_param)
        request_to.set_output_file(self.topo_file)

        for period in self.chunks:
            year = period["year"]
            month = period["month"]
            day = period["day"]

            pl = Era5RequestParameters(product_type=self.product_type,
                                       variable=cf_to_cds_pressure(variables) + ['geopotential'],
                                       year=year,
                                       month=month,
                                       day=day,
                                       time=Era5RequestParameters.all_times(),
                                       area=area,
                                       pressure_level=self.levels,
                                       format='netcdf',
                                       data_format= 'netcdf',
                                       download_format = 'unarchived')
            rpl = Era5Request('reanalysis-era5-pressure-levels', self.directory, pl)
            requests_pl.append(rpl)

            sl = Era5RequestParameters(product_type=self.product_type,
                                       variable=cf_to_cds_single(variables),
                                       year=year,
                                       month=month,
                                       day=day,
                                       time=Era5RequestParameters.all_times(),
                                       area=area,
                                       format='netcdf',
                                       data_format= 'netcdf',
                                       download_format = 'unarchived')

            rsl = Era5Request('reanalysis-era5-single-levels', self.directory, sl)
            requests_sl.append(rsl)

        return [request_to] + requests_pl + requests_sl

    def rename_files(self):
        rename_pl_dir(self.directory)
        rename_sl_dir(self.directory)

    def get_dotrc(self):
        credentials_dir = self.par.get('credentials_directory', None)
        if credentials_dir is not None:
            dotrc = Path(credentials_dir, ".cdsapirc") 
            if dotrc.is_file():
                return dotrc

    def download_threadded(self, cds_requests, workers=6):
        logger.info(f"Starting ERA5 multi-threaded download with {workers} workers")
        incomplete_requests = list()
        for request in cds_requests:
            if request.is_downloaded():
                logger.warning(f"Found {[file.name for file in request.renamed_files]}. Will not re-download.")
            elif request.exists():
                logger.warning(f"Found {request.output_file.name}. Will not re-download.")
            else:
                incomplete_requests.append(request)

        dotrc = self.get_dotrc()
        failed_downloads = download_threadded(incomplete_requests, workers, dotrc=dotrc)
        logger.info(f"Finished ERA5 multi-threaded download with {workers} workers")
        
        logger.info(f"Renaming files in {self.directory}")
        rename_pl_dir(self.directory)
        rename_sl_dir(self.directory)

        failfile_path = Path(self.directory, "failed_downloads.txt")
        
        if len(failed_downloads) > 0:
            logger.warning(f"Failed to download {len(failed_downloads)} files. See {failfile_path} for more information")
        
            with open(failfile_path, 'w') as f:
                for failed in failed_downloads:
                    f.write(f"{failed.output_file}\n")
        
        else:
            logger.info("All files downloaded successfully")
            if failfile_path.is_file():
                failfile_path.unlink()

    def retrieve(self, workers=6):
        requests = self.list_requests()
        failed_downloads = self.download_threadded(requests, workers)
        return failed_downloads


def _download_request(request: Era5Request, dotrc=None):
    max_retries = 5
    current_retries = 0
    
    while (not request.exists()) and (current_retries < max_retries):
        try:
            request.download(dotrc=dotrc)
        except ChunkedEncodingError as e:
            logger.warning(f"ChunkedEncodingError on download attempt #{current_retries}")
            request.purge()
        except Exception as e:
            logger.warning(f"Unknown error on download attempt #{current_retries}: {e}")
            request.purge()
        finally:
            current_retries += 1
    
    if not request.exists():
        logger.error(f"Failed to download {request.output_file}")
        return request
    else:
        return None


def download_threadded(cds_requests, workers=6, dotrc=None):
    failed_downloads = []
    
    with ThreadPoolExecutor(max_workers=workers) as executor:
        for result in executor.map(partial(_download_request,dotrc=dotrc), cds_requests):
            if result is not None:
                failed_downloads.append(result)
        executor.shutdown()
    
    return failed_downloads

def download_serial(cds_requests, dotrc=None):
    failed_downloads = []
    
    for r in cds_requests:
        result = _download_request(r, dotrc=dotrc)
        if result is not None:
                failed_downloads.append(result)
    
    return failed_downloads


def rename_pl_dir(dir):
    pl_pattern = re.compile(r"era5_re_repl_(\d{8}_to_\d{8}).nc")
    files = [str(f) for f in Path(dir).iterdir() if pl_pattern.search(str(f))]
    for f in files:
        rename_pl_file(f)


def rename_pl_file(f, overwrite=True):
    pl_pattern = re.compile(r"era5_re_repl_(\d{8}_to_\d{8}).nc")
    new_file = pl_pattern.sub(r"era5_pl_\1.nc", f)
    if not overwrite and Path(new_file).exists():
        print(f"Skipping {new_file}")
        return
    logger.debug(f"Renaming {Path(f).name}")
    print(new_file)
    rename(f, new_file)


def rename_sl_dir(dir):
    orig = re.compile(r"era5_re_resl_(\d{8}_to_\d{8}).nc")
    files = [str(f) for f in Path(dir).iterdir()]
    matched_files = [f for f in files if orig.search(f)]
    for f in matched_files:
        split_sl(f)


def split_sl(f, overwrite=True, time_var='valid_time'):
    orig = re.compile(r"era5_re_resl_(\d{8}_to_\d{8}).nc")
    sf = orig.sub(r'era5_sf_\1.nc', f)
    sa = orig.sub(r'era5_sa_\1.nc', f)
    if not overwrite and Path(sf).exists():
        print(f"Skipping {sf}")
        return
    if not overwrite and Path(sa).exists():
        print(f"Skipping {sa}")
        return
    logger.debug(f"Splitting {Path(f).name}")
    if zipped:
        split_resl_zip(zipf, sa, sf)
    else:
        split_resl_nc(f, sa, sf, time_var)


def split_resl_zip(zipf:zipfile.ZipFile, sa:str, sf:str):
    zname = zipf.filename
    parent_dir = Path(zname).parent
    zipf.extractall(parent_dir)
    f_acc = Path(zname).with_name("data_stream-oper_stepType-accum.nc")
    f_ins = Path(zname).with_name("data_stream-oper_stepType-instant.nc")
    cmd1 = f'ncks -C -O -x -v number,expver {f_acc} {sf}'
    cmd2 = f'ncks -C -O -x -v number,expver {f_ins} {sa}'
    logger.debug(cmd1)
    p1 = subprocess.Popen(cmd1.split(" "))
    p1.wait()
    if not wait_for_file_scandir(parent_dir, f"{sf}$", 360, 1):
        logger.error(f"Timed out waiting for {sf}")
    logger.debug(cmd2)
    subprocess.Popen(cmd2.split(" "))
    if not wait_for_file_scandir(parent_dir, f"{sa}$", 360, 1):
        logger.error(f"Timed out waiting for {sa}")
    if Path(sf).exists() and Path(sa).exists():
        logger.debug(f"Removing {zipf.filename}")
        # Path(f).unlink()
    for f in [f_acc, f_ins]:
        f.unlink()

def split_resl_nc(f, sa, sf, time_var):
    cmd1 = f"nccopy -V {time_var},latitude,longitude,ssrd,strd,tp {f} {sf}"
    cmd2 = f"nccopy -V {time_var},latitude,longitude,d2m,t2m,tco3,tcwv,u10,v10 {f} {sa}"
    logger.debug(cmd1)
    p1 = subprocess.Popen(cmd1.split(" "))
    p1.wait()
    logger.debug(cmd2)
    subprocess.Popen(cmd2.split(" "))
    if Path(sf).exists() and Path(sa).exists():
        logger.debug(f"Removing {f}")
        # Path(f).unlink()


def wait_for_file_scandir(directory, pattern, timeout=60, check_interval=1):
    regex = re.compile(pattern)
    elapsed_time = 0

    while elapsed_time < timeout:
        with scandir(directory) as entries:
            if any(regex.match(entry.name) for entry in entries if entry.is_file()):
                return True  #  matching files exist
        time.sleep(check_interval)
        elapsed_time += check_interval

    return False  # Timed out waiting for file 

if __name__ == "__main__":
    import argparse
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    logger.addHandler(console_handler)

    parser = argparse.ArgumentParser(description="Download ERA5 data")
    
    # add subparser for splitting
    subparsers = parser.add_subparsers(dest='command')

    # add parser for download
    split = subparsers.add_parser('split', help="Download ERA5 data")
    split.add_argument('-s', '--sl', type=str, default=None, help="Single-level file")
    split.add_argument('-S', '--sl-dir', type=str, dest="sldir", default=None, help="Single-level file")
    split.add_argument('-p', '--pl', type=str, default=None, help="Pressure-level file")
    split.add_argument('-o', '--overwrite', action='store_true', help="Overwrite existing files")

    args = parser.parse_args()
    
    if args.command == 'split':
        if args.sl:
            split_sl(args.sl, overwrite=args.overwrite)
        if args.pl:
            rename_pl_file(args.pl, overwrite=args.overwrite)
        if args.sldir:
            rename_sl_dir(args.sldir)
