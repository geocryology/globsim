import logging
import re
import subprocess

from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from os import rename

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

    def retrieve(self):
        requests = self.list_requests()
        self.download_threadded(requests, 12)
        
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

    def download_threadded(self, cds_requests, workers=6):
        logger.info(f"Starting ERA5 multi-threaded download with {workers} workers")
        incomplete_requests = list()
        for request in cds_requests:
            if request.exists():
                logger.warning(f"Found {request.output_file}. Will not re-download.")
            elif request.is_downloaded():
                logger.warning(f"Found {request.renamed_files}. Will not re-download.")
            else:
                incomplete_requests.append(request)

        download_threadded(incomplete_requests, workers)
        logger.info(f"Finished ERA5 multi-threaded download with {workers} workers")
        
        rename_pl_dir(self.directory)
        rename_sl_dir(self.directory)

    def retrieve(self, workers=6):
        requests = self.list_requests()
        self.download_threadded(requests, workers)


def download_threadded(cds_requests, workers=6):
    def download_request(request):
        if not request.exists():
            request.download()
        else:
            print("skipping request")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(download_request, cds_requests)
        executor.shutdown()


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
    print(new_file)
    rename(f, new_file)


def rename_sl_dir(dir):
    orig = re.compile(r"era5_re_resl_(\d{8}_to_\d{8}).nc")
    files = [str(f) for f in Path(dir).iterdir()]
    matched_files = [f for f in files if orig.search(f)]
    for f in matched_files:
        split_sl(f)


def split_sl(f, overwrite=False):
    orig = re.compile(r"era5_re_resl_(\d{8}_to_\d{8}).nc")
    sf = orig.sub(r'era5_sf_\1.nc', f)
    sa = orig.sub(r'era5_sa_\1.nc', f)
    if not overwrite and Path(sf).exists():
        print(f"Skipping {sf}")
        return
    if not overwrite and Path(sa).exists():
        print(f"Skipping {sa}")
        return
    cmd1 = f"nccopy -V time,latitude,longitude,ssrd,strd,tp {f} {sf}"
    cmd2 = f"nccopy -V time,latitude,longitude,d2m,t2m,tco3,tcwv,u10,v10 {f} {sa}"
    logger.debug(cmd1)
    p1 = subprocess.Popen(cmd1.split(" "))
    p1.wait()
    logger.debug(cmd2)
    subprocess.Popen(cmd2.split(" "))


if __name__ == "__main__":
    import argparse
    logger.setLevel(logging.DEBUG)
    parser = argparse.ArgumentParser(description="Download ERA5 data")
    
    # add subparser for splitting
    subparsers = parser.add_subparsers(dest='command')

    # add parser for download
    split = subparsers.add_parser('split', help="Download ERA5 data")
    split.add_argument('-S', '--sl', type=str, default=None, help="Single-level file")
    split.add_argument('-P', '--pl', type=str, default=None, help="Pressure-level file")
    split.add_argument('-o', '--overwrite', action='store_true', help="Overwrite existing files")

    args = parser.parse_args()
    
    if args.command == 'split':
        if args.sl:
            split_sl(args.sl, overwrite=args.overwrite)
        if args.pl:
            rename_pl_file(args.pl, overwrite=args.overwrite)
