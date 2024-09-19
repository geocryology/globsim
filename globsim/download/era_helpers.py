import cdsapi
import re

import numpy as np

from datetime import datetime
from pathlib import Path
from collections.abc import MutableMapping
from globsim.meteorology import pressure_from_elevation
from pandas import DataFrame
from typing import Optional

from globsim.download.ERA5download import ERA5generic


class Era5RequestParameters(MutableMapping):
    """ Request dictionary """
    VALID_KEYS = ['product_type','format','year',
                  'month','day','time','area','variable','data_format','download_format',
                  'pressure_level']

    REQUIRED_KEYS = ['product_type','data_format','year',
                     'month','day','time','variable']

    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __repr__(self):
        return str(self.store)

    def __str__(self):
        return str(self.store)

    def __getitem__(self, key):
        return self.store[key]

    def __setitem__(self, key, value):
        if key in self.VALID_KEYS:
            self.store[key] = value
        else:
            raise KeyError(f"{key}:{value} is not a valid ERA5 request parameter")

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def validate(self):
        if self.missing():
            print("missing parameters")
            return False

        return True

    def missing(self):
        missing = []
        for key in self.REQUIRED_KEYS:
            if key not in self.store.keys():
                missing.append(key)

        return(missing)

    @property
    def start(self):
        return self.__date(min)

    @property
    def end(self):
        return self.__date(max)

    def as_dict(self):
        return self.store

    def __date(self, f):
        if isinstance(self['year'], list):
            y = f([int(year) for year in self['year']])
        else:
            y = int(self['year'])

        if isinstance(self['month'], list):
            m = f([int(month) for month in self['month']])
        else:
            m = int(self['month'])

        if isinstance(self['day'], list):
            d = f([int(day) for day in self['day']])
        else:
            d = int(self['day'])

        return "{:04d}{:02d}{:02d}".format(y,m,d)

    @staticmethod
    def all_times():
        return ["{:02d}:00".format(H) for H in range(0, 24)]

    @staticmethod
    def all_days():
        return ["{:02d}".format(d) for d in range(1, 32)]

    @staticmethod
    def all_months():
        return ["{:02d}".format(m) for m in range(1, 13)]


class Era5Request(ERA5generic):
    DATASETS = {"reanalysis-era5-pressure-levels": "repl",
                'reanalysis-era5-single-levels': "resl"}

    PRODUCTTYPES = {
        'ensemble_members': 'ens',
        'reanalysis': 're'}

    def __init__(self, dataset: str, directory: str, request_params: Era5RequestParameters):
        self.params = request_params
        self.directory  = directory
        self.dataset = dataset
        self._output_file = None

    def __repr__(self):
        return f"ERA5 Request: {self.dataset} {self.params}"

    def __str__(self):
        return f"{self.dataset} request for {len(self.params)} variables at {self.params['area']}"
    
    def download(self, target=None):
        ''' download '''
        server = cdsapi.Client()

        query = self.params.as_dict()

        if not target:
            target = self.output_file

        server.retrieve(self.dataset, query, target)

    @property
    def dataset(self):
        return self._dataset

    def set_output_file(self, filename):
        if Path(self.directory) in Path(filename).parents:  # full path specified
            self._output_file = Path(filename)
        else:                                               # just the filename specified
            self._output_file = Path(self.directory, filename)

    @dataset.setter
    def dataset(self, value):
        if value in self.DATASETS:
            self._dataset = value
        else:
            raise KeyError(f"Not a valid dataset. Must be in {self.DATASETS.keys()}")

    def exists(self) -> bool:
        """ check whether the file originally downloaded from the server is present"""
        return self.output_file.is_file()

    def is_downloaded(self) -> bool:
        """ checks for presence of either the original or renamed file """
        exists = self.exists()
        renamed = all([p.is_file() for p in self.renamed_files])

        return (exists or renamed)

    @property
    def product_type_alias(self) -> "Optional[str]":
        try:
            p_type = self.PRODUCTTYPES[self.params["product_type"]]
        except KeyError:
            p_type = None

        return p_type

    @property
    def datset_alias(self) -> "Optional[str]":
        try:
            dataset = self.DATASETS[self.dataset]
        except KeyError:
            dataset = None

        return dataset

    @property
    def renamed_files(self) -> "list[Path]":
        """ How is output_file renamed """
        time = f"{self.params.start}_to_{self.params.end}"
        files = []

        if self.product_type_alias == "re":
            if self.dataset == "reanalysis-era5-single-levels":
                files = [Path(self.directory, f"era5_{x}_{time}.nc")
                         for x in ["sa", "sf"]]

            elif self.dataset == "reanalysis-era5-pressure-levels":
                files = [Path(self.directory, f"era5_pl_{time}.nc")]

        elif self.product_type_alias == "ens":
            pass

        return files

    @property
    def output_file(self) -> Path:
        if self._output_file:
            file = self._output_file
        else:
            time = f"{self.params.start}_to_{self.params.end}"
            era_type = self.PRODUCTTYPES[self.params["product_type"]]
            dataset = self.DATASETS[self.dataset]
            file = Path(self.directory, f"era5_{era_type}_{dataset}_{time}.nc")

        return file
    
    @property
    def globsim_outputs(self) -> "list[Path]":
        if self.dataset == "reanalysis-era5-single-levels": 
            orig = re.compile(r"era5_re_resl_(\d{8}_to_\d{8}).nc")
            sf = orig.sub(r'era5_sf_\1.nc', self.output_file.name)
            sa = orig.sub(r'era5_sa_\1.nc', self.output_file.name)
            
            files = [self.output_file.with_name(n) for n in [sa, sf]]

        elif self.dataset == "reanalysis-era5-pressure-levels":
            pl_pattern = re.compile(r"era5_re_repl_(\d{8}_to_\d{8}).nc")
            pl = pl_pattern.sub(r"era5_pl_\1.nc", self.output_file.name)
            
            files = [self.output_file.with_name(pl)]

        return files


def era5_pressure_levels(elev_min: float, elev_max: float) -> "list[str]":
    """Restrict list of ERA5 pressure levels to be downloaded"""
    Pmax = pressure_from_elevation(elev_min) + 55  # Pressure inversely correlated
    Pmin = pressure_from_elevation(elev_max) - 55
    levs = np.array([300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775,
                     800, 825, 850, 875, 900, 925, 950, 975, 1000])
    mask = (levs >= Pmin) * (levs <= Pmax)  # select
    levs = [str(levi) for levi in levs[mask]]
    return levs


def era5_area_string(north, west, south, east):
    """Converts numerical coordinates into string: North/West/South/East"""
    res = str(round(north,2)) + "/"
    res += str(round(west, 2)) + "/"
    res += str(round(south,2)) + "/"
    res += str(round(east, 2))
    return(res)


def era5_area_list(north: float, west:float, south:float, east:float):
    return [north, west, south, east]


def make_monthly_chunks(start: datetime, end: datetime) -> "list[dict]":
    chunks = []

    for year in range(start.year, end.year + 1):

        for month in range(1, 13):

            chunk = {'year': str(year),
                     'month': "{:02d}".format(month),
                     'time': ["{:02d}:00".format(H) for H in range(0, 24)]}
            if year == start.year and month < start.month:
                continue

            elif year == start.year and month == start.month:
                day = ["{:02d}".format(d) for d in range(start.day, 32)]

            elif year == end.year and month == end.month:
                day = ["{:02d}".format(d) for d in range(1, end.day + 1)]

            elif year == end.year and month > end.month:
                continue

            else:
                day = ["{:02d}".format(d) for d in range(1, 32)]

            chunk['day'] = day

            if len(chunk['day'])  == 1:
                chunk['day'] = chunk['day'][0]

            chunks.append(chunk)

    return chunks

# https://www.ecmwf.int/en/forecasts/datasets/list-parameters-incorrectly-calculated-era-20cmv0-edmo


def cf_to_cds_single(vars: "list[str]"):
    sl = DataFrame.from_dict(
        {0: {"mars": '228.128', "cf":'precipitation_amount', "cds":"total_precipitation"},
         1: {"mars": '169.128', "cf":"downwelling_shortwave_flux_in_air", "cds": "surface_solar_radiation_downwards"},
         2: {"mars": "175.128", "cf":"downwelling_longwave_flux_in_air", "cds": "surface_thermal_radiation_downwards"},
         3: {"mars": "167.128", "cf":"air_temperature", "cds": "2m_temperature"},
         4: {"mars": "168.128", "cf": "relative_humidity", "cds": "2m_dewpoint_temperature"},
         5: {"mars": "165.128", "cf": "wind_speed", "cds": "10m_u_component_of_wind"},
         6: {"mars": "166.128", "cf": "wind_speed", "cds": "10m_v_component_of_wind"},
         7: {"mars": "206.128", "cf": "downwelling_shortwave_flux_in_air_assuming_clear_sky", "cds": "total_column_ozone"},
         8: {"mars": "137.128", "cf": "downwelling_shortwave_flux_in_air_assuming_clear_sky", "cds": "total_column_water_vapour"}},
        orient='index')
    index = sl['cf'].isin(vars)
    return sl['cds'][index].tolist()


def cf_to_cds_pressure(vars: "list[str]"):

    pl = DataFrame.from_dict(
        {1: {"mars": "130.128", "cf": "air_temperature", "cds": "temperature"},
         2: {"mars": "157.128", "cf": "relative_humidity", "cds": "relative_humidity"},
         3: {"mars": "131.128", "cf": "wind_speed", "cds": "u_component_of_wind"},
         4: {"mars": "132.128", "cf": "wind_speed", "cds": "v_component_of_wind"}},
        orient='index')

    index = pl['cf'].isin(vars)
    return pl['cds'][index].tolist()
