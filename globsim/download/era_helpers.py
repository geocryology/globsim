import cdsapi

from datetime import datetime
from pathlib import Path
from collections.abc import MutableMapping


from globsim.download.ERA5download import ERA5generic


class Era5RequestParameters(MutableMapping):
    """ Request dictionary """
    VALID_KEYS = ['product_type','format','year',
                 'month','day','time','area','variable',
                 'pressure_level']

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
            raise KeyError("Not a valid ERA5 request parameter")

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)
    
    def __len__(self):
        return len(self.store)

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

    def download(self):
        ''' download '''
        server = cdsapi.Client()

        query = self.params.as_dict()

        target = self.output_file
       
        server.retrieve(self.dataset, query, target)

    @property
    def dataset(self):
        return self._dataset

    @dataset.setter
    def dataset(self, value):
        if value in self.DATASETS:
            self._dataset = value
        else:
            raise KeyError(f"Not a valid dataset. Must be in {self.DATASETS.keys()}")

    def exists(self):
        return self.output_file.is_file()
    
    @property
    def output_file(self) -> Path:
        time = f"{self.params.start}_to_{self.params.end}"
        era_type = self.PRODUCTTYPES[self.params["product_type"]]
        dataset = self.DATASETS[self.dataset]
        file = Path(self.directory, f"era5_{era_type}_{dataset}_{time}.nc")

        return file

def make_monthly_chunks(start: datetime, end: datetime) -> "list[dict]":
    chunks = []
    
    for year in range(start.year, end.year):
        
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