import re

from typing import Optional

from globsim.meteorology import pressure_from_elevation
from globsim.download.RDA import Rdams


def getPressureLevels(levels: list, min_elev: float, max_elev: float) -> "list[float]":
    # flip max and min because 1000 is the bottom and 0 is the top
    elevationMax = pressure_from_elevation(min_elev)
    elevationMin = pressure_from_elevation(max_elev)

    minNum = min(levels, key=lambda x:abs(x - elevationMin))
    maxNum = min(levels, key=lambda x:abs(x - elevationMax))

    if (minNum > elevationMin and levels.index(minNum) > 0):
        elevationMinRange = levels.index(minNum) - 1
    else:
        elevationMinRange = levels.index(minNum)

    if (maxNum < elevationMin and levels.index(maxNum) < 36):
        elevationMaxRange = levels.index(maxNum) - 1
    else:
        elevationMaxRange = levels.index(maxNum)

    elevation = []
    for e in range(elevationMinRange, elevationMaxRange + 1):
        elevation.append(levels[e])

    return elevation


class JRAformatter:

    VALID_LEVELS = [100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550,
                    600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
                    925, 950, 975, 1000]
    
    def __init__(self, date, area, elevation={'min': 0, 'max': 2000}, variables=[]):
        '''Returns an object for JRA55 data that has methods for querying the
        NCAR server for pressure level variables (prec, swin, lwin). '''

        self.date = date
        self.area = area
        self.elevation = elevation
        self.variables = variables

    def makeDate(self):
        '''convert data format to NCAR RDA request'''

        beg = self.date['beg'].strftime('%Y%m%d%H%M')
        end = self.date['end'].strftime('%Y%m%d%H%M')
        dateRange = beg + '/to/' + end

        return dateRange

    def getParam(self, dpar, variables, add=None):

        varlist = []
        for var in variables:
            varlist.append(dpar.get(var))
        
        varlist = [item for item in varlist if item is not None]
        varlist = [item for sublist in varlist for item in sublist]

        if add is not None:
            for var in add:
                varlist.append(var)

        return varlist

    @staticmethod
    def param_descriptions_to_name(params: str, param_dict: dict) -> str:
        if isinstance(params, list):
            param_list = params
        elif isinstance(params, str):
            param_list = params.split('/')
        
        param_codes = {v['param_description']:v['param'] for v in param_dict}
        
        shortpar = "/".join([i for i in param_list])
        
        return shortpar
    

class J55DictFormatter(JRAformatter):
    
    DATASET = 'ds628.0'

    pl_dict = {'air_temperature'   : ['Temperature'],
                'relative_humidity' : ['Relative humidity'],
                'wind_speed'        : ['u-component of wind',
                                       'v-component of wind']}
    
    sf_dict = {'precipitation_amount':
            ['Total precipitation'],
        'downwelling_shortwave_flux_in_air':
            ['Downward solar radiation flux'],
        'downwelling_longwave_flux_in_air':
            ['Downward longwave radiation flux'],
        'downwelling_shortwave_flux_in_air_assuming_clear_sky':
            ['Clear sky downward solar radiation flux'],
        'downwelling_longwave_flux_in_air_assuming_clear_sky':
            ['Clear sky downward longwave radiation flux']}

    sa_dict = {'air_temperature'  : ['Temperature'],
                'relative_humidity': ['Relative humidity'],
                'specific_humidity': ['Specific humidity'],
                'wind_speed'       : ['u-component of wind',
                                      'v-component of wind']}

    def get_dict_template(self):
        dictionary = {
                'dataset': self.DATASET,
                'date': self.makeDate(),
                'param': None,
                'level': None,
                'oformat': 'netCDF',
                'nlat': str(self.area['north']),
                'slat': str(self.area['south']),
                'wlon': str(self.area['west']),
                'elon': str(self.area['east']),
                'product': 'Analysis',
                'compression': 'NN',
                'gridproj': 'latLon',
                'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'}
        
        return dictionary

    def get_pl_dict(self, variables):
        levels = getPressureLevels(self.VALID_LEVELS, self.elevation['min'], self.elevation['max'])
        elevation = [str(ele) for ele in levels]
        temp = self.get_dict_template()
        param = self.getParam(self.pl_dict, variables, add=['Geopotential height', 'Pressure'])
        temp['param'] = '/'.join(param)
        temp['level'] = 'Isobaric surface:' + '/'.join(elevation)
        return temp

    def get_sf_dict(self, variables):
        temp = self.get_dict_template()
        param = self.getParam(self.sf_dict, variables, add=['Pressure'])
        temp['param'] = '/'.join(param)
        temp['level'] = 'Ground or water surface:0'
        temp['product'] = '3-hour Average (initial+0 to initial+3)'

        return temp
    
    def get_sa_dict(self, variables):
        temp = self.get_dict_template()
        param = self.getParam(self.sa_dict, variables)
        temp['param'] = '/'.join(param)
        temp['level'] = 'Specified height above ground:2/10'

        return temp

    def get_to_dict(self):
        temp = self.get_dict_template()
        temp['date'] = "194709010000/to/194709010000"
        temp['param'] = 'Geopotential'
        temp['level'] = 'Ground or water surface:0'
        
        return temp
    

class J3QDictFormatter(JRAformatter):
    DATASET = 'd640000'
    REANALYSIS = 'jra3q'

    _pl_dict = {'air_temperature'   : ['tmp-pres-an-{grid}'],
                'relative_humidity' : ['rh-pres-an-{grid}'],
                'specific_humidity': ['spfh-pres-an-{grid}'],
                'wind_speed'        : ['ugrd-pres-an-{grid}',
                                       'vgrd-pres-an-{grid}']}

    _sa_dict = {'air_temperature'  : ['tmp2m-hgt-an-{grid}'],
                'relative_humidity': ['rh2m-hgt-an-{grid}'],
                'specific_humidity': ['spfh2m-hgt-an-{grid}'],
                'wind_speed'       : ['ugrd10m-hgt-an-{grid}',
                                      'vgrd10m-hgt-an-{grid}']}

    _to_dict = {'geopotential': ['gp-sfc-cn-{grid}']}

    _sf_dict = {'precipitation_amount':
            ['tprate1have-sfc-fc-{grid}'],
        'downwelling_shortwave_flux_in_air':
            ['dswrf1have-sfc-fc-{grid}'],
        'downwelling_longwave_flux_in_air':
            ['dlwrf1have-sfc-fc-{grid}'],
        'downwelling_shortwave_flux_in_air_assuming_clear_sky':
            ['dswrfcs1have-sfc-fc-{grid}'],
        'downwelling_longwave_flux_in_air_assuming_clear_sky':
            ['dlwrfcs1have-sfc-fc-{grid}']}
    
    def __init__(self, grid='ll125', metadata=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_grid(grid)
        rda = Rdams(metadata)
        self.metadata = rda.get_metadata(self.DATASET)['data']['data']
        
    def lookup_param(self, param, returns='param_description'):
        return lookup_param(self.metadata, param, returns)
    
    @property
    def pl_dict(self):
        return {k:[i.format(grid=self.grid) for i in v] 
                for k,v in self._pl_dict.items()}
    
    @property
    def sa_dict(self):
        return {k:[i.format(grid=self.grid) for i in v] 
                for k,v in self._sa_dict.items()}

    @property
    def sf_dict(self):
        return {k:[i.format(grid=self.grid) for i in v] 
                for k,v in self._sf_dict.items()}
    
    @property
    def to_dict(self):
        return {k:[i.format(grid=self.grid) for i in v] 
                for k,v in self._to_dict.items()}

    def get_pl_dict(self, variables):
        levels = getPressureLevels(self.VALID_LEVELS, self.elevation['min'], self.elevation['max'])
        elevation = [str(ele) for ele in levels]
        temp = self.get_dict_template()
        param = self.getParam(self.pl_dict, variables, add=[f'hgt-pres-an-{self.grid}'])
        temp['param'] = '/'.join(param)
        temp['level'] = 'pressure (isobaric) level:' + '/'.join(elevation)

        return temp

    def get_sf_dict(self, variables):
        temp = self.get_dict_template()
        param = self.getParam(self.sf_dict, variables, add=[f'pres-sfc-fcst-{self.grid}'])
        temp['param'] = '/'.join(param)
        temp['level'] = 'Surface:0'
        

        return temp
    
    def get_sa_dict(self, variables):
        temp = self.get_dict_template()
        param = self.getParam(self.sa_dict, variables)
        temp['param'] = '/'.join(param)
        temp['level'] = 'Surface:0'

        return temp

    def get_to_dict(self, *args, **kwargs):
        temp = self.get_dict_template()
        param = self.getParam(self.to_dict, ['geopotential'])
        metadata = find_param(self.metadata, level_description_pattern=None, 
                              name_pattern=None, 
                              param_pattern=f'.*-sfc-cn-{self.grid}',
                              names_only=False, drop_levels=False)

        temp['date'] = f"{metadata[0]['start_date']}/to/{metadata[0]['end_date']}"
        temp['param'] = '/'.join(param)
        temp['level'] = 'Surface:0'
        
        return temp
    
    def set_grid(self, grid):
        if grid not in ['ll125', 'gauss']:
            raise ValueError(f'Invalid grid ({grid}). Must be "ll125" or "gauss"')
        self.griddef = {'ll125': '288:145:90N:0E:90S:358.75E:1.25:1.25',
                        'gauss': '960:480:89.713N:0E:89.713S:359.625E:0.375:0.375'}[grid]
        self.grid = grid

    def get_dict_template(self):
        dictionary = {
                'dataset': self.DATASET,
                'date': self.makeDate(),
                'param': None,
                'level': None,
                'oformat': 'netCDF',
                'nlat': str(self.area['north']),
                'slat': str(self.area['south']),
                'wlon': str(self.area['west']),
                'elon': str(self.area['east']),
                'product': 'Analysis',
                'compression': 'NN',
                'gridproj': 'latLon',
                'griddef': self.griddef}
        
        return dictionary
    
    @staticmethod
    def param_descriptions_to_name(params: str, param_dict: dict) -> str:
        if isinstance(params, list):
            param_list = params
        elif isinstance(params, str):
            param_list = params.split('/')
        
        shortpar = "/".join(param_list)
        
        return shortpar


def find_param(md:dict, name_pattern:Optional[str]='wind', 
               level_description_pattern:Optional[str]='isobaric',
               param_pattern:Optional[str]='ll125', 
               names_only:bool=False,
               drop_levels:bool=True) -> list:
    param = md.copy()
    if name_pattern is not None:
        param = filter(lambda x: re.search(name_pattern, x['param_description']), md)

    if level_description_pattern is not None:
        param = filter(lambda x: re.search(level_description_pattern, x['levels'][0]['level_description']), param)
    
    if param_pattern is not None:
        param = filter(lambda x: re.search(param_pattern, x['param']), param)
    
    if names_only:
        return set([p['param'] for p in param])
    else:
        if drop_levels:
            param = [{k:v for k,v in p.items() if k != 'levels'} for p in param]
        
        return list(param)
    

def lookup_param(param, metadata=None, dsid=None, returns='param_description'):
    if metadata is None:
        if dsid is not None:
            rda = Rdams(None)
            metadata = rda.get_metadata(dsid)['data']['data']
        else:
            raise ValueError('metadata or dsid must be provided')
    
    return list(filter(lambda x: x['param'] == param, metadata))[0][returns]


if __name__ =="__main__":
    
    rda = Rdams(None)
    md = rda.get_metadata('d640000')['data']['data']
    a = find_param(md, 
                level_description_pattern=None, 
                name_pattern=None, 
                param_pattern=f'.*-sfc-cn',
                names_only=False, drop_levels=False)
    print(a)
    # [print(p['param_description'], p['param']) for p in q]
