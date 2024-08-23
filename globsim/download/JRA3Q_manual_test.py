
from rdams import *
import traceback

dsid = 'd640000' 

get_summary(dsid)

# parms  = get_all_params(dsid)
md = get_metadata(dsid)['data']['data']
param = list(filter(lambda x: x['param']=='tmp-pres-an-ll125', md['data']['data']))[0]

print(param['param_description'])
print(param['levels'][5])

control = { 
         'dataset' : dsid,
         'date':'201609200000/to/201609210000',
         'datetype': 'init',
         'param': 'tmp-pres-an-ll125',  
         'level': 'pressure (isobaric) level:200', 
         'oformat':'netCDF4',
         'nlat':53,
         'slat': 50,
         'elon':-100,
         'wlon':-120,
         'product':'Analysis'
         } 
         
control = {'dataset': 'd640000',
           'date': '201909010000/to/202008312300',
           'param': 'tmp-pres-an-ll125/rh-pres-an-ll125/ugrd-pres-an-ll125/vgrd-pres-an-ll125/hgt-pres-an-ll125',
           'level': 'pressure (isobaric) level:750/775/800/825/850/875/900/925/950/975/1000',
           'oformat': 'netCDF',
           'nlat': '66', 'slat': '62', 'wlon': '-112', 'elon': '-108',
           'product': 'Analysis', 'compression': 'NN',
           'gridproj': 'latLon',
           'griddef': '288:145:90N:0E:90S:358.75E:1.25:1.25'}

response = submit_json(control)
get_status(response['data']['request_id'])
print(response)
