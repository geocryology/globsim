
d = {'dataset': 'ds628.0', 
     'date': '198201010000/to/198202010000', 
     'param': 'GP', 
     'level': 'Ground or water surface:0', 
     'oformat': 'netCDF', 
     'nlat': '45', 
     'slat': '40', 
     'wlon': '39',
     'elon': '44', 
     'product': 'Analysis', 
     'compression': 'NN', 
     'gridproj': 'latLon',
     'griddef': '288:145:90N:0E:90S:1.25W:1.25:1.25'}
r = submit_json(d)

dsid ='ds084.1'
dsid = 'ds628.0'  # JRA 55
 
dsid = 'ds640.0'  # jra 3Q  
# Get dataset info 
get_summary('ds084.1')
get_summary('ds640.0')

get_summary(dsid)
parms  = get_all_params(dsid)
md = get_metadata(dsid)['data']['data']


for i, param in enumerate(md):
    if param['param'] == 'R H':
        print(i, param['levels'][0])

for i, param in enumerate(md['data']['data']):
    if 'tmp' in param['param']:
        print(f"{i:03d}: [{len(param['levels']):03d} levels] {param['param']:<25s} {param['param_description']} , " )
        

# Let's get a template. 
response = get_control_file_template(dsid)
template = response['data']['template'] # Template string

# Parse the string
template_dict = read_control_file(template)
from pprint import pprint
pprint(template_dict)
# Insert our TMP param
template_dict['param'] = 'TMP'
template_dict['level'] = 'ISBL:' + '/'.join(ISBL_levels)
template_dict

s = get_status()
if __name__ == "__main__":
    """Calls main method"""
    query(sys.argv[1:])

[print(v['param'], v['param_description']) for v in par['data']['data'] if 'eopot' in v['param_description']]

'''
A working example

param = md['data']['data'][84]
param = md['data']['data'][97]

param = 
temps = list(filter(lambda x: ('an-ll125' in x['param']) and ('tmp' in x['param']) and (len(x['levels']) > 1), md['data']['data']))
[f"{p['param_description']}:{p['param']}" for p in temps]

ps = list(filter(lambda x: ('tmp-pres-an-ll125' in x['param']), md['data']['data']))

for key in ps[0].keys():
    print(key, ps[0][key], ps[1][key])
    
control = { 
         'dataset' : 'ds084.1',
         'date':'201609200000/to/201609200000',
         'datetype':'init',
         'param':'V GRD',
         'level':'HTGL:100',
         'oformat':'csv',
         'nlat':-10,
         'slat':-10,
         'elon':45,
         'wlon':45,
         'product':'Analysis'
         } 
         
response = submit_json(control)
assert response['http_response'] == 200
rqst_id = response['data']['request_id']

def check_ready(rqst_id, wait_interval=120):
    """Checks if a request is ready."""
    for i in range(100): # 100 is arbitrary. This would wait 200 minutes for request to complete
        res = get_status(rqst_id)
        request_status = res['data']['status']
        if request_status == 'Completed':
            return True
        print(request_status)
        print('Not yet available. Waiting ' + str(wait_interval) + ' seconds.' )
        time.sleep(wait_interval)
    return False
    
check_ready(rqst_id)
download(rqst_id)
purge_request(rqst_id)

print(response)
         
'''

'''
A JRA3Q example

A broken example

control = { 
         'dataset' : 'd640000',
         'date':'20160920/to/20160921',
         'datetype':'init',
         'param': 'tmp-pres-an-ll125',  # param['param'],
         'level': "pressure (isobaric) level:650",  # param['levels'][0]['level_description'],
         'oformat':'netCDF4',
         'nlat':"60",
         'slat':"50",
         'elon':"-100",
         'wlon':"-120",
         'product':'Analysis'
         } 
         
         
response = submit_json(control)
print(response)
assert response['http_response'] == 200
rqst_id = response['data']['request_id']

def check_ready(rqst_id, wait_interval=120):
    """Checks if a request is ready."""
    for i in range(100): # 100 is arbitrary. This would wait 200 minutes for request to complete
        res = get_status(rqst_id)
        request_status = res['data']['status']
        if request_status == 'Completed':
            return True
        elif request_status == 'Error':
            raise ValueError('Request failed')
        print(request_status)
        print('Not yet available. Waiting ' + str(wait_interval) + ' seconds.' )
        time.sleep(wait_interval)
    return False
    
check_ready(rqst_id)
download(rqst_id)
purge_request(rqst_id)

print(response)



A JRA55 example

A working example

param = list(filter(lambda x: ('emperat' in x['param_description']) and (len(x['levels'])  > 10),  md['data']['data']))[1]
control = { 
         'dataset' : 'ds628.0',
         'date':'19810101/to/19810102',
         # 'datetype':'init',
         'param': param['param'],
         'level': "Isobaric surface:850/875/900",  # param['levels'][0]['level_description'],
         'oformat':'netCDF4',
         'nlat':"55",
         'slat':"50",
         'elon':"-120",
         'wlon':"-115",
         'product':'Analysis'
         } 
response = submit_json(control)
print(response)
rqst_id = response['data']['request_id']

check_ready(rqst_id)



control = { 
         'dataset' : 'ds628.0',
         'date':'19810101/to/19810101',
         'datetype':'init',
         'param': 'tmp-pres-an-ll125',  # param['param'],
         'level': "Ground or water surface:0",  # param['levels'][0]['level_description'],
         'oformat':'netCDF4',
         'nlat':"55",
         'slat':"50",
         'elon':"-120",
         'wlon':"-115",
         'product':'Analysis'
         } 
response = submit_json(control)
print(response)
rqst_id = response['data']['request_id']

check_ready(rqst_id)
'''