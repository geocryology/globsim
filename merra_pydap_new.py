from pydap.client import open_url
from pydap.cas.urs import setup_session

username = "quanxj17"
password = "Qxj17carleton"

url = ('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/'
       'M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4')

session = setup_session(username, password, check_url=url)
ds = open_url(url, session=session)

#get variable keys
print ds.keys

#get latitudes
lat = ds.lat[:]

#find shape of one variable: T2M
ds.T2M.shape

#get subset of one variable
ds.T2M[0:2,0:2,1]