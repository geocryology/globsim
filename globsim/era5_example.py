# -*- coding: utf-8 -*-
# topofeature
# ==== DOWNLOAD ===============================================================
import cdsapi

from datetime import timedelta
from math import floor
from os import path
from ecmwfapi.api import ECMWFDataServer

from .era5 import ERA5download, ERA5scale, ERA5interpolate

pfile = "/path/to/examples.globsim_download"
ifile = "/path/to/examples.globsim_interpolate"
sfile = "/path/to/examples.globsim_scale"

# ======================================================

# download
ERAd = ERA5download(pfile, 'ensemble_members')


'''
# prepare time loop
date_i = {}
slices = floor(float((ERAd.date['end'] - ERAd.date['beg']).days)/
               ERAd.chunk_size)+1

# topography
#if path.isfile(path.join(ERAd.directory,'era5_to.nc')):
#    print("WARNING: File 'era5_to.nc' already exists. Skipping.")
#else:
#    top = ERA5to(ERAd.area, ERAd.directory, ERAd.product)
#    top.download(ERAd.api, ERAd.storage)
ind = 0
date_i['beg'] = ERAd.date['beg'] + timedelta(days =
                ERAd.chunk_size * ind)
date_i['end'] = ERAd.date['beg'] + timedelta(days =
                ERAd.chunk_size * (ind+1) - 1)
if ind == (slices-1):
    date_i['end'] = ERAd.date['end']
from era5 import ERA5pl, ERA5sa, ERA5sf, ERA5to
#actual functions
to = ERA5to(ERAd.product, ERAd.area, ERAd.directory)
sa = ERA5sa(ERAd.product, date_i, ERAd.area, ERAd.variables, ERAd.directory)
sf = ERA5sf(ERAd.product, date_i, ERAd.area, ERAd.variables, ERAd.directory)
pl = ERA5pl(ERAd.product, date_i, ERAd.area, ERAd.elevation, ERAd.variables, ERAd.directory)
to.dictionary = {
    'product_type': to.productType,
    'variable'    : ['orography', 'land_sea_mask'],
       }
to.dictionary.update(sf.getDictionaryGen(sf.area, sf.date))
era = to
to.download(api, storage)
# ======================================================
# era.download(ERAd.api, ERAd.storage)
api = ERAd.api
storage = ERAd.api
server = cdsapi.Client()

query   = era.getDictionary()
target  = era.file_ncdf
levtype = era.levelType
#query = era.ECM2CDS(query)
server.retrieve(levtype, query, target)
# report inventory
self.inventory()
'''

ERAd.retrieve()
ERAd.inventory()

# ERAsc = ERAscale(sfile)
# ERAsc.process()
