# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 01:48:46 2018

@author: BinCao
"""

# ==== INTERPOLATION ==========================================================

from era5 import ERAinterpolate
from os import path

ifile = 'C:/OneDrive/Bitbucket/era5/examples.globsim_interpolate'

ERAinterp = ERAinterpolate(ifile)
#ERAinterp.process()

     
ERAinterp.ERA2station(path.join(ERAinterp.dir_inp,'era_to.nc'), 
                      path.join(ERAinterp.dir_out,'era_to_' + 
                                ERAinterp.list_name + '.nc'), 
                                ERAinterp.stations,
                           ['z', 'lsm'], date = ERAinterp.date)  
                          
# === 2D Interpolation for Surface Analysis Data ===    
# dictionary to translate CF Standard Names into ERA5
# pressure level variable keys. 
dpar = {'air_temperature'   : ['t2m'],  # [K] 2m values
        'relative_humidity' : ['d2m'],  # [K] 2m values
        'downwelling_shortwave_flux_in_air_assuming_clear_sky' : 
            ['tco3', 'tcwv'],   # [kg m-2] Total column ozone 
                                # [kg m-2] Total column W vapor                                                             
        'wind_speed' : ['u10', 'v10']}   # [m s-1] 10m values   
varlist = ERAinterp.TranslateCF2short(dpar)
ERAinterp.ERA2station(path.join(ERAinterp.dir_inp,'era_sa_*.nc'), 
                 path.join(ERAinterp.dir_out,'era_sa_' + 
                           ERAinterp.list_name + '.nc'), ERAinterp.stations,
                           varlist, date = ERAinterp.date)          

# 2D Interpolation for Surface Forecast Data    'tp', 'strd', 'ssrd' 
# dictionary to translate CF Standard Names into ERA5
# pressure level variable keys.       
# [m] total precipitation
# [J m-2] short-wave downward
# [J m-2] long-wave downward
dpar = {'precipitation_amount'              : ['tp'],   
        'downwelling_shortwave_flux_in_air' : ['ssrd'], 
        'downwelling_longwave_flux_in_air'  : ['strd']} 
varlist = ERAinterp.TranslateCF2short(dpar)                           
ERAinterp.ERA2station(path.join(ERAinterp.dir_inp,'era_sf_*.nc'), 
                 path.join(ERAinterp.dir_out,'era_sf_' + 
                           ERAinterp.list_name + '.nc'), 
                           ERAinterp.stations,
                           varlist, date = ERAinterp.date)          
                 

# === 2D Interpolation for Pressure Level Data ===
# dictionary to translate CF Standard Names into ERA5
# pressure level variable keys. 
dpar = {'air_temperature'   : ['t'],           # [K]
        'relative_humidity' : ['r'],           # [%]
        'wind_speed'        : ['u', 'v']}      # [m s-1]
varlist = self.TranslateCF2short(dpar).append('z')
self.ERA2station(path.join(self.dir_inp,'era_pl_*.nc'), 
                 path.join(self.dir_out,'era_pl_' + 
                           self.list_name + '.nc'), self.stations,
                           varlist, date = self.date)  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
# 1D Interpolation for Pressure Level Data
self.levels2elevation(path.join(self.dir_out,'era_pl_' + 
                                self.list_name + '.nc'), 
                      path.join(self.dir_out,'era_pl_' + 
                                self.list_name + '_surface.nc'))
                 
from os import listdir


import numpy as np
import netCDF4 as nc

import glob

ncfile_in  = path.join(ERAinterp.dir_inp,'era_sa_*.nc')

#ncsingle = filter((listdir(ERAinterp.dir_inp), path.basename(ncfile_in)))[0]
#ncsingle = path.join(ERAinterp.dir_inp, ncsingle)
ncsingle = glob.glob(path.join(ERAinterp.dir_inp, 
                                     path.basename(ncfile_in)))[0]


sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)