#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic classes, methods, functions used for more than one reanalysis.
#
#
# (C) Copyright 2017-2019 Stephan Gruber
#				2017-2018 Xiaojing Quan
#				2018-2019 Nicholas Brown
#				2018-2019 Bin Cao
#         
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# === CONTRIBUTIONS ============================================================
# Intial code for globsim is designed by Stephan Gruber. Xiaojing Quan developed
# MERRA-2, download part of ERA-Interim, and intial interpolation and scale  
# parts of globsim. Bin Cao wrote the download parts of ERA5, JRA-55, and
# finilized the interpolation and scale part. Nicholas Brown improved ERA5, 
# wrote the CF conventions part, and tested globsim.
#===============================================================================


from globsim.era_interim import ERAIdownload, ERAIinterpolate, ERAIscale
from globsim.era5 import ERA5download, ERA5interpolate, ERA5scale
from globsim.merra_2 import MERRAdownload, MERRAinterpolate, MERRAscale
from globsim.JRA import JRAdownload, JRAinterpolate, JRAscale
from multiprocessing.dummy import Pool as ThreadPool


def GlobsimDownload(pfile, multithread = True, 
                    ERAI=True, ERA5=True, 
                    ERA5ENS=True, MERRA=True, JRA=True):
    """
    Download data from multiple reanalyses. Each reanalysis is run as one 
    separate thread if 'multithread = True'. If 'multithread = False', each
    download is run sequentially. This is easier for interpreting the output.
    """
    # make list of objects to execute
    objects = [] 
    
    # === ERA-Interim ===
    if ERAI:
        ERAIdownl = ERAIdownload(pfile)
        objects.append(ERAIdownl)
    
    # === ERA5 ===
    if ERA5:
        ERA5REAdownl = ERA5download(pfile, 'reanalysis')
        objects.append(ERA5REAdownl)
    
    # === ERA5 10-member ensemble ===
    if ERA5ENS:
        ERA5ENSdownl = ERA5download(pfile, 'ensemble_members')
        objects.append(ERA5ENSdownl)
    
    # === MERRA-2 ===
    if MERRA:
        MERRAdownl = MERRAdownload(pfile)
        objects.append(MERRAdownl)

    # === JRA-55 ===
    if JRA:
        JRAdownl = JRAdownload(pfile)
        objects.append(JRAdownl)

    # serial of parallel execution
    if multithread == True:
        # Make the Pool of workers and run as lambda functions
        pool = ThreadPool(len(objects)) 
        results = pool.map(lambda ob: ob.retrieve(), objects)
        #close the pool and wait for the work to finish 
        pool.close() 
        pool.join() 
        print('Multithreaded download finished')
    else:
        for result in (ob.retrieve() for ob in objects):
            print(result)
        print('Serial download finished')
    
def GlobsimInterpolateStation(ifile, ERAI=True, ERA5=True, ERA5ENS=True, 
                              MERRA=True, JRA=True):
    """
    Interpolate re-analysis data to individual stations (points: lat, lon, ele).
    The temporal granularity and variables of each re-analysis are preserved. 
    The resulting data is saved in netCDF format. THis is not parallelised to 
    differing processors as memory may be a limiting factor.
    """
    
    # === ERA-Interim ===
    if ERAI:
        ERAIinterp = ERAIinterpolate(ifile)
        ERAIinterp.process()
    
    # === ERA5 ===
    if ERA5:
        ERA5interp = ERA5interpolate(ifile)
        ERA5interp.process()
        
    # === ERA5ENS ===
    if ERA5ENS:
        ERA5interp = ERA5interpolate(ifile, 'ensemble_members')
        ERA5interp.process()
    
    # === MERRA-2 ===
    if MERRA:
        MERRAinterp = MERRAinterpolate(ifile)
        MERRAinterp.process()

    # === JRA-55 ===
    if JRA:
        JRAinterp = JRAinterpolate(ifile)
        JRAinterp.process()
  
            
def GlobsimScale(sfile, ERAI=True, ERA5=True, ERA5ENS=True, MERRA=True, JRA=True):
    """
    Use re-analysis data that has been interpolated to station locations to 
    derive fluxes scaled / converted / adjusted to drive point-scale 
    near-surface models. The resulting data has coherent variables for all 
    reanalyses, optionally coherent temporal granularity, and is saved as netCDF.
    """
    # === ERA-Interim ===
    if ERAI:
        ERAIsc = ERAIscale(sfile)
        ERAIsc.process()
    
    # === ERA5 ===
    if ERA5:
        ERA5sc = ERA5scale(sfile)
        ERA5sc.process()
        
    # === ERA5ENS ===
    if ERA5ENS:
        ERA5sc = ERA5scale(sfile, 'ensemble_members')
        ERA5sc.process()
    
    # # === MERRA-2 ===
    if MERRA:
        MERRAsc = MERRAscale(sfile)
        MERRAsc.process()
    #     
    # # === JRA-55 ===
    if JRA:
        JRAsc = JRAscale(sfile)
        JRAsc.process()
                  
