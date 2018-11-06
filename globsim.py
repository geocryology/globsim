#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Generic classes, methods, functions used for more than one reanalysis.
#
#
# (C) Copyright Stephan Gruber (2017)
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
#===============================================================================
from era_interim import ERAdownload, ERAinterpolate, ERAscale
from merra_2 import MERRAdownload, MERRAinterpolate, MERRAscale 
from JRA import JRAdownload, JRAinterpolate, JRAscale
from multiprocessing.dummy import Pool as ThreadPool
from era5 import ERAdownload as ERA5download

def GlobsimDownload(pfile, multithread = True, ERAI=True, ERA5=True, MERRA=True, JRA=True):
    """
    Download data from multiple re-analyses. Each re-analysis is run as one 
    separate thread if 'multithread = True'. If 'multithread = False', each
    download is run sequentially. This is easier for interpreting the output.
    """
    # make list of objects to execute
    objects = [] 
    
    # === ERA-Interim ===
    if ERAI:
        ERAdownl = ERAdownload(pfile)
        objects.append(ERAdownl)
    
    # === ERA-5 ===
    if ERA5:
        ERA5downl = ERA5download(pfile)
        objects.append(ERA5downl)
    
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
    
def GlobsimInterpolateStation(ifile):
    """
    Interpolate re-analysis data to individual stations (points: lat, lon, ele).
    The temporal granularity and variables of each re-analysis are preserved. 
    The resulting data is saved in netCDF format. THis is not parallelised to 
    differing processors as memory may be a limiting factor.
    """
    
    # === ERA-Interim ===
    ERAinterp = ERAinterpolate(ifile)
    ERAinterp.process()
    
    # === MERRA-2 ===
    MERRAinterp = MERRAinterpolate(ifile)
    MERRAinterp.process()

    # === JRA-55 ===
    JRAinterp = JRAinterpolate(ifile)
    JRAinterp.process()
  
            
def GlobsimScale(sfile):
    """
    Use re-analysis data that has been interpolated to station locations to 
    derive fluxes scaled / converted / adjusted to drive point-scale 
    near-surface models. The resulting data has coherent variables for all 
    reanalyses, optionally coherent temporal granularity, and is saved as netCDF.
    """
    # === ERA-Interim ===
    ERAsc = ERAscale(sfile)
    ERAsc.process()
    
    # # === MERRA-2 ===
    # MERRAsc = MERRAscale(sfile)
    # MERRAsc.process()
    #     
    # # === JRA-55 ===
    # JRAsc = JRAscale(sfile)
    # JRAsc.process()
                  