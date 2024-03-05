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

from multiprocessing.dummy import Pool as ThreadPool
from globsim.LazyLoader import LazyLoader
from globsim.download.era5_monthly import download_threadded

download = LazyLoader('globsim.download')
interpolate = LazyLoader('globsim.interpolate')
scale = LazyLoader('globsim.scale')

# from globsim.download import *
# from globsim.scale import *
# from globsim.interpolate import *

# from globsim.JRA import JRAdownload, JRAinterpolate, JRAscale


def GlobsimDownload(pfile, multithread=True, 
                    ERAI=True, ERA5=True, 
                    ERA5ENS=True, MERRA=True, JRA=True, JRA3Q=True):
    """
    Download data from multiple reanalyses. Each reanalysis is run as one 
    separate thread if 'multithread = True'. If 'multithread = False', each
    download is run sequentially. This is easier for interpreting the output.
    """
    # make list of objects to execute
    objects = [] 
    
    # === ERA-Interim ===
    if ERAI:
        ERAIdownl = download.ERAIdownload(pfile)
        objects.append(ERAIdownl)
    
    # === ERA5 ===
    if ERA5:
        ERA5REAdownl = download.ERA5MonthlyDownload(pfile, False)
        objects.append(ERA5REAdownl)
    
    # === ERA5 10-member ensemble ===
    if ERA5ENS:
        ERA5ENSdownl = download.ERA5Edownload(pfile)
        objects.append(ERA5ENSdownl)
    
    # === MERRA-2 ===
    if MERRA:
        MERRAdownl = download.MERRAdownload(pfile)
        objects.append(MERRAdownl)

    # === JRA-55 ===
    if JRA:
        JRAdownl = download.JRAdownload(pfile)
        objects.append(JRAdownl)

    # === JRA-3Q ===
    if JRA3Q:
        JRA3Qdownl = download.J3QD(pfile)
        objects.append(JRA3Qdownl)

    # serial of parallel execution
    if multithread:
        # Make the Pool of workers and run as lambda functions
        pool = ThreadPool(len(objects)) 
        results = pool.map(lambda ob: ob.retrieve(), objects)
        # close the pool and wait for the work to finish 
        pool.close()
        pool.join()
        print('Multithreaded download finished')
    else:
        for result in (ob.retrieve() for ob in objects):
            print(result)
        print('Serial download finished')
    

def GlobsimInterpolateStation(ifile, ERAI=True, ERA5=True, ERA5ENS=True,
                              MERRA=True, JRA=True, JRA3Q=True, **kwargs):
    """
    Interpolate re-analysis data to individual stations (points: lat, lon, ele).
    The temporal granularity and variables of each re-analysis are preserved.
    The resulting data is saved in netCDF format. THis is not parallelised to
    differing processors as memory may be a limiting factor.
    """
    
    # === ERA-Interim ===
    if ERAI:
        ERAIinterp = interpolate.ERAIinterpolate(ifile, **kwargs)
        ERAIinterp.process()
    
    # === ERA5 ===
    if ERA5:
        ERA5interp = interpolate.ERA5interpolate(ifile, **kwargs)
        ERA5interp.process()
        
    # === ERA5ENS ===
    if ERA5ENS:
        ERA5interp = interpolate.ERA5interpolate(ifile, 'ensemble_members', **kwargs)
        ERA5interp.process()
    
    # === MERRA-2 ===
    if MERRA:
        MERRAinterp = interpolate.MERRAinterpolate(ifile, **kwargs)
        MERRAinterp.process()

    # === JRA-55 ===
    if JRA:
        JRAinterp = interpolate.JRAinterpolate(ifile)
        JRAinterp.process()
  
    if JRA3Q:
        JRA3Qinterp = interpolate.J3QI(ifile)
        JRA3Qinterp.process()
            
def GlobsimScale(sfile, ERAI=True, ERA5=True, ERA5ENS=True, MERRA=True, JRA=True, JRA3Q=True):
    """
    Use re-analysis data that has been interpolated to station locations to 
    derive fluxes scaled / converted / adjusted to drive point-scale 
    near-surface models. The resulting data has coherent variables for all 
    reanalyses, optionally coherent temporal granularity, and is saved as netCDF.
    """
    # === ERA-Interim ===
    if ERAI:
        ERAIsc = scale.ERAIscale(sfile)
        ERAIsc.process()
    
    # === ERA5 ===
    if ERA5:
        ERA5sc = scale.ERA5scale(sfile)
        ERA5sc.process()
        
    # === ERA5ENS ===
    if ERA5ENS:
        ERA5sc = scale.ERA5Escale(sfile)
        ERA5sc.process()
    
    # # === MERRA-2 ===
    if MERRA:
        MERRAsc = scale.MERRAscale(sfile)
        MERRAsc.process()
    #     
    # # === JRA-55 ===
    if JRA:
        JRAsc = scale.JRAscale(sfile)
        JRAsc.process()

    # # === JRA-3Q ===
    if JRA3Q:
        JRA3Qsc = scale.J3QS(sfile)
        JRA3Qsc.process()
                  
