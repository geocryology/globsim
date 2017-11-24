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
from merra_2 import MERRAdownload, MERRAinterpolate 
from JRA import JRAdownload

def GlobsimDownload(pfile):
    """
    Download data from multiple re-analyses.
    """
    #TODO: make each re-analysis one sub-process
    
    # # === ERA-Interim ===
    ERAdownl = ERAdownload(pfile)
    ERAdownl.retrieve()
    
    # === MERRA-2 ===
    MERRAdownl = MERRAdownload(pfile)
    MERRAdownl.retrieve()

    # === JRA-55 ===
    JRAdownl = JRAdownload(pfile)
    JRAdownl.retrieve()

def GlobsimInterpolateStation(ifile):
    """
    Interpolate re-analysis data to individual stations (points: lat, lon, ele).
    The temporal granularity and variables of each re-analysis are preserved. 
    The resulting data is saved in netCDF format.
    """
    
    # === ERA-Interim ===
    ERAinterp = ERAinterpolate(ifile)
    ERAinterp.process()
    
    # === MERRA-2 ===
    # MERRAinterp = MERRAinterpolate(ifile)
    # MERRAinterp.process()

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
    
    # === MERRA-2 ===
    #TODO
    
    # === JRA-55 ===
    JRAsc = JRAscale(sfile)
    JRAsc.process()
                  