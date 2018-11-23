import re

import numpy as np
import pandas as pd

from geotop import geotop_interaction as GI
from scipy.interpolate import RectBivariateSpline
from os import path
from .CLASS_interaction import classp, INI


class SoilTemperature:
    '''
    A representation of soil temperature arrays with methods for
    calculating associated parameters such as MAGST. Relies on a bivariate
    spline interpolation to calculate averages for any desired depth.
    
    Instances can be generated from several different output file types (e.g. geotop or 
    class output files)
    '''
    def __init__(self):
        pass
        
    @classmethod
    def from_OF5_G(cls, OF5_G):
        ''' load temperature grid from class output file'''

        # get soil layer midpoint depths from INI file
        ini = re.sub("OF5_G", "INI", OF5_G) 
        ini = INI(ini)
        ini.readINI()
        DELZ = ini.params['DELZ']
        nodes = np.array([float(x) * 1000 for x in DELZ])
        nodes = np.cumsum(nodes) - nodes / 2
        n_layer = len(DELZ)
        
        # open OF5_G
        data = classp.ts_read(OF5_G)
        
        # also open OF6_G if there are more than 3 soil layers
        if n_layer > 3:
            OF6_G = re.sub("OF5_G$", "OF5_G", OF5_G) 
            OF6_G = classp.ts_read(OF6_G)
            OF6_G.drop(['Date'], inplace = True)
            
            # we assume that the dates are the same and combine the arrays
            data = pd.concat(data, OF6_G)
        

        # subset columns to include only temperature
        for column in data.columns[1:]:
            if not re.match("TG\\d{1,2}", column):
                data.drop(column, 1, inplace=True)

        # rename column headers to cell midpoint depth (mm)
        data.columns = ['Date'] + [str(n) for n in nodes]
        
        # create soil temperature object
        ST = cls()
       
        # interpolate
        x = [float(d) for d in data.columns[1:]] # depths
        y = [i for i in data.index]              # time index
        z = np.array(data.drop('Date', 1))
        
        if n_layer <= 3:
            # use a 2-point spline instead of 3 point (not enough data for 3)
            S = RectBivariateSpline(x=y, y=x, z=z, ky=2) 
        else:
            S = RectBivariateSpline(x=y, y=x, z=z)
        # add reference info
        ST.dates = data['Date']
        ST.depths = x
        ST.data = z
        ST.grid = S
        
        ST.model = 'class'
        return(ST)
    
    @classmethod
    def from_geotop_file(cls, gt_file):
        ''' load temperature from geotop file'''
        reader = GI.DirectoryReader(directory='', inpts='')
        data = reader.ascii_read(gt_file)
        drop_columns = ['Simulation_Period', 'Run', 'IDpoint']
        for col in drop_columns:
            if col in data.columns:
                data.drop(col, 1, inplace=True)

        # interpolate
        x = [float(d) for d in data.columns[1:]] # depths
        y = [i for i in data.index]              # time index
        z = np.array(data.drop('Date', 1))
        S = RectBivariateSpline(x=y, y=x, z=z)
    
        # create object
        ST = cls()
        
        # add reference info
        ST.dates = data['Date']
        ST.depths = x
        ST.data = z
        ST.grid = S
        
        ST.model = 'geotop'
        return(ST)
    
    def ALT(self, date=None, prec=10):
        ''' return ALT at given date. Defaults to end-of-simulation date'''
        
        # default to last available date
        if date is None:
            date = max(self.dates)
        
        # minimum and maximum depth to check
        min_d = np.min(self.depths)
        max_d = np.max(self.depths)
        
        ALT = min_d
        # from the lowest depth, check each lower depth to see if permafrost
        for dpth in np.arange(min_d, max_d, prec):
            maxT = self.Tstat_last_n_years(depth=dpth, date=date, years=2, stat=np.max)

            if maxT < 0: # permafrost condition
                break
            else: 
                ALT += prec  # ALT depth
        
        if ALT == min_d:
            raise Exception("[WARNING]: ALT calculation not possible due to deep minimum observation depth")
            return(-999)
        else:
            return(ALT)

    def Tstat_last_n_years(self, depth, date=None, years=1, stat=None):
        '''
        return ground temperature statistic for the last 'n' years
        prior to target date at specified depth. Date must be a 
        pandas.Timestamp object e.g. pd.Timestamp(2014, 1, 1, 12)
        '''
        
        # default to last available date
        if date is None:
            date = max(self.dates)
            
        # default statistic is mean
        if stat is None:
            stat = np.mean
        
        # establish date range 
        start = date - pd.DateOffset(years=years)
        ix = (self.dates[(self.dates > start) & (self.dates <= date)]).index
        
        # compute statistic at desired depth using interpolation mesh
        T = self.grid(ix, depth)
        MT = stat(T)
        
        return(MT)
    
    def MAGST(self, date=None, years=1):
        ''' a wrapper for Tstat_last_n_years to produce MAGST (5cm depth)''' 
        d_MAGST = 50
        min_dep = np.min(self.depths)
        cutoff = 50

        if min_dep > cutoff:
            raise Exception("Shallowest output depth is {} mm. A shallower depth"
            "is required to calculate MAGST ({} mm)".format(min_dep, cutoff))

        MAGST = self.Tstat_last_n_years(depth=d_MAGST, date=date, years=years)       
        
        return(MAGST)
        
    def ZAA_depth(self, date=None, cutoff=0.1):
        '''return depth of  zero-annual amplitude'''
        pass

    def surface_offset(self, meteo_path, date, years=1):
        ''' return surface offset, given a meteo file '''
        # read air_temp file 
        self.read_meteofile(meteo_path, model=self.model)
        
        # make sure dates are equivalent! TODO[NB]
        
        # calculate MAAT 
        MAAT = self.MAAT(date=date, years=years)
        
        # calculate MAGST
        MAGST = self.MAGST(date=date, years=years)
        
        # calculate offset
        SO = MAAT - MAGST
        
        return(SO)
    
    def read_meteofile(self, meteo_path, model='geotop'):
        ''' generic handle to read different kinds of meteo files from various models'''

        if model == 'geotop':
            met = self.read_geotop_forcing(meteo_path)
        else:
            raise NotImplementedError("no such model type")
        
        self.met = met
    
    def read_geotop_forcing(self, directory):
        '''reads geotop forcing file. relies on inpts file within directory'''
        # read data
        if path.isfile(directory):
            raise Exception("must provide a geotop directory, not a file")

        G = GI.Gtpy(directory)
        forc = G.read_forcing()
        
        # change header names so that they match header keywords
        sdct = {G.inpts[key][1:-1]:key for key in G.inpts if "Header" in key}
        forc.columns = ['Date']+[sdct[key] for key in forc.columns[1:]]
        
        return(forc)
    
    def MAAT(self, date=None, years=1):
        # default to last available date
        if date is None:
            date = max(self.dates)
        
        # establish date range 
        start = date - pd.DateOffset(years=years)
        
        # caclculate mean
        metdat = self.met['Date']
        data = self.met['HeaderAirTemperature'][(metdat > start) & (metdat <= date)]
        MAAT = np.mean(data)
        
        return(MAAT)


# R = SoilTemperature.from_OF5_G(r"Q:\Projects\GLOBSIM\jobs3\test4\1009_ST03_test4.OF5_G")   
# 
# R = SoilTemperature.from_geotop_file(r"E:\Users\Nick\Desktop\temp\E1\soil.txt")    
#  
# R.read_meteofile(r"E:\Users\Nick\Desktop\temp\E1")
# 
# R.Tstat_last_n_years(6500, date = pd.Timestamp(2009, 1, 1, 12), years=1)
# R.ALT()
# R.MT_last_n_years(1500, date = pd.Timestamp(2009, 1, 1, 12), years=1, stat=np.max)
