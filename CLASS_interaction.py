# -*- coding: utf-8 -*-
import subprocess
import pandas as pd
from datetime import datetime
from os import path, remove
import numpy as np
import re

class classp():
    """
    This class provides methods for running the Canadian Land Surface Scheme
    (CLASS) from python
         
    """
    SUCCESSFUL_RUN_INDICATOR = "OF1_G"
    
    def __init__(self, jobopt, site_name, executable):
        """
        joboptions file, job site_name and path to CLASS executable
        """
        self.executable = executable
        self.site_name = site_name
        self.jobopt = jobopt
        
    def run(self, site_name=None):
        """
        Run CLASS based on current settings and return True/False on whether
        it completed successfully.
        """
        if site_name is None:
            site_name = self.site_name
              
        # Write inpts
        self.inpts_write()
        
        # Delete successful run file TODO
        if self.successful_run():
            sfile = "{}.{}".format(self.site_name, self.SUCCESSFUL_RUN_INDICATOR)
            remove(sfile) 

        # run model
        term = subprocess.check_output([self.executable, self.jobopt, site_name])

        # Return True / False if completed successfully
        return self.successful_run()
    
    def successful_run(self):
        """
        Check for existence of output file indicating a successful run.
        """   
        f =  "{}.{}".format(self.site_name, self.SUCCESSFUL_RUN_INDICATOR)
        return path.isfile(f)
        
    def readINI(self, file):
        inifile = "{}.INI".format(self.site_name)
        return INI(inifile).readINI().params
        
    def inpts_write(self):
        pass
    
    @staticmethod
    def ts_read(ascii):
        ''' reads ascii files corresponding to timeseries outputs of class model '''
        # open file
        raw = pd.read_csv(ascii, skiprows=2, delim_whitespace=True)
        
        # create datetime column
        datestr = raw['YEAR'].map(lambda x: "{:04}".format(x)) + "-" + \
                    raw['DAY'].map(lambda x: "{:02}".format(x)) + " " + \
                    raw['HOUR'].map(lambda x: "{:02}".format(x)) + ":" + \
                    raw['MIN'].map(lambda x: "{:02}".format(x)) 
        
        raw["Date"] = datestr.map(lambda x: datetime.strptime(x, "%Y-%j %H:%M"))
        raw.drop(['HOUR', 'MIN', 'DAY', 'YEAR'], 1, inplace=True)
        
        # reorder and return
        output = raw[['Date'] + list(raw.columns[:-1])]

        return(output)

class INI():
    ''' reads a class INI file into a dictionary ''' 
    
    def __init__(self, inifile):
        self.file = inifile
        self.readINI()
    
    def __repr__(self):
        x = self.columndefs()
        defs = list(x.keys() )
        defs.sort()
        
        print('{:*^79}'.format(' GLOBAL PARAMETERS '))
        for key in defs:
            if key in self.params.keys():
                print("{}: {}".format(key, x[key]))
                print(self.params[key])
        
        print('{:*^79}'.format(' MOSAIC PARAMETERS '))
        
        M = -1 
        
        while True:
            M += 1
            if not "M{}".format(M) in self.params.keys():
                print('{:-^79}'.format(' End of Mosaics '))
                break
                
            else:
                mos = "M{}".format(M)
                print('{:-^79}'.format(' Mosaic ' + mos + ' '))
                
                for key in defs:
                    if key in self.params[mos].keys():
                        print("{}: {}".format(key, x[key]))
                        print(self.params[mos][key])
        return('')

    def nlayers(self, raw_ini): 
        ''' determines how many soil layers are represented in an INI file'''
        ix = -3    # starting index (start at the end and work backwards)
        n_layers = 0
        
        # check how many lines before 3rd from last have length == 2
        while True:
            if len(raw_ini[ix]) == 2:
                n_layers += 1
                ix -= 1    
            
            else:
                break
                
        return(n_layers)
        
    def inpts_read(self, inpts_file):
        ''' reads keyword inputs into INI dictionary.'''
        pass
    
    def readINI(self):
        ''' 
        read CLASS *.INI file into a dictionary, accounting for the possibility
        of more than 3 layers
        '''
        with open(self.file) as f:
            raw = f.readlines()
            raw = [x.strip() for x in raw]
            raw[3:] = [x.split() for x in raw[3:]]
            raw[3:] = [[float(x) for x in line]  for line in raw[3:]]
        
        self.raw = raw 
        
        # determine number of soil layers
        n_lyr = self.nlayers(raw)
        
        ## Assume ICAN is 4 (default)
        ICAN = 4

        # determine the number of mosaicks
        # 6 non-mosaic lines (excluding soil thicknesses), 15 rows per mosaic
        n_mos = int((len(raw) - (6 + n_lyr)) / 15) 

        x = dict()

        # first 3 lines are general information
        x['info1'] = raw[0]
        x['info2'] = raw[1]
        x['info3'] = raw[2]
        
        # first line of actual data
        x['DEGLAT'] = [raw[3][0]]
        x['DEGLON'] = [raw[3][1]]
        x['ZRFM']   = [raw[3][2]]
        x['ZRFH']   = [raw[3][3]]
        x['ZBLD']   = [raw[3][4]]
        x['GC']     = [raw[3][5]]
        x['NLTEST'] = [raw[3][6]]
        x['NMTEST'] = [raw[3][7]]

        # begin mosaics        
        for M in range(n_mos):
            
            mos = {} # create empty dictionary for each separate mosaic
            d = M * 15 # offset
            
            mos['FCAN'] = raw[4 + d][0:(ICAN + 1)]
            mos['PAMX'] = raw[4 + d][(ICAN + 1):(2 * ICAN + 1)]
            
            mos['LNZ0'] = raw[5 + d][0:(ICAN + 1)]
            mos['PAMN'] = raw[5 + d][(ICAN + 1):(2 * ICAN + 1)]
            
            mos['ALVC'] = raw[6 + d][0:(ICAN + 1)]
            mos['CMAS'] = raw[6 + d][(ICAN + 1):(2 * ICAN + 1)]
            
            mos['ALIC'] = raw[7 + d][0:(ICAN + 1)]
            mos['ROOT'] = raw[7 + d][(ICAN + 1):(2 * ICAN + 1)]
            
            mos['RSMN'] = raw[8 + d][0:ICAN]
            mos['QA50'] = raw[8 + d][ICAN:(2 * ICAN)]
            
            mos['VPDA'] = raw[9 + d][0:ICAN]
            mos['VPDB'] = raw[9 + d][ICAN:(2 * ICAN)]
            
            mos['PSGA'] = raw[10 + d][0:ICAN]
            mos['PSGB'] = raw[10 + d][ICAN:(2 * ICAN)]
            
            mos['DRN']  = [raw[11 + d][0]]
            mos['SDEP'] = [raw[11 + d][1]]
            mos['FARE'] = [raw[11 + d][2]]
            
            mos['XSLP'] = [raw[12 + d][0]]         
            mos['GRKF'] = [raw[12 + d][1]] # unused
            mos['WFSF'] = [raw[12 + d][2]] # unused
            mos['WFCI'] = [raw[12 + d][3]] # unused         
            mos['MID']  = [raw[12 + d][4]]
            
            mos['SAND'] = raw[13 + d]
            mos['CLAY'] = raw[14 + d] 
            mos['ORGM'] = raw[15 + d]
            
            mos['TBAR'] = raw[16 + d][0:n_lyr]
            mos['TCAN'] = raw[16 + d][n_lyr:(n_lyr + 1)]
            mos['TSNO'] = raw[16 + d][(n_lyr + 1):(n_lyr + 2)]
            mos['TPND'] = raw[16 + d][(n_lyr + 2):(n_lyr + 3)]
            
            mos['THLQ'] = raw[17 + d][0:n_lyr]
            mos['THIC'] = raw[17 + d][n_lyr:(2 * n_lyr)]
            mos['ZPND'] = raw[17 + d][(2 * n_lyr):(2 * n_lyr + 1)]
           
            mos['RCAN'] = [raw[18 + d][0]]
            mos['SCAN'] = [raw[18 + d][1]]
            mos['SNO']  = [raw[18 + d][2]]
            mos['ALBS'] = [raw[18 + d][3]]
            mos['RHOS'] = [raw[18 + d][4]]
            mos['GRO']  = [raw[18 + d][5]]
            
            x['M{}'.format(M)] = mos            
        
        # end mosaics
        
        mos = (n_mos - 1) * 15 # extra mosaic index offset
        los = (n_lyr - 3)      # extra layer index offset
        
        x['DELZ'] = [z[0] for z in raw[(19 + mos):(19 + n_lyr + mos)]]
        x['ZBOT'] = [z[1] for z in raw[(19 + mos):(19 + n_lyr + mos)]]
        
        x['hrly_out_days']   = raw[22 + mos + los][0:2]
        x['hrly_out_years']  = raw[23 + mos + los][0:2]
        x['dly_out_days']    = raw[22 + mos + los][2:4]
        x['dly_out_years']   = raw[23 + mos + los][2:4]
        
        self.params = x
        
    def writeINI(self, out_file):
        '''
        write stored dictionary to an *.INI file for CLASS, using the
        proper formatting.  Automatically detects the number of layers and
        the number of mosaics
        '''
        ini = self.params
        
        # get number of soil layers and mosaics
        n_lyr = len(ini['DELZ'])
        n_mos = sum([1 if re.match("M\\d", x) else 0 for x in ini.keys()])
        
        # Write lines to output file
        with open(out_file, 'w') as f:
            
            # define writer to automatically add line return
            writenew = lambda line: f.writelines(line + "\n")
                 
            # write top 3 info lines
            writenew("  {}".format(ini['info1']))
            writenew("  {}".format(ini['info2']))
            writenew("  {}".format(ini['info3']))
            
            # write first data line
            L1 = "{:10.2f}{:10.2f}{:9.2f}{:9.2f}{:10.2f}{:7.1f}{:5.0f}{:5.0f}"
            writenew(L1.format( ini['DEGLAT'][0], ini['DEGLON'][0], 
                                    ini['ZRFM'][0],   ini['ZRFH'][0], 
                                    ini['ZBLD'][0] ,  ini['GC'][0],  
                                    ini['NLTEST'][0], ini['NMTEST'] [0]))

            # BEGIN MOSAICS
            
            # define formatting strings for mosaics
            F1 = "{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}"
            F2 = "{:8.3f}{:8.3f}{:8.3f}{:8.3f}{:16.3f}{:8.3f}{:8.3f}{:8.3f}"
            F3 = "{:8.3f}{:8.3f}{:8.3f}"
            F4 = "{:8.1E}{:8.1E}{:8.1E}{:8.1E}{:8.0f}"
            F5 = "{:10.1f}" * n_lyr
            F6 = "{:10.2f}" * (n_lyr + 3)
            F7 = "{:10.3f}" * (2 * n_lyr + 1)
            F8 = "{:10.4f}{:10.4f}{:10.1f}{:10.3f}{:10.4f}{:10.3f}"

            # write mosaics

            for M in range(n_mos):
                
                # subset main dictionary for this mosaic
                D = ini['M{}'.format(M)]

                writenew(F1.format(*(D['FCAN'] + D['PAMX'])))
                writenew(F1.format(*(D['LNZ0'] + D['PAMN'])))
                writenew(F1.format(*(D['ALVC'] + D['CMAS'])))
                writenew(F1.format(*(D['ALIC'] + D['ROOT'])))

                writenew(F2.format(*(D['RSMN'] + D['QA50'])))
                writenew(F2.format(*(D['VPDA'] + D['VPDB'])))
                writenew(F2.format(*(D['PSGA'] + D['PSGB'])))

                writenew(F3.format(*(D['DRN'] + D['SDEP'] + D['FARE'])))

                LF4 = D['XSLP'] + D['GRKF'] + D['WFSF'] + D['WFCI'] + D['MID']
                writenew(F4.format(*LF4))
             
                # write sand / silt / clay thicknesses
                writenew(F5.format(*D['SAND']))
                writenew(F5.format(*D['CLAY']))
                writenew(F5.format(*D['ORGM']))

                # write starting temperatures
                LF6 = D['TBAR'] + D['TCAN'] + D['TSNO'] + D['TPND']
                writenew(F6.format(*LF6))

                # write starting moisture
                LF7 = D['THLQ'] + D['THIC'] + D['ZPND']
                writenew(F7.format(*LF7))

                # write snow & vegetation
                LF8 = D['RCAN'] + D['SCAN'] + D['SNO'] + D['ALBS'] + D['RHOS'] + D['GRO']
                writenew(F8.format(*LF8))
                
                ## END MOSAIC

       
            # write soil layers thickness / bottom depths

            for (dz, zb) in zip(ini['DELZ'], ini['ZBOT']):
                writenew("  {:8.2f}{:8.2f}".format(dz, zb))
            
            # write dates (last and second last lines)

            Ls2 = "{:10.0f}{:10.0f}{:10.0f}{:10.0f}"
            Ls2_vals = ini['hrly_out_days'] + ini['dly_out_days']
            writenew(Ls2.format(*Ls2_vals))
            
            Ls1 = "{:10.0f}{:10.0f}{:10.0f}{:10.0f}"
            Ls1_vals = ini['hrly_out_years'] + ini['dly_out_years']
            writenew(Ls1.format(*Ls1_vals))

    
    @staticmethod
    def columndefs():
        defs = {
        'DEGLAT' : 'Latitude',
        'DEGLON' : 'Longitude',
        'ZRFM'   : 'Ref.height',
        'ZRFH'   : 'Ref.height',
        'ZBLD'   : 'Blend ht.',
        'NLTEST' : '???',
        'NMTEST' : '???',
        'GC'     : 'Land/sea cover',
        'FCAN'   : 'Fractional coverage for CLASS’ 4 PFTs+Urban',
        'PAMX'   : 'Maximum LAI for CLASS’ 4 PFTs',
        'LNZ0'   : 'Log of roughness length(m) for 4 PFTs+Urban',
        'PAMN'   : 'Minimum LAI for 4 PFTs',
        'ALVC'   : 'Visible albedo for 4 PFTs + Urban',
        'CMAS'   : 'Aboveground canopy mass for 4 PFTs',
        'ALIC'   : 'Near-infrared albedo for 4 PFTs + Urban',
        'ROOT'   : 'Rooting depth for 4 PFTs',
        'RSMN'   : 'Minimum stomatal resistance for 4 PFTs  ',
        'QA50'   : 'Reference value of SW radiation for 4 PFTs',
        'VPDA'   : 'First of two vapor pressure deficit coefficients for 4 PFTs for CLASS’ stomatal resistance formulation',
        'VPDB'   : 'Second of two vapor pressure deficit coefficients for 4 PFTs for CLASS’ stomatal resistance formulation',
        'PSGA'   : 'First of two soil moisture suction coefficients for 4 PFTs',
        'PSGB'   : 'Second of two soil moisture suction coefficients for 4 PFTs',
        'DRN'    : 'Drainage',
        'SDEP'   : 'Soil permeable depth (m)  ',
        'FARE'   : 'Fractional coverage of mosaic tile on modeled area',
        'XSLP'   : 'Surface  slope  ',
        'GRKF'   : 'parameter not used',
        'WFSF'   : 'parameter not used',
        'WFCI'   : 'parameter not used',
        'MID'    : 'Mosaic tile type identifier (1 for land surface, 0 for inland lake)',
        'SAND'   : 'Percentage SAND in each of the soil layers',
        'CLAY'   : 'Percentage CLAY in each of the soil layers',
        'ORGM'   : 'Percentage ORGANIC MATTER in each of the soil layers',
        'TBAR'   : 'Temperatures of each soil layer (ºC)',
        'TCAN'   : 'Temperature of vegetation canopy  (ºC)',
        'TSNO'   : 'Temperature of snow (ºC)',
        'TPND'   : 'Temperature of ponded water (ºC)',
        'THLQ'   : 'Soil moisture  of each soil layer  (m3 m-3)',
        'THIC'   : 'Frozen soil moisture of each soil layer  (m3 m-3)',
        'ZPND'   : 'Depth of ponded water (m)',
        'RCAN'   : 'Water on canopy',
        'SCAN'   : 'Snow on canopy',
        'SNO'    : 'Snow on ground',
        'ALBS'   : 'Snow albedo',
        'RHOS'   : 'Snow density',
        'GRO'    : 'Plant growth index',
        'DELZ'   : 'Thickness of soil layers (m)',
        'ZBOT'   : 'Depth at bottom of soil layers (m)',
        }
        return(defs)


def Gonzalo(n): 
    ''' Soil layer thicknesses defined by Sapriza-Azuri et al. (2018) ''' 
    
    x = [0.1, 0.35, 0.86, 1.77, 3.09, 4.81, 6.9, 9.5, 12.4, 15.8, 
                19.5, 23.7, 28.2, 33.2, 38.6, 44.4, 50.6, 57.2, 64.2, 71.6]
    x = x[0:n]
    ZBOT = np.array(x)            
    DELZ = np.diff([0] + x)
    out = dict(ZBOT = ZBOT, DELZ = DELZ)
    return(out)
    



# Testing
# 
# if __name__ == "__main__":
#     CL = classp("/home/nbrown/class/test1/job_options_LOOP.txt" ,
#     "/home/nbrown/class/test1/2005_ST02_test1", executable="/home/himani/src/bin/CLASS36CTEM" )
#     x = CL.run()
#     print("the result is: {}".format(x))