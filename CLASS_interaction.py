import subprocess
import pandas as pd
from datetime import datetime
from os import path, remove
import numpy as np

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
    
    def readINI(self):
        with open(self.file) as f:
            raw = f.readlines()
            raw = [x.strip() for x in raw]
            raw[3:] = [x.split() for x in raw[3:]]
        
        self.raw = raw 
        
        # determine number of soil layers
        n_lyr = self.nlayers(raw)
        
        ## Assume ICAN is 4 (default)
        ICAN = 4

        # determine the number of mosaicks
        # 6 non-mosaic lines (excluding soil thicknesses), 15 rows per mosaic
        n_mos = int((len(raw) - (6 + n_lyr)) / 15) 

        x = dict()
        
        # first line
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
            # skip 3 unused parameters
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

        os = (n_mos - 1) * 15 # index offset
        
        x['DELZ'] = [z[0] for z in raw[(19 + os):(22 + os)]]
        x['ZBOT'] = [z[1] for z in raw[(19 + os):(22 + os)]]
        
        x['hrly_out_days']   = raw[22 + os][0:2]
        x['hrly_out_years']  = raw[23 + os][0:2]
        x['dly_out_days']    = raw[22 + os][2:4]
        x['dly_out_years']   = raw[23 + os][2:4]
        
        self.params = x

## Soil layer thicknesses defined by Sapriza-Azuri et al. (2018)
def Gonzalo(n):  
    x = [0.1, 0.35, 0.86, 1.77, 3.09, 4.10, 4.81, 6.9, 9.5, 12.4, 15.8, 
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