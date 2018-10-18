import subprocess

from os import path, remove

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


class INI():
    
    def __init__(self, inifile):
        self.file = inifile
        
    def readINI(self):
        with open(self.file) as f:
            raw = f.readlines()
            raw = [x.strip() for x in raw]
            raw[3:] = [x.split() for x in raw[3:]]
        
        self.raw = raw 
        
        x = dict()
        
        # first line
        x['DEGLAT'] = raw[3][0]
        x['DEGLON'] = raw[3][1]
        x['ZRFM']   = raw[3][2]
        x['ZRFH']   = raw[3][3]
        x['ZBLD']   = raw[3][4]
        x['GC']     = raw[3][5]
        x['NLTEST'] = raw[3][6]
        x['NMTEST'] = raw[3][7]

        # begin mosaics 
        n_mos = int((len(raw) - 9) / 15)
        x['mosaic'] = int(n_mos == 1)
        
        for M in range(n_mos):
            mos = {}

            mos['FCAN'] = None
            mos['PAMX'] = None
            mos['LNZ0'] = None
            mos['PAMN'] = None
            mos['ALVC'] = None
            mos['CMAS'] = None
            mos['ALIC'] = None
            mos['ROOT'] = None
            mos['RSMN'] = None
            mos['QA50'] = None
            mos['VPDA'] = None
            mos['VPDB'] = None
            mos['PSGA'] = None
            mos['PSGB'] = None
            mos['DRN']  = None
            mos['SDEP'] = None
            mos['FARE'] = None
            mos['XSLP'] = None
            # skip 3 unused parameters
            mos['MID']  = None
            mos['SAND'] = None
            mos['CLAY'] = None
            mos['ORGM'] = None
            mos['TBAR'] = None
            mos['TCAN'] = None
            mos['TSNO'] = None
            mos['TPND'] = None
            mos['THLQ'] = None
            mos['THIC'] = None
            mos['ZPND'] = None
            mos['RCAN'] = None
            mos['SCAN'] = None
            mos['SNO']  = None
            mos['ALBS'] = None
            mos['RHOS'] = None
            mos['GRO']  = None
            
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
    
    
# Testing

if __name__ == "__main__":
    CL = classp("/home/nbrown/class/test1/job_options_LOOP.txt" ,
    "/home/nbrown/class/test1/2005_ST02_test1", executable="/home/himani/src/bin/CLASS36CTEM" )
    x = CL.run()
    print("the result is: {}".format(x))