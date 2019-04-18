#
# Provides methods to create directories which can be called by models.
# such as GEOtop or CLASS
#
# Also provides methods to batch-create directories from a csv
#

import shutil
import os 
import time
import glob

import pandas as pd

from globsim.generic import ParameterIO
from globsim.CLASS_interaction import classp, INI
from shutil import rmtree, copyfile

try:
    from geotop.geotop_interaction import DirectoryReader
except:
    print("Missing module: geotop. Clone it from github.com/geocryology/geotop")


class Jobs():
    
    def __init__(self, model, forcing):
        pass
   
class BaseModelDir(object):
    ''' shared code for creating CLASS and GeoTOP directories'''
    def __init__(self, rootdir, jobname, overwrite=True):
        self.rootdir = rootdir
        self.jobname = jobname
        
        # create new job directory, overwriting if appropriate
        self.jobdir = os.path.join(rootdir, jobname)
        if os.path.isdir(self.jobdir) and overwrite:
            rmtree(self.jobdir)
            time.sleep(0.1)
        os.mkdir(self.jobdir)

class CLASSDirectory(BaseModelDir):
    """ Create a directory containing files to run the CLASS model"""
    def __init__(self, rootdir, jobname, overwrite=True):
        super(CLASSDirectory, self).__init__(rootdir, jobname, overwrite=True)
    
    @classmethod
    def from_file_paths(cls, rootdir, jobname, MET, INI, CTM, LOOP, **kwargs):
        ''' initialize using direct paths to data files '''
        CLD = cls(rootdir, jobname, **kwargs)
        CLD.MET  = MET
        CLD.INI  = INI
        CLD.CTM  = CTM
        CLD.LOOP = LOOP
        return(CLD)
        
    def write_parameters(self):
        out_ini = os.path.join(self.jobdir, "{}.INI".format(self.jobname))
        copyfile(self.INI, out_ini)
        
        out_ctm = os.path.join(self.jobdir, "{}.CTM".format(self.jobname))
        copyfile(self.CTM, out_ctm)

        out_loop = os.path.join(self.jobdir, "job_options_LOOP.txt".format(self.jobname))
        copyfile(self.LOOP, out_loop)

    def write_forcing(self):
        out_MET = os.path.join(self.jobdir, "{}.MET".format(self.jobname))
        copyfile(self.MET, out_MET)

    def build(self):
        # write parameters files (CTM, LOOP, INI)
        self.write_parameters()
        
        # write forcing file (MET)
        self.write_forcing()
        
    @staticmethod
    def is_valid(directory):
        ''' test whether directory is a valid class job'''
        
        # are basic required file types available
        MET  = glob.glob(os.path.join(directory, "*.MET"))
        INI  = glob.glob(os.path.join(directory, "*.INI"))
        CTM  = glob.glob(os.path.join(directory, "*.CTM"))
        LOOP = glob.glob(os.path.join(directory, "job_options_LOOP.txt"))
        
        if all([len(x) == 1 for x in [MET, INI, CTM, LOOP]]):
            # test if names are consistent
            names = [os.path.basename(x[0]).split('.')[0] for x in [MET, INI, CTM]]
            if all(x==names[0] for x in names):
                return True
            
        return False
        
class GTDirectory(BaseModelDir):
    """
    usage 
    
    to create a directory using paths to files, use the .from_file_paths constructor.
    The file names will be constructed appropriately in the new directory.
    
    GTDirectory.from_file_paths(rootdir="~/jobs", jobname="job1", 
                                inpts_file="", forcing_file="").build()
    """
    def __init__(self, rootdir, jobname, overwrite=True):
        super(GTDirectory, self).__init__(rootdir, jobname, overwrite=True)
    
    @classmethod
    def from_file_paths(cls, rootdir, jobname, inpts_file, forcing_file, **kwargs):
        ''' initialize using direct paths to data files '''
        GTD = cls(rootdir, jobname, **kwargs)
        GTD.inpts_file = inpts_file
        GTD.forcing_file = forcing_file
        
        # create geotop object
        GTD.gtpy = DirectoryReader(GTD.jobdir, 'geotop')
        GTD.gtpy.inpts_read(inpts_file)
        GTD.raw_forcing = GTD.gtpy.ascii_read(forcing_file, convert_geotop2=False) # read in forcing data
        
        # ensure start and end dates are appropriate
        if not ('InitDateDDMMYYYYhhmm' in GTD.gtpy.inpts.keys()) or ('InitDateDDMMYYYYhhmm' in GTD.gtpy.inpts.keys()) :
            print("No dates specified - using date range from forcing file")
            GTD.gtpy.inpts['InitDateDDMMYYYYhhmm'] = min(GTD.raw_forcing['Date'])
            GTD.gtpy.inpts['EndDateDDMMYYYYhhmm'] = max(GTD.raw_forcing['Date'])
        else:
            print("Using start and end dates from *.inpts file")
            
        return(GTD)
        
    @classmethod
    def from_abstraction(cls, rootdir, jobname):
        ''' initialize using abstracted notions of forcing and initial conditions'''
        GTD = cls(rootdir, jobname)
        return(GTD)

    def write_parameters(self):
        #out_params = os.path.join(self.jobdir, "geotop.inpts")
        #copyfile(self.inpts_file, out_params)
        self.gtpy.inpts_write()
        
    def write_forcing(self):
        ''' write forcing data from gtpy object '''
        #par = ParameterIO(self.inpts_file)
        #out_forcing = os.path.join(self.jobdir, par.MeteoFile.strip("\"") + "0001.txt")
        #copyfile(self.forcing_file, out_forcing)
        meteofile = os.path.join(self.jobdir, self.gtpy.inpts['MeteoFile'].strip("\"") + "0001.txt")
        self.gtpy.ascii_write(meteofile, self.raw_forcing)
        
    def build(self):
        ''' build directory from currently stored parameters '''
        # write parameters file
        self.write_parameters()
        
        # write forcing file
        self.write_forcing()
    
    @staticmethod
    def is_valid(directory):
        ''' test whether directory is a valid geotop job'''
        try:
            # check if directory has a geotop.inpts file
            inpts_file = os.path.join(directory, 'geotop.inpts')
        
            # check if directory has an appropriately named forcing file
            meteo = ParameterIO(inpts_file).MeteoFile.strip("\"")
            meteo_file = glob.glob(os.path.join(directory, meteo  + "*"))
            if os.path.isfile(meteo_file[0]):
                return True
        except:
            pass
        
        return False
        
class jobs_from_csv:
    ''' 
    Create a directory containing CLASS and GEOtop job directories from a csv file.
    The csv file should have the following columns:
    
    ID: integer index, 
    SITE: name of 
    MODEL: one of  {geotop, CLASS}
    FORCING: path to getotop forcing file or CLASS MET file
    INPTS: path to geotop.inpts file
    INI: path to CLASS INI file
    CTM: path to CLASS CTM file
    LOOP: path to CLASS job options file
    
    The directory is named according to the ID, SITE and MODEL chosen.
    The files in the directory are checked for consistency, for example, 
    geotop forcing files are cross-referenced with the geotop.inpts to
    ensure that the 
    '''
    
    def __init__(self, csv, rootdir, overwrite=True):
        ''' 
        @args
        csv: path to csv file specifying model directories
        rootdir: path to directory in which model directories will be created
        '''
        self.rootdir   = rootdir
        self.overwrite = overwrite 
        
        # make root directory
        if os.path.isdir(rootdir) and overwrite:
            shutil.rmtree(rootdir)
        os.mkdir(rootdir)
    
        # read csv 
        self.df = pd.read_csv(csv)
    
    def build(self):
        ''' build directories based on parameters in csv file'''
        for index, row in self.df.iterrows():
            model = row['MODEL'].upper()
            
            if model == 'GEOTOP':
                self.build_GT(row)
                
            elif model == 'CLASS':
                self.build_CLASS(row)
                
            else:
                raise NotImplementedError("Model type not implemented: {}".format(model))
            
    def build_GT(self, line):
        '''build geotop directory from a row in the csv file'''
        inpts   = line['INPTS']
        forcing = line['FORCING']
        site    = line['SITE']
        ID      = line['ID']
        jobname = "{}-{}-{}".format(ID, site, "GT")  
        
        D = GTDirectory.from_file_paths(rootdir=self.rootdir, 
                                        jobname=jobname,  
                                        inpts_file=inpts,
                                        forcing_file=forcing,
                                        overwrite=self.overwrite)
        D.build()
        
    def build_CLASS(self, line):
        ''' build CLASS directory from a row in the csv file '''
        ini   = line['INI']
        ctm   = line['CTM']
        loop  = line['LOOP']
        met   = line['FORCING']
        site    = line['SITE']
        ID      = line['ID']
        jobname = "{}-{}-{}".format(ID, site, "CLS")  
        
        D = CLASSDirectory.from_file_paths(rootdir=self.rootdir, 
                                           jobname=jobname,  
                                           MET=met, 
                                           INI=ini, 
                                           CTM=ctm, 
                                           LOOP=loop, 
                                           overwrite=self.overwrite)
        D.build()



#jobs_from_csv("/home/nbrown/storage/Projects/GLOBSIM/temp/todo.csv", "/home/nbrown/storage/Projects/GLOBSIM/jobs3" ).build()

# GTDirectory.from_file_paths("/home/nbrown/storage/Projects/GLOBSIM/temp/jobs2", "gtjob3",
# "/home/nbrown/storage/Projects/GLOBSIM/temp/data/geotop.inpts",
# "/home/nbrown/storage/Projects/GLOBSIM/temp/data/meteo1.txt").build()

# print(GTDirectory.is_valid("/home/nbrown/storage/Projects/GLOBSIM/temp/jobs2/gtjob3"))

# GTDirectory.from_file_paths("E:/NB/jobs", "cats2",
# r"E:\Users\Nick\Desktop\temp\geotop.inpts",
# r"E:\Users\Nick\Desktop\temp\Forcing_0001.txt").build()
#         
# CLASSDirectory.from_file_paths("E:/NB/jobs", "CLASS1_site",
# MET=r"E:\Users\Nick\Desktop\temp\metfiles\station_1.txt",
# CTM=r"E:\Users\Nick\Desktop\temp\class\2020_ST01_test1.CTM",
# INI=r"E:\Users\Nick\Desktop\temp\class\2020_ST01_test1.INI",
# LOOP=r"E:\Users\Nick\Desktop\temp\class\job_options_LOOP.txt").build()
#                            

        

class CollectJobs():
    def __init__(self):
        pass
        
if __name__ == "__main__":
    import sys
    jobs_from_csv(sys.argv[1], sys.argv[2]).build()
    
