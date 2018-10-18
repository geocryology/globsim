from multiprocessing import Pool, cpu_count
from RunningJobs import CLASSDirectory, GTDirectory
from CLASS_interaction import classp
from generic import ParameterIO
from geotop_interaction import Gtpy

import glob
import os


class RunJobs():
    def __init__(self, rootdir, params):
        """
        Launches a directory of model jobs, can handle CLASS and GEOTOP jobs
        @args
            rootdir: where are the jobs
            params: text file that provides, at the minimum, paths to executables
                    (GT_PATH and CLS_PATH) and how many cpus to use (CORES)
        """
        self.rootdir = rootdir
        self.jobs = []
        self.par = ParameterIO(params)

        if 'CORES' in list(self.par.keys()):
            self.cores = self.par['CORES']
        else:
            self.cores = cpu_count() // 4 
    
    def create_jobs_list(self):
        '''create jobs list from directory'''
        # reset jobs list
        self.jobs = []

        # list files in directory
        dir_list = glob.glob(os.path.join(self.rootdir, "*"))

        for D in dir_list:
            # check what kind of job it is
            if GTDirectory.is_valid(D):
                self.add_geotop_job(D)
            
            elif CLASSDirectory.is_valid(D):
                self.add_class_job(D)
            
            else:
                pass
    
    def add_geotop_job(self, directory):
        '''add geotop_interaction object to jobs list'''
        # instantiate geotop wrapper and add to list
        GT = Gtpy(directory=directory, executable=self.par['GT_PATH'])
        self.jobs.append(GT)
    
    def add_class_job(self, directory):
        '''add class_interaction object to jobs list'''
        # get paths required to instantiate CLASS wrapper
        jobopt = glob.glob(os.path.join(directory, "job_options_LOOP.txt"))
        ctm = glob.glob(os.path.join(directory, "*.CTM"))
        site_name = os.path.basename(ctm[0]).split('.')[0] 
        
        # instantiate CLASS wrapper and add to list
        CL = classp(jobopt=jobopt, site_name=site_name, executable=self.par['CLS_PATH'])
        self.jobs.append(CL)

    @staticmethod
    def run_job(model_wrapper):
        ''' launch job.  model_wrapper must have a run method'''
        model_wrapper.run()
        
    def run_jobs(self):
        poolparty = Pool(self.ncores)

        # run all jobs from joblist in parallel
        poolparty.wrap(self.jobs, self.run_job)
