from multiprocessing import Pool, cpu_count
from model_directories import CLASSDirectory, GTDirectory
from CLASS_interaction import classp
from generic import ParameterIO
from geotop import geotop_interaction as GI

import subprocess
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

        if hasattr(self.par, "CORES"):
            self.cores = int(self.par.CORES)
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
                print("[WARNING] {} is not a valid directory".format(os.path.basename(D)))
    
    def add_geotop_job(self, directory):
        '''add geotop_interaction object to jobs list'''
        # instantiate geotop wrapper and add to list
        GT = GI.Gtpy(directory=directory, executable=self.par.GT_PATH)
        GT.inpts_read()
        
        #assign job id so it can be identified if it dies
        GT.JOBID = os.path.basename(directory)
        self.jobs.append(GT)
    
    def add_class_job(self, directory):
        '''add class_interaction object to jobs list'''
        # get paths required to instantiate CLASS wrapper
        jobopt = glob.glob(os.path.join(directory, "job_options_LOOP.txt"))[0]
        ctm = glob.glob(os.path.join(directory, "*.CTM"))
        site_name = ctm[0].split('.')[0] 
        
        # instantiate CLASS wrapper and add to list
        CL = classp(jobopt=jobopt, site_name=site_name, executable=self.par.CLS_PATH)
        
        #assign job id so it can be identified if it dies
        CL.JOBID = os.path.basename(directory)
        self.jobs.append(CL)
        
    def run_jobs(self):
        n_jobs = len(self.jobs)
        cores = min(n_jobs, self.cores)
        print("Found {} valid directories. Running jobs with {} cores".format(n_jobs, cores))
        
        # Allocate resources
        poolparty = Pool(cores)

        # run all jobs from joblist in parallel
        results = poolparty.map(global_run_job, self.jobs)
        jobids = [x.JOBID for x in self.jobs]
        success = list(zip(jobids, results))
        print(success)

def global_run_job(job):
    ''' 
    Runs a model with a .run() method.  If the model is successful ( run() returns
    True), then returns True.  Returns false otherwise and prints any errors to the terminal
    '''
    try:
        result = job.run()
    except subprocess.CalledProcessError as e:
        print("Job '{}' crashed with the following error: \n {}".format(job.JOBID, e))
        result = False
    return(result)
    
if __name__ == "__main__":
    import sys
    M = RunJobs(rootdir = sys.argv[1], params = sys.argv[2])
    M.create_jobs_list()
    M.run_jobs()