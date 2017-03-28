#
# DOWNLOAD AND MANIPULATE DATA FOR THE EKATI AREA
# 
# Change in header
#
#===============================================================================
from datetime import datetime
from os import path




#settings
dir_src  = '/Users/stgruber/src/globsim'

execfile(path.join(dir_src, 'ERA_Interim_CF.py'))

# set parameters
date      = {'beg' : datetime(1979, 1, 1), 'end' : datetime(2016,12, 31)}
area      = {'north':  65.0, 'south': 64.0, 'west': -111.0, 'east': -110.0}
elevation = {'min': 0, 'max': 2000}           
directory = '/Users/stgruber/Desktop/data'             

# get it
ts = ERAbatch(date, area, elevation, directory, 15000) 
ts.retrieve()