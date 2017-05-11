#
# DOWNLOAD AND MANIPULATE DATA FOR THE EKATI AREA
# 
# Change in header
#
#===============================================================================
from datetime import datetime
import era_interim as ei


# set parameters
date      = {'beg' : datetime(2015, 7, 1), 'end' : datetime(2015,7, 2)}
area      = {'north':  65.0, 'south': 64.0, 'west': -111.0, 'east': -110.0}
elevation = {'min': 0, 'max': 2000}           
directory = '/Users/stgruber/Desktop/data'             

# get it
ts = ei.ERAbatch(date, area, elevation, directory, 15) 
ts.retrieve()