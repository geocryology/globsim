#
# DOWNLOAD AND MANIPULATE DATA FOR THE EKATI AREA
#
#===============================================================================

import ERA_Interim

# set parameters
date      = {'beg' : datetime(1979, 1, 1), 'end' : datetime(2016, 9, 1)}
area      = {'north':  65.0, 'south': 64.0, 'west': -111.0, 'east': -110.0}
elevation = {'min': 0, 'max': 2000}           
directory = '/home/sgruber/data'             

# get it
ts = ERAbatch(date, area, elevation, directory, 15) 
ts.retrieve()
ts.deleteGrib()