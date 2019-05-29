Examples
=========

Introduction
-------------
This is a detailed tutorial 

Requirements
^^^^^^^^^^^^^^
To complete this example, you will need python3, a working version of GLOBSIM and :ref:`credentials`  for MERRA-2 and ERA5 (make sure they're in your home directory!). You also need the GLOBSIM package to be discoverable by python. To do this, we recommend (adding a *\*.pth* file to your site-packages folder)[https://stackoverflow.com/questions/700375/how-to-add-a-python-import-path-using-a-pth-file] that gives the globsim folder's parent directory.

Tutorial
---------

Create directory structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can skip this step if you're using the example directory provided in the GLOBSIM package. 

Create the following directory structure::

    ./
      /par
      /merra2
      /era5
    

    
Then create the following three control files in the */par* directory::

    examples.globsim_download
    examples.globsim_interpolate
    examples.globsim_scale
    
The contents of each file should be as follows. 

.. warning:: Make sure to replace any file paths marked with << >> to match your file structure:

*examples.globsim_download*::

    # This is an example of a Globsim download parameter file
    # Use extension like this: 'project_a.globsim_download'

    # Path to the parent directory of /par - It should match the interpolate and scale files
    project_directory = <</home/user/globsim/example>>

    # where to find login credentials for servers (this should be your home directory)
    credentials_directory = <</home/user01>>

    # chunk size for splitting files and download [days]
    chunk_size = 2

    # area bounding box [decimal degrees]
    bbN = 66
    bbS = 62
    bbW = -112
    bbE = -109

    # ground elevation range within area [m]
    ele_min = 0
    ele_max = 2500

    # time slice [YYYY/MM/DD]
    beg = 2017/07/01
    end = 2017/07/10

    # variables to download [CF Standard Name Table]
    variables = precipitation_amount,downwelling_longwave_flux_in_air, downwelling_longwave_flux_in_air_assuming_clear_sky,downwelling_shortwave_flux_in_air, downwelling_shortwave_flux_in_air_assuming_clear_sky, relative_humidity,2-metre_air_temperature,10-metre_eastward_wind,10-metre_northward_wind
    air_temperature, relative_humidity, precipitation_amount,
    wind_from_direction, wind_speed

    
    
**examples.globsim_interpolate**::

    # Path to the parent directory of /par - It should match the download and scale files
    project_directory = <</home/user/globsim/example>>

    # Filename (without path) of csv containing site information (must match scaling control file)
    station_list = siteslist.csv

    # How many time steps to interpolate at once? This helps memory management.
    # Keep small for large area files and small memory computer, make larger to get 
    # speed on big machines and when working with small area files.
    # for a small area, we suggest values up to 2000, but consider the memory limit of your computer
    chunk_size = 2000

    # time slice [YYYY/MM/DD] assuming 00:00 hours
    beg = 1980/01/01
    end = 1980/01/30

    # variables to interpolate [CF Standard Name Table]
    variables = air_temperature, relative_humidity, precipitation_amount, downwelling_shortwave_flux_in_air,downwelling_longwave_flux_in_air, downwelling_shortwave_flux_in_air_assuming_clear_sky, downwelling_longwave_flux_in_air_assuming_clear_sky, wind_from_direction, wind_speed


    
**examples.globsim_scale**::

    # Path to the parent directory of /par - It should match the download and interpolate files
    project_directory = <</home/user/globsim/example>>

    # Full path (with .nc extension) to output netcdf file
    output_file = <</home/user/siteslist_globsim_out.nc>>

    # Filename (without path) of csv containing site information (must match interpolation control file)
    station_list = siteslist.csv

    # processing kernels to be used.  Unavailable kernels will be ignored
    kernels = PRESS_ERA_Pa_pl, AIRT_ERA_C_pl, AIRT_ERA_C_sur, PREC_ERA_mm_sur, RH_ERA_per_sur, WIND_ERA_sur, SW_ERA_Wm2_sur, LW_ERA_Wm2_sur, SH_ERA_kgkg_sur

    # desired time step for output data [hours]
    time_step = 1

    # Should the output file be overwritten if it exists?
    overwrite = False
    
Finally, create a text file called siteslist.csv inside the */par* directory with the following contents::
    
    station_number, station_name, longitude_dd, latitude_dd, elevation_m ,sky_view
    1,Site_1,-111.5929,64.8672,415,1
    2,Site_2,-110.431,64.7026,460,1
    3,Site_3,-109.5314,62.4563,197,1
    

This file is used to specify the location of the sites for which you want data.   


Running GLOBSIM from the command line
--------------------------------------

Navigate to the folder with the globsim code

If you have your paths set up properly, you should be able to run the help document:: 
    
    python3 globsim_download.py -h 
 
Now, to download data, run the following code, replacing the text in <<>> with the appropriate path:: 
   
    python3 globsim_download.py -f <<full path to examples.globsim_download>> -d ERA5
    python3 globsim_download.py -f <<full path to examples.globsim_download>> -d MERRA

Although it is possible to download both reanalyses simultaneously by specifying them both after the -d flag, it is recommended to run them separately. To save time, you can run them on different screens. 



