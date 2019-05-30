Operation
=========

Set-up
------

Parameter files
^^^^^^^^^^^^^^^
Three parameter files are used to control GLOBSIM. There is one for each step in the procedure (download, interpolate, scale). The parameter files should all be in the /par subdirectory of the project directory. 


Downloading
--------------


=========================         =============
   **Keyword**                    **Description** 
-------------------------         ------------- 

**project_directory**             This is the full path to the project directory which stores the downloaded files and the control files. It should include a subdirectory called /par which contains parameter files (control files) as well as a csv describing the sites to which data are scaled.
**credentials_directory**         The location of the credential files (e.g. `.merrarc` and `.jrarc`).  Does not apply to credential file .ecmwfapi which defaults to your home directory. It is recommended to set this parameter to your home directory
**chunk_size**                    How many days to include in each download file.  Larger chunk size values mean that a smaller number of files will be downloaded, each with a larger size
**bbN**                           Coordinates for northern boundary of bounding box describing the area for which data will be downloaded.  Coordinates must be in decimal degrees.
**bbS**                           Coordinates for southern boundary of bounding box describing the area for which data will be downloaded. Coordinates must be in decimal degrees.
**bbW**                           Coordinates for western boundary of bounding box describing the area for which data will be downloaded.  Coordinates must be in decimal degrees with negative values for locations west of 0.
**bbE**                           Coordinates for eastern boundary of bounding box describing the area for which data will be downloaded. Coordinates must be in decimal degrees with negative values for locations west of 0.    
**ele_min**                       Minimum elevation that will be downloaded. Recommended to leave at 0.
**ele_max**                       Maximum elevation that will be downloaded. Should be at least 2500.
**beg**                           First date for which data is downloaded YYYY/MM/DD
**end**                           Last date for which data is downloaded YYYY/MM/DD
**variables**                     Which variables should be downloaded from the server. The variables names come from the `CF Standard Names table <http://cfconventions.org/Data/cf-standard-names/59/build/cf-standard-name-table.html>`_.  It is recommended that the variables parameter be left to include all relevant variables: air_temperature, relative_humidity, precipitation_amount, downwelling_longwave_flux_in_air, downwelling_longwave_flux_in_air_assuming_clear_sky, downwelling_shortwave_flux_in_air, downwelling_shortwave_flux_in_air_assuming_clear_sky,  wind_from_direction, wind_speed
=========================         =============

.. note:: To check download progress, you can use your credentials to log onto the website for `JRA <https://rda.ucar.edu/#ckrqst>`_ and `ERA5 (CDS API) <https://cds.climate.copernicus.eu/cdsapp#!/yourrequests>`_

Interpolating
-------------

=========================         ===============
   **Keyword**                    **Description** 
-------------------------         ---------------
**project_directory**             This is the full path to the project directory which stores the downloaded files and the control files. It should include a subdirectory called /par which contains parameter files (control files) as well as a csv describing the sites to which data are scaled. 
**station_list**                  The filename (without path) of csv containing site information such as *sitelist.csv* (note that this must match the scaling parameter file)
**chunk_size**                    How many time-steps to interpolate at once. This helps memory management. Keep small for large area files and/or computers with little memory. Make larger to get performance improvements on computers with lots of memory.
**beg**                           Beginning of date range for which data will be interpolated in YYYY/MM/DD format.  Note that this date range must include dates that are represented in the downloaded data.
**end**                           End of date range for which data will be interpolated in YYYY/MM/DD format.  Note that this date range must include dates that are represented in the downloaded data.
**variables**                     Which variables should be downloaded from the server. The variables names come from the `CF Standard Names table <http://cfconventions.org/Data/cf-standard-names/59/build/cf-standard-name-table.html>`_.  It is recommended that the variables parameter be left to include all relevant variables
=========================         ===============


Rescaling
---------

=========================         ===============
   **Keyword**                    **Description** 
-------------------------         ---------------
**project_directory**             This is the full path to the project directory which stores the downloaded files and the control files. It should include a subdirectory called /par which contains parameter files (control files) as well as a csv describing the sites to which data are scaled.
**station_list**                  The filename (without path) of csv containing site information such as *sitelist.csv* (note that this must match the interpolation parameter file)
**output_file**                   Path to output netCDF to be created. 
**overwrite**                     Either *True* or *False*. Whether or not to overwrite the `output_file` if it exists.
**time_step**                     The desired output time step in hours
**kernels**                       Which processing kernels should be used. Missing or misspelled kernels will be ignored by globsim.
=========================         ===============

