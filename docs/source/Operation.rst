Operation
=========


Parameter files
---------------
Parameter files are used to control GLOBSIM. The parameter files follow the TOML standard (https://toml.io/en/). Each step in the procedure (download, interpolate, scale) can be split into its own file, or combined. The important thing is that they always have the appropriate headings (e.g. [download], [interpolate], or [scale]). The parameter files should be in the /par subdirectory of the project directory. 


Downloading
^^^^^^^^^^^


=========================         =============
   **Keyword**                    **Description** 
-------------------------         ------------- 

``project_directory``             This is the full path to the project directory which stores the downloaded files and the control files. It should include a subdirectory called /par which contains parameter files (control files) as well as a csv describing the sites to which data are scaled.
``credentials_directory``         The location of the credential files (e.g. `.merrarc` and `.jrarc`).  Does not apply to credential file .ecmwfapi which defaults to your home directory. It is recommended to set this parameter to your home directory
``chunk_size``                    How many days to include in each download file.  Larger chunk size values mean that a smaller number of files will be downloaded, each with a larger size
``bbN``                           Coordinates for northern boundary of bounding box describing the area for which data will be downloaded.  Coordinates must be in decimal degrees.
``bbS``                           Coordinates for southern boundary of bounding box describing the area for which data will be downloaded. Coordinates must be in decimal degrees.
``bbW``                           Coordinates for western boundary of bounding box describing the area for which data will be downloaded.  Coordinates must be in decimal degrees with negative values for locations west of 0.
``bbE``                           Coordinates for eastern boundary of bounding box describing the area for which data will be downloaded. Coordinates must be in decimal degrees with negative values for locations west of 0.    
``ele_min``                       Minimum elevation that will be downloaded. Recommended to leave at 0.
``ele_max``                       Maximum elevation that will be downloaded. Should be at least 2500.
``beg``                           First date for which data is downloaded YYYY/MM/DD
``end``                           Last date for which data is downloaded YYYY/MM/DD
``variables``                     Which variables should be downloaded from the server. The variables names come from the `CF Standard Names table <http://cfconventions.org/Data/cf-standard-names/59/build/cf-standard-name-table.html>`_.  It is recommended that the variables parameter be left to include all relevant variables: air_temperature, relative_humidity, precipitation_amount, downwelling_longwave_flux_in_air, downwelling_longwave_flux_in_air_assuming_clear_sky, downwelling_shortwave_flux_in_air, downwelling_shortwave_flux_in_air_assuming_clear_sky,  wind_from_direction, wind_speed
=========================         =============

.. note:: To check download progress, you can use your credentials to log onto the website for `JRA <https://rda.ucar.edu/#ckrqst>`_ and `ERA5 (CDS API) <https://cds.climate.copernicus.eu/cdsapp#!/yourrequests>`_

Interpolating
^^^^^^^^^^^^^

=========================         ===============
   **Keyword**                    **Description** 
-------------------------         ---------------
``project_directory``             This is the full path to the project directory which stores the downloaded files and the control files. It should include a subdirectory called /par which contains parameter files (control files) as well as a csv describing the sites to which data are scaled. 
``output_directory``              The directory to which interpolated files will be saved. If not provided, or if the directory does not exist, globsim will write to the ``project_directory`` by default. 
``station_list``                  The filename of a csv containing site names and coordinates. If just the filename is specified without a path, globsim will look in the ``par`` folder within the project directory. If a full filepath is used, globsim will use that file instead. Typically the same ``station_list`` file will also be used in the scaling step.
``chunk_size``                    How many time-steps to interpolate at once. This helps memory management. Keep small for large area files and/or computers with little memory. Make larger to get performance improvements on computers with lots of memory.
``beg``                           Beginning of date range for which data will be interpolated in ``YYYY/MM/DD`` format.  Note that this date range must include dates that are represented in the downloaded data.
``end``                           End of date range for which data will be interpolated in ``YYYY/MM/DD`` format.  Note that this date range must include dates that are represented in the downloaded data.
``variables``                     Which variables should be downloaded from the server. The variables names come from the `CF Standard Names table <http://cfconventions.org/Data/cf-standard-names/59/build/cf-standard-name-table.html>`_.  It is recommended that the variables parameter be left to include all relevant variables
=========================         ===============


Rescaling
^^^^^^^^^

=========================         ===============
   **Keyword**                    **Description** 
-------------------------         ---------------
``project_directory``             This is the full path to the project directory that contains the interpolated files. By default, it contains a subdirectory called /par where site list csv files are kept.
``output_directory``              The directory to which scaled files will be saved. If not provided, or if the directory does not exist, globsim will write to the ``project_directory`` by default. 
``station_list``                  The filename (without path) of csv containing site information such as *sitelist.csv* (note that this must match the interpolation parameter file)
``output_file``                   Path to output netCDF to be created. 
``overwrite``                     Either ``true`` or ``false``. Whether or not to overwrite the ``output_file`` if it exists.
``time_step``                     The desired output time step in hours
``kernels``                       Which processing kernels should be used. Missing or misspelled kernels will be ignored by globsim.
``scf``                           (optional) snow correction factor, a positive real number used to manually scale the precipitation for all sites.
``rh_approximation``              (optional) which relative humidity approximation to use. Choose from  'rh_liston' (default) or 'rh_lawrence'
=========================         ===============

Scaling kernels
^^^^^^^^^^^^^^^

=========================         ===============
   **Keyword**                    **Description** 
-------------------------         ---------------
``PRESS_Pa_pl``                   ``PRESS_pl``: Surface air pressure from pressure levels (long_name='air_pressure pressure levels only', standard_name='surface_air_pressure', units='Pa').
``AIRT_C_pl``                     ``AIRT_pl``: Air temperature derived from pressure levels, exclusively (long_name='air_temperature pressure levels only', standard_name='air_temperature', units='degrees_C').
``AIRT_C_sur``                    ``AIRT_sur``: Air temperature derived from surface data, exclusively (long_name='2_metre_temperature surface only', standard_name='air_temperature', units='degrees_C').
``AIRT_DReaMIT``                  Air temperature derived from surface data, pressure level data, and dynamically-computed inversion metrics as shown by the method DReaMIT. Produces the following metrics:

                                  - ``z_top_inversion_m`` (long_name='Top of inversion elevation', units='m')
                                  - ``T_lapse_grid_C`` (long_name='Temperature at grid level if solely determined by dynamically-computed linear lapse rate', units='degree_C')
                                  - ``T_lapse_station_C`` (long_name='Temperature at station level if solely determined by dynamically-computed linear lapse rate', units='degree_C')
                                  - ``lapse_Cperm`` (long_name='Lower atmosphere dynamically-computed linear lapse rate', units='degree_C m-1')
                                  - ``AIRT_DReaMIT_C`` (long_name='DReaMIT model air temperature corrected for surface-based inversion', units='degree_C')
                                  - ``beta_t_C`` (long_name='reanalysis bias independent of stations', units='degree_C')
``AIRT_redcapp``                  ``AIRT_redcapp_DeltaT``: (long_name='Land-surface influences on surface air temperature at elevation of coarse scale', units='degree_C').
``RH_per_sur``                    ``RH_sur``: Relative Humidity derived from surface data, exclusively (long_name='relative humidity surface only', standard_name='relative_humidity', units='percent').
``RH_per_pl``                     ``RH_pl``: Relative Humidity derived from pressure-level data, exclusively (long_name='relative humidity pressure-level only', standard_name='relative_humidity', units='percent').
``WIND_sur``                      Wind at 10 metre derived from surface data, exclusively. Produces the following metrics:

                                  - ``10 metre U wind component`` (long_name='10 metre U wind component', units='m s-1')
                                  - ``10 metre V wind component`` (long_name='10 metre V wind component', units='m s-1')
                                  - ``WSPD_sur`` (long_name='10 metre wind speed surface only', standard_name='wind_speed', units='m s-1')
                                  - ``WSPD_dir`` (long_name='10 metre wind direction surface only', standard_name='wind_from_direction', units='degree')
``SW_Wm2_sur``                    ``SW_sur``: solar radiation downwards derived from surface data, exclusively (long_name='Surface solar radiation downwards surface only', standard_name='surface_downwelling_shortwave_flux', units='W m-2').
``SW_Wm2_topo``                   Short-wave downwelling radiation corrected using a modified version of TOPOscale. Partitions into direct and diffuse . Produces the following metrics:

                                  - ``SW_topo_diffuse`` (long_name='TOPOscale-corrected diffuse solar radiation', standard_name='surface_diffuse_downwelling_shortwave_flux_in_air', units='W m-2')
                                  - ``SW_topo_direct`` (long_name='TOPOscale-corrected direct solar radiation', standard_name='surface_direct_downwelling_shortwave_flux_in_air', units='W m-2')
                                  - ``SW_topo_global`` (long_name='TOPOscale-corrected global solar radiation', standard_name='surface_downwelling_shortwave_flux_in_air', units='W m-2')
``LW_Wm2_sur``                    ``LW_sur``: Long-wave radiation downwards derived from surface data, exclusively (long_name='Surface thermal radiation downwards surface only', standard_name='surface_downwelling_longwave_flux', units='W m-2').
``LW_Wm2_topo``                   ``LW_topo``: Long-wave downwelling scaled using TOPOscale with surface- and pressure-level data (long_name='TOPOscale-corrected thermal radiation downwards', standard_name='surface_downwelling_longwave_flux', units='W m-2').
``PREC_mm_sur``                   ``PREC_sur``: Precipitation derived from surface data, exclusively. Convert unit: to mm/s (kg m-2 s-1) (long_name='Total precipitation surface only', standard_name='precipitation_flux', units='kg m-2 s-1').
``SH_kgkg_sur``                   ``SH_sur``: Specific humidity [kg/kg] derived from surface data, exclusively (long_name='specific humidity surface only', standard_name='specific_humidity', units='kg kg-1').
=========================         ===============

Example Parameter File
^^^^^^^^^^^^^^^^^^^^^^
Here is an example of a TOML parameter file with all three sections (download, interpolate, scale) combined into one section.

::

    title = "Globsim Control File"

    [download]
    # logistics
    project_directory = "/opt/globsim/examples/Example1"
    credentials_directory = "/root"

    # chunk size for splitting files and download [days]
    chunk_size = 2

    # area bounding box [decimal degrees]
    bbN = 66
    bbS = 62
    bbW = -112
    bbE = -108

    # ground elevation range within area [m]
    ele_min = 0
    ele_max = 2500

    # time slice [YYYY/MM/DD]
    beg = "2017/07/01"
    end = "2017/07/05"

    # variables to download [CF Standard Name Table]
    variables = ["air_temperature", "relative_humidity", "wind_speed", "wind_from_direction", "precipitation_amount", "downwelling_shortwave_flux_in_air", "downwelling_longwave_flux_in_air", "downwelling_shortwave_flux_in_air_assuming_clear_sky", "downwelling_longwave_flux_in_air_assuming_clear_sky"]

    [interpolate]
    # Path to the parent directory of /par - It should match the download and scale files
    project_directory = "/opt/globsim/examples/Example1"

    # Filename (without path) of csv containing site information (must match scaling control file)
    station_list = "siteslist.csv"

    # How many time steps to interpolate at once? This helps memory management.
    # Keep small for large area files and small memory computer, make larger to get 
    # speed on big machines and when working with small area files.
    # for a small area, we suggest values up to 2000, but consider the memory limit of your computer
    chunk_size = 2000

    # time slice [YYYY/MM/DD] assuming 00:00 hours
    beg = "2017/07/01"
    end = "2017/07/05"

    # variables to interpolate [CF Standard Name Table]
    variables = ["air_temperature", "relative_humidity", "wind_speed", "wind_from_direction", "precipitation_amount", "downwelling_shortwave_flux_in_air", "downwelling_longwave_flux_in_air", "downwelling_shortwave_flux_in_air_assuming_clear_sky", "downwelling_longwave_flux_in_air_assuming_clear_sky"]

    [scale]
    # Path to the parent directory of /par - It should match the download and interpolate files
    project_directory = "/opt/globsim/examples/Example1"

    # Filename (without path) of csv containing site information (must match interpolation control file)
    station_list = "siteslist.csv"

    # processing kernels to be used.  Unavailable kernels will be ignored
    kernels = ["PRESS_Pa_pl", "AIRT_C_pl", "AIRT_C_sur", "AIRT_DReaMIT", "AIRT_redcapp", "RH_per_sur", "RH_per_pl", "WIND_sur", "SW_Wm2_sur", "SW_Wm2_topo", "LW_Wm2_sur", "LW_Wm2_topo", "PREC_mm_sur", "SH_kgkg_sur"]

    # desired time step for output data [hours]
    time_step = 1

    # Should the output file be overwritten if it exists?
    overwrite = true

Station list for interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is an example of a Globsim station list file. The resulting netCDF file will use the station numbers as identifiers. ::

    station_number,station_name,longitude_dd,latitude_dd,elevation_m,sky_view,slope,aspect,hypsometry
    1,'hilltop',-129.1,64.1,2500,0.9,70,90,0
    2,'swamp_2',-129.4,61.1,1000,0.5,0,0,0.5
    3,'research_station',-126.1,64.1,1500,1.0,10,180,1
   
More information about the station list can be found on the :ref:`stationlist` page
     
Project directory
^^^^^^^^^^^^^^^^^     
The ``project directory`` is the location to which data is downloaded and where processed data is found. The project directory is subdivided by re-analysis type and by the type of derived product::

     project_a/              (project directory)
     project_a/par/          (parameter files for data download and interpolation)
     project_a/jra-55/       (JRA-55 data)
     project_a/eraint/       (ERA-Interim data)
     project_a/era5/         (ERA5 data)
     project_a/merra2/       (MERRA 2 data)
     project_a/station/      (data interpolated to stations)
     project_a/scale/        (final scaled files)