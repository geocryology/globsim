.. _example:

Examples
=========

Introduction
-------------
This is a brief tutorial that assumes you are :ref:`docker`.

Requirements
^^^^^^^^^^^^^^
To complete this example, you will need to have set up GlobSim in a Docker container using the instructions provided. You will also need :ref:`credentials` files for ERA5 and JRA-55

Tutorial
--------
Navigate to the GlobSim folder::

    cd /opt/GlobSim/globsim
   
Run the download script::

    python3 globsim.globsim_download.py -d ERA5 -f ./examples/Example1/par/examples.globsim_download
    
This downloads meteorological variables from ERA5 for a small region in the NWT. It should take a few minutes to complete (although depending on how busy the servers are it may take longer), and will download files into the Example1/era5 directory. Once it has finished, you can interpolate the gridded data to points::

    python3 globsim.globsim_interpolate.py -d ERA5 -f ./examples/Example1/par/examples.globsim_interpolate
    
This creates three new files per reanalysis in the /Example1/stations files. Finally, run the final scaling operation::

    python3 globsim.globsim_scale.py -d ERA5 -f ./examples/Example1/par/examples.globsim_scale

This creates a single file in the /Example1/scale folder with scaled meteorological variables at each of the stations in the /Example1/par/siteslist.csv.

You can now try with another reanalysis::

    python3 globsim.globsim_download.py -d JRA -f ./examples/Example1/par/examples.globsim_download
    python3 globsim.globsim_interpolate.py -d JRA -f ./examples/Example1/par/examples.globsim_interpolate
    python3 globsim.globsim_scale.py -d JRA -f ./examples/Example1/par/examples.globsim_scale

