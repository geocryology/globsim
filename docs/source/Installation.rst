Installation
============

This describes the steps required to get GlobSim up and running.  You can either install all the packages and dependencies yourself, or you can learn about :ref:`docker`. 

Requirements
------------

GlobSim source files can be obtained in one of two ways

Using pip (does not work yet!!)::

    pip3 install globsim 

From the github source repository::

    git clone https://github.com/geocryology/globsim
    checkout working

NetCDF
------
NetCDF files are used to store data in a standard format. The NCO libraries must be installed and built on your computer for GlobSim to work.  Instructions can be found on the `Unidata website <https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html>`_. 

ECMWF Client libraries
----------------------
ECMWF is used to access ERA files. Python libraries (supporting python 2.7 and 3) to access the API are available from the `ECMWF website <https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch>`_

Grib API and pygrib
--------------------
GlobSim uses the `GRIB API <https://confluence.ecmwf.int/display/GRIB/What+is+GRIB-API>`_. 
Python bindings to the GRIB API are also necessary `(pygrib) <https://jswhit.github.io/pygrib/docs/>`_. The code was tested using pygrib version 2.0.1. 

ESMF
----
GLOBSIM uses ESMF libraries to do efficient regridding. These libraries must be built on your machine and have additional dependencies.  ESMP versions 7.0.1 and 7.1.0r are supported. To download ESMF, consult the `ESMF Users Guide <http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/ESMF_usrdoc/>`_, particularly sections 5 and 8.

During installation, several environment variables are set::

    ESMF_NETCDF
    ESMF_NETCDF_LIBPATH
    ESMF_NETCDF_LIBS
    ESMF_NETCDF_INCLUDE
    ESMF_COMPILER
    ESMF_COMM

ESMPy
^^^^^^
ESMPy provides python bindings for the ESMF libraries.  If you have successfully installed the ESMF libraries, you can follow the instructions `here <http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmpy_doc/html/install.html#installing-esmpy>`_ to extract the python bindings.  More info is available at the `ESMPy main page <https://www.earthsystemcog.org/projects/esmpy/>`_.

Example
-------

The following setup was used to install on Ubuntu 16.04::


    # GRIB API
    apt install libgrib-api-dev

    #python libraries
    apt install python3-pip
    pip3 install numpy
    pip3 install nco
    pip3 install netCDF4
    pip3 install scipy
    pip3 install pandas
    pip3 install pydap
    pip3 install python-dateutil
    pip3 install pyproj
    pip3 install pygrib==2.0.1  # requires grib API
    pip3 install lxml

    # ECMWFAPI
    pip3 install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz 

    ## Netcdf - C
    apt install netcdf-bin
    apt install libnetcdf-dev

    ## Netcdf - Fortran
    apt install libnetcdff-dev

    ## fortran compiler
    apt install gfortran

    ## open mpi
    apt install libopenmpi-dev

    # ESMF
    wget \"http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_7_1_0r/esmf_7_1_0r_src.tar.gz\"

To install ESMF, the following script was then used (again, tested on Ubuntu 16.04)::
   
    #!/bin/bash

    set -e

    export BASE_DIR=$(pwd)

    export ESMF_DIR=${BASE_DIR}/esmf
    export ESMF_INSTALL_PREFIX=${BASE_DIR}/esmf-install
    export ESMF_NETCDF=split
    export ESMF_NETCDF_LIBPATH=/usr/lib/x86_64-linux-gnu/
    export ESMF_NETCDF_LIBS="-lnetcdff -lnetcdf"
    export ESMF_NETCDF_INCLUDE=/usr/include

    export ESMF_COMPILER=gfortran
    export ESMF_COMM=openmpi


    tar xvf ~/esmf_7_1_0r_src.tar

    cd esmf
    make -j 12
    # make check # (optional)

    make install
    # make installcheck # (optional)

    cd src/addon/ESMPy
   
    python setup.py  build --ESMFMKFILE=${ESMF_DIR}/lib/libO/Linux.gfortran.64.openmpi.default/esmf.mk install


    echo "To use this version of ESMPy, run:"
    echo "  export PYTHONPATH='$BASE_DIR/python/lib/python2.7/site-packages'"


