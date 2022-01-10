FROM ubuntu:18.04
RUN apt-get update

# Basics
RUN apt-get install -y dialog
RUN apt-get install -y apt-utils

# Linux packages
RUN apt-get install -y gfortran
RUN apt-get install -y libgrib-api-dev
RUN apt-get install -y libopenmpi-dev
RUN apt-get install -y netcdf-bin
RUN apt-get install -y libnetcdf-dev
RUN apt-get install -y libnetcdff-dev
RUN apt-get install -y nano
RUN apt-get install -y git
RUN apt-get install -y git-lfs

# Python
RUN apt-get install -y python3.6
RUN apt-get install -y python3-pip
RUN python3.6 -m pip install -U pip
RUN cd /usr/local/bin/ && ln -s /usr/bin/python3.6 python


# Python packages
RUN python3.6 -m pip install 'numpy>=1.18.4'
RUN python3.6 -m pip install 'nco'
RUN python3.6 -m pip install 'netCDF4'
RUN python3.6 -m pip install 'scipy'
RUN python3.6 -m pip install 'pandas'
RUN python3.6 -m pip install 'pydap'
RUN python3.6 -m pip install 'python-dateutil'
RUN python3.6 -m pip install 'pyproj'
RUN python3.6 -m pip install 'cdsapi'
RUN python3.6 -m pip install 'requests'
RUN python3.6 -m pip install 'lxml'
# RUN python3.6 -m pip install 'cftime>=1.0.3.4'
RUN python3.6 -m pip install 'pygrib>=2.0.1'
RUN python3.6 -m pip install 'tomlkit'

# RUN python3.6 -m pip install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz

# Install GlobSim
RUN git lfs clone https://github.com/geocryology/globsim --depth 1 /opt/globsim

RUN cd /opt/globsim && python3.6 setup.py install

# Copy pre-built ESMF package with python bidings
RUN tar xvfz /opt/globsim/lib/esmf-python.tar.gz -C /opt

# Add to python path so ESMF can be found
#ENV PYTHONPATH=/opt/python/lib/python3.6/site-packages
#RUN echo /opt > /usr/lib/python3.6/dist-packages/esmf.pth
#RUN ln -s /opt/esmf-install/ /opt/esmf

# Add metadata
LABEL description="GlobSim"

# create entrypoint
RUN echo "#!/bin/bash" > /opt/globsim.sh && echo 'globsim "$@"' >> /opt/globsim.sh && chmod +x /opt/globsim.sh
ENTRYPOINT [ "/opt/globsim.sh" ]

