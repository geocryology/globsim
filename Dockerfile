FROM ubuntu:18.04
RUN apt-get update

# ENV PYTHONPATH=/opt/python/lib/python3.6/site-packages

RUN apt-get install -y dialog
RUN apt-get install -y apt-utils
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN cd /usr/local/bin/ && ln -s /usr/bin/python3 python
RUN apt-get install -y gfortran
RUN apt-get install -y libgrib-api-dev
RUN apt-get install -y libopenmpi-dev
RUN apt-get install -y git

RUN pip3 install 'numpy>=1.18.4'
RUN pip3 install nco
RUN pip3 install netCDF4
RUN pip3 install scipy
RUN pip3 install pandas
RUN pip3 install pydap
RUN pip3 install python-dateutil
RUN pip3 install pyproj
RUN pip3 install cdsapi
RUN pip3 install requests
RUN pip3 install lxml
RUN pip3 install 'cftime==1.0.3.4'

RUN pip3 install pygrib==2.0.1
RUN pip3 install tomlkit

RUN pip3 install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz

RUN apt-get install -y netcdf-bin
RUN apt-get install -y libnetcdf-dev
RUN apt-get install -y libnetcdff-dev
RUN apt-get install -y nano
RUN apt-get install -y git-lfs

# Install GlobSim
RUN git lfs clone https://github.com/geocryology/globsim --depth 1 /opt/globsim

RUN cd /opt/globsim && python3 setup.py install

# Copy pre-built ESMF package with python bidings
RUN tar xvfz /opt/globsim/lib/esmf-python.tar.gz -C /opt

# Add to python path so ESMF can be found
RUN echo /opt >/usr/lib/python3/dist-packages/globsim.pth
RUN ln -s /opt/esmf-install/ /opt/esmf

# Add metadata
LABEL description="A container for globsim"

ENTRYPOINT [ "globsim" ]
