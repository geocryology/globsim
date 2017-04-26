#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# (C) Copyright Xiaojing Quan & Stephan Gruber
#
# For variable codes and units of MERRA AND MERRA-2, see: 
#  MERRA: https://gmao.gsfc.nasa.gov/products/documents/MERRA_File_Specification.pdf
#  MERRA-2: https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2C0NXLND.5.12.4/doc/MERRA2.README.pdf
#
#==============================================================================
# A script for downloading MERRA data 
# Saved as netCDF 4
#
#==============================================================================
# REFERENCE:
#
# The code is referenced from the sample code provided at NASA Earthdata website:
# HOW TO Access Data With PyDAP:
# https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+PyDAP
#
#====================HOW TO RUN THIS ==========================================
#
# (1) Register a New User in Earthdata Login: 
#  https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login
#
# (2) Authorize NASA GESDISC DATA ARCHIVE in Earthdata Login
# https://disc.gsfc.nasa.gov/registration/authorizing-gesdisc-data-access-in-earthdata_login
#
#(3) Configuring USERNAME and PASSWORD for authentication using a .netrc file by following commands:
#   > cd ~
#   > touch .netrc
#   > echo "machine urs.earthdata.nasa.gov login <USERNAME> password <PASSWORD>" > .netrc
#   > chmod 0600 .netrc
# 
#(4) Adapt the script below with authrized USERNAME and PASSWROD, then run it
#
#==============================================================================

from pydap.client import open_url
from pydap.cas.urs import setup_session

username = "<quanxj17>"
password = "<Qxj17carleton>"

session=setup_session(username, password)

dataset = open_url('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4/2016/06/MERRA2_400.tavg1_2d_slv_Nx.20160601.nc4')
 
 






