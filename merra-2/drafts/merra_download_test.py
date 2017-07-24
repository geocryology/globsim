#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# (C) Copyright 2007-2008 Xiaojing Quan & Stephan Gruber
#
#==============================================================================
#CONTRIBUTIONS:
# 1. Original is generated: 01/05/2017
# 2. The First draft for downloading single dataset is tested: 05/05/2017-08/05/2017
# 3. The Second draft for downloading specific dataset is continued (updated by 08/05/2017)
#
#==============================================================================
# For variable codes and units of MERRA AND MERRA-2, see: 
#  MERRA: https://gmao.gsfc.nasa.gov/products/documents/MERRA_File_Specification.pdf
#  MERRA-2: https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/M2C0NXLND.5.12.4/doc/MERRA2.README.pdf
#
#==============================================================================
# A script for downloading MERRA data 
# Saved as netCDF 4
#==============================================================================
# REFERENCE:
#
# The code is referenced from the sample code provided at NASA Earthdata website:
# HOW TO Access Data With Python:
# https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
#
#====================HOW TO RUN THIS ==========================================
#
# (1) Register a New User in Earthdata Login: 
#  https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login
#
# (2) Authorize NASA GESDISC DATA ARCHIVE in Earthdata Login
# https://disc.gsfc.nasa.gov/registration/authorizing-gesdisc-data-access-in-earthdata_login
#
# (3) Adapt the script below with authrized USERNAME and PASSWROD
# 
# (4) Obtaining the URL address of the objected single dataset at:
#     https://disc.sci.gsfc.nasa.gov/daac-bin/FTPSubset2.pl
# 
# (5) Obtianing the mutiple dataset with spefici spacial and temporal(Continued to work on it) 
# 
# 
# IMPORTANT Notes: 
# For downloading dataset as .nc4 file, extend the related URL address with double .nc4.nc4 in the end as following examples:
# url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4" TO
# url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4.nc4" 
#
#
#==============================================================================

from cookielib import CookieJar

from urllib import urlencode

import urllib2

#import ssl

# The user credentials that will be used to authenticate access to the data

username = "quanxj17"
password = "Qxj17carleton"

# The url of the file we wish to retrieve

#url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4"
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4.nc4"
#url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2_DIURNAL/M2IUNXASM.5.12.4/2016/MERRA2_400.instU_2d_asm_Nx.201601.nc4.ascii"


# Create a password manager to deal with the 401 reponse that is returned from
# Earthdata Login

password_manager = urllib2.HTTPPasswordMgrWithDefaultRealm()
password_manager.add_password(None, "https://urs.earthdata.nasa.gov", username, password)
 

# Create a cookie jar for storing cookies. This is used to store and return
# the session cookie given to use by the data server (otherwise it will just
# keep sending us back to Earthdata Login to authenticate).  Ideally, we
# should use a file based cookie jar to preserve cookies between runs. This
# will make it much more efficient.

cookie_jar = CookieJar()

# Install all the handlers.
 
opener = urllib2.build_opener(
urllib2.HTTPBasicAuthHandler(password_manager),
    #urllib2.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
    #urllib2.HTTPSHandler(debuglevel=1),   # details of the requests/responses
     urllib2.HTTPCookieProcessor(cookie_jar))
urllib2.install_opener(opener)

# Create and submit the request. There are a wide range of exceptions that
# can be thrown here, including HTTPError and URLError. These should be
# caught and handled.
 
request = urllib2.Request(url)

response = urllib2.urlopen(request)
 
 
# Print out the result (not a good idea with binary data!)
 
dataset = response.read()
#print dataset

# Save the result 

file = open("/Users/xquan/data/output.nc4", "w")
file.write(dataset)
file.close()
