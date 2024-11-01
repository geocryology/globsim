.. _credentials:

Credentials
===========

The credential files are usually stored under the directory on PC or user's home directory on virtual machine (~/home/<username>). It is necessary to configure the directory containing credential files in the GLOBSIM download parameter file (e.g. examples.globsim_download).

JRA-55
^^^^^^

To `access JRA <https://rda.ucar.edu/datasets/ds628.4/#!access>`_ using the NCAR/UCAR servers you need an RDA account. Fill out the form at `The RDA user registration page <https://rda.ucar.edu/index.html?hash=data_user&action=register>`_ and activate the account by following the email link provided.

Create a token file called rdams_token.txt and put it in your home/credentials directory.  It should contain one line with the token obtained from your `account profile <https://rda.ucar.edu/accounts/profile/>`_. 

ERA5
^^^^
Downloaded from the CDS server. 

1. Login with existing ECMWF account (or create a new one) at https://www.ecmwf.int/
2. Complete form
3. Activate your profile
4. Navigate to https://cds-beta.climate.copernicus.eu/how-to-api
5. Create your personal .cdsapirc file under /fs/yedoma/home/<user_name>/ and use the url and key provided
6. Go to https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download
scroll down a lot and accept the Licence to use Copernicus Products under the 'Terms of use'

MERRA-2
^^^^^^^
The MERRA credential should be named ``.merrarc``. It can be obtained from the `Earthdata website <https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login>`_.  Once you have signed up for an account, log onto your profile `here <https://urs.earthdata.nasa.gov/home>`_. and click on applications > Approved Applications.  Click on *Approve more applications* and approve the following:

- NASA GESDISC DATA ARCHIVE
- GES DISC
- LP DAAC Data Pool

Your credential should be a text file named ``.merrarc`` with two lines as follows::

<username>
<password>

