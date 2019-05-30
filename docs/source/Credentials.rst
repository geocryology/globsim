.. _credentials:

Credentials
===========

The credential files for ERA-I, MERRA2, JRA55 are usually stored under the directory on PC or user’s home directory on virtual machine (~/home/<username>). It is necessary to configure the directory containing credential files in the GLOBSIM download parameter file (e.g. examples.globsim_download).

JRA-55
^^^^^^
"The JRA-55 credential should be named ``.jra55rc``.  Globsim accesses data from the NCAR/UCAR servers rather than the JMA servers.

NCAR / UCAR
-----------

To `access JRA <https://rda.ucar.edu/datasets/ds628.4/#!access>`_ using the NCAR/UCAR servers you need an RDA account. Fill out the form at `The RDA user registration page <https://rda.ucar.edu/index.html?hash=data_user&action=register>`_ and activate the account by following the email link provided.
You will need an eligible email to access restricted data collections such as the JRA55.  If you did not sign up with an email from an eligible institution, you will need to edit your profile to include your institutional email (usually a university account will work) and follow the verification link that is sent to your institutional account.
Next, in your profile, authorize the JRA-55 reanalysis and click 'update profile'

The old way to access the data is to get it from the Japan Meteorological Agency servers. Credentials for these servers can be obtained from `The JRA Website <http://jra.kishou.go.jp/JRA-55/index_en.html#application>`_. by signing up for an account.  Once you have signed up, you will recieve an email with your UID and password"
"Your credential should be a text file named ``.jrarc`` with two lines as follows::

<UID>
<password>

where the *UID* and *password* values are either the same ones as the RDA account (NCAR-UCAR Servers), or were sent to you in an email (JMA servers)."

ERA Interim (ERA-I)
^^^^^^^^^^^^^^^^^^^
The ERA-I credential should be named ``.erarc``. It can be obtained by following the instructions on `The EMCWF website <https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch>`_. under the subheading *Install ECMWF Key*. Once obtained, the file structure should look like::

    {
        "url"   : "https://api.ecmwf.int/v1",
        "key"   : "XXXXXXXXXXXXXXXXXXXXXX",
        "email" : "john.smith@example.com"
    }

Where the *key* and *email* fields are specified when you sign up to the service. 

ERA5
^^^^
As of February 2019, the ECMWF servers can no longer be used to download ERA5 data. Instead, the copernicus climate data store (CDS) servers are used. 


The same credential will be used as for ERAI but you will also need to accept the `copernicus user agreement <http://apps.ecmwf.int/datasets/licences/copernicus/>`_.

MERRA-2
^^^^^^^
The MERRA credential should be named ``.merrarc``. It can be obtained from the `Earthdata website <https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login>`_.  Once you have signed up for an account, log onto your profile `here <https://urs.earthdata.nasa.gov/home>`_. and click on applications > Approved Applications.  Click on *Approve more applications* and approve the following:

- NASA GESDISC DATA ARCHIVE
- GES DISC
- LP DAAC Data Pool

Your credential should be a text file named ``.merrarc`` with two lines as follows::

<username>
<password>

