![Tests](https://github.com/geocryology/globsim/workflows/Tests/badge.svg)
# Globsim + DReaMIT

last revision November 2025

Welcome to the **DReaMIT_demo** branch of Globsim, a tool for automating the downloading, interpolation, and scaling of different reanalyses to produce meteorological time series for user-defined point locations.

The new **DReaMIT_demo** branch is the result of the integration of the **DReaMIT** (**D**ynamical **Rea**nalysis **M**odel for **I**nversions of **T**emperature) model into GlobSim. This allows the user to directly produce surface-based inversion (SBI) metrics and air temperature for any point location. Full desciption available in a coming GMD publication, authored by Victor Pozsgay, Nick Noad, Philip Bonnaventure, and Stephan Gruber.

We provide a demonstration folder: **globsim/DReaMIT_demo/**

The main files for the  integration are:
1. core file for the model [dreamit.py](/globsim/dreamit.py)
2. integration into the GlobSim scaling process for [ERA5](/globsim/scale/ERA5scale.py) and [JRA-3Q](/globsim/scale/JRAscale.py).

## Installation

Start by cloning this branch
```bash
git clone --branch DReaMIT_demo --single-branch https://github.com/geocryology/globsim.git
```
Use a conda environment to handle most dependencies to avoid building ESMF yourself. Navigate to the **globsim/DReaMIT_demo/** folder, install the **globsim_DReaMIT** conda environment from the YAML file, and activate it, before navigating back to the root of the **globsim/** folder and installing it:
```bash
cd globsim/DReaMIT_demo/
conda env create -f environment.yml
conda activate globsim_DReaMIT
cd ..
python -m pip install -e .
cd DReaMIT_demo/
```

For the remaining of the demonstration use, users will need to work from the folder **globsim/DReaMIT_demo/**.

Details for usage (and outdated installation instructions) for Globsim can be found on our [ReadTheDocs page](https://globsim.readthedocs.io/en/latest/?).

## Credentials for ERA5 and JRA-3Q

_All paths in this section are given relative to the root of /globsim/DReaMIT_demo/._

See the [Credentials](/docs/source/Credentials.rst) page for details about credentials.

Get personal tokens (credentials) to be able to download ERA5 and JRA-3Q data.
- ERA5:
    1. Create an ECMWF account, complete the form, activate your profile, and go to https://cds-beta.climate.copernicus.eu/how-to-api.
    2. Fill the credential file **./user_input/.cdsapirc** with your info.
    3. Install the CDS API client (pip install 'cdsapi>=0.7.0').
    4. Go to: https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download and scroll down a lot and accept the Licence to use Copernicus Products under the 'Terms of use'  
- JRA-3Q:
    1. Create a GDEX account at https://gdex.ucar.edu/
    2. Go to your Profile->API Token, and copy it into **./user_input/rdams_token.txt**

## Run the Notebook
_All paths in this section are given relative to the root of /globsim/DReaMIT_demo/._

_Depending on your Jupyter notebook interpreter, you might be prompted to install ipykernel too._

_The downloaded DEM data is not included in GitHub due to large file size limits. However, the user should feel free to use the notebook to download the data locally._

_The folder provides the reanalysis data, downloaded for a short time window. If the user wishes to download it for a longer time period, or different locations, they should feel free to modify the GlobSim configuration files in the ./user_input/ folder._


The main notebook, **DReaMIT_demo.ipynb** is an interactive Jupyter notebook. There, you can
1. Download ArcticDEM
2. Compute hypsometry of point locations from DEM data
3. Download, scale, and interpolate reanalysis data at point locations, including **DReaMIT** model metrics, such as air temperature, top of inversion, lapse rate, etc.
4. Produce plots of the hourly **DReaMIT** time-series.

### Initial structure of the notebook
```bash
DReaMIT_demo/
├── DReaMIT_demo.ipynb               # main Jupyter notebook for the demonstration
├── environment.yml                  # file from which to build the conda environment
├── dem_to_hypso/                    # where all the ArcticDEM data will be stored
│   ├── dem_download.py              # python script to download DEM data
│   └── __pycache__
├── reanalysis/                      # where all the reanalysis data will be stored
│   ├── reanalysis.py                # python script to download reanalysis data
│   └── __pycache__
└── user_input/                      # where the user should modify location, period, etc.
    ├── config_globsim_era5.toml     # GlobSim config file for ERA5
    ├── config_globsim_jra3qg.toml   # GlobSim config file for JRA-3Q
    ├── config_globsim_pre_hypso.csv # csv list of stations
    ├── .cdsapirc                    # ECMWF credentials for ERA5
    └── rdams_token.txt              # GDEX credentials for JRA-3Q

```

### Final structure of the notebook
```bash
DReaMIT_demo/
├── DReaMIT_demo.ipynb
├── environment.yml
├── dem_to_hypso/
│   ├── dem_download.py
│   ├── config_globsim_with_hypso.csv
│   ├── __pycache__
│   ├── DMP_WS01                              # All the DEM data for station 1: DMP_WS01
|   |   ├── df_grid_DMP_WS01.pkl              # pickled DEM data as a panda dataframe
|   |   ├── arcticdem_clipped.tif
|   |   ├── arcticdem_merged.tif
|   |   └── arcticdem_data/
|   |       └── lots of files and folders
│   └── DMP_WS02                              # All the DEM data for station 2: DMP_WS02
|       ├── df_grid_DMP_WS02.pkl              # pickled DEM data as a panda dataframe
|       ├── arcticdem_clipped.tif
|       ├── arcticdem_merged.tif
|       └── arcticdem_data/
|           └── lots of files and folders
├── reanalysis/
│   ├── reanalysis.py
│   ├── __pycache__
│   ├── era5/
│   │   ├── par/
│   │   │   └── config_globsim_with_hypso.csv # automatically-generated file with 'hypsometry' column
│   │   ├── some .nc files
│   │   ├── era5_to.grib
│   │   ├── grib_files/
│   │   │   └── some .grib files
│   │   ├── interpolated/
│   │   │   └── some .nc files
│   │   └── scaled/
│   │       └── scaled_era5.nc                # final scaled ERA5 netCDF file
│   └── jra3qg/
│       ├── par/
│       │   └── config_globsim_with_hypso.csv # automatically-generated file with 'hypsometry' column
│       ├── some .nc files
│       ├── interpolated/
│       │   └── some .nc files
│       └── scaled/
│           └── scaled_jra3qg.nc              # final scaled JRA-3Q netCDF file
├── user_input/
│   ├── config_globsim_era5.toml
│   ├── config_globsim_jra3qg.toml
│   ├── config_globsim_pre_hypso.csv
│   ├── .cdsapirc
│   └── rdams_token.txt
└── plots/                                    # all the produced plots are found here
    └── some .pdf files
```

## Adapting the code

### Change locations
1. Make sure to update the list of stations in **./user_input/config_globsim_pre_hypso.csv**
2. Update the area bounding box of the ERA5 and JRA-3Q TOML configuration files (MAKE SURE TO HAVE AT LEAST 1.5 DECIMAL DEGREES IN BOTH LATITUDE AND LATITUDE)

See the [Siteslist](/docs/source/Siteslist.rst) page for details about sites list.

### Change period
Change 'beg' and 'end' fields under [scale], [interpolate], and [scale] sections of both TOML configuration files.

### Change kernels
See the [Operation](/docs/source/Operation.rst) page for details about available scaling kernels and associated netCDF variables.

### Change model parameters  (e.g., for new regional calibration)
The file that contains the calibrated model parameters (alpha, beta) is found at **globsim/globsim/data/DReaMIT_params.csv**. If you would like to update the parameters with your own, modify this csv file.

## Observational data
Observational data used in the study to calibrate and test the model is available in the **./observations/** folder. There, you will find hourly air temperature data for all WS01 and WS02 sites, together with all Dawson sites. Note that data from ECCC/NAVCAN sites is accessible for download with the R-package weathercan (Lazerte, 2018).

## Disclaimer
GlobSim is made available for use under the GNU GPL-3 license. We do not guarantee that this software will work with your particular hardware or software. We also make no claim of offering technical support or continued development. However, any issues or bugs should be reported using the github issue tracking tool.

# Have fun
Please customise the code and use it for your project. Then let us know how things work. We hope this is useful for you.