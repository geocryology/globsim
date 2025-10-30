"""Module downloads, interpolates, scales, and plots reanalysis data,
   including the new DReaMIT metrics."""

import subprocess
import os
from netCDF4 import Dataset, num2date
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def format_data(scaled_files, path_csv, list_reanalysis):
    """Formats scaled netCDF data into panda dataframes"""
    bboxes = pd.read_csv(path_csv).set_index('station_name')
    list_stations = list(bboxes.index)

    ds_final = {ra: Dataset(scaled_files[ra]) for ra in list_reanalysis}

    df = {ra: {station: [] for station in list_stations} for ra in list_reanalysis}

    list_vars = ['AIRT_pl', 'AIRT_sur', 'AIRT_DReaMIT_C', 'z_top_inversion_m',
                 'T_lapse_grid_C', 'T_lapse_station_C', 'lapse_Cperm']

    for ra in list_reanalysis:
        for i,station in enumerate(list_stations):
            time = num2date(ds_final[ra]['time'][:], ds_final[ra]['time'].units)
            full_array = np.array([time]+[ds_final[ra][var][:,i] for var in list_vars]).T
            df[ra][station] = pd.DataFrame(full_array,
                                           columns=['time']+list_vars)
            df[ra][station] = df[ra][station].astype({'time':'str'})
            df[ra][station] = df[ra][station].astype({'time':'datetime64[ns]'})
            df[ra][station] = df[ra][station].set_index('time')

    return df

def plotting(df, station, ra, formatted_reanalysis):
    """Plots the panda dataframe version of the scaled netCDF time series"""
    fig,ax=plt.subplots(2,2,sharex=True,squeeze=False,constrained_layout=True,figsize=(12,8))

    df_temp = df[ra][station]

    df_temp[['AIRT_pl', 'AIRT_sur', 'AIRT_DReaMIT_C']].plot(ax=ax[0,0])
    df_temp[['T_lapse_grid_C', 'T_lapse_station_C']].plot(ax=ax[0,1])
    df_temp[['z_top_inversion_m']].plot(ax=ax[1,0])
    (df_temp[['lapse_Cperm']]*1000).plot(ax=ax[1,1])

    ax[0,0].set_ylabel('Air temperature [°C]')
    ax[0,1].set_ylabel('Lapse air temperature [°C]')
    ax[1,0].set_ylabel('Top of inversion elevation [m]')
    ax[1,0].set_xlabel('')
    ax[1,1].set_ylabel('Lapse rate [°C km-1]')
    ax[1,1].set_xlabel('')

    fig.align_ylabels()
    title = f'Temperature and inversion metrics for {station}, with {formatted_reanalysis[ra]}'
    fig.suptitle(title)

    plt.close()

    fig.savefig(f'./plots/plot_{ra}_{station}.pdf', bbox_inches='tight')

    return fig

def plotting_from_reanalysis(list_reanalysis, formatted_reanalysis):
    """Plots the scaled netCDF time series directly from netCDF"""
    scaled_files = {ra: f'./reanalysis/{ra}/scaled/scaled_{ra}.nc' for ra in list_reanalysis}
    df = format_data(scaled_files, './user_input/config_globsim_pre_hypso.csv', list_reanalysis)
    figs = {ra: {station: plotting(df, station, ra, formatted_reanalysis)
                 for station in df[list_reanalysis[0]].keys()}
            for ra in list_reanalysis}
    return figs

def reanalysis_download(list_reanalysis, formatted_reanalysis):
    """Download reanaltsis data (ERA5 and JRA3Q)"""
    for ra in list_reanalysis:
        f = f'./user_input/config_globsim_{ra}.toml'
        ra_form = formatted_reanalysis[ra]
        print(f'Starting the {ra_form} DOWNLOAD')
        cmd = f"globsim download -d {ra_form} -f {f}"
        p = subprocess.Popen(cmd.split(" "))
        p.wait()
        print(f'{ra_form} DOWNLOAD concluded\n\n')

def reanalysis_interpolate(list_reanalysis, formatted_reanalysis):
    """Interpolates reanaltsis data (ERA5 and JRA3Q)"""
    for ra in list_reanalysis:
        f = f'./user_input/config_globsim_{ra}.toml'
        ra_form = formatted_reanalysis[ra]
        print(f'Starting the {ra_form} INTERPOLATION')
        cmd = f"globsim interpolate -d {ra_form} -f {f} --skip-checks"
        p = subprocess.Popen(cmd.split(" "))
        p.wait()
        print(f'{ra_form} INTERPOLATION concluded\n\n')

def reanalysis_scale(list_reanalysis, formatted_reanalysis):
    """Scales reanaltsis data (ERA5 and JRA3Q).
       This step includes the DReaMIT model."""
    for ra in list_reanalysis:
        f = f'./user_input/config_globsim_{ra}.toml'
        ra_form = formatted_reanalysis[ra]
        print(f'Starting the {ra_form} SCALING')
        cmd = f"globsim scale -d {ra_form} -f {f}"
        p = subprocess.Popen(cmd.split(" "))
        p.wait()
        os.rename(f'./reanalysis/{ra}/scaled/scaled_{ra}_1h_scf1.nc',
                  f'./reanalysis/{ra}/scaled/scaled_{ra}.nc')
        print(f'{ra_form} SCALING concluded\n\n')
