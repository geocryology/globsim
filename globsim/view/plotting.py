import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

import math


DATA_VARS = ['PRESS_pl', 'AIRT_pl', 'AIRT_sur', 'AIRT_redcapp_DeltaT',
             'PREC_sur', 'RH_sur', 'RH_pl', 'WSPD_sur', 'WDIR_sur',
             'SW_sur', 'LW_sur', 'SH_sur', 'SW_topo_diffuse', 
             'SW_topo_direct', 'LW_topo']


""" make plots of interpolated variables """


def plot_var(df: pd.DataFrame, ylab:str=""):
    """ plot variable 
    
    """
    variable = df.columns[0]
    fig, ax = plt.subplots()
    for station, group in df.groupby(level="station"):
        group.reset_index(level="station", drop=True).plot(
            y=variable, ax=ax, label=station, title=variable)
    ax.set_ylabel(ylab)
    
    return fig, ax


def format_fig(fig, ax, dataset:xr.Dataset, variable:str):
    array = dataset[variable]
    ax.set_ylabel(f"{array.units}")
    ax.set_title(f"{array.long_name} (GlobSim-Interpolated)")
    ax.legend()
    # get file path from dataset
    watermark(fig, dataset)


def watermark(fig, dataset):
    _ = fig.text(0.95, 0.05, dataset.encoding['source'], 
         horizontalalignment='right', 
         verticalalignment='bottom', 
         transform=plt.gca().transAxes,  # Use axes coordinates
         fontsize=10, 
         color='blue')

def get_df(ds: xr.Dataset, variable:str, aggregate:str="D"):
    arr = ds[variable]
    df = arr.to_dataframe()
    if aggregate:
        df = df.groupby('station').apply(
                    lambda x: x.droplevel('station').resample(aggregate).mean()
                    )
    return df


def optimal_layout(n):
    """
    Given an integer n, returns the optimal (rows, cols) for arranging n axes in a grid.
    
    The function finds the closest balanced layout to a square while minimizing empty spaces.
    """
    if n <= 0:
        raise ValueError("Number of axes must be a positive integer")

    # Compute the square root and round down to get the base number of rows
    rows = math.floor(math.sqrt(n))
    
    # Compute columns as the ceiling of (n / rows) to fit all plots
    cols = math.ceil(n / rows)

    return rows, cols


def plot_scaled(dataset, station, agg='ME'):
    included_vars = [v for v in dataset.variables if v in DATA_VARS]
    layout = optimal_layout(len(included_vars))
    fig, axs = plt.subplots(layout[1], layout[0], figsize=(15, 15))
    
    for i, var in enumerate(included_vars):
        ax = axs.flat[i]
        if agg:
            dataset[var].sel(station=station).resample({'time':agg}).mean().plot(ax=ax)
        else:
            dataset[var].sel(station=station).plot(ax=ax)
        #except Exception:
        ax.set_title(var)
        ax.set_ylabel(dataset[var].units)
        ax.set_xlabel('Time')
        
    plt.tight_layout()
    watermark(fig, dataset)
    return fig, axs
