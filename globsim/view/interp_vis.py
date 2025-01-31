import xarray as xr
import tomlkit as tk
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

from pathlib import Path
""" make plots of interpolated variables """

DATA_VARS = ['PRESS_pl', 'AIRT_pl', 'AIRT_sur', 'AIRT_redcapp_DeltaT',
              'PREC_sur', 'RH_sur', 'RH_pl', 'WSPD_sur', 'WDIR_sur',
                'SW_sur', 'LW_sur', 'SH_sur', 'SW_topo_diffuse', 
                'SW_topo_direct', 'LW_topo']


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
    fig.text(0.95, 0.05, dataset.encoding['source'], 
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

def main(file, variable, aggregate="D", output_directory=None, reanalysis=None, ftype=None):
    """ plot interpolated variable """
    if Path(file).suffix.lower() == ".toml":
        if ftype is None:
            print("Please specify the type of file to plot")
            sys.exit(1)
        with open(file, 'r') as f:
            toml = tk.parse(f.read())
        stationlist = str(toml['interpolate']['station_list'])
        name_suffix = Path(stationlist).stem
        output_directory = toml.get('interpolate').get('output_directory', 
                                                       toml.get('interpolate').get('project_directory', None))
        if ftype.lower()=='pl':
            file = Path(output_directory, 'interpolated', f"{reanalysis}_{ftype}_{name_suffix}_surface.nc")
        else:
            file = Path(output_directory, 'interpolated', f"{reanalysis}_{ftype}_{name_suffix}.nc")


    dat = xr.open_dataset(file)
    if (variable is None) or (variable not in dat.data_vars):
        print(f"Variable ({variable}) not in file. Choose from:{[var for var in dat.data_vars]}")
        sys.exit(1)

    df  = get_df(dat, variable, aggregate)
    fig, ax = plot_var(df)

    format_fig(fig, ax, dat, variable)

    figname = f"{os.path.basename(file)}_{variable}_{aggregate}.png"
    if output_directory is None:
        output_directory = os.getcwd()
    figpath = os.path.join(output_directory, figname)
    
    fig.savefig(figpath)
    print(figpath)
    

def main_args(args):
    main(args.file, args.variable, args.aggregate, args.output, args.reanalysis, args.ftype)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument("file", type=str, help="file to plot, or TOML file with one of ()")
    parser.add_argument("reanalysis",  type=str, choices=('era5','jra3qg','merra'), nargs='?', default=None,
                        help="if file is a TOML file, specify the reanalysis to plot.")
    parser.add_argument("ftype",  type=str, choices=('sa','pl','sf', 'pls'), nargs='?',
                        help="if file is a TOML file, specify the type of file to plot.")
    parser.add_argument("-v", "--var", type=str, dest='variable', help="variable to plot")
    parser.add_argument("-a", "--agg", choices=["1h", "6h", "D", "ME", "YE"], dest='aggregate', default="ME", help="aggregate data")
    parser.add_argument("-o", "--output", type=str, dest='output', help="output directory")
    
    args = parser.parse_args()
    main(args.file, args.variable, args.aggregate, args.ftype, args.output)
