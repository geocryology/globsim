import xarray as xr
import tomlkit as tk
import sys
import os

from pathlib import Path

from globsim.view.plotting import plot_var, format_fig, get_df, plot_scaled


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


    # Scaled data
    if Path(file).name.startswith("scaled"):
        if variable is None:
            station = dat['station'].to_numpy()[0]
            fig, ax = plot_scaled(dat, station, aggregate)
        else:
            pass

    # Interpolated data
    else:
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
    parser.add_argument("-v", "--var", type=str, dest='variable', help="variable to plot", default=None)
    parser.add_argument("-a", "--agg", choices=["1h", "6h", "D", "ME", "YE"], dest='aggregate', default="ME", help="aggregate data")
    parser.add_argument("-o", "--output", type=str, dest='output', help="output directory")
    
    args = parser.parse_args()
    main(args.file, args.variable, args.aggregate, args.ftype, args.output)
