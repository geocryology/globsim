"""Module downloads ArcticDEM and ofrmatis it into a panda dataframe,
   before computing point location hypsometry."""

from io import BytesIO
import tarfile
import zipfile
import os
import shutil
import pickle
import requests
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

from geopy.distance import distance
from geopy.point import Point

import rasterio
from rasterio.merge import merge
from rasterio.mask import mask

from pyproj import Transformer


def bboxes(csv_in, d=60):
    """Creates bounding boxes from point locations"""
    df = pd.read_csv(csv_in)

    dict_bboxes = {k: [] for k in df['station_name']}

    for i in df.index:
        sta, lat, lon = df.loc[i,'station_name'], df.loc[i,'latitude_dd'], df.loc[i,'longitude_dd']

        origin = Point(lat, lon)

        north = distance(kilometers=d).destination(origin, 0)   # bearing 0°
        south = distance(kilometers=d).destination(origin, 180) # bearing 180°
        east  = distance(kilometers=d).destination(origin, 90)  # bearing 90°
        west  = distance(kilometers=d).destination(origin, 270) # bearing 270°

        bbox_text = {
            "min_lon": west.longitude,
            "min_lat": south.latitude,
            "max_lon": east.longitude,
            "max_lat": north.latitude
        }

        bbox = tuple(list(bbox_text.values()))

        dict_bboxes[sta] = bbox

    pd_bbox = pd.DataFrame(dict_bboxes, index=list(bbox_text)).T
    pd_bbox.index.name = 'station_name'

    pd_bbox.to_csv('./dem_to_hypso/bboxes.csv')

def download_dem_files(bbox, indir, outdir):
    """Downloads ArcticDEM files"""
    # DEM resolution (for filtering URLs)
    dem_res = "2m"  # use "2m" or "10m" or "32m", whatever is available

    # --------------------------
    # STEP 1: Download and load shapefile index
    # --------------------------
    index_url = "https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Mosaic_Index_latest_shp.zip"
    print("Downloading ArcticDEM index...")
    r = requests.get(index_url)
    r.raise_for_status()

    with zipfile.ZipFile(BytesIO(r.content)) as z:
        shapefile_name = [f for f in z.namelist() if f.endswith(".shp")][0]
        z.extractall(outdir)

    index_path = os.path.join(outdir, shapefile_name)
    gdf = gpd.read_file(index_path)
    print("Index loaded, total tiles:", len(gdf))

    # --------------------------
    # STEP 2: Reproject bounding box into EPSG:3413
    # --------------------------
    bbox_geom = box(*bbox)
    bbox_gdf = gpd.GeoDataFrame(index=[0], crs="EPSG:4326", geometry=[bbox_geom])
    bbox_proj = bbox_gdf.to_crs(gdf.crs)  # reproject to match DEM index

    # Find intersecting tiles
    subset = gdf[gdf.intersects(bbox_proj.iloc[0].geometry)]
    print(f"Tiles intersecting bbox: {len(subset)}")

    if subset.empty:
        raise ValueError("No tiles found — check if your bounding box overlaps ArcticDEM coverage.")

    # --------------------------
    # STEP 3: Download DEM tiles
    # --------------------------
    dem_files = []

    # Filter tiles by requested resolution first
    urls = [url for url in subset["fileurl"] if dem_res in url]

    # If none match, pick the first available resolution automatically
    if not urls:
        print(f"No tiles found for dem_res='{dem_res}', using available resolution instead.")
        urls = subset["fileurl"].tolist()

    for url in urls:
        fname = os.path.join(outdir, os.path.basename(url))
        if not os.path.exists(fname):
            print(f"Downloading {fname}")
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            with open(fname, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        dem_files.append(fname)

    if not dem_files:
        raise ValueError("No DEM files found. Check bounding box intersects ArcticDEM coverage.")

    # --------------------------
    # STEP 4: Extract DEM TIFFs
    # --------------------------
    extracted_files = []

    for archive in dem_files:
        if archive.endswith(".tar.gz"):
            with tarfile.open(archive, "r:gz") as tar:
                # Extract to same directory or a subfolder
                tar.extractall(path=indir)
                # Collect the paths of all .tif files
                extracted_files += [os.path.join(indir, f) for f in tar.getnames() if f.endswith(".tif")]

    print(f"Extracted {len(extracted_files)} DEM TIFFs")

def dem_clip_merge(bbox, indir, out_mosaic, out_clipped):
    """Merges and clips ArcticDEM files to bounding box"""
    # -----------------------
    # Step 1: Merge all DEM tiles
    # -----------------------
    tif_files = [os.path.join(indir, f) for f in os.listdir(indir) if f.endswith("_dem.tif")]
    srcs = [rasterio.open(f) for f in tif_files]

    mosaic, trans = merge(srcs)
    meta = srcs[0].meta.copy()
    meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": trans,
        "crs": srcs[0].crs  # should be EPSG:3413 already
    })
    for s in srcs: s.close()

    with rasterio.open(out_mosaic, "w", **meta) as dst:
        dst.write(mosaic)

    print(f"Merged DEM saved: {out_mosaic}")

    # -----------------------
    # Step 2: Clip to bbox (in 3413, meters)
    # -----------------------

    # Make a polygon in EPSG:4326
    bbox_gdf = gpd.GeoDataFrame(index=[0], crs="EPSG:4326", geometry=[box(*bbox)])
    # Reproject to EPSG:3413
    bbox_3413 = bbox_gdf.to_crs(meta["crs"])

    with rasterio.open(out_mosaic) as src:
        out_image, out_transform = mask(src, bbox_3413.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "crs": src.crs
        })

        with rasterio.open(out_clipped, "w", **out_meta) as dst:
            dst.write(out_image)

    print(f"Clipped DEM saved: {out_clipped}")

def save_to_panda(out_clipped, pickle_path, station):
    """Saves clipped DEM data to panda dataframe and pickles it"""
    # Open DEM
    with rasterio.open(out_clipped) as src:
        data = src.read(1)
        transform = src.transform
        height, width = src.height, src.width

        # Generate 1D coordinates for pixel centers
        cols = np.arange(width)
        rows = np.arange(height)
        x_coords_1d, _ = rasterio.transform.xy(transform, np.zeros(width), cols, offset='center')
        _, y_coords_1d = rasterio.transform.xy(transform, rows, np.zeros(height), offset='center')

    # Convert to numpy arrays
    x_coords_1d = np.array(x_coords_1d).astype(int)
    y_coords_1d = np.array(y_coords_1d).astype(int)

    # Replace nodata
    data[data == -9999] = np.nan

    # Build DataFrame
    df_grid = pd.DataFrame(data, index=y_coords_1d, columns=x_coords_1d)
    df_grid = df_grid.T.iloc[::-1]

    # Save df_grid
    with open(pickle_path, "wb") as f:
        pickle.dump(df_grid, f)

    print(f"df_grid_{station} saved to {pickle_path}")

def dem_full_download():
    """Full workflow for ArcticDEM data downloading
       and saving to panda dataframe (and pickle)"""
    pd_bboxes = pd.read_csv('./dem_to_hypso/bboxes.csv').set_index('station_name')
    list_stations = list(pd_bboxes.index)
    dict_bbox = {sta: tuple(pd_bboxes.loc[sta]) for sta in list_stations}

    for station in list_stations:
        bbox = dict_bbox[station]

        dir_sta = f'./dem_to_hypso/{station}'
        os.makedirs(dir_sta, exist_ok=True)
        outdir = f'{dir_sta}/arcticdem_data'
        os.makedirs(outdir, exist_ok=True)

        # Input directory with extracted DEM tiles
        indir = f'{outdir}/extracted'
        out_mosaic = f'{dir_sta}/arcticdem_merged.tif'
        out_clipped = f'{dir_sta}/arcticdem_clipped.tif'

        # Path to save pickle
        pickle_path = f'{dir_sta}/df_grid_{station}.pkl'

        download_dem_files(bbox, indir, outdir)
        dem_clip_merge(bbox, indir, out_mosaic, out_clipped)
        save_to_panda(out_clipped, pickle_path, station)

def load_df_grid(list_stations):
    """Loads the ArcticDEM panda dataframe from pickle"""
    df_grid = {station: [] for station in list_stations}
    for station in list_stations:
        pickle_path = f'./dem_to_hypso/{station}/df_grid_{station}.pkl'
        with open(pickle_path, "rb") as f:
            df_grid[station] = pickle.load(f)

    return df_grid

def create_df_sites(list_stations, df_grid):
    """Creates a table of point locations with their elevation from DEM"""
    cols = ['station_name', 'longitude_dd', 'latitude_dd', 'elevation_m']
    df_sites = pd.read_csv('./user_input/config_globsim_pre_hypso.csv')[cols]
    df_sites = df_sites.rename(columns={'station_name': 'Site_Name',
                                        'longitude_dd': 'Longitude',
                                        'latitude_dd': 'Latitude'}).set_index('Site_Name')

    transformer = Transformer.from_crs('EPSG:4326', 'EPSG:3413', always_xy=True)

    vec = []

    for station in list_stations:
        x_array = np.array(df_grid[station].columns)
        y_array = np.array(df_grid[station].index)

        y, x = transformer.transform(df_sites.loc[station, 'Longitude'],
                                     df_sites.loc[station, 'Latitude'])
        sub_x = [xi for xi in x_array if np.abs(x-xi)<10]
        sub_y = [yi for yi in y_array if np.abs(y-yi)<10]
        df_test = df_grid[station].loc[sub_y, sub_x]
        elev_interp = (1/((sub_x[1]-sub_x[0])*(sub_y[1]-sub_y[0])) * np.array([sub_x[1]-x, x-sub_x[0]]) @ df_test.to_numpy() @ np.array([[sub_y[1]-y], [y-sub_y[0]]]))[0]
        vec.append(elev_interp)

    df_sites['Elevat_interp'] = vec

    return df_sites

def points_dem(df_sites):
    """Creates a table of point locations data into northing/easting"""
    list_stations = list(df_sites.index)
    points = {sta: {"lat": df_sites.loc[sta, 'Latitude'], "lon": df_sites.loc[sta, 'Longitude']}
            for sta in list_stations}

    transformer = Transformer.from_crs('EPSG:4326', 'EPSG:3413', always_xy=True)

    # Convert points to x/y in projection
    for _, val in points.items():
        y, x = transformer.transform(val["lon"], val["lat"])
        val["x_proj"] = x
        val["y_proj"] = y

    return points

def radius(x, y, x0, y0, r):
    """Creates a radius mask"""
    mask_r = (x - x0)**2 + (y - y0)**2 <= r**2
    return np.where(mask_r, 1, np.nan)

def format_panda_dem_radius(list_stations, df_grid, points):
    """Applies radius mask to ArcticDEM panda dataframe"""
    df_grid_radius = {station: [] for station in list_stations}

    for station in list_stations:
        x_array = np.array(df_grid[station].columns)
        y_array = np.array(df_grid[station].index)

        x_mesh, y_mesh = np.meshgrid(x_array, y_array)
        # Apply f(x,y) on full grid
        p = points[station]
        z_r = np.where(radius(x_mesh, y_mesh, p["x_proj"], p["y_proj"], 50000) == 1,
                       df_grid[station].values, np.nan)
        # Wrap back into DataFrame
        df_grid_radius[station]=pd.DataFrame(z_r, index=df_grid[station].index,
                                             columns=df_grid[station].columns)
        df_grid_radius[station]=df_grid_radius[station].dropna(how='all').dropna(axis=1, how='all')

    return df_grid_radius

def plot_dem_with_station(list_stations, df_sites, df_grid_radius, points):
    """For each station, plots the station's location on the DEM data"""
    figures_dem = {}

    for station in list_stations:
        threshold = df_sites.loc[station, 'Elevat_interp']
        step = 10
        data = df_grid_radius[station].values
        x = np.array(df_grid_radius[station].columns)
        y = np.array(df_grid_radius[station].index)

        fig, ax = plt.subplots(figsize=(10, 10))

        # Base DEM
        ax.imshow(data[::step, ::step], extent=[x[0], x[-1], y[-1], y[0]], cmap='gray')

        # Overlay contour line where elevation = threshold
        cs = plt.contour(x[::step], y[::step], data[::step, ::step],
                         levels=[threshold], colors='brown', linewidths=2, zorder=4)

        ax.clabel(cs, inline=True, fontsize=8, fmt='%1.0f m')

        # Plot point
        ax.scatter(points[station]["x_proj"], points[station]["y_proj"],
                   label=station, s=100, edgecolors='black', zorder=5)

        ax.set_xlabel("Easting [m]")
        ax.set_ylabel("Northing [m]")
        ax.legend()
        ax.set_title(f'ArcticDEM around {station} — elevation contour at {int(threshold)} m')
        fig.savefig(f"./plots/ArcticDEM_{station}.pdf", format='pdf', bbox_inches='tight')
        figures_dem[station] = fig
        plt.close()

    return figures_dem

def compute_hypso(list_stations, df_sites, df_grid_radius):
    """Computes the hypsometric position of a point location given its elevation and DEM data"""
    h_vec = []
    for station in list_stations:
        threshold = df_sites.loc[station, 'Elevat_interp']
        all_cells = df_grid_radius[station].count().sum()
        cells_above_threshold = (df_grid_radius[station] > threshold).sum().sum()
        h_vec.append(cells_above_threshold/all_cells)
    return h_vec

def add_hypso_to_config(h_vec, list_reanalysis):
    """Add the computed hypsometry to the GlobSim configuration file
       under the new 'hyposmetry' column needed for the DReaMIT model"""
    data = pd.read_csv('./user_input/config_globsim_pre_hypso.csv')
    data['hypsometry']=h_vec
    data.to_csv('./dem_to_hypso/config_globsim_with_hypso.csv',index=False)
    for ra in list_reanalysis:
        os.makedirs(f'./reanalysis/{ra}/par', exist_ok=True)
        shutil.copy('./dem_to_hypso/config_globsim_with_hypso.csv', f'./reanalysis/{ra}/par')

def hypso_compute_and_plot(list_stations, list_reanalysis):
    """Full workflow to compute hypsometry and plot DEM"""
    df_grid = load_df_grid(list_stations)
    df_sites = create_df_sites(list_stations, df_grid)
    points = points_dem(df_sites)
    df_grid_radius = format_panda_dem_radius(list_stations, df_grid, points)
    figures_dem = plot_dem_with_station(list_stations, df_sites, df_grid_radius, points)
    h_vec = compute_hypso(list_stations, df_sites, df_grid_radius)
    add_hypso_to_config(h_vec, list_reanalysis)

    return figures_dem
