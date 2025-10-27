import pandas as pd
import netCDF4 as nc
# from netCDF4 import Dataset, num2date
import numpy as np
import os
from lmfit import Model
import globsim.constants as const

ref_dreamit = 'DReaMIT: A Dynamical Reanalysis Framework for Modelling Surface-Based Temperature Inversions in Northern Mountain Environments'

def add_var_z_top_inversion(group: nc.Dataset):
    """ Add DReaMIT top of inversion elevation the netCDF file"""
    vn  = 'z_top_inversion_m'  # variable name
    var = group.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Top of inversion elevation'
    var.comment   = ref_dreamit
    var.units     = 'm'
    
    return group[vn]

def add_var_T_lapse_grid(group: nc.Dataset):
    """ Add DReaMIT lapse temperature at grid level to the netCDF file"""
    vn  = 'T_lapse_grid_C'  # variable name
    var = group.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Temperature at grid level if solely determined by dynamically-computed linear lapse rate'
    var.comment   = ref_dreamit
    var.units     = 'degree_C'
    
    return group[vn]

def add_var_T_lapse_station(group: nc.Dataset):
    """ Add DReaMIT lapse temperature at station level to the netCDF file"""
    vn  = 'T_lapse_station_C'  # variable name
    var = group.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Temperature at station level if solely determined by dynamically-computed linear lapse rate'
    var.comment   = ref_dreamit
    var.units     = 'degree_C'
    
    return group[vn]

def add_var_lapse(group: nc.Dataset):
    """ Add DReaMIT lapse rate to the netCDF file"""
    vn  = 'lapse_Cperm'  # variable name
    var = group.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Lower atmosphere dynamically-computed linear lapse rate'
    var.comment   = ref_dreamit
    var.units     = 'degree_C m-1'
    
    return group[vn]

def add_var_AIRT_DReaMIT(group: nc.Dataset):
    """ Add DReaMIT air T to the netCDF file"""
    vn  = 'AIRT_DReaMIT_C'  # variable name
    var = group.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'DReaMIT model air temperature corrected for surface-based inversions'
    var.comment   = ref_dreamit
    var.units     = 'degree_C'
    
    return group[vn]

def add_var_beta_t(group: nc.Dataset):
    """ Add DReaMIT reanalysis bias air T to the netCDF file"""
    vn  = 'beta_t_C'  # variable name
    var = group.createVariable(vn,'f4',('time',))
    var.long_name = 'reanalysis bias independent of stations'
    var.comment   = ref_dreamit
    var.units     = 'degree_C'
    
    return group[vn]

def format_temp_elev_inversions(reanalysis: str, T_pl_in: np.ndarray, h_pl_in: np.ndarray, T_pl_surface_in: np.ndarray, h_pl_surface_in: np.ndarray, grid_elev_in: np.ndarray):
    """
    Formatting of the interpolated GlobSim values to only account for elevations between 
    maximal elevation (5000m) and minimal elevation (lowest of grid or station)
    
    Parameters
    ----------
    reanalysis : str
        Name of the reanalysis, 'era5' or 'jra3qg
    T_pl_in : np.ndarray
        Pressure-level air temperature [K] with dimensions (time, level, station)
    h_pl_in : np.ndarray
        Pressure-level elevation [m] with dimensions (time, level, station)
    T_pl_surface_in : np.ndarray
        Pressure-level air temperature at surface (station) [K] with dimensions (time, station)
    h_pl_surface_in : np.ndarray
        Pressure-level elevation at surface (station) [m] with dimensions (station,)
    grid_elev_in : np.ndarray
        Time-invariant grid elevation [m] with dimensions (time, station)

    Returns
    ------- 
    T_pl : np.ndarray
        Pressure-level air temperature [K] with dimensions (time, level, station)
        all values for elevations outside the prescribed range are nans
    h_pl : np.ndarray
        Pressure-level elevation [m] with dimensions (time, level, station)
        all values for elevations outside the prescribed range are nans
        converted to meters for ERA5
    T_pl_surface : np.ndarray
        Pressure-level air temperature at surface (station) [K] with dimensions (time, station)
        all values for elevations outside the prescribed range are nans
    h_pl_surface : np.ndarray
        Pressure-level elevation at surface (station) [m] with dimensions (station,)
        all values for elevations outside the prescribed range are nans
    grid_elev : np.ndarray
        Time-invariant grid elevation [m] with dimensions (time, station)
        all values for elevations outside the prescribed range are nans
    """
    num_station = grid_elev_in.shape[-1]
    grid_elev = float(grid_elev_in[0,0])/const.G
    zmax = 5000 ############# !!!!!!!!!!!!!!!!!!! #############
    h_pl_in = h_pl_in[:]/(const.G if reanalysis=='era5' else 1)
    zmin = {station: np.min([h_pl_surface_in[station],grid_elev]) for station in range(num_station)}

    T_pl = {station: [] for station in range(num_station)}
    h_pl = {station: [] for station in range(num_station)}
    T_pl_surface = {station: [] for station in range(num_station)}
    h_pl_surface = {station: [] for station in range(num_station)}

    for station in range(num_station):
        mask_z = (h_pl_in[:,:,station] < zmin[station]) | (h_pl_in[:,:,station] > zmax)
        T_pl[station] = np.ma.array(T_pl_in[:,:,station], mask=mask_z).filled(np.nan)
        h_pl[station] = np.ma.array(h_pl_in[:,:,station], mask=mask_z).filled(np.nan)

        T_pl_surface[station] = T_pl_surface_in[:,station]
        h_pl_surface[station] = h_pl_surface_in[station]

    T_pl = np.stack(list(T_pl.values()), axis=-1) # dimensions (time, level, station)
    h_pl = np.stack(list(h_pl.values()), axis=-1) # dimensions (time, level, station)
    T_pl_surface = np.stack(list(T_pl_surface.values()), axis=-1) # dimensions (time, station)
    h_pl_surface = np.stack(list(h_pl_surface.values()), axis=-1) # dimensions (station)

    return T_pl, h_pl, T_pl_surface, h_pl_surface, grid_elev

def lin(X, a=0, b=0):
    T_lin = a * X + b
    return T_lin 

def quad(X, a=0, b=0, c=0):
    T_lin = a * X**2 + b * X + c
    return T_lin 

def compute_dreamit_metrics(reanalysis, T_pl, h_pl, T_pl_surface, h_pl_surface, grid_elev):
    """
    Calculate inversion metrics dynamically (top of inversion, lapse rate, etc)
    from formatted GlobSim interpolate input
    
    Parameters
    ----------
    reanalysis : str
        Name of the reanalysis, 'era5' or 'jra3qg
    T_pl : np.ndarray
        Pressure-level air temperature [K] with dimensions (time, level, station)
        all values for elevations outside the prescribed range are nans
    h_pl : np.ndarray
        Pressure-level elevation [m] with dimensions (time, level, station)
        all values for elevations outside the prescribed range are nans
        converted to meters for ERA5
    T_pl_surface : np.ndarray
        Pressure-level air temperature at surface (station) [K] with dimensions (time, station)
        all values for elevations outside the prescribed range are nans
    h_pl_surface : np.ndarray
        Pressure-level elevation at surface (station) [m] with dimensions (station,)
        all values for elevations outside the prescribed range are nans
    grid_elev : np.ndarray
        Time-invariant grid elevation [m] with dimensions (time, station)
        all values for elevations outside the prescribed range are nans

    Returns
    ------- 
    z_top_inversion_m : np.ndarray
        Top of inversion elevation [m] with dimensions (time, station)
    T_lapse_grid_C : np.ndarray
        Temperature at grid level if solely determined by dynamically-computed linear lapse rate [K] with dimensions (time, station)
    T_lapse_station_C : np.ndarray
        Temperature at station level if solely determined by dynamically-computed linear lapse rate [K] with dimensions (time, station)
    lapse_Cperm : np.ndarray
        Lower atmosphere dynamically-computed linear lapse rate [K/m] with dimensions (time, station)
    """
    num_station = len(h_pl_surface)

    df_inversion = {station: [] for station in range(num_station)}

    for station in range(num_station):
        # initializing empty lists to be filled iteratively
        list_z_top, list_T_lapse_grid, list_T_lapse_station, list_lapse, list_sea_level_temp = [], [], [], [], []
        
        # iteration through every time step (hourly here)
        for loc in range(T_pl.shape[0]):
            y = h_pl[loc,:,station]
            y = y[~np.isnan(y)]
            # x: list of temperatures (pressure-level) at each of the y elevations
            x = T_pl[loc,:,station]
            x = x[~np.isnan(x)]

            if reanalysis == 'era5': # flip so that we have decreasing elevations (already True for jra3qg)
                y = y[::-1]
                x = x[::-1]
            
            y = np.append(y,h_pl_surface[station])
            x = np.append(x,T_pl_surface[loc,station])

            #####################################
            # DETERMINGING THE TOP OF INVERSION #
            #####################################

            # we find where the evolution of x changes signs (e.g., goes from increasing to decreasing with elevation)
            # this will lossely correspond to the top of the inversion event, if it exists, but then we need to refine it
            loc_ind = np.argmax(np.diff(x)<0)

            if loc_ind == 0: # no inversion
                top_inv = np.nan
            else: # inversion 
                # we identify a 'band' around the turning point identified by y[loc_ind]
                # we look at this particular elevation, one above, and one beyond
                # the y spacing is variable with 100m to 300m
                # here the integration into GlobSim will need to involve some thinking
                band = np.array([x[loc_ind-1:loc_ind+2], y[loc_ind-1:loc_ind+2]])
                # then, we 'model' the turning point as a quadratic function of elevation and find
                # the actual turning point
                gmodel_quad = Model(quad)
                params_quad = gmodel_quad.make_params()
                a_quad = gmodel_quad.fit(band[0], params_quad, X=band[1])
                # top of inversion z_top [m]
                top_inv = -a_quad.params['b']/(2*a_quad.params['a'])
                if top_inv > np.max(y):
                    top_inv = np.nan
                if top_inv < np.min(y):
                    top_inv = np.nan

            list_z_top.append(top_inv)

            ############################
            # FITTING THE LINEAR LAPSE #
            ############################

            # here we will try to find the optimal range to fit the linear lapse rate
            # this is dne by including more and more points into the linear fitting
            # routine and stopping when AIC is the best   

            # iteratively adding elevation points (starting from top of column)
            # to fit the linear lapse rate, while staying above the top of inversion z_top
            # we stop when AIC is the best
            index=3 # we need at least 3 points to fit the lapse rate
            while y[index] > np.nanmax([np.min(y), top_inv]):
                gmodel_lin = Model(lin)
                params_lin = gmodel_lin.make_params()
                a_lin = gmodel_lin.fit(x[:index], params_lin, X=y[:index])
                best_aic_lin = 1e10 # ridiculously high AIC to begin with, mark to beat (easy)
                if a_lin.aic > best_aic_lin:
                    break
                else:
                    best_aic_lin = a_lin.aic
                    best_params_lin = a_lin.params.valuesdict().values()
                index+=1

            # lapse rate [C/m]
            list_lapse.append(list(best_params_lin)[0])
            # linear continuation of the laspe rate at sea level: T_lapse_sea_level [C]
            list_sea_level_temp.append(list(best_params_lin)[1])

        # surface temperature from reanalysis data T_sur [C]
        list_T_lapse_grid = lin(grid_elev, np.array(list_lapse), np.array(list_sea_level_temp))
        list_T_lapse_station = lin(np.min(y), np.array(list_lapse), np.array(list_sea_level_temp))

        # we create the panda dataframe putting everything together
        df_inversion[station] = pd.DataFrame([list_z_top, list_T_lapse_grid, list_T_lapse_station, list_lapse]).T
        list_new_cols = ['z_top_inversion_m', 'T_lapse_grid_C', 'T_lapse_station_C', 'lapse_Cperm']
        df_inversion[station] = df_inversion[station].rename(columns={k:v for k,v in zip(list(df_inversion[station].columns), list_new_cols)})

        # we make sure we are using good variable types
        df_inversion[station] = df_inversion[station].apply(pd.to_numeric, errors='coerce')

    z_top_inversion_m = np.stack([i['z_top_inversion_m'] for i in df_inversion.values()], axis=-1)
    T_lapse_grid_C    = np.stack([i['T_lapse_grid_C'] for i in df_inversion.values()], axis=-1)
    T_lapse_station_C = np.stack([i['T_lapse_station_C'] for i in df_inversion.values()], axis=-1)
    lapse_Cperm       = np.stack([i['lapse_Cperm'] for i in df_inversion.values()], axis=-1)
    
    return z_top_inversion_m, T_lapse_grid_C, T_lapse_station_C, lapse_Cperm

def dreamit_metrics(reanalysis:str, T_pl_in: np.ndarray, h_pl_in: np.ndarray, T_pl_surface_in: np.ndarray, h_pl_surface_in: np.ndarray, grid_elev_in: np.ndarray):
    """
    Calculate inversion metrics dynamically (top of inversion, lapse rate, etc)
    from raw GlobSim interpolate input
    
    Parameters
    ----------
    reanalysis : str
        Name of the reanalysis, 'era5' or 'jra3qg
    T_pl_in : np.ndarray
        Pressure-level air temperature [K] with dimensions (time, level, station)
    h_pl_in : np.ndarray
        Pressure-level elevation [m] with dimensions (time, level, station)
    T_pl_surface_in : np.ndarray
        Pressure-level air temperature at surface (station) [K] with dimensions (time, station)
    h_pl_surface_in : np.ndarray
        Pressure-level elevation at surface (station) [m] with dimensions (station,)
    grid_elev_in : np.ndarray
        Time-invariant grid elevation [m] with dimensions (time, station)

    Returns
    ------- 
    z_top_inversion_m : np.ndarray
        Top of inversion elevation [m] with dimensions (time, station)
    T_lapse_grid_C : np.ndarray
        Temperature at grid level if solely determined by dynamically-computed linear lapse rate [K] with dimensions (time, station)
    T_lapse_station_C : np.ndarray
        Temperature at station level if solely determined by dynamically-computed linear lapse rate [K] with dimensions (time, station)
    lapse_Cperm : np.ndarray
        Lower atmosphere dynamically-computed linear lapse rate [K/m] with dimensions (time, station)
    """
    T_pl, h_pl, T_pl_surface, h_pl_surface, grid_elev = format_temp_elev_inversions(reanalysis, T_pl_in, h_pl_in, T_pl_surface_in, h_pl_surface_in, grid_elev_in)
    z_top_inversion_m, T_lapse_grid_C, T_lapse_station_C, lapse_Cperm = compute_dreamit_metrics(reanalysis, T_pl, h_pl, T_pl_surface, h_pl_surface, grid_elev)

    return z_top_inversion_m, T_lapse_grid_C, T_lapse_station_C, lapse_Cperm

def get_model_params(reanalysis:str):
    """
    Get DReaMIT model parameters (found in globsim/globsim/data/DReaMIT_params.csv)
    
    Parameters
    ----------
    reanalysis : str
        Name of the reanalysis, 'era5' or 'jra3qg

    Returns
    ------- 
    list_params : list
        list of parameter values for DReaMIT model and given reanalysis
    """
    current_dir = os.path.dirname(__file__)
    csv_path = os.path.join(current_dir, 'data', 'DReaMIT_params.csv')
    list_name_params = ['alpha_slope', 'alpha_intercept', 'beta_amp', 't_star', 'beta_bias']
    mod_pars = pd.read_csv(csv_path).set_index('reanalysis')
    list_params = mod_pars.loc[reanalysis, list_name_params] # here we extract the calibrated model parameters (list)

    return list_params

def time_frac_year(nc_time: nc.Variable):
    """
    Compute time fraction within a year
    (0=first second on January 1st, 1=last second on December 31st)
    
    Parameters
    ----------
    nc_time : netCDF4.Variable
        Name of the reanalysis, 'era5' or 'jra3qg

    Returns
    ------- 
    time_frac_year : np.ndarray
        time fraction within a year, no unit [-], with dimensions (time,)
    """
    df_time = pd.DataFrame(nc.num2date(nc_time[:], nc_time.units, nc_time.calendar, only_use_cftime_datetimes=False),columns=['time'])
    # df_time = pd.DataFrame(nc.num2date(time_in[:], time_in.units),columns=['time'])
    df_time = df_time.astype({'time': 'str'}).astype({'time': 'datetime64[ns]'})
    time_frac_year = np.array((pd.to_datetime(df_time['time'])-pd.to_datetime(df_time['time'].dt.year, format='%Y')).dt.total_seconds() / (pd.to_datetime(df_time['time'].dt.year+1, format='%Y')-pd.to_datetime(df_time['time'].dt.year, format='%Y')).dt.total_seconds())

    return time_frac_year

def dreamit_air_T(T_lapse_grid: np.ndarray, T_lapse_station: np.ndarray, T_sur: np.ndarray, time_frac_year: np.ndarray, hyps: np.ndarray, params: list):
    """
    Compute the DReaMIT model air temperature corrected for surface-based inversions
    
    Parameters
    ----------
    T_lapse_grid : np.ndarray
        Temperature at grid level if solely determined by dynamically-computed linear lapse rate [K] with dimensions (time, station)
    T_lapse_station : np.ndarray
        Temperature at station level if solely determined by dynamically-computed linear lapse rate [K] with dimensions (time, station)
    T_sur : np.ndarray
        Surface temperature [K] with dimensions (time, station)
    time_frac_year : np.ndarray
        time fraction within a year, no unit [-], with dimensions (time,)
    hyps : np.ndarray
        hypsometry (from GlobSim csv configuration file), no unit [-], with dimensions (station,)
    params : list
        list of parameter values for DReaMIT model and given reanalysis (found in globsim/globsim/data/DReaMIT_params.csv)

    Returns
    ------- 
    AIRT_DReaMIT_C : np.ndarray
        DReaMIT model air temperature corrected for surface-based inversions [C], with dimensions (time,station)
    beta_t_C : np.ndarray
        reanalysis bias independent of stations [C], with dimensions (time,)
    """
    alpha_slope, alpha_intercept, beta_amp, t_star, beta_bias = params
    beta_t_C  = (beta_amp * np.cos(2*np.pi*(time_frac_year-t_star)) + beta_bias)
    num_station = len(hyps)
    T_mod = {n: [] for n in range(num_station)}
    for n in range(num_station):
        DT_grid = T_sur[:,n] - T_lapse_grid[:,n]
    
        alpha_h = (alpha_intercept + (np.exp(alpha_slope*hyps[n])-1))
    
        T_mod[n] = (alpha_h * DT_grid + beta_t_C + T_lapse_station[:,n])

    AIRT_DReaMIT_C = np.stack(list(T_mod.values()), axis=-1)
    
    return AIRT_DReaMIT_C, beta_t_C
