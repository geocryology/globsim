import netCDF4 as nc

import globsim.dreamit as dreamit
from globsim.redcapp import add_var_delta_T


def PRESS_Pa_pl(rg: nc.Dataset, reanalysis: str):
    vn = 'PRESS_pl'
    var           = rg.createVariable(vn, 'f4',('time','station'))
    var.long_name = 'air_pressure {} pressure levels only'.format(reanalysis)
    var.units     = 'Pa'
    var.standard_name = 'surface_air_pressure'
    
    return vn


def AIRT_C_pl(rg: nc.Dataset, reanalysis: str):
    vn = 'AIRT_pl'
    var           = rg.createVariable(vn,'f4',('time','station'))
    var.long_name = 'air_temperature {} pressure levels only'.format(reanalysis)
    var.units     = 'degrees_C'
    var.standard_name = 'air_temperature'

    return vn


def AIRT_C_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'AIRT_sur'
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = '2_metre_temperature {} surface only'.format(reanalysis)
    var.units     = 'degrees_C'
    var.standard_name = 'air_temperature'

    return vn


def PREC_mm_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'PREC_sur'  
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Total precipitation {} surface only'.format(reanalysis)
    var.units     = 'kg m-2 s-1'
    var.comment = "units [kg m-2 s-1] corresponds to [mm s-1] for water (density 1000 [kg m-3])"
    var.standard_name = 'precipitation_flux'

    return vn


def RH_per_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'RH_sur' 
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'relative humidity {} surface only'.format(reanalysis)
    var.units     = 'percent'
    var.standard_name = 'relative_humidity'

    return vn


def RH_per_pl(rg: nc.Dataset, reanalysis: str):
    vn = 'RH_pl' 
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'relative humidity {} pressure-level only'.format(reanalysis)
    var.units     = 'percent'
    var.standard_name = 'relative_humidity'

    return vn

def WIND_sur(rg: nc.Dataset, reanalysis: str):
    vn_u = '10 metre U wind component' 
    var           = rg.createVariable(vn_u,'f4',('time', 'station'))
    var.long_name = '10 metre U wind component'
    var.units     = 'm s-1' 
       
    vn_v = '10 metre V wind component'  
    var           = rg.createVariable(vn_v,'f4',('time', 'station'))
    var.long_name = '10 metre V wind component'
    var.units     = 'm s-1'

    vn_spd = 'WSPD_sur'
    var           = rg.createVariable(vn_spd,'f4',('time', 'station'))
    var.long_name = '10 metre wind speed {} surface only'.format(reanalysis)
    var.units     = 'm s-1'
    var.standard_name = 'wind_speed'

    vn_dir = 'WDIR_sur'  
    var           = rg.createVariable(vn_dir,'f4',('time', 'station'))
    var.long_name = '10 metre wind direction {} surface only'.format(reanalysis)
    var.units     = 'degree'
    var.standard_name = 'wind_from_direction'

    return vn_u, vn_v, vn_spd, vn_dir


def SW_Wm2_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'SW_sur' 
    var           = rg.createVariable(vn, 'f4', ('time', 'station'))
    var.long_name = 'Surface solar radiation downwards {} surface only'.format(reanalysis)
    var.units     = 'W m-2'
    var.standard_name = 'surface_downwelling_shortwave_flux'

    return vn


def LW_Wm2_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'LW_sur'  # variable name
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Surface thermal radiation downwards {} surface only'.format(reanalysis)
    var.units     = 'W m-2'
    var.standard_name = 'surface_downwelling_longwave_flux'

    return vn


def SH_kgkg_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'SH_sur' 
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Specific humidity {} surface only'.format(reanalysis)
    var.units     = 'kg kg-1'
    var.standard_name = 'specific_humidity'

    return vn


def LW_Wm2_topo(rg: nc.Dataset, reanalysis: str):
    vn = 'LW_topo' 
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'TOPOscale-corrected thermal radiation downwards {}'.format(reanalysis)
    var.standard_name = 'surface_downwelling_longwave_flux'
    var.units     = 'W m-2'

    return vn


def SW_Wm2_topo(rg: nc.Dataset, reanalysis: str):
    vn_diff = 'SW_topo_diffuse'  # variable name
    var           = rg.createVariable(vn_diff,'f4',('time', 'station'))
    var.long_name = 'TOPOscale-corrected diffuse solar radiation {}'.format(reanalysis)
    var.units     = 'W m-2'
    var.standard_name = 'surface_diffuse_downwelling_shortwave_flux_in_air'

    vn_dir = 'SW_topo_direct'  # variable name
    var           = rg.createVariable(vn_dir,'f4',('time', 'station'))
    var.long_name = 'TOPOscale-corrected direct solar radiation {}'.format(reanalysis)
    var.units     = 'W m-2'
    var.standard_name = 'surface_direct_downwelling_shortwave_flux_in_air'

    vn_glob = 'SW_topo_global'  # variable name
    var           = rg.createVariable(vn_glob,'f4',('time', 'station'))
    var.long_name = 'TOPOscale-corrected global solar radiation {}'.format(reanalysis)
    var.units     = 'W m-2'
    var.standard_name = 'surface_downwelling_shortwave_flux_in_air'

    return vn_dir, vn_diff, vn_glob


def PRECCORR_mm_sur(rg: nc.Dataset, reanalysis: str):
    vn = 'PRECCORR_sur'  # variable name
    var           = rg.createVariable(vn,'f4',('time', 'station'))
    var.long_name = 'Corrected Total precipitation {} surface only'.format(reanalysis)
    var.units     = 'kg m-2 s-1'
    var.comment = "units [kg m-2 s-1] corresponds to [mm/s] for water (density 1000 [kg m-3])"
    var.standard_name = 'precipitation_flux'

    return vn


def AIRT_DReaMIT(rg: nc.Dataset):
    dreamit.add_var_z_top_inversion(rg)
    dreamit.add_var_T_lapse_grid(rg)
    dreamit.add_var_T_lapse_station(rg)
    dreamit.add_var_lapse(rg)
    dreamit.add_var_AIRT_DReaMIT(rg)
    dreamit.add_var_beta_t(rg)

def AIRT_redcapp_DeltaT(rg: nc.Dataset):
    _ = add_var_delta_T(rg)
    return 'AIRT_redcapp_DeltaT'

__all__ = ['PRESS_Pa_pl',
           'AIRT_C_pl', 
           'AIRT_C_sur',
           'PREC_mm_sur',
           'RH_per_sur',
           'RH_per_pl',
           'WIND_sur',
           'SW_Wm2_sur',
           'LW_Wm2_sur',
           'SH_kgkg_sur',
           'LW_Wm2_topo',
           'SW_Wm2_topo',
           'PRECCORR_mm_sur',
           'AIRT_redcapp_DeltaT',
           'AIRT_DReaMIT']
