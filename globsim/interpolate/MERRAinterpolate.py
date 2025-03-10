import netCDF4 as nc
import numpy as np
import xarray as xr
import sys
import logging
import gc

from datetime          import datetime
from os                import path, makedirs
from pathlib import Path

from globsim.common_utils import str_encode, variables_skip
from globsim.interpolate.GenericInterpolate import GenericInterpolate
from globsim.nc_elements import netcdf_base, new_interpolated_netcdf
from globsim.interp import calculate_weights, ele_interpolate, extrapolate_below_grid

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='netCDF4')

logger = logging.getLogger('globsim.interpolate')


class MERRAinterpolate(GenericInterpolate):
    """
    Algorithms to interpolate MERRA-2 netCDF files to station coordinates.
    All variables retains their original units and time-steps.

    Referenced from era_interim.py (Dr.Stephan Gruber): Class ERAinterpolate()
    """
    REANALYSIS = 'merra2'
    
    def __init__(self, ifile, **kwargs):
        super().__init__(ifile, **kwargs)
        par = self.par

        self.input_dir = path.join(par['project_directory'],'merra2')

        # Load MF Datasets
        self.mf_sc = xr.open_mfdataset(path.join(self.input_dir,'merra_sc.nc'), decode_times=False)
        self.mf_sa = xr.open_mfdataset(self.prefilter_mf_paths(path.join(self.input_dir,'merra_sa_*.nc')), decode_times=False)
        self.mf_sf = xr.open_mfdataset(self.prefilter_mf_paths(path.join(self.input_dir,'merra_sf_*.nc')), decode_times=False)
        self.mf_pl = xr.open_mfdataset(self.prefilter_mf_paths(path.join(self.input_dir,'merra_pl_*.nc')), decode_times=False)

        # Check dataset integrity
        if not self.skip_checks:
            logger.info("Check data integrity (sa)")
            self.ensure_datset_integrity(self.mf_sa['time'], 1)
            logger.info("Check data integrity (sf)")
            self.ensure_datset_integrity(self.mf_sf['time'], 1)
            logger.info("Check data integrity (pl)")
            self.ensure_datset_integrity(self.mf_pl['time'], 6)
            logger.info("Data integrity ok")

    def MERRA2station(self, ncf_in, ncfile_out, points,
                      variables=None, date=None):
        """
        Given the type of variables to interpoalted from MERRA2 downloaded diretory
        Create the empty of netCDF file to hold the interpolated results, by calling
        self.netCDF_empty
        Get the interpolated results from MERRA2station
        Append all variables into the empty netCDF file
        Close all files

        Args:
        ncf_in:  MFDatast of an MERRA-2 derived netCDF file. 

        ncfile_out: Full path to the output netCDF file to write.

        points: A dataframe of locations. See method StationListRead in
                common_utils.py for more details.

        variables:  List of variable(s) to interpolate such as
                    ['T','RH','U','V',' T2M', 'U2M', 'V2M', 'U10M', 'V10M',
                    'PRECTOT', 'SWGDN','SWGDNCLR','LWGDN', 'LWGDNCLR'].
                    Defaults to using all variables available.

        date: Directory to specify begin and end time for the derived time
                series. Defaluts to using all times available in ncf_in.

        """
        # Check station bounds
        self.validate_stations_extent(ncf_in)

        # is it a file with pressure levels?
        pl = self.vn_level in ncf_in.sizes.keys()
        
        level_var = self.vn_level if pl else None
        
        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]
        t_unit = "hours since 1980-01-01 00:00:00"  # ncf_in.variables['time'].units
        try :
            t_cal = ncf_in.variables['time'].calendar
        except AttributeError :  # Attribute doesn't exist
            t_cal = u"gregorian"  # or standard
        # TODO: rm time = [nc.num2date(timei, units=t_unit, calendar=t_cal) for timei in nctime]
        # TODO: rm time = np.asarray(time)

        # detect invariant files (topography etc.)
        invariant = True if len(np.unique(nctime)) <= 2 else False

        # restrict to date/time range if given
        if date is None:
            tmask = nctime < nc.date2num(datetime(3000, 1, 1), units=t_unit, calendar=t_cal)
        else:
            beg_num = nc.date2num(date['beg'], units=t_unit, calendar=t_cal)
            end_num = nc.date2num(date['end'], units=t_unit, calendar=t_cal)
            tmask = (nctime < end_num) * (nctime >= beg_num)

        if not any(tmask):
            sys.exit('''\n ERROR: No downloaded data exist within date range specified by interpolation control file.
                    Download new data or change 'beg' / 'end' in interpolation control file''')

        # get time indices
        time_in = nctime[tmask]
        
        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter  = len(time_in) // self.cs
        niter += ((len(time_in) % self.cs) > 0)

        # Create source grid
        sgrid = self.create_source_grid(ncf_in)
        subset_grid, lon_slice, lat_slice = self.create_subset_source_grid(sgrid, self.stations_bbox)
        # build the output of empty netCDF file
        if self.resume:
            self.require_file_can_be_resumed(ncfile_out)
            
        if not (self.resume and Path(ncfile_out).exists()):
            #rootgrp = self.netCDF_empty(ncfile_out, self.stations, ncf_in)
            rootgrp = new_interpolated_netcdf(ncfile_out, self.stations, ncf_in,
                                              time_units=t_unit,
                                              level_var=level_var,
                                              n_time = len(time_in))
            rootgrp.source = f'{self.REANALYSIS}, interpolated bilinearly to stations'
            rootgrp.globsim_interpolate_start = self.par['beg']
            rootgrp.globsim_interpolate_end = self.par['end']
            rootgrp.globsim_chunk_size = self.cs
            rootgrp.globsim_interpolate_success = 0
            rootgrp.globsim_last_chunk_written = -1
            rootgrp.close()

        # open the output netCDF file, set it to be appendable ('a')
        with nc.Dataset(ncfile_out, 'a') as ncf_out:
            # write time
            ncf_out.variables['time'][:] = time_in

            # loop in chunk size cs
            for n in range(niter):
                if self.resume and n <= ncf_out.globsim_last_chunk_written:
                    if n == ncf_out.globsim_last_chunk_written:
                        logger.info(f"Resuming interpolation at chunk {n+1}")
                    continue
                
                self.require_safe_mem_usage()

                # indices
                beg = n * self.cs
                # restrict last chunk to lenght of tmask plus one (to get last time)
                if invariant:
                    end = beg
                else:
                    end = min(n * self.cs + self.cs, len(time_in)) - 1

                # make tmask for chunk
                beg_time = time_in[beg]
                if invariant:
                    # allow topography to work in same code, len(nctime) = 1
                    # TODO: rm end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                    end_time = nctime[0]
                    # end = 1
                else:
                    # TODO: rm end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)
                    end_time = time_in[end]

                # !! CAN'T HAVE '<= end_time', would damage appeding
                tmask_chunk = (nctime <= end_time) * (nctime >= beg_time)
                if invariant:
                    # allow topography to work in same code
                    tmask_chunk = np.array([True])

                # get the interpolated variables
                dfield, variables = self.interp2D(ncf_in,
                                                self.stations, tmask_chunk,
                                                subset_grid, lon_slice, lat_slice,
                                                variables=None, date=None)
                tbeg = nc.num2date(ncf_in[self.vn_time][np.where(tmask_chunk)[0][0]], ncf_in[self.vn_time].units, ncf_in[self.vn_time].calendar)
                tend = nc.num2date(ncf_in[self.vn_time][np.where(tmask_chunk)[0][-1]],  ncf_in[self.vn_time].units, ncf_in[self.vn_time].calendar)
                logger.info(f"{Path(ncf_out.filepath()).name} -- {tbeg} to {tend}")
                
                # append variables
                self.write_dfield_to_file(dfield, variables, ncf_out, beg, end, pl)
                ncf_out.globsim_last_chunk_written = n

                del dfield, tmask_chunk
                gc.collect()
                
            # Write success flag
            ncf_out.globsim_interpolate_success = 1
        
        ncf_in.close()

    def levels2elevation(self, ncfile_in, ncfile_out):
        """
        Linear 1D interpolation of pressure level data available for individual
        stations to station elevation. Where and when stations are below the
        lowest pressure level, they are assigned the value of the lowest
        pressure level.

        """
        # open file

        ncf = nc.Dataset(ncfile_in, 'r')
        height = ncf.variables['height'][:]
        nt = len(ncf.variables['time'][:])
        nl = len(ncf.variables['level'][:])

        # list variables
        varlist = [str_encode(x) for x in ncf.variables.keys()]
        for V in ['time', 'station', 'latitude', 'longitude', 'level', 'height', 'H']:
            varlist.remove(V)

        # === open and prepare output netCDF file =============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        if not (self.resume and Path(ncfile_out).exists()):
            rootgrp = netcdf_base(ncfile_out, 
                                  len(height),
                                  nt,
                                  'hours since 1980-01-01 00:00:00')
            rootgrp.source  = 'MERRA-2, interpolated (bi)linearly to stations'

            time = rootgrp["time"]
            station = rootgrp["station"]
            latitude = rootgrp["latitude"]
            longitude = rootgrp["longitude"]
            height = rootgrp["height"]

            # assign base variables
            time[:]      = ncf.variables['time'][:]
            station[:]   = ncf.variables['station'][:]
            latitude[:]  = ncf.variables['latitude'][:]
            longitude[:] = ncf.variables['longitude'][:]
            height[:]    = ncf.variables['height'][:]

            rootgrp.globsim_interpolate_success = 0
            rootgrp.last_station_written = -1
            rootgrp.vars_written = ""

            # create and assign variables from input file
            for var in varlist:
                tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))
                tmp.long_name = str_encode(ncf.variables[var].long_name)
                tmp.units     = str_encode(ncf.variables[var].units)

            # add air pressure as new variable
            var = 'air_pressure'
            varlist.append(var)
            tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))
            tmp.long_name = str_encode(var)
            tmp.units     = str_encode('hPa')
            rootgrp.close()
        # end file prepation ==================================================
        
        with nc.Dataset(ncfile_out, 'a') as rootgrp:            
            # loop over stations
            for n, h in enumerate(height):
                if self.resume and n <= rootgrp.last_station_written:
                    if n == rootgrp.last_station_written:
                        logger.info(f"Resuming interpolation at station {n+1}")
                    continue

                self.require_safe_mem_usage()
                logger.debug(f"Interpolating station {n+1} to station elevation using pressure-level data")
                # geopotential unit: height [m]
                # shape: (time, level)
                elevation = ncf.variables['H'][:,:,n]
                # TODO: check if height of stations in data range (+50m at top,
                # lapse r.)

                elev_diff, va, vb = ele_interpolate(elevation, h, nl)
                wa, wb = calculate_weights(elev_diff, va, vb)

                # loop over variables and apply interpolation weights
                for v, var in enumerate(varlist):
                    if var in str(rootgrp.vars_written).split(" "):
                            logger.debug(f"Skipping {var}")
                            continue
                    
                    if var == 'air_pressure':
                        # pressure [hPa] variable from levels, shape: (time, level)
                        data = np.repeat([ncf.variables['level'][:]],
                                        len(time),axis=0).ravel()
    
                        # 2025-01-28 [NB]: I don't think this block of code is necessary
                        """level_highest = ncf.variables['level'][:][-1]
                        level_lowest = ncf.variables['level'][:][0]

                        for j, value in enumerate(ipol):
                            if value == level_highest:
                                ipol[j] = level_lowest"""
                    else:
                        data = ncf.variables[var][:,:,n].ravel()
                    
                    ipol = data[va] * wa + data[vb] * wb   # interpolated value

                    if self.extrapolate_below_grid:
                            extrapolated_values = extrapolate_below_grid(elevation, data, h)
                            ipol = np.where(~extrapolated_values.mask, extrapolated_values, ipol)
                            
                    rootgrp.variables[var][:,n] = ipol  # assign to file
                    rootgrp.vars_written = " ".join(set(str(rootgrp.vars_written).split(" ") + [var]))
                
                rootgrp.vars_written = ""
                rootgrp.last_station_written = n
                
            rootgrp.globsim_interpolate_success = 1

        ncf.close()
        # closed file ==========================================================
    
    def _preprocess(self):
        if not path.isdir(self.output_dir):
            makedirs(self.output_dir)
        
        self.MERRA2station(self.mf_sc,
                           path.join(self.output_dir,'merra2_sc_' + self.list_name + '.nc'),
                           self.stations, ['PHIS','FRLAND'], date=None)
    
    def _process_sa(self):
        # === 2D Interpolation for Surface Analysis Data ===
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['T2M'],  # [K] 2m values
                'wind_speed' : ['U2M', 'V2M', 'U10M','V10M'],   # [m s-1] 2m & 10m values
                'relative_humidity' : ['QV2M']}  # 2m value
        varlist = self.TranslateCF2short(dpar)
        if self.resume and self.completed_successfully(self.getOutFile('sa')):
            logger.info("Skipping surface analysis interpolation")
        else:
            self.MERRA2station(self.mf_sa,
                            self.getOutFile('sa'),
                            self.stations, varlist, date=self.date)
    
    def _process_sf(self):
        # 2D Interpolation for Single-level Radiation Diagnostics Data 'SWGDN',
        # 'LWGDN', 'SWGDNCLR'. 'LWGDNCLR'
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['T2MDEW'],  # [K] 2m values
                'precipitation_amount' : ['PRECTOT','PRECTOTCORR'],  # [kg/m2/s] total precipitation
                'downwelling_shortwave_flux_in_air' : ['SWGDN'],  # [W/m2] short-wave downward
                'downwelling_longwave_flux_in_air'  : ['LWGDN'],  # [W/m2] long-wave downward
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['SWGDNCLR'],  # [W/m2] short-wave downward assuming clear sky
                'downwelling_longwave_flux_in_air_assuming_clear_sky': ['LWGDNCLR']}  # [W/m2] long-wave downward assuming clear sky
        varlist = self.TranslateCF2short(dpar)
        if self.resume and self.completed_successfully(self.getOutFile('sf')):
            logger.info("Skipping surface forecast interpolation")
        else:
            self.MERRA2station(self.mf_sf,
                            self.getOutFile('sf'),
                            self.stations, varlist, date=self.date)
            
    def _process_pl(self):
        # NEED ADD 'H' in it!
        # === 2D Interpolation for Pressure-Level, Analyzed Meteorological DATA ===
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['T'],           # [K]
                'wind_speed'        : ['U', 'V'],      # [m s-1]
                'relative_humidity' : ['RH']}          # [1]
        varlist = self.TranslateCF2short(dpar).append('H')
        if self.resume and self.completed_successfully(self.getOutFile('pl')):
            logger.info("Skipping pressure level interpolation")
        else:
            self.MERRA2station(self.mf_pl,
                            self.getOutFile('pl'),
                            self.stations, varlist, date=self.date)

        # 1D Interpolation for Pressure Level Analyzed Meteorological Data
        self.levels2elevation(path.join(self.output_dir,'merra2_pl_' + self.list_name + '.nc'),
                              path.join(self.output_dir,'merra2_pl_' + self.list_name + '_surface.nc'))

        

        
