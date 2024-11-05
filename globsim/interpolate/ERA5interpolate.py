import netCDF4 as nc
import gc
import xarray as xr
import logging
import numpy as np
import urllib3

from datetime import datetime
from os import path
from pathlib import Path

from globsim.common_utils import variables_skip, str_encode
from globsim.interpolate.GenericInterpolate import GenericInterpolate
from globsim.nc_elements import netcdf_base, new_interpolated_netcdf
from globsim.interp import ele_interpolate, calculate_weights
import globsim.constants as const

logger = logging.getLogger('globsim.interpolate')


urllib3.disable_warnings()


class ERA5interpolate(GenericInterpolate):
    """
    Collection of methods to interpolate ERA5 netCDF files to station
    coordinates. All variables retain their original units and time stepping.
    """

    def __init__(self, ifile, era5type='reanalysis', **kwargs):
        super().__init__(ifile=ifile, **kwargs)
        self.era5type = era5type

        if self.era5type == 'reanalysis':
            self._set_input_directory("era5")
            self.ens = False
        elif self.era5type == 'ensemble_members':
            self._set_input_directory("era5ens")
            self.ens = True
            
        # convert longitude to ERA notation if using negative numbers
        self.stations['longitude_dd'] = self.stations['longitude_dd'] % 360

        self.mf_to = xr.open_mfdataset(self.get_input_file_glob('to'), decode_times=False)
        
        # Check dataset integrity
        if not self.skip_checks:
            with xr.open_mfdataset(self.get_input_file_paths('sa'), decode_times=False) as sa:
                logger.info("Check data integrity (sa)")
                self.ensure_datset_integrity(sa[self.vn_time], 3600)
            
            with xr.open_mfdataset(self.get_input_file_paths('sf'), decode_times=False) as sf:
                logger.info("Check data integrity (sf)")
                self.ensure_datset_integrity(sf[self.vn_time], 3600)

            with xr.open_mfdataset(self.get_input_file_paths('pl'), decode_times=False) as pl:
                logger.info("Check data integrity (pl)")
                self.ensure_datset_integrity(pl[self.vn_time], 3600)

            logger.info("Data integrity ok")

    @property
    def vn_time(self):
        return 'valid_time'
    
    @property
    def vn_level(self):
        return 'pressure_level'
    
    def get_input_file_glob(self, levStr):
        """ Return the input file(s) for the interpolation. """
        # edited naming conventions for simplicity and to avoid errors
        if self.ens:
            typeStr = '_ens'
        else:
            typeStr = ''

        f_glob = 'era5{}_{}_*.nc'

        if levStr == 'to':
            infile = path.join(self.input_dir, 'era5{}_to.nc'.format(typeStr))
        else:
            infile = path.join(self.input_dir, f_glob.format(typeStr, levStr))

        return infile

    def get_input_file_paths(self, levStr):
        glob = self.get_input_file_glob(levStr)
        infiles = self.prefilter_mf_paths(glob)
        return infiles
    
    def getOutFile(self, levStr):
        if self.ens:
            nome = 'era5_ens_{}_'.format(levStr) + self.list_name + '.nc'
        else:
            nome = 'era5_{}_'.format(levStr) + self.list_name + '.nc'
        outfile = path.join(self.output_dir, nome)
        return outfile

    def ERA2station(self, ncf_in, ncfile_out, points,
                    variables=None, date=None):

        """
        Biliner interpolation from fields on regular grid (latitude, longitude)
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all ERA5 files). The type
        of variable and file structure are determined from the input.

        This function creates an empty of netCDF file to hold the interpolated
        results, by calling new_interpolated_netcdf(). Then, data is
        interpolated in temporal chunks and appended. The temporal chunking can
        be set in the interpolation parameter file.

        Args:
        ncf_in: Globsim-derived ERA5 dataset (single file or MFDataset)

        ncfile_out: Full path to the output netCDF file to write.

        points: A dataframe of locations. See method StationListRead in
                common_utils.py for more details.

        variables:  List of variable(s) to interpolate such as
                    ['r', 't', 'u','v', 't2m', 'u10', 'v10', 'ssrd', 'strd', 'tp'].
                    Defaults to using all variables available.

        date: Dictionary to specify begin and end time for the derived time
                series. Defaluts to using all times available in ncfile_in.
        """

        # Check station bounds
        self.validate_stations_extent(ncf_in)

        # is it a file with pressure levels?
        pl = self.vn_level in ncf_in.sizes.keys()
        ens = 'number' in ncf_in.sizes.keys()

        # reduce chunk size for pressure-level interpolation
        cs = self.cs
        if pl:
            cs = cs // len(ncf_in.variables[self.vn_level][:])  # get actual number of levels

        # get time and convert to datetime object
        nctime = ncf_in.variables[self.vn_time][:]
        
        t_unit = ncf_in.variables[self.vn_time].attrs['units']
        try:
            t_cal = ncf_in.variables[self.vn_time].attrs['calendar']
        except AttributeError:  # attribute doesn't exist
            t_cal = u"gregorian"  # standard
        time = nc.num2date(nctime, units=t_unit, calendar=t_cal)

        # detect invariant files (topography etc.)
        invariant = True if len(time) == 1 else False

        # restrict to date/time range if given
        if date is None:
            tmask = np.array([True])
        else:
            tmask = (time < date['end']) * (time >= date['beg'])

        # get time vector for output
        time_in = nctime[tmask]
        
        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter = len(time_in) // cs
        niter += ((len(time_in) % cs) > 0)
        
        # Create source grid
        sgrid = self.create_source_grid(ncf_in)
        subset_grid, lon_slice, lat_slice = self.create_subset_source_grid(sgrid, self.stations_bbox)

        # build the output of empty netCDF file
        level_var = self.vn_level if pl else None
        
        if self.resume:
            self.require_file_can_be_resumed(ncfile_out)
        
        if not (self.resume and Path(ncfile_out).exists()):
            rootgrp = new_interpolated_netcdf(ncfile_out, self.stations, ncf_in,
                                            time_units=ncf_in[self.vn_time].units,  # 'hours since 1900-01-01 00:00:0.0'
                                            calendar=ncf_in[self.vn_time].calendar,
                                            level_var=level_var,
                                            n_time = len(time_in))
            
            rootgrp.globsim_interpolate_start = self.par['beg']
            rootgrp.globsim_interpolate_end = self.par['end']
            rootgrp.globsim_chunk_size = self.cs
            rootgrp.globsim_interpolate_success = 0
            rootgrp.globsim_last_chunk_written = -1

            if self.ens:
                rootgrp.source = 'ERA5 10-member ensemble, interpolated bilinearly to stations'
            else:
                rootgrp.source = 'ERA5, interpolated bilinearly to stations'

            rootgrp.close()  # close the file to write the header (needed for appending?)

        # open the output netCDF file, set it to be appendable ('a')
        with nc.Dataset(ncfile_out, 'a') as ncf_out:
            # write time
            ncf_out.variables['time'][:] = time_in

            # loop over chunks
            for n in range(niter):
                if self.resume and n <= ncf_out.globsim_last_chunk_written:
                    if n == ncf_out.globsim_last_chunk_written:
                        logger.info(f"Resuming interpolation at chunk {n+1}")
                    continue

                self.require_safe_mem_usage()

                # indices (relative to index of the output file)
                beg = n * cs
                # restrict last chunk to length of tmask plus one (to get last time)
                end = min(n * cs + cs, len(time_in)) - 1

                # time to make tmask for chunk
                beg_time = nc.num2date(time_in[beg], units=t_unit, calendar=t_cal)
                if invariant:
                    # allow topography to work in same code, len(nctime) = 1
                    end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                    # end = 1
                else:
                    end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)

                # '<= end_time', would damage appending
                tmask_chunk = (time <= end_time) * (time >= beg_time)
                if invariant:
                    # allow topography to work in same code
                    tmask_chunk = np.array([True])

                # get the interpolated variables
                dfield, variables = self.interp2D(ncf_in,
                                                self.stations, tmask_chunk,
                                                subset_grid, lon_slice, lat_slice,
                                                variables=None, date=None)
                
                # append variables
                for i, var in enumerate(variables):
                    if variables_skip(var):
                        continue

                    if ens:
                        num = ncf_in.variables['number'][:]
                        for ni in num:
                            if pl:
                                # dfield [station, variables, time, levels, number]
                                vi = dfield[ni].data[:,i,:,:].transpose((1,2,0))
                                ncf_out.variables[var][beg:end+1,ni,:,:] = vi
                            else:
                                vi = dfield[ni].data[:,i,:].transpose((1,0))
                                ncf_out.variables[var][beg:end+1,ni,:] = vi
                    else:
                        if pl:
                            vi = dfield.data[:,i,:,:].transpose((1,2,0))
                            ncf_out.variables[var][beg:end+1,:,:] = vi
                        else:
                            vi = dfield.data[:,i,:].transpose((1,0))
                            ncf_out.variables[var][beg:end+1,:] = vi

                # delete objects and free memory
                del dfield, vi, tmask_chunk
                gc.collect()

                ncf_out.globsim_last_chunk_written = n

            # Write success flag
            ncf_out.globsim_interpolate_success = 1
        ncf_in.close()

    def levels2elevation(self, ncfile_in, ncfile_out):
        """
        ncfile_in : nc.Dataset
            2-d Interpolated ERA5 file on pressure-levels
        ncfile_out : str
            output file name
        Linear 1D interpolation of pressure level data available for individual
        stations to station elevation. Where and when stations are below the
        lowest pressure level, they are assigned the value of the lowest
        pressure level.

        """
        # open file
        ncf = nc.Dataset(ncfile_in)
        height = ncf.variables['height'][:]

        nt = len(ncf.variables['time'][:])  # these are reading from a globsim file, not an ERA5 file
        nl = len(ncf.variables['level'][:])

        # list variables
        varlist = [key for key in ncf.variables.keys()]
        for V in ['time', 'station', 'latitude', 'longitude',
                  'level','height','z','expver']:
            if V in varlist:
                varlist.remove(V)
        if self.ens:
            varlist.remove('number')

        # === open and prepare output netCDF file ==============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        if not (self.resume and Path(ncfile_out).exists()):
            rootgrp = netcdf_base(ncfile_out, len(height), nt,
                                time_units=ncf['time'].units,
                                nc_in=ncf,
                                calendar=ncf['time'].calendar)
            if self.ens:
                rootgrp.source = 'ERA5 10-member ensemble, interpolated (bi)linearly to stations'
            else:
                rootgrp.source = 'ERA5, interpolated (bi)linearly to stations'

            time = rootgrp['time']
            station = rootgrp['station']
            latitude = rootgrp['latitude']
            longitude = rootgrp['longitude']
            height = rootgrp['height']

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
                if self.ens:
                    tmp = rootgrp.createVariable(var,
                                                'f4',('time','number','station'))
                else:
                    tmp = rootgrp.createVariable(var,'f4',('time', 'station'))

                tmp.long_name = ncf[var].long_name
                tmp.units     = ncf[var].units

            # add air pressure as new variable
            var = 'air_pressure'
            varlist.append(var)
            if self.ens:
                tmp = rootgrp.createVariable(var,'f4',('time','number','station'))
            else:
                tmp = rootgrp.createVariable(var,'f4',('time','station'))
            tmp.long_name = var.encode('UTF8')
            tmp.units = 'hPa'.encode('UTF8')

            rootgrp.close()
            # end file prepation ===================================================

        with nc.Dataset(ncfile_out, 'a') as rootgrp:
            # loop over stations
            for n, h in enumerate(height):
                if self.resume and n <= rootgrp.last_station_written:
                        if n == rootgrp.last_station_written:
                            logger.info(f"Resuming interpolation at station {n+1}")
                        continue

                self.require_safe_mem_usage()
                logger.debug(f"Interpolating station {n+1} to station elevation using pressure-level data")

                if self.ens:
                    num = ncf.variables['number'][:]
                    for ni in num:
                        elevation = ncf.variables['z'][:,ni,:,n] / const.G
                        elev_diff, va, vb = ele_interpolate(elevation, h, nl)
                        wa, wb = calculate_weights(elev_diff, va, vb)
                        for v, var in enumerate(varlist):
                            if var == 'air_pressure':
                                # pressure [Pa] variable from levels, shape: (time, level)
                                data = np.repeat([ncf.variables['level'][:]],
                                                len(time),axis=0).ravel()
                            else:
                                # read data from netCDF
                                data = ncf.variables[var][:,ni,:,n].ravel()

                            ipol = data[va] * wa + data[vb] * wb   # interpolated value
                            rootgrp.variables[var][:,ni,n] = ipol  # assign to file
                else:
                    # convert geopotential [mbar] to height [m], shape: (time, level)
                    elevation = ncf.variables['z'][:,:,n] / const.G
                    elev_diff, va, vb = ele_interpolate(elevation, h, nl)
                    wa, wb = calculate_weights(elev_diff, va, vb)

                    # loop over variables and apply interpolation weights
                    for v, var in enumerate(varlist):
                        if var in str(rootgrp.vars_written).split(" "):
                            logger.debug(f"Skipping {var}")
                            continue
                        
                        if var == 'air_pressure':
                            # pressure [Pa] variable from levels, shape: (time, level)
                            data = np.repeat([ncf.variables['level'][:]],
                                            len(time),axis=0).ravel()
                        else:
                            # read data from netCDF
                            logger.debug(f"Reading {var}")
                            data = ncf.variables[var][:,:,n].ravel()

                        ipol = data[va] * wa + data[vb] * wb   # interpolated value
                        rootgrp.variables[var][:,n] = ipol  # write to file
                    
                        rootgrp.vars_written = " ".join(set(str(rootgrp.vars_written).split(" ") + [var]))
                
                rootgrp.vars_written = ""
                rootgrp.last_station_written = n
                
            
            rootgrp.globsim_interpolate_success = 1

        ncf.close()
        # closed file ==========================================================

    def _preprocess(self):
        # 2D Interpolation for Invariant Data
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.
        self.ERA2station(self.mf_to, self.getOutFile('to'),
                         self.stations, ['z', 'lsm'], date=None)
        
    def _process_sa(self):
        # === 2D Interpolation for Surface Analysis Data ===
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['t2m'],  # [K] 2m values
                'relative_humidity' : ['d2m'],  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky':
                    ['tco3', 'tcwv'],   # [kg m-2] Total column ozone
                                        # [kg m-2] Total column W vapor
                'wind_speed': ['u10', 'v10']}   # [m s-1] 10m values
        varlist = self.TranslateCF2short(dpar)
        
        if self.resume and self.completed_successfully(self.getOutFile('sa')):
            logger.info("Skipping surface analysis interpolation")
        else:
            with xr.open_mfdataset(self.get_input_file_paths('sa'), decode_times=False) as sa:
                self.ERA2station(sa, self.getOutFile('sa'),
                                 self.stations, varlist, date=self.date)
    
    def _process_pl(self):
        # === 2D Interpolation for Pressure Level Data ===
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['t'],           # [K]
                'relative_humidity' : ['r'],           # [%]
                'wind_speed'        : ['u', 'v']}      # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('z')
        if self.resume and self.completed_successfully(self.getOutFile('pl')):
            logger.info("Skipping pressure level interpolation")
        else:
            with xr.open_mfdataset(self.get_input_file_paths('pl'), decode_times=False) as pl:
                self.ERA2station(pl, self.getOutFile('pl'),
                                self.stations, varlist, date=self.date)

        # 1D Interpolation for Pressure Level Data
        if self.ens:
            outf = 'era5_ens_pl_'
        else:
            outf = 'era5_pl_'
        outf = path.join(self.output_dir, outf + self.list_name + '_surface.nc')
        
        if self.resume and self.completed_successfully(outf):
            logger.info("Skipping pl surface interpolation")
        else:
            self.levels2elevation(self.getOutFile('pl'), outf)
    
    def _process_sf(self):
        # 2D Interpolation for Surface Forecast Data    'tp', 'strd', 'ssrd'
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.
        # [m] total precipitation
        # [J m-2] short-wave downward
        # [J m-2] long-wave downward
        dpar = {'precipitation_amount'              : ['tp'],
                'downwelling_shortwave_flux_in_air' : ['ssrd'],
                'downwelling_longwave_flux_in_air'  : ['strd']}
        varlist = self.TranslateCF2short(dpar)
        if self.resume and self.completed_successfully(self.getOutFile('sa')):
            logger.info("Skipping surface forecast interpolation")
        else:
            with xr.open_mfdataset(self.get_input_file_paths('sf'), decode_times=False) as sf:
                self.ERA2station(sf, self.getOutFile('sf'),
                                self.stations, varlist, date=self.date)


# vv = nc.Dataset(self.getOutFile('pl'))
# xx = xr.open_dataset(self.getOutFile('sa'))