import netCDF4 as nc
import gc
import xarray as xr
import logging
import numpy as np
import urllib3

from os import path
from pathlib import Path

from globsim.interp import calculate_weights, ele_interpolate, extrapolate_below_grid
from globsim.common_utils import variables_skip
from globsim.interpolate.GenericInterpolate import GenericInterpolate
from globsim.nc_elements import new_interpolated_netcdf
import globsim.constants as const

logger = logging.getLogger('globsim.interpolate')


urllib3.disable_warnings()


class ERA5interpolate(GenericInterpolate):
    """
    Collection of methods to interpolate ERA5 netCDF files to station
    coordinates. All variables retain their original units and time stepping.
    """
    REANALYSIS = 'era5'
    PL_SKIP_VARS = {'time', 'station', 'latitude', 'longitude',
                    'level', 'height', 'z', 'expver', 'number',
                    'station_name'}

    def __init__(self, ifile,  **kwargs):
        super().__init__(ifile=ifile, **kwargs)
        self._set_input_directory("era5")
        self.ens = False

        # convert longitude to ERA notation if using negative numbers
        self.stations['longitude_dd'] = self.stations['longitude_dd'] % 360

        try:
            self.mf_to = xr.open_mfdataset(self.get_input_file_glob('to'), decode_times=False)
        except OSError as e:
            logger.error("No files found for invariant (era_to) data")
        
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

        f_glob = 'era5_{}_*.nc'

        if levStr == 'to':
            infile = path.join(self.input_dir, 'era5_to.nc')
        else:
            infile = path.join(self.input_dir, f_glob.format(levStr))

        return infile

    def get_input_file_paths(self, levStr):
        glob = self.get_input_file_glob(levStr)
        infiles = self.prefilter_mf_paths(glob)
        return infiles
    
    def get_output_file(self, levStr):
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
                                              n_time = len(time_in),
                                              station_names=self.stations['station_name'])
            
            rootgrp.globsim_interpolate_start = self.par['beg']
            rootgrp.globsim_interpolate_end = self.par['end']
            rootgrp.globsim_chunk_size = self.cs
            rootgrp.globsim_interpolate_success = 0
            rootgrp.globsim_last_chunk_written = -1
            rootgrp.source = f'{self.REANALYSIS}, interpolated bilinearly to stations'

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

    def get_elevation(self, nc_pl_interp, station_index):
        # geopotential → meters
        return nc_pl_interp.variables['z'][:, :, station_index] / const.G
    
    def _preprocess(self):
        # 2D Interpolation for Invariant Data
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.
        if self._skip_invariant or (self.resume and self.completed_successfully(self.get_output_file('to'))):
            logger.info("Skipping invariant interpolation")
        else:
            self.ERA2station(self.mf_to, self.get_output_file('to'),
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
        
        if self.resume and self.completed_successfully(self.get_output_file('sa')):
            logger.info("Skipping surface analysis interpolation")
        else:
            with xr.open_mfdataset(self.get_input_file_paths('sa'), decode_times=False) as sa:
                self.ERA2station(sa, self.get_output_file('sa'),
                                 self.stations, varlist, date=self.date)
    
    def _process_pl(self):
        # === 2D Interpolation for Pressure Level Data ===
        # dictionary to translate CF Standard Names into ERA5
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['t'],           # [K]
                'relative_humidity' : ['r'],           # [%]
                'wind_speed'        : ['u', 'v']}      # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('z')
        if self.resume and self.completed_successfully(self.get_output_file('pl')):
            logger.info("Skipping pressure level interpolation")
        else:
            with xr.open_mfdataset(self.get_input_file_paths('pl'), decode_times=False) as pl:
                self.ERA2station(pl, self.get_output_file('pl'),
                                self.stations, varlist, date=self.date)
    
    def _process_pl_sur(self):
        # 1D Interpolation for Pressure Level Data
        outf = self.get_output_file('pl')[:-3] + "_surface.nc"
  
        if self.resume and self.completed_successfully(outf):
            logger.info("Skipping pl surface interpolation")
        else:
            self.levels2elevation(self.get_output_file('pl'), outf)
    
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
        if self.resume and self.completed_successfully(self.get_output_file('sf')):
            logger.info("Skipping surface forecast interpolation")
        else:
            with xr.open_mfdataset(self.get_input_file_paths('sf'), decode_times=False) as sf:
                self.ERA2station(sf, self.get_output_file('sf'),
                                self.stations, varlist, date=self.date)


class ERA5EnsembleInterpolate(ERA5interpolate):

    REANALYSIS = 'era5ens'

    def __init__(self, ifile, **kwargs):
        super().__init__(ifile=ifile, **kwargs)

        self._set_input_directory("era5ens")
        self.ens = True

    def get_input_file_glob(self, levStr):
        """ Return the input file(s) for the interpolation. """
        # edited naming conventions for simplicity and to avoid errors
        base = super(ERA5EnsembleInterpolate, self).get_input_file_glob(levStr)
        return base.replace('era5', 'era5_ens')
    
    def get_output_file(self, levStr):
        nome = 'era5_ens_{}_'.format(levStr) + self.list_name + '.nc'
        outfile = path.join(self.output_dir, nome)
        return outfile

    def write_dfield_to_file(self, 
                            dfield, 
                            variables:list,
                            ncf_out,
                            beg:int,
                            end:int,
                            pl:bool):

        for i, var in enumerate(variables):
            if variables_skip(var):
                continue
                        
            else:
                for j, _dfield in enumerate(dfield):
                    if pl:
                        # dfield [station, variables, time, levels, number]
                        vi = _dfield.data[:,i,:,:].transpose((1,2,0))
                        ncf_out.variables[var][beg:end+1,j,:,:] = vi
                    else:
                        vi = _dfield.data[:,i,:].transpose((1,0))
                        ncf_out.variables[var][beg:end+1,j,:] = vi

    def create_source_field(self,
                            sgrid, 
                            variables, 
                            nt,
                            ncf_in, 
                            pl:bool):
            
        num = ncf_in.variables['number'][:]
        sfield = []
        for ni in num:
            f = super(ERA5EnsembleInterpolate, self).create_source_field(sgrid, variables, nt, ncf_in, pl)
            sfield.append(f)

        return sfield
    
    @staticmethod
    def nc_data_subset_to_source_field(variables, sfield_list: list, ncf_in,
                                         tmask_chunk, pl: bool, lon_subset, lat_subset):
        for n, var in enumerate(variables):
            for ni in ncf_in['number'][:]:
                if pl:
                    # sfield.data [lon, lat, var, time, number, level]
                    # vi [time, number, level, lat, lon]
                    vi = ncf_in[var][tmask_chunk,ni,:,:,:]
                    sfield_list[ni].data[:,:,n,:,:] = vi.transpose((3,2,0,1))

                else:

                    vi = ncf_in[var][tmask_chunk,ni,:,:]
                    sfield_list[ni].data[:,:,n,:] = vi.transpose((2,1,0))

            logger.debug(f"Wrote {var} data to source field for regridding")

    def create_pl_output_variables(self, rootgrp, ncf, varlist):
        """
        Same as base, but output dimensions are (time, number, station)
        instead of (time, station).
        """
        # create 'number' dimension if it doesn't exist
        if 'number' not in rootgrp.dimensions:
            num = ncf.variables['number'][:]
            rootgrp.createDimension('number', len(num))
            rootgrp.createVariable('number', 'i4', ('number',))
            rootgrp['number'][:] = num

        for var in varlist:
            tmp = rootgrp.createVariable(var, 'f4',
                                         ('time', 'number', 'station'))
            for attr in ['long_name', 'units']:
                if attr in ncf.variables[var].ncattrs():
                    tmp.setncattr(attr, ncf.variables[var].getncattr(attr))

        var = 'air_pressure'
        varlist.append(var)
        tmp = rootgrp.createVariable(var, 'f4',
                                     ('time', 'number', 'station'))
        tmp.long_name = 'Air pressure'
        tmp.units = 'hPa'

    def get_elevation(self, ncf, station_index):
        """
        For ensemble data, returns a dict keyed by ensemble member number,
        each value being a (time, level) elevation array in meters.
        """
        import globsim.constants as const
        num = ncf.variables['number'][:]
        return {ni: ncf.variables['z'][:, ni, :, station_index] / const.G
                for ni in num}

    def interpolate_and_write_station(self, ncf, rootgrp, n, h,
                                       elevation_dict, varlist, nl, time):
        """
        Loop over ensemble members, interpolate each independently,
        and write to the (time, number, station) output.
        """
        for ni, elevation in elevation_dict.items():
            elev_diff, va, vb = ele_interpolate(elevation, h, nl)
            wa, wb = calculate_weights(elev_diff, va, vb)

            for var in varlist:
                if var == 'air_pressure':
                    data = np.repeat([ncf.variables['level'][:]],
                                     len(time), axis=0).ravel()
                else:
                    data = ncf.variables[var][:, ni, :, n].ravel()

                ipol = data[va] * wa + data[vb] * wb

                if self.extrapolate_below_grid:
                    extrapolated_values = extrapolate_below_grid(
                        elevation, data, h)
                    ipol = np.where(~extrapolated_values.mask,
                                    extrapolated_values, ipol)

                rootgrp.variables[var][:, ni, n] = ipol

        rootgrp.vars_written = ""
        rootgrp.last_station_written = n
    
    