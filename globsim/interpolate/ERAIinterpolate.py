from __future__ import print_function

import re
import netCDF4 as nc
import numpy as np

from datetime import datetime
from os import path

from globsim.common_utils import variables_skip, str_encode
from globsim.interpolate.GenericInterpolate import GenericInterpolate
from globsim.nc_elements import netcdf_base, new_interpolated_netcdf

try:
    import ESMF
    # Check ESMF version.  7.0.1 behaves differently than 7.1.0r
    ESMFv = int(re.sub("[^0-9]", "", ESMF.__version__))
    ESMFnew = ESMFv > 701

except ImportError:
    print("*** ESMF not imported, interpolation not possible. ***")
    pass


class ERAIinterpolate(GenericInterpolate):
    """
    Collection of methods to interpolate ERA-Interim netCDF files to station
    coordinates. All variables retain theit original units and time stepping.
    """

    def __init__(self, ifile):
        super().__init__(ifile)
        par = self.par
        self.input_dir = path.join(par['project_directory'], 'erai')

        #convert longitude to ERA notation if using negative numbers
        self.stations['longitude_dd'] = self.stations['longitude_dd'] % 360

    def ERA2station(self, ncfile_in, ncfile_out, points,
                    variables = None, date = None):

        """
        Biliner interpolation from fields on regular grid (latitude, longitude)
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all Era_Interim files). The type
        of variable and file structure are determined from the input.

        This function creates an empty of netCDF file to hold the interpolated
        results, by calling ERAIgeneric().netCDF_empty. Then, data is
        interpolated in temporal chunks and appended. The temporal chunking can
        be set in the interpolation parameter file.

        Args:
        ncfile_in: Full path to an Era-Interim derived netCDF file. This can
                   contain wildcards to point to multiple files if temporal
                  chunking was used.

        ncfile_out: Full path to the output netCDF file to write.

        points: A dataframe of locations. See method StationListRead in
                common_utils.py for more details.

        variables:  List of variable(s) to interpolate such as
                    ['r', 't', 'u','v', 't2m', 'u10', 'v10', 'ssrd', 'strd', 'tp'].
                    Defaults to using all variables available.

        date: Directory to specify begin and end time for the derived time
                series. Defaluts to using all times available in ncfile_in.

        cs: chunk size, i.e. how many time steps to interpolate at once. This
            helps to manage overall memory usage (small cs is slower but less
            memory intense).
        """

        # read in one type of mutiple netcdf files
        ncf_in = nc.MFDataset(ncfile_in, 'r', aggdim ='time')

        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # build the output of empty netCDF file
        rootgrp = new_interpolated_netcdf(ncfile_out, self.stations, ncf_in,
                                          time_units= 'hours since 1900-01-01 00:00:0.0')
        rootgrp.source = 'ERA_Interim, interpolated bilinearly to stations'
        rootgrp.close()

        # open the output netCDF file, set it to be appendable ('a')
        ncf_out = nc.Dataset(ncfile_out, 'a')

        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]
        #"hours since 1900-01-01 00:00:0.0"
        t_unit = ncf_in.variables['time'].units
        try :
            t_cal = ncf_in.variables['time'].calendar
        except AttributeError :  # attribute doesn't exist
            t_cal = u"gregorian" # standard
        time = nc.num2date(nctime, units = t_unit, calendar = t_cal)

        # detect invariant files (topography etc.)
        invariant = True if len(time) == 1 else False

        # restrict to date/time range if given
        if date is None:
            tmask = time < datetime(3000, 1, 1)
        else:
            tmask = (time < date['end']) * (time >= date['beg'])

        # get time vector for output
        time_in = nctime[tmask]

        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter  = len(time_in) // self.cs
        niter += ((len(time_in) % self.cs) > 0)

        # loop over chunks
        for n in range(niter):
            # indices (relative to index of the output file)
            beg = n * self.cs
            # restrict last chunk to lenght of tmask plus one (to get last time)
            end = min(n*self.cs + self.cs, len(time_in))-1

            # time to make tmask for chunk
            beg_time = nc.num2date(time_in[beg], units=t_unit, calendar=t_cal)
            if invariant:
                # allow topography to work in same code, len(nctime) = 1
                end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                #end = 1
            else:
                end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)

            #'<= end_time', would damage appending
            tmask_chunk = (time <= end_time) * (time >= beg_time)
            if invariant:
                # allow topography to work in same code
                tmask_chunk = [True]

            # get the interpolated variables
            dfield, variables = self.interp2D(ncfile_in, ncf_in,
                                                     self.stations, tmask_chunk,
                                                     variables=None, date=None)

            # append time
            ncf_out.variables['time'][:] = np.append(ncf_out.variables['time'][:],
                                                     time_in[beg:end+1])

            #append variables
            for i, var in enumerate(variables):
                if variables_skip(var):
                    continue

                if pl:
                    if ESMFnew:
                        # dfield has dimensions (station, variables, time, pressure levels)
                        ncf_out.variables[var][beg:end+1,:,:] = dfield.data[:,i,:,:].transpose((1,2,0))
                    else:
                        # dimension: time, level, station (pressure level files)
                        ncf_out.variables[var][beg:end+1,:,:] = dfield.data[i,:,:,:]
                else:
                    if ESMFnew:
                        # dfield has dimensions (station, variables, time)
                        ncf_out.variables[var][beg:end+1,:] = dfield.data[:,i,:].transpose((1,0))
                    else:
                        # dfield has dimensions time, station (2D files)
                        ncf_out.variables[var][beg:end+1,:] = dfield.data[i,:,:]

        #close the file
        ncf_in.close()
        ncf_out.close()


    def levels2elevation(self, ncfile_in, ncfile_out):
        """
        Linear 1D interpolation of pressure level data available for individual
        stations to station elevation. Where and when stations are below the
        lowest pressure level, they are assigned the value of the lowest
        pressure level.

        """
        # open file

        ncf = nc.MFDataset(ncfile_in, 'r', aggdim='time')
        height = ncf.variables['height'][:]
        nt = len(ncf.variables['time'][:])
        nl = len(ncf.variables['level'][:])

        # list variables
        varlist = [str_encode(x) for x in ncf.variables.keys()]
        for V in ['time', 'station', 'latitude', 'longitude', 'level','height','z']:
            varlist.remove(V)

        # === open and prepare output netCDF file ==============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = netcdf_base(ncfile_out=ncfile_out, n_stations=len(height), n_time=nt, time_units='hours since 1900-01-01 00:00:0.0')
        rootgrp.source = 'ERA-Interim, interpolated (bi)linearly to stations'

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

        # create and assign variables from input file
        for var in varlist:
            tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))
            tmp.long_name = str_encode(ncf.variables[var].long_name)
            tmp.units     = str_encode(ncf.variables[var].units)

        # add air pressure as new variable
        var = 'air_pressure'
        varlist.append(var)
        tmp   = rootgrp.createVariable(var,'f4',('time', 'station'))
        tmp.long_name = var.encode('UTF8')
        tmp.units     = 'hPa'.encode('UTF8')
        # end file prepation ===================================================

        # loop over stations
        for n, h in enumerate(height):
            # convert geopotential [mbar] to height [m], shape: (time, level)
            elevation = ncf.variables['z'][:,:,n] / 9.80665
            # TODO: check if height of stations in data range

            # difference in elevation, level directly above will be >= 0
            elev_diff = elevation - h
            # vector of level indices that fall directly above station.
            # Apply after ravel() of data.
            va = np.argmin(elev_diff + (elev_diff < 0) * 100000, axis=1)
            # mask for situations where station is below lowest level
            mask = va < (nl-1)
            va += np.arange(elevation.shape[0]) * elevation.shape[1]

            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask # +1 when OK, +0 when below lowest level

            wa, wb = self.calculate_weights(elev_diff, va, vb)

            #loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                if var == 'air_pressure':
                    # pressure [Pa] variable from levels, shape: (time, level)
                    data = np.repeat([ncf.variables['level'][:]],
                                      len(time),axis=0).ravel()
                else:
                    #read data from netCDF
                    data = ncf.variables[var][:,:,n].ravel()
                ipol = data[va]*wa + data[vb]*wb   # interpolated value
                rootgrp.variables[var][:,n] = ipol # assign to file

        rootgrp.close()
        ncf.close()
        # closed file ==========================================================


    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to
        the more generically ERA-like interpolation functions.
        """

        # 2D Interpolation for Invariant Data
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys.
        dummy_date  = {'beg' : datetime(1979, 1, 1, 12, 0),
                       'end' : datetime(1979, 1, 1, 12, 0)}
        self.ERA2station(path.join(self.input_dir,'erai_to.nc'),
                         path.join(self.output_dir,'erai_to_' +
                                   self.list_name + '.nc'), self.stations,
                                   ['z', 'lsm'], date = dummy_date)

        # === 2D Interpolation for Surface Analysis Data ===
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['t2m'],  # [K] 2m values
                'relative_humidity' : ['d2m'],  # [K] 2m values
                'downwelling_shortwave_flux_in_air_assuming_clear_sky' :
                    ['tco3', 'tcwv'],   # [kg m-2] Total column ozone
                                        # [kg m-2] Total column W vapor
                'wind_speed' : ['u10', 'v10']}   # [m s-1] 10m values
        varlist = self.TranslateCF2short(dpar)
        self.ERA2station(path.join(self.input_dir,'erai_sa_*.nc'),
                         path.join(self.output_dir,'erai_sa_' +
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)

        # 2D Interpolation for Surface Forecast Data    'tp', 'strd', 'ssrd'
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys.
        # [m] total precipitation
        # [J m-2] short-wave downward
        # [J m-2] long-wave downward
        dpar = {'precipitation_amount'              : ['tp'],
                'downwelling_shortwave_flux_in_air' : ['ssrd'],
                'downwelling_longwave_flux_in_air'  : ['strd']}
        varlist = self.TranslateCF2short(dpar)
        self.ERA2station(path.join(self.input_dir,'erai_sf_*.nc'),
                         path.join(self.output_dir,'erai_sf_' +
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)


        # === 2D Interpolation for Pressure Level Data ===
        # dictionary to translate CF Standard Names into ERA-Interim
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['t'],           # [K]
                'relative_humidity' : ['r'],           # [%]
                'wind_speed'        : ['u', 'v']}      # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('z')
        self.ERA2station(path.join(self.input_dir,'erai_pl_*.nc'),
                         path.join(self.output_dir,'erai_pl_' +
                                   self.list_name + '.nc'), self.stations,
                                   varlist, date = self.date)

        # 1D Interpolation for Pressure Level Data
        self.levels2elevation(path.join(self.output_dir,'erai_pl_' +
                                        self.list_name + '.nc'),
                              path.join(self.output_dir,'erai_pl_' +
                                        self.list_name + '_surface.nc'))

