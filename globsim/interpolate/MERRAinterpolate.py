from __future__        import print_function

import netCDF4 as nc
import numpy as np
import sys
import logging

from datetime          import datetime
from os                import path, makedirs

from globsim.common_utils import str_encode, variables_skip
from globsim.interpolate.GenericInterpolate import GenericInterpolate
from globsim.nc_elements import netcdf_base

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module='netCDF4')

logger = logging.getLogger('globsim.interpolate')

class MERRAinterpolate(GenericInterpolate):
    """
    Algorithms to interpolate MERRA-2 netCDF files to station coordinates.
    All variables retains their original units and time-steps.

    Referenced from era_interim.py (Dr.Stephan Gruber): Class ERAinterpolate()
    """

    def __init__(self, ifile):
        super().__init__(ifile)
        par = self.par

        self.input_dir = path.join(par['project_directory'],'merra2')

    def netCDF_empty(self, ncfile_out, stations, nc_in):
        # TODO: change date type from f4 to f8 for lat and lon
        '''
        Creates an empty station file to hold interpolated reults. The number of
        stations is defined by the variable stations, variables are determined by
        the variable list passed from the gridded original netCDF.

        ncfile_out: full name of the file to be created
        stations:   station list read with common_utils.StationListRead()
        variables:  variables read from netCDF handle
        lev:        list of pressure levels, empty is [] (default)
        '''
        rootgrp = netcdf_base(ncfile_out, len(stations), None,
                              'hours since 1980-01-01 00:00:00')

        station = rootgrp["station"]
        latitude = rootgrp["latitude"]
        longitude = rootgrp["longitude"]
        height = rootgrp["height"]

        # assign station characteristics
        station[:]   = list(stations['station_number'])
        latitude[:]  = list(stations['latitude_dd'])
        longitude[:] = list(stations['longitude_dd'])
        height[:]    = list(stations['elevation_m'])

        # extra treatment for pressure level files
        try:
            lev = nc_in.variables['level'][:]
            print("== 3D: file has pressure levels")
            level           = rootgrp.createDimension('level', len(lev))
            level           = rootgrp.createVariable('level','i4',('level'))
            level.long_name = 'pressure_level'
            level.units     = 'hPa'
            level[:] = lev
        except Exception:
            print("== 2D: file without pressure levels")
            lev = []

        # remove extra variables
        varlist_merra = [str_encode(x) for x in nc_in.variables.keys()]

        # create and assign variables based on input file
        for n, var in enumerate(varlist_merra):
            if variables_skip(var):
                continue
            print("VAR: ", str_encode(var))
            # extra treatment for pressure level files
            if len(lev):
                tmp = rootgrp.createVariable(var,'f4', ('time', 'level', 'station'))
            else:
                tmp = rootgrp.createVariable(var,'f4', ('time', 'station'))
            tmp.long_name = str_encode(nc_in.variables[var].long_name)  # for merra2
            tmp.units     = str_encode(nc_in.variables[var].units)

        # close the file
        rootgrp.close()

    def MERRA2station(self, ncfile_in, ncfile_out, points,
                      variables=None, date=None):
        """
        Given the type of variables to interpoalted from MERRA2 downloaded diretory
        Create the empty of netCDF file to hold the interpolated results, by calling
        self.netCDF_empty
        Get the interpolated results from MERRA2station
        Append all variables into the empty netCDF file
        Close all files

        Args:
        ncfile_in: Full path to an MERRA-2 derived netCDF file. This can
                    contain wildcards to point to multiple files if temporal
                    chunking was used.

        ncfile_out: Full path to the output netCDF file to write.

        points: A dataframe of locations. See method StationListRead in
                common_utils.py for more details.

        variables:  List of variable(s) to interpolate such as
                    ['T','RH','U','V',' T2M', 'U2M', 'V2M', 'U10M', 'V10M',
                    'PRECTOT', 'SWGDN','SWGDNCLR','LWGDN', 'LWGDNCLR'].
                    Defaults to using all variables available.

        date: Directory to specify begin and end time for the derived time
                series. Defaluts to using all times available in ncfile_in.

        """

        # read in one type of mutiple netcdf files
        ncf_in = nc.MFDataset(ncfile_in, 'r', aggdim='time')

        # is it a file with pressure levels?
        pl = 'level' in ncf_in.dimensions.keys()

        # build the output of empty netCDF file
        self.netCDF_empty(ncfile_out, self.stations, ncf_in)

        # open the output netCDF file, set it to be appendable ('a')
        ncf_out = nc.Dataset(ncfile_out, 'a')

        # get time and convert to datetime object
        nctime = ncf_in.variables['time'][:]
        # "hours since 1980-01-01 00:00:00"
        t_unit = "hours since 1980-01-01 00:00:00"  # ncf_in.variables['time'].units
        try :
            t_cal = ncf_in.variables['time'].calendar
        except AttributeError :  # Attribute doesn't exist
            t_cal = u"gregorian"  # or standard
        time = [nc.num2date(timei, units=t_unit, calendar=t_cal) for timei in nctime]
        time = np.asarray(time)

        # detect invariant files (topography etc.)
        invariant = True if len(time) == 1 else False

        # restrict to date/time range if given
        if date is None:
            tmask = time < datetime(3000, 1, 1)
        else:
            tmask = (time < date['end']) * (time >= date['beg'])

        if not any(tmask):
            sys.exit('''\n ERROR: No downloaded data exist within date range specified by interpolation control file.
                     Download new data or change 'beg' / 'end' in interpolation control file''')

        # get time indices
        time_in = nctime[tmask]

        # ensure that chunk sizes cover entire period even if
        # len(time_in) is not an integer multiple of cs
        niter  = len(time_in) // self.cs
        niter += ((len(time_in) % self.cs) > 0)

        # loop in chunk size cs
        for n in range(niter):
            # indices
            beg = n * self.cs
            # restrict last chunk to lenght of tmask plus one (to get last time)
            end = min(n * self.cs + self.cs, len(time_in)) - 1

            # time to make tmask for chunk
            beg_time = nc.num2date(time_in[beg], units=t_unit, calendar=t_cal)
            if invariant:
                # allow topography to work in same code, len(nctime) = 1
                end_time = nc.num2date(nctime[0], units=t_unit, calendar=t_cal)
                # end = 1
            else:
                end_time = nc.num2date(time_in[end], units=t_unit, calendar=t_cal)

            # !! CAN'T HAVE '<= end_time', would damage appeding
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
                                                     time_in[beg:end + 1])
            # append variables
            for i, var in enumerate(variables):
                if variables_skip(var):
                    continue

                # extra treatment for pressure level files
                if pl:
                    lev = ncf_in.variables['level'][:]
                    ncf_out.variables[var][beg:end + 1, :, :] = dfield.data[:, i, :, :].transpose(1,2,0)
                else:
                    ncf_out.variables[var][beg:end + 1, :] = dfield.data[:, i, :].transpose(1,0)

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
        for V in ['time', 'station', 'latitude', 'longitude', 'level', 'height', 'H']:
            varlist.remove(V)

        # === open and prepare output netCDF file =============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = netcdf_base(ncfile_out, len(height), nt,
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
        # end file prepation ==================================================

        # loop over stations
        for n, h in enumerate(height):
            # geopotential unit: height [m]
            # shape: (time, level)
            elevation = ncf.variables['H'][:,:,n]
            # TODO: check if height of stations in data range (+50m at top,
            # lapse r.)

            # difference in elevation.
            # level directly above will be >= 0
            elev_diff = -(elevation - h)
            # vector of level indices that fall directly above station.
            # Apply after ravel() of data.
            va = np.argmin(elev_diff + (elev_diff < 0) * 100000, axis=1)
            # mask for situations where station is below lowest level
            mask = va < (nl - 1)
            va += np.arange(elevation.shape[0]) * elevation.shape[1]

            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask  # +1 when OK, +0 when below lowest level

            wa, wb = self.calculate_weights(elev_diff, va, vb)

            # loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                if var == 'air_pressure':
                    # pressure [hPa] variable from levels, shape: (time, level)
                    data = np.repeat([ncf.variables['level'][:]],
                                     len(time),axis=0).ravel()
                    ipol = data[va] * wa + data[vb] * wb   # interpolated value

                    # if mask[pixel] == false, pass the maximum of pressure level to pixles
                    level_highest = ncf.variables['level'][:][-1]
                    level_lowest = ncf.variables['level'][:][0]

                    for j, value in enumerate(ipol):
                        if value == level_highest:
                            ipol[j] = level_lowest

                else:
                    # read data from netCDF
                    data = ncf.variables[var][:,:,n].ravel()
                    ipol = data[va] * wa + data[vb] * wb   # interpolated value
                rootgrp.variables[var][:,n] = ipol  # assign to file

        rootgrp.close()
        ncf.close()
        # closed file ==========================================================

    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to
        the more generically MERRA-like interpolation functions.
        """

        # 2D Interpolation for Constant Model Parameters
        # dictionary to translate CF Standard Names into MERRA
        # pressure level variable keys.
        dummy_date = {'beg' : datetime(1992, 1, 2, 3, 0),
                      'end' : datetime(1992, 1, 2, 4, 0)}

        if not path.isdir(self.output_dir):
            makedirs(self.output_dir)

        # === 2D Interpolation for Surface Analysis Data ===
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['T2M'],  # [K] 2m values
                'wind_speed' : ['U2M', 'V2M', 'U10M','V10M'],   # [m s-1] 2m & 10m values
                'relative_humidity' : ['QV2M']}  # 2m value
        varlist = self.TranslateCF2short(dpar)
        self.MERRA2station(path.join(self.input_dir,'merra_sa_*.nc'),
                           path.join(self.output_dir,'merra2_sa_' + self.list_name + '.nc'),
                           self.stations, varlist, date=self.date)

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
        self.MERRA2station(path.join(self.input_dir,'merra_sf_*.nc'),
                           path.join(self.output_dir,'merra2_sf_' + self.list_name + '.nc'),
                           self.stations, varlist, date=self.date)

        # NEED ADD 'H' in it!
        # === 2D Interpolation for Pressure-Level, Analyzed Meteorological DATA ===
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys.
        dpar = {'air_temperature'   : ['T'],           # [K]
                'wind_speed'        : ['U', 'V'],      # [m s-1]
                'relative_humidity' : ['RH']}          # [1]
        varlist = self.TranslateCF2short(dpar).append('H')
        self.MERRA2station(path.join(self.input_dir,'merra_pl_*.nc'),
                           path.join(self.output_dir,'merra2_pl_' + self.list_name + '.nc'),
                           self.stations, varlist, date=self.date)

        # 1D Interpolation for Pressure Level Analyzed Meteorological Data
        self.levels2elevation(path.join(self.output_dir,'merra2_pl_' + self.list_name + '.nc'),
                              path.join(self.output_dir,'merra2_pl_' + self.list_name + '_surface.nc'))
