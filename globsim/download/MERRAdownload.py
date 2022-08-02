
from __future__        import print_function

import numpy as np
import netCDF4 as nc
import pandas as pd
import sys
import re
import warnings
import logging

from datetime          import datetime, timedelta
from netCDF4           import Dataset
from os                import path
from pathlib           import Path
from pydap.client      import open_url, open_dods
from pydap.cas.urs     import setup_session
from scipy.interpolate import interp1d
from typing            import Optional

from globsim.download.GenericDownload import GenericDownload
from globsim.meteorology import pressure_from_elevation

logger = logging.getLogger("globsim.download")

warnings.filterwarnings("ignore", category=UserWarning, module='netCDF4')


class SaveNCDF_pl_3dm():
    """ write output netCDF file for abstracted variables from original
    meteorological data at pressure levels
    """

    def safety_checks(self, data_3dmasm, data_3dmana):
        assert(np.array_equal(data_3dmasm[0]['lat'].data, data_3dmana[0]['lat'].data))
        assert(np.array_equal(data_3dmasm[0]['lev'].data, data_3dmana[0]['lev'].data))
        assert(np.array_equal(data_3dmasm[0]['lon'].data, data_3dmana[0]['lon'].data))

    def saveData(self, data_3dmasm, data_3dmana, dir_data, elevation):
        """
        build NetCDF file for saving output variables (T, U, V, H, RH, lat, lon, levels, time)
        """
        self.safety_checks(data_3dmasm, data_3dmana)

        # set up File
        file_ncdf  = path.join(dir_data, f"merra_pl_{data_3dmana[0]['time'].begin_date}_to_{data_3dmana[-1]['time'].begin_date}.nc")
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        rootgrp.source      = 'MERRA-2, meteorological variables from metadata at pressure levels'
        logger.debug(f"Creating netcdf file {file_ncdf} for pl_3dm data")
        
        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('level', len(data_3dmana[0]['lev']))
        _  = rootgrp.createDimension('lat', len(data_3dmana[0]['lat']))
        _  = rootgrp.createDimension('lon', len(data_3dmana[0]['lon']))

        # Create coordinate variables
        Latitudes               = rootgrp.createVariable('latitude', "f8",('lat'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = data_3dmana[0]['lat'].data    # pass the values of latitude

        Longitudes               = rootgrp.createVariable('longitude', "f8",('lon'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = data_3dmana[0]['lon'].data

        Level = rootgrp.createVariable('level','i4', ('level'))
        Level.standard_name = "air_pressure"
        Level.long_name = "vertical level"
        Level.units = "hPa"
        Level.positive = "down"
        Level.axis = "Z"
        Level[:] = data_3dmana[0]['lev'].data

        Time  = rootgrp.createVariable('time', 'i4', ('time'))
        Time.standard_name = "time"
        Time.units  = "hours since 1980-1-1 00:00:00"
        Time.calendar      = "gregorian"

        netCDFTime = []
        for x in data_3dmana:
            time_datetime = nc.num2date(x['time'].data, units=x['time'].units, calendar=Time.calendar)  # Convert units
            time_nc = nc.date2num(time_datetime, units=Time.units, calendar=Time.calendar)
            netCDFTime.extend(time_nc)

        Time[:] = netCDFTime

        data_variables_3dmasm = [v for v in list(data_3dmasm[0]) if v not in ['time', 'lat', 'lon', 'lev']]
        data_variables_3dmana = [v for v in list(data_3dmana[0]) if v not in ['time', 'lat', 'lon', 'lev']]

        # Treat 6-hourly data
        for x in data_variables_3dmana:
            out_var = rootgrp.createVariable(x, 'f4', ('time', 'level','lat','lon'), fill_value=9.9999999E14)
            out_var.long_name = data_3dmana[0][x].long_name
            out_var.units         = data_3dmana[0][x].units
            out_var.missing_value = np.array(data_3dmana[0][x].missing_value, dtype=out_var.dtype)

            # stack data arrays
            all_data = np.concatenate([dataset[x].data for dataset in data_3dmana], axis=0)
            if x in ["U", "V"]:
                all_data = MERRAdownload.constant_extrapolation(all_data)

            elif x in ["T"]:
                h_data = np.concatenate([dataset["H"].data for dataset in data_3dmana], axis=0)  # could also get from rootgrp["H"]
                all_data = MERRAdownload.constant_extrapolation(all_data)

            out_var[:] = all_data

        # Treat 3-hourly data
        for x in data_variables_3dmasm:
            out_var = rootgrp.createVariable(x, 'f4', ('time', 'level','lat','lon'), fill_value=9.9999999E14)
            out_var.long_name = data_3dmasm[0][x].long_name
            out_var.units         = data_3dmasm[0][x].units
            out_var.missing_value = np.array(data_3dmasm[0][x].missing_value, dtype=out_var.dtype)

            # stack data arrays
            all_data = np.concatenate([dataset[x].data for dataset in data_3dmasm], axis=0)

            if x in ["RH"]:
                all_data = MERRAdownload.constant_extrapolation(all_data)

            out_var[:] = all_data[::2,:,:,:]  # downsample to 6-hourly

        rootgrp.close()


class SaveNCDF_sa():
    """ write output netCDF file for abstracted variables from original 2D
        meteorological Diagnostics dataset and suface flux Diagnostics datasets
    """
    def saveData(self, data_2dm, dir_data):
        """build a NetCDF file for saving output variables
        (T2M,U2M, V2M, U10M,V10M, PRECTOT, PRECTOTCORR,T2MDEW,QV2M,
        lat, lon, time)

        Parameters
        ----------
        data_2m : list
            list of DatasetTypes retrieved from
        dir_data : str
            Path to directory in which files will be saved
        """

        # get the varable list and time indices

        # set up File
        file_ncdf  = path.join(dir_data, f"merra_sa_{data_2dm[0]['time'].begin_date}_to_{data_2dm[-1]['time'].begin_date}.nc")
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        logger.debug(f"Creating netcdf file {file_ncdf} for sa data")
        rootgrp.source      = 'MERRA-2, meteorological variables from metadata at surface level'

        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('lat', len(data_2dm[0]['lat']))
        _  = rootgrp.createDimension('lon', len(data_2dm[0]['lon']))

        # Create coordinate variables
        Latitudes               = rootgrp.createVariable('latitude', "f8",('lat'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = data_2dm[0]['lat'].data    # pass the values of latitude

        Longitudes               = rootgrp.createVariable('longitude', "f8",('lon'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = data_2dm[0]['lon'].data

        data_variables = [v for v in list(data_2dm[0]) if v not in ['time', 'lat', 'lon']]

        # Fill in data variables
        for x in data_variables:
            out_var = rootgrp.createVariable(x, 'f4', ('time','lat','lon'), fill_value=9.9999999E14)
            out_var.long_name = data_2dm[0][x].long_name
            out_var.units         = data_2dm[0][x].units
            out_var.missing_value = np.array(data_2dm[0][x].missing_value, dtype=out_var.dtype)

            # stack data arrays
            all_data = np.concatenate([dataset[x].data for dataset in data_2dm], axis=0)
            out_var[:] = all_data

        # Fill in the time
        Time  = rootgrp.createVariable('time', 'i4', ('time'))
        Time.standard_name = "time"
        Time.units  = "hours since 1980-1-1 00:00:00"
        Time.calendar      = "gregorian"

        netCDFTime = []
        for x in data_2dm:
            time_datetime = nc.num2date(x['time'].data, units=x['time'].units, calendar=Time.calendar)  # Convert units
            time_nc = nc.date2num(time_datetime, units=Time.units, calendar=Time.calendar)
            netCDFTime.extend(time_nc)

        Time[:] = netCDFTime

        rootgrp.close()


class SaveNCDF_sf():
    """ write output netCDF file for abstracted variables from original 2D
        radiation Diagnostics datasets  datasets
    """
    def saveData(self, data_2dr, data_2ds, data_2dv, dir_data):
        """Build a NetCDF file for saving output variables (SWGDN, SWGDNCLR,
           LWRNT, LWGEM, LWGNTCLR, LWGAB, LWGABCLR, lat, lon, time)
           """
        self.safety_checks(data_2dr, data_2ds, data_2dv)

        # set up File
        file_ncdf  = path.join(dir_data, f"merra_sf_{data_2dr[0]['time'].begin_date}_to_{data_2dr[-1]['time'].begin_date}.nc")
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        logger.debug(f"Creating netcdf file {file_ncdf} for sf data")
        rootgrp.source      = 'MERRA-2, radiation variables from metadata at surface level'

        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('lat', len(data_2dr[0]['lat']))
        _  = rootgrp.createDimension('lon', len(data_2dr[0]['lon']))

        # Create coordinate variables
        Latitudes               = rootgrp.createVariable('latitude', "f8",('lat'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = data_2dr[0]['lat'].data    # pass the values of latitude

        Longitudes               = rootgrp.createVariable('longitude', "f8",('lon'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = data_2dr[0]['lon'].data

        data_variables_2dr = [v for v in list(data_2dr[0]) if v not in ['time', 'lat', 'lon']]
        data_variables_2ds = [v for v in list(data_2ds[0]) if v not in ['time', 'lat', 'lon']]
        data_variables_2dv = [v for v in list(data_2dv[0]) if v not in ['time', 'lat', 'lon']]

        # Fill in data variables from the three sources
        for source, namelist in zip([data_2dr, data_2ds, data_2dv],
                                    [data_variables_2dr, data_variables_2ds, data_variables_2dv]):

            for x in namelist:
                out_var = rootgrp.createVariable(x, 'f4', ('time','lat','lon'), fill_value=9.9999999E14)
                out_var.long_name = source[0][x].long_name
                out_var.units         = source[0][x].units
                out_var.missing_value = np.array(source[0][x].missing_value, dtype=out_var.dtype)

                # stack data arrays
                all_data = np.concatenate([dataset[x].data for dataset in source], axis=0)
                out_var[:] = all_data

        # Fill in the time
        Time  = rootgrp.createVariable('time', 'i4', ('time'))
        Time.standard_name = "time"
        Time.units  = "hours since 1980-1-1 00:00:00"
        Time.calendar      = "gregorian"

        netCDFTime = []
        for x in data_2dr:
            time_datetime = nc.num2date(x['time'].data, units=x['time'].units, calendar=Time.calendar)  # Convert units
            time_nc = nc.date2num(time_datetime, units=Time.units, calendar=Time.calendar)
            netCDFTime.extend(time_nc)

        Time[:] = netCDFTime

        # Calculate variables
        all_variables = list(rootgrp.variables)

        if 'LWGEM' in all_variables and 'LWGNT' in all_variables:
            self.derive_lwgdn(rootgrp)

        if 'LWGEM' in all_variables and 'LWGNTCLR' in all_variables:
            self.derive_lwgdnclr(rootgrp)

        rootgrp.close()

    def derive_lwgdn(self, rootgrp):
        out_var = rootgrp.createVariable("LWGDN", 'f4', ('time','lat','lon'), fill_value=9.9999999E14)
        out_var.long_name = "downwelling_longwave_flux_in_air"
        out_var.units = 'W/m2'
        out_var[:] = rootgrp['LWGEM'][:].data + rootgrp['LWGNT'][:].data

    def derive_lwgdnclr(self, rootgrp):
        out_var = rootgrp.createVariable("LWGDNCLR", 'f4', ('time','lat','lon'), fill_value=9.9999999E14)
        out_var.long_name = "downwelling_longwave_flux_in_air_assuming_clear_sky"
        out_var.units = 'W/m2'
        out_var[:] = rootgrp['LWGEM'][:].data + rootgrp['LWGNTCLR'][:].data

    def safety_checks(self, data_2dr, data_2ds, data_2dv):
        """ ensure assumptions are true of datasets """
        assert(np.array_equal(data_2dr[0]['time'].data, data_2ds[0]['time'].data))
        assert(np.array_equal(data_2dr[0]['time'].data, data_2dv[0]['time'].data))
        assert(data_2dr[0]['time'].units == data_2dv[0]['time'].units)


class SaveNCDF_sc():
    """ write output netCDF file for abstracted variables from original 2D
        Constant Model Parameters
    """

    def saveData(self, just_the_data, dir_data):
        """ Create a NetCDF file for saving output variables """

        # set up file path and names
        file_ncdf  = path.join(dir_data,("merra_sc" + ".nc"))
        rootgrp = Dataset(file_ncdf, 'w', format='NETCDF4_CLASSIC')
        logger.debug(f"Creating netcdf file {file_ncdf} for sc data")
        rootgrp.source      = 'MERRA2 constant model parameters'

        # Create dimensions
        _  = rootgrp.createDimension('time', None)
        _  = rootgrp.createDimension('lat', len(just_the_data['lat']))
        _  = rootgrp.createDimension('lon', len(just_the_data['lon']))

        # Output the results of extracted variables
        data_variables = [v for v in list(just_the_data) if v not in ['time', 'lat', 'lon']]

        for x in list(data_variables):
            out_var = rootgrp.createVariable(x, 'f4',
                                                ('time','lat','lon'),
                                                fill_value=9.9999999E14)

            out_var.long_name = just_the_data[x].long_name
            out_var.units = just_the_data[x].units
            out_var.missing_value = np.array(just_the_data[x].missing_value, dtype=out_var.dtype)

            out_var[:] = just_the_data[x].data[0]

        # Create dimension variables
        Time  = rootgrp.createVariable('time', 'i4', ('time'))

        Time.standard_name = "time"
        Time.units  = "hours since 1992-01-02 03:00:00"
        Time.calendar      = "gregorian"
        Time[:] = just_the_data['time'].data

        Latitudes  = rootgrp.createVariable('latitude', "f8", ('lat'))
        Latitudes.standard_name = "latitude"
        Latitudes.units         = "degrees_north"
        Latitudes.axis          = "Y"
        Latitudes[:]  = just_the_data['lat'].data

        Longitudes  = rootgrp.createVariable('longitude', "f8", ('lon'))
        Longitudes.standard_name = "longitude"
        Longitudes.units         = "degrees_east"
        Longitudes.axis          = "X"
        Longitudes[:] = just_the_data['lon'].data

        rootgrp.close()


# TODO: update functions of checking available files
class MERRAdownload(GenericDownload):
    """
    Class for MERRA-2 data that has methods for querying NASA GES DISC server,
    and returning data and coordinate variables.
    """
    # Translate between config file variable names (CF standard_names) and MERRA-2

    full_variables_dic = {
            'air_temperature': ['air_temperature',
                                '2-meter_air_temperature'],

            'relative_humidity' : ['relative_humidity',
                                   '2-metre_dewpoint_temperature',
                                   '2-metre_specific_humidity'],

            'precipitation_amount': ['total_precipitation',
                                     'total_precipitation_corrected'],

            'wind_from_direction': ['eastward_wind',
                                    'northward_wind',
                                    '2-meter_eastward_wind',
                                    '2-meter_northward_wind',
                                    '10-meter_eastward_wind',
                                    '10-meter_northward_wind'],

            'wind_speed':['eastward_wind',
                          'northward_wind',
                          '2-meter_eastward_wind',
                          '2-meter_northward_wind',
                          '10-meter_eastward_wind',
                          '10-meter_northward_wind'],

            'downwelling_shortwave_flux_in_air': ['surface_incoming_shortwave_flux'],

            'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['surface_incoming_shortwave_flux_assuming_clear_sky'],

            'downwelling_longwave_flux_in_air': ['surface_net_downward_longwave_flux',
                                                 'longwave_flux_emitted_from_surface'],

            'downwelling_longwave_flux_in_air_assuming_clear_sky':['surface_net_downward_longwave_flux_assuming_clear_sky',
                                                                   'longwave_flux_emitted_from_surface']}

    # build variables Standards Names and referenced Names for downloading
    # from orginal MERRA-2 datasets

    # 3D Analyzed Meteorological fields data
    full_variables_pl_ana = {'geopotential_height':'H',
                             'air_temperature':'T',
                             'eastward_wind':'U',
                             'northward_wind': 'V'}
    # 3D Assimilated Meteorological fields data
    full_variables_pl_asm = {'relative_humidity': 'RH'}

    # 2D Single-Level Diagnostics data (instantaneous)
    full_variables_sm = {'2-meter_air_temperature': 'T2M',
                         '2-meter_eastward_wind': 'U2M',
                         '2-meter_northward_wind':'V2M',
                         '10-meter_eastward_wind':'U10M',
                         '10-meter_northward_wind':'V10M',
                         '2-metre_specific_humidity':'QV2M'}

    # 2D surface flux diagnostics data
    full_variables_sf = {'total_precipitation': 'PRECTOT',
                         'total_precipitation_corrected': 'PRECTOTCORR'}

    # 2D single-level diagnostics (time-averageed)
    full_variables_sv = {'2-metre_dewpoint_temperature': 'T2MDEW'}

    # 2D radiation diagnostics data
    full_variables_sr = {'surface_incoming_shortwave_flux': 'SWGDN',
                         'surface_incoming_shortwave_flux_assuming_clear_sky': 'SWGDNCLR',
                         'surface_net_downward_longwave_flux':'LWGNT',
                         'longwave_flux_emitted_from_surface': 'LWGEM',
                         'surface_net_downward_longwave_flux_assuming_clear_sky': 'LWGNTCLR'}

    def __init__(self, pfile):
        super().__init__(pfile)
        par = self.par

        self._set_input_directory("merra2")

        # chunk size for downloading and storing data [days]
        self.chunk_size = par['chunk_size']

        try:
            self.mode = par['merramode'].lower()
            if self.mode not in ("download", "links", "combine"):
                raise KeyError
        except KeyError:
            self.mode = 'download'

        # time bounds
        init_date  = {'beg': datetime.strptime(par['beg'], '%Y/%m/%d'),
                      'end': datetime.strptime(par['end'], '%Y/%m/%d')}

        self.date = self.update_time_bounds(init_date)

        if self.mode != 'combine':  # could possibly exclude links as well?
            # start connection session
            self.credential = path.join(par['credentials_directory'], ".merrarc")
            self.account = open(self.credential, "r")
            self.inf = self.account.readlines()
            self.username = ''.join(self.inf[0].split())
            self.password = ''.join(self.inf[1].split())
            self.start_session()

            # Create Subsetter objects
            self.build_subsetters()

    @staticmethod
    def constant_extrapolation(all_data):
        """ arr:  array with dimensions [time, lev, lat, lon] """
        
        if np.any(all_data[:, -1, :, :] > 99999):
            raise ValueError(("Missing data for an entire atmospheric column (masked at the highest level for one or more (lat,lon,time) coordinates)."
                            "This might mean the terrain mask includes all the data. Try increasing the max elevation"))

        n_lev = all_data.shape[1]
        
        for lev in range(n_lev - 2, -1, -1):  # start from the second-highest pressure level
            level_slice = all_data[:, lev, :, :]
            missing_index = np.where(level_slice > 99999)

            # replace with data from higher level
            level_slice[missing_index] = all_data[:, lev + 1, :, :][missing_index]  
            all_data[:, lev, :, :] = level_slice
        
        return(all_data)

    def update_time_bounds(self, date):
        """ Update time bounds for downloading so that files aren't re-downloaded
        date : dict
            Dictionary with keys 'beg' and 'end' and values as datetime objects
        Returns
        -------
        date : dict
            Updated dictionary with keys 'beg' and 'end' and values as datetime objects
        """
        import glob
        import re
        files = glob.glob(f"{self.directory}/*.nc")

        pl = re.compile(r"merra_pl_(\d{8})_to_(\d{8})")
        sa = re.compile(r"merra_sa_(\d{8})_to_(\d{8})")
        sf = re.compile(r"merra_sf_(\d{8})_to_(\d{8})")

        last_dates = [[F.search(x).group(2) for x in files if F.search(x)] for F in [pl, sa, sf]]

        if self.mode != "download":
            return date

        if not all(last_dates):  # No matching files
            return date

        else:
            last_completed = [datetime.strptime(max(D), "%Y%m%d") for D in last_dates]

            if min(last_completed) == date['end']:
                logger.error("All data files have already been downloaded. Exiting.")
                sys.exit(0)

            if min(last_completed) != max(last_completed):  # Partially finished chunk
                date['beg'] = min(last_completed) - timedelta(days=self.chunk_size - 1)

            if max(last_completed) < date['beg']:
                # e.g. Download begins at 2000, but there are 1980 downloads in directory.
                # This line ensures the download starts at date['beg'], not in the 80s.
                logger.warning("There are other download files in this directory, beyond the scope of this project.")

            else:  # Completely finished chunk
                date['beg'] = max(last_completed) + timedelta(days=1)

            logger.info(f"Previously downloaded files detected in {self.directory}. Beginning download at {date['beg']}")

        return date

    def getURLs(self, date: "dict[str, datetime]"):
        """ Set up urls by given range of date and type of data to get objected
        url address

        Args:
        date : dict

        Returns
        -------
        list
        """
        # Setup the based url strings
        baseurl_2d = ('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/')
        baseurl_3d = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/')

        baseurl_3dn = (baseurl_3d + 'M2I6NPANA.5.12.4/{YM}/MERRA2_{FN}.inst6_3d_ana_Np.{YMD}.nc4')
        baseurl_3da = (baseurl_3d + 'M2I3NPASM.5.12.4/{YM}/MERRA2_{FN}.inst3_3d_asm_Np.{YMD}.nc4')
        baseurl_2dm = (baseurl_2d + 'M2I1NXASM.5.12.4/{YM}/MERRA2_{FN}.inst1_2d_asm_Nx.{YMD}.nc4')
        baseurl_2dr = (baseurl_2d + 'M2T1NXRAD.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_rad_Nx.{YMD}.nc4')
        baseurl_2ds = (baseurl_2d + 'M2T1NXFLX.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_flx_Nx.{YMD}.nc4')
        baseurl_2dv = (baseurl_2d + 'M2T1NXSLV.5.12.4/{YM}/MERRA2_{FN}.tavg1_2d_slv_Nx.{YMD}.nc4')

        # build the urls list
        urls_3dmana = []
        urls_3dmasm = []
        urls_2dm = []
        urls_2dr = []
        urls_2ds = []
        urls_2dv = []

        for d in pd.date_range(date['beg'], date['end']):
            ym = d.strftime("%Y/%m")
            ymd = d.strftime("%Y%m%d")

            fn = self.get_file_number(d.year, d.month)

            urls_3dmana.append(baseurl_3dn.format(YM=ym, FN=fn, YMD=ymd))
            urls_3dmasm.append(baseurl_3da.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dm.append(baseurl_2dm.format(YM=ym, FN=fn, YMD=ymd))
            urls_2ds.append(baseurl_2ds.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dr.append(baseurl_2dr.format(YM=ym, FN=fn, YMD=ymd))
            urls_2dv.append(baseurl_2dv.format(YM=ym, FN=fn, YMD=ymd))

        # Setup URL for getting constant model parameters (2D, single-level, full horizontal resolution)
        url_2dc = ['https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4']

        return urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv

    @staticmethod
    def get_file_number(year, month=-1):
        """
        The file names consist of a number and a meta data string.
        The number changes over the years. 1980 until 1991 it is 100,
        1992 until 2000 it is 200, 2001 until 2010 it is  300
        and from 2011 until now it is 400.

        401 - reprocessed files
        https://disc.gsfc.nasa.gov/information/documents?title=Records%20of%20MERRA-2%20Data%20Reprocessing%20and%20Service%20Changes
        """
        file_number = ''

        if year >= 1980 and year < 1992:
            file_number = '100'
        elif year >= 1992 and year < 2001:
            file_number = '200'
        elif year >= 2001 and year < 2011:
            file_number = '300'
        elif year >= 2011:
            if (year == 2020 and month == 9):
                file_number = '401'
            elif (year == 2021 and (6 <= month <= 9)):
                file_number = '401'
            else:
                file_number = '400'
        else:
            raise Exception('The specified year is out of range.')
        return file_number

    def start_session(self):
        self.session = setup_session(self.username, self.password,
                                     check_url="https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2_MONTHLY/M2C0NXASM.5.12.4/1980/MERRA2_101.const_2d_asm_Nx.00000000.nc4")

    def build_subsetters(self):
        self.subsetters = {"3dmana": MERRASubsetter('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I6NPANA.5.12.4/2016/01/MERRA2_400.inst6_3d_ana_Np.20160101.nc4', self.session),
                           "3dmasm": MERRASubsetter('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2016/01/MERRA2_400.inst3_3d_asm_Np.20160101.nc4', self.session),
                           "2dm": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I1NXASM.5.12.4/2016/01/MERRA2_400.inst1_2d_asm_Nx.20160102.nc4', self.session),
                           "2ds": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXFLX.5.12.4/2016/01/MERRA2_400.tavg1_2d_flx_Nx.20160101.nc4', self.session),
                           "2dr": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXRAD.5.12.4/1981/01/MERRA2_100.tavg1_2d_rad_Nx.19810101.nc4', self.session),
                           "2dv": MERRASubsetter('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2T1NXSLV.5.12.4/2016/01/MERRA2_400.tavg1_2d_slv_Nx.20160101.nc4', self.session)}

        for s in self.subsetters.values():
            s.set_lon_range(self.area['west'], self.area['east'])
            s.set_lat_range(self.area['south'], self.area['north'])
            s.set_elev_range(self.elevation['min'], self.elevation['max'])

        # set target variables for subsetter objects
        self.subsetters['3dmasm'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_pl_asm))  # Add H manually to match old code
        self.subsetters['3dmana'].set_variables(["H"] + self.getVariables(self.full_variables_dic, self.full_variables_pl_ana))

        self.subsetters['2dm'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sm))
        self.subsetters['2ds'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sf))
        self.subsetters['2dr'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sr))
        self.subsetters['2dv'].set_variables(self.getVariables(self.full_variables_dic, self.full_variables_sv))

    def getVariables(self, full_variables_dic, full_variables_type):
        """
        build the major variable list for retrieving between lists from
        download parameter file and the product
        """
        get_variables = []

        for i in range(0, len(self.variables)):
            for var in full_variables_dic.keys():
                if var == self.variables[i]:
                    var_names = full_variables_dic[var]
                    #  Set up the variables list for accessing original type
                    # of MERRA-2 datasets (3D and 2D)
                    for j in range(0, len(var_names)):
                        for var1 in full_variables_type.keys():
                            if var1 == var_names[j]:
                                get_variables.append(full_variables_type[var1])
        # set the list
        get_variables = list(set(get_variables))
        """
        #add the basic variables into the list
        if 'relative_humidity' in list(full_variables_type.keys()) :
            get_variables.extend(['lat','lon','lev','time'])
        elif 'air_temperature' in list(full_variables_type.keys()):
            # !ADD Geopotential Height in the first element of downloading
            # list. Must be the first one
            get_variables.insert(0,'H')
            # add the variables names of latitude, longitude, levels and time
            get_variables.extend(['lat','lon','lev','time'])
        else:
            # add the variables names of latitude, longitude and time
            get_variables.extend(['lat','lon','time'])
        """
        return get_variables

    def retrieve(self):
        """
        Retrive all required MERRA-2 data from NASA Goddard Earth Sciences Data
        and Information Services Center
        """

        # Chunk size for spliting files and download [days], Format:Integer
        chunk_size = int(self.chunk_size)

        # Get merra-2 2d Constant Model Parameters (outside of time & date looping!)
        logger.info("Download variables From MERRA2 2d Time-Invariant (Single-level, Constant Model Parameters)")
        self.download_merra_to()
        logger.info("MERRA2 2d Time-Invariant Complete")

        # build chunks
        all_dates = pd.date_range(self.date['beg'], self.date['end'])
        chunks = [{'beg': x, 'end': y} for (x,y) in zip(all_dates[0::chunk_size], all_dates[chunk_size - 1::chunk_size])]

        if not chunks[-1]['end'] == self.date['end']:  # if chunks not an even multiple
            chunks.append({'beg': chunks[-1]['beg'] + timedelta(days=chunk_size),
                           'end': self.date['end']})

        for date_range in chunks:
            if self.mode == "download":
                self.download(date_range)
            
            elif self.mode == "links":
                url_file = Path(self.directory, "merra-wishlist.txt")
                self.download_links(date_range, str(url_file))
                
                logger.info(f"Created OPeNDAP links file: {url_file}")
                logger.info(f"To download the files, use the command 'cat {url_file} | xargs -n 1 -P 16 wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition'")
                logger.info("For more information, visit https://disc.gsfc.nasa.gov/data-access#mac_linux_wget")
            
            elif self.mode == "combine":
                M = MerraAggregator(self.directory)
                M.combine(self.elevation, date_range)

    def download_links(self, date_range: "dict[str, datetime]", url_file: str):
        """ Save the links to the datasets (which can be downloaded in paralell with wget + xargs)"""

        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = self.getURLs(date_range)

        _ = [self.subsetters['2dm'].get_download_link(url, url_file=url_file) for url in urls_2dm]
        _ = [self.subsetters['2ds'].get_download_link(url, url_file=url_file) for url in urls_2ds]
        _ = [self.subsetters['2dr'].get_download_link(url, url_file=url_file) for url in urls_2dr]
        _ = [self.subsetters['2dv'].get_download_link(url, url_file=url_file) for url in urls_2dv]
        _ = [self.subsetters['3dmana'].get_download_link(url, url_file=url_file) for url in urls_3dmana]
        _ = [self.subsetters['3dmasm'].get_download_link(url, url_file=url_file) for url in urls_3dmasm]

    def download(self, date_range):
        logger.info(f"Downloading chunk {date_range['beg']} to {date_range['end']}")

        # Build the urls for the chunk range
        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = self.getURLs(date_range)

        # Access, subset (server-side) and download the data (in memory)
        logger.debug(f"Downloading chunk {date_range['beg']} to {date_range['end']} from dataset M2I1NXASM")
        download_sa_start = datetime.now()

        # Download 2d surface analysis

        data_2dm = [self.subsetters['2dm'].subset_dataset(url, metadata=True) for url in urls_2dm]

        logger.info(f"Downloaded surface analysis data in {datetime.now() - download_sa_start}")

        SaveNCDF_sa().saveData(data_2dm, self.directory)

        # Download 2d surface forecast
        download_sf_start = datetime.now()
        logger.info(f"Downloading chunk {date_range['beg']} to {date_range['end']} from dataset M2T1NXFLX")
        data_2ds = [self.subsetters['2ds'].subset_dataset(url, metadata=True) for url in urls_2ds]

        logger.info(f"Downloading chunk {date_range['beg']} to {date_range['end']} from dataset M2T1NXRAD")
        data_2dr = [self.subsetters['2dr'].subset_dataset(url, metadata=True) for url in urls_2dr]

        logger.info(f"Downloading chunk {date_range['beg']} to {date_range['end']} from dataset M2T1NXSLV")
        data_2dv = [self.subsetters['2dv'].subset_dataset(url, metadata=True) for url in urls_2dv]

        logger.info(f"Downloaded surface forecast data in {datetime.now() - download_sf_start}")

        SaveNCDF_sf().saveData(data_2dr, data_2ds, data_2dv, self.directory)

        # Download 3d pressure-level data
        download_pl_start = datetime.now()
        logger.info(f"Downloading chunk {date_range['beg']} to {date_range['end']} from dataset M2I6NPANA")
        data_3dmana = [self.subsetters['3dmana'].subset_dataset(url, metadata=True) for url in urls_3dmana]

        logger.info(f"Downloading chunk {date_range['beg']} to {date_range['end']} from dataset M2I3NPASM")
        data_3dmasm = [self.subsetters['3dmasm'].subset_dataset(url, metadata=True) for url in urls_3dmasm]

        logger.info(f"Downloaded pressure level data in {datetime.now() - download_pl_start}")

        SaveNCDF_pl_3dm().saveData(data_3dmasm, data_3dmana, self.directory, self.elevation)

        logger.info(f"Downloaded chunk in {datetime.now() - download_sa_start}")

    def download_merra_to(self):
        """ Download time-invariant data for the specified lat-lon boundaries """
        # check if it already exists
        if path.isfile(path.join(self.directory, ("merra_sc" + ".nc"))):
            logger.warning("file 'merra_to.nc' already exists and is being skipped")
            return()

        urls_3dmana, urls_3dmasm, urls_2dm, urls_2ds, urls_2dr, url_2dc, urls_2dv = self.getURLs(self.date)
        S = MERRASubsetter(url_2dc[0], self.session)

        get_variables_2dc = ['PHIS','FRLAKE','FRLAND','FRLANDICE',
                             'FROCEAN','SGH','lat','lon','time']

        S.set_elev_range(self.elevation['min'], self.elevation['max'])
        S.set_lat_range(self.area['south'], self.area['north'])
        S.set_lon_range(self.area['west'], self.area['east'])
        S.set_variables(get_variables_2dc[:-3])

        data = S.subset_dataset(url_2dc[0], metadata=True)

        SaveNCDF_sc().saveData(data, self.directory)


class MERRASubsetter:

    LEVS = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750,
                     725, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250,
                     200, 150, 100, 70, 50, 40, 30, 20, 10, 7.0, 5.0, 4.0,
                     3.0, 2.0, 1.0, 0.7, 0.5, 0.4, 0.3, 0.1])

    lat_values = np.arange(-90, 90.5, 0.5)
    lon_values = np.arange(-180, 180, 0.625)

    def __init__(self, url, session):
        """ For subsetting datasets on the MERRA OpenDAP server

        url : str
            path to pydap-accessible URL endpoint for example dataset (same variables and dimensions)

        session : pydap access session object created from setup_session

        2d VAR[lat:lat][lon:lon][time:time] """
        self.session = session
        logger.info(f"Requesting dataset information from {url}")
        self.dataset = open_url(url, session=session)
        self.set_time_step(self.dataset)

    def set_time_step(self, dataset_type):
        try:
            self.n_timesteps = dataset_type['time'].shape[0]
        except KeyError:
            raise KeyError("The selected dataset does not have a time dimension, or the time dimension is not named 'time'.")

    def subset_time(self):
        """ Get the full day """
        return "[0:{}]".format(self.n_timesteps - 1)

    def subset_lat(self, lat_min, lat_max):
        indices = np.where((self.lat_values >= lat_min) & (self.lat_values <= lat_max))
        return "[{}:{}]".format(np.min(indices), np.max(indices))

    def subset_lon(self, lon_min, lon_max):
        indices = np.where((self.lon_values >= lon_min) & (self.lon_values <= lon_max))
        return "[{}:{}]".format(np.min(indices), np.max(indices))

    def subset_lev(self, elev_min, elev_max):
        Pmax = pressure_from_elevation(elev_min) + 55
        Pmin = pressure_from_elevation(elev_max) - 55

        indices = np.where((self.LEVS >= Pmin) & (self.LEVS <= Pmax))
        return "[{}:{}]".format(np.min(indices), np.max(indices))

    def subset_2d_variable(self, variable_name, lat_min, lat_max, lon_min, lon_max):
        var_string = "".join([variable_name, self.subset_time(),
                              self.subset_lat(lat_min, lat_max), self.subset_lon(lon_min, lon_max)])

        return var_string

    def subset_3d_variable(self, variable_name, elev_min, elev_max, lat_min, lat_max, lon_min, lon_max):
        var_string = "".join([variable_name, self.subset_time(), self.subset_lev(elev_min, elev_max),
                              self.subset_lat(lat_min, lat_max), self.subset_lon(lon_min, lon_max)])

        return var_string

    def set_lat_range(self, lat_min, lat_max):
        self.lat_min = lat_min
        self.lat_max = lat_max

    def set_lon_range(self, lon_min, lon_max):
        self.lon_min = lon_min
        self.lon_max = lon_max

    def set_elev_range(self, elev_min, elev_max):
        self.elev_min = elev_min
        self.elev_max = elev_max

    def set_variables(self, variables):
        self.variables = variables

    def subset_variable(self, variable):
        """ Generate subset slice string for a variable

        Returns
        -------
        str : variable name with slices, e.g. "H[0:3][4:59][23:44]"
        """
        is_3d = ('lev' in self.dataset[variable].dimensions)

        if is_3d:
            slice_string = self.subset_3d_variable(variable, self.elev_min, self.elev_max,
                                                   self.lat_min, self.lat_max,
                                                   self.lon_min, self.lon_max)
        else:
            slice_string = self.subset_2d_variable(variable,
                                                   self.lat_min, self.lat_max,
                                                   self.lon_min, self.lon_max)

        return slice_string

    def create_dods_url(self, dataset_url):
        """ Create a DODS uri link to access a subset of the data by adding variables and slices """
        url = self.create_request_url(dataset_url, 'dods')
        return url

    def create_request_url(self, dataset_url, type='nc4'):
        """ Create a uri link to access a subset of the data by adding variables and slices """
        uri_parameters = []

        for var in self.variables:
            uri_parameters.append(self.subset_variable(var))

        uri_parameters.append('lon' + self.subset_lon(self.lon_min, self.lon_max))
        uri_parameters.append('lat' + self.subset_lat(self.lat_min, self.lat_max))
        uri_parameters.append('time' + self.subset_time())

        if 'lev' in self.dataset:
            uri_parameters.append('lev' + self.subset_lev(self.elev_min, self.elev_max))

        dods_url = dataset_url + f".{type}?" + ",".join(uri_parameters)

        return dods_url
    
    def get_download_link(self, url, url_file):
        """ download a link to the dataset that can be downloaded with wget
        
        Parameters
        ----------
        url : [type]
            base URL for a dataset that can be opened with pydap.open_url
        
        """
        ncurl = self.create_request_url(url)
        with open(url_file, 'a') as f:
            f.write(ncurl + "\n")
        
        return None

    def subset_dataset(self, url, metadata=False):
        """ Return a subset dataset as

        Parameters
        ----------
        url : [type]
            base URL for a dataset that can be opened with pydap.open_url

        Returns
        -------
        DatasetType : pydap DatasetType
        """
        dods_url = self.create_dods_url(url)

        logger.debug(f"Downloading subset from {dods_url}")
        ds = open_dods(dods_url, session=self.session, metadata=metadata)

        return ds


class MerraAggregator():
    """ Turn downloaded Merra files into *_pl, *_sf, and *_sa files """
    def __init__(self, directory: str):
        self.directory = directory
        self.inst6_3d_ana_Np = map_dates(directory,"MERRA2_?0?.inst6_3d_ana_Np*")
        self.inst3_3d_asm_Np = map_dates(directory,"MERRA2_?0?.inst3_3d_asm_Np*")
        self.inst1_2d_asm_Nx = map_dates(directory,"MERRA2_?0?.inst1_2d_asm_Nx*")
        self.tavg1_2d_flx_Nx = map_dates(directory,"MERRA2_?0?.tavg1_2d_flx_Nx*")
        self.tavg1_2d_rad_Nx = map_dates(directory,"MERRA2_?0?.tavg1_2d_rad_Nx*")
        self.tavg1_2d_slv_Nx = map_dates(directory,"MERRA2_?0?.tavg1_2d_slv_Nx*")

    def combine(self, elevation: dict, date_range: "Optional[dict[str, datetime]]"):
        if date_range is None:
            date_range = {'beg': datetime(year=1900, month=1, day=1),
                          'end': datetime(year=2100, month=1, day=1)}

        logger.info(f"Creating surface analysis (sa) file for {date_range['beg']} to {date_range['end']}")

        data_2dm = [VirtualDataset(x) for (d, x) in self.inst1_2d_asm_Nx if self.in_date_range(d, date_range)]

        SaveNCDF_sa().saveData(data_2dm, self.directory)
        del data_2dm  # close datasets to avoid "too many files open error"

        logger.info(f"Creating surface forecast (sf) file for {date_range['beg']} to {date_range['end']}")

        data_2dr = [VirtualDataset(x) for (d, x) in self.tavg1_2d_rad_Nx if self.in_date_range(d, date_range)]
        data_2ds = [VirtualDataset(x) for (d, x) in self.tavg1_2d_flx_Nx if self.in_date_range(d, date_range)]
        data_2dv = [VirtualDataset(x) for (d, x) in self.tavg1_2d_slv_Nx if self.in_date_range(d, date_range)]

        SaveNCDF_sf().saveData(data_2dr, data_2ds, data_2dv, self.directory)

        del data_2dr, data_2ds, data_2dv  # close datasets to avoid "too many files open error"

        logger.info(f"Creating pressure level (pl) file for {date_range['beg']} to {date_range['end']}")

        data_3dmasm = [VirtualDataset(x) for (d, x) in self.inst3_3d_asm_Np if self.in_date_range(d, date_range)]
        data_3dmana = [VirtualDataset(x) for (d, x) in self.inst6_3d_ana_Np if self.in_date_range(d, date_range)]

        SaveNCDF_pl_3dm().saveData(data_3dmasm, data_3dmana, self.directory, elevation)

        del data_3dmana, data_3dmasm  # close datasets to avoid "too many files open error"

    def in_date_range(self, test_date: datetime, date_range: "dict[str, datetime]") -> bool:
        return ((test_date >= date_range['beg']) and (test_date <= date_range['end']))
        

def map_dates(directory: str, globtext: str) -> "zip[tuple[datetime, Path]]":
    files = Path(directory).glob(globtext)
    files_sorted = sorted(files, key=lambda x: x.name[11:])  # don't distinguish MERRA2_?00* and MERRA2_?01*
    merra_date = re.compile(r"(\d{8})\.nc4")
    dates = [datetime.strptime(merra_date.search(str(x)).group(1), r"%Y%m%d") for x in files_sorted]
    # dates = [x for x in files_sorted]
    return zip(dates, files_sorted)
    

class VirtualDataset:
    """" treat netcdf4 files as pydap DODS dataset """
    def __init__(self, dataset):
        if isinstance(dataset, str) or isinstance(dataset, Path):
            self.dataset = nc.Dataset(str(dataset))
        elif isinstance(dataset, nc.Dataset):
            self.dataset = dataset
        else: 
            raise TypeError("dataset wrong type")
        
        self.dataset.set_auto_mask(False)
        self.variables = dict()

        for variable in self.dataset.variables:
            self.variables[variable] = VirtualVariable(self.dataset[variable])
    
    def __getitem__(self, value):
        return self.variables[value]

    def __getattr__(self, value):
        return self.dataset.getnccattr(value)

    def __iter__(self):
        varnames = [v for v in self.dataset.variables]
        return iter(varnames)

    def __del__(self):
        self.close()
        del self.dataset

    def close(self):
        self.dataset.close()


class VirtualVariable:
    """ treat netcdf4 variables as pydap DODS dataset """ 
    def __init__(self, variable:"nc.Variable"):
        self.variable = variable

    def __len__(self):
        return len(self.variable)

    def __getattr__(self, value):
        if value == 'data':
            return self.variable[:]
        else: 
            try:
                return self.variable.getncattr(value)
            except AttributeError:
                return getattr(self.variable, value)


