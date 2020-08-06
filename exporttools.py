"""
Export functions to convert globsim output to other file types
"""


def globsimScaled2Pandas(ncdf_in, station_nr):
    """
    Read a scaled (or interpolated) globsim netCDF file and return all values
    for one station as a Pandas data frame.

    ncdf_in: full path to a globsim netCDF (by station)

    station_nr: station_number, as given in the stations .csv file to identify
                the station.

    """
    # open file
    ncf = nc.Dataset(ncdf_in, 'r')

    # station mask
    sm = ncf.variables['station'][:] == int(station_nr)
    # list variables
    varlist = [x.encode('UTF8') for x in ncf.variables.keys()]

    # get and convert time
    time = ncf.variables['time'][:]
    t_unit = ncf.variables['time'].units
    t_cal = ncf.variables['time'].calendar
    time = nc.num2date(time, units = t_unit, calendar = t_cal)

    # make data frame with time
    df = pd.DataFrame(data=time, columns=['time'])
    # add variables
    for var in varlist:
        if variables_skip(var):
            continue
        data = ncf.variables[var][:,sm]
        df = pd.concat([df, pd.DataFrame(data=data, columns=[var])], axis=1)

    return df


def globsim2GEOtop(ncdf_globsim, txt_geotop):
    """
    Convert globsim scaled netCDF to GEOtop meteo file.
    """

    outfile = '/Users/stgruber/Supervision/MSc/Mary_Pascale_Laurentian/reanalysis/station/Meteo_0001.txt'

    # read time object
    time = self.rg.variables['time']
    date = self.rg.variables['time'][:]

    # read all other values
    columns = ['Date', 'AIRT_ERA_C_pl', 'AIRT_ERA_C_sur', 'PREC_ERA_mm_sur', 'RH_ERA_per_sur', 'SW_ERA_Wm2_sur', 'LW_ERA_Wm2_sur', 'WSPD_ERA_ms_sur', 'WDIR_ERA_deg_sur']
    metdata = np.zeros((len(date), len(columns)))
    metdata[:,0] = date
    for n, vn in enumerate(columns[1:]):
        metdata[:, n+1] = self.rg.variables[vn][:, 0]

    # make data frame
    data = pd.DataFrame(metdata, columns=columns)
    data[['Date']] = nc.num2date(date, time.units, calendar=time.calendar)

    # round
    decimals = pd.Series([2, 1, 1, 1, 1, 1, 1, 1], index=columns[1:])
    data.round(decimals)

    # export to file
    fmt_date = "%d/%m/%Y %H:%M"
    data.to_csv(outfile, date_format=fmt_date, index=False, float_format='%.2f')


def globsim2CLASS(ncdf_globsim, met_class, station_nr):
    """
    Convert globsim scaled netCDF to CLASS-CTEM .met file.

    ncdf_globsim: full path to a globsim scaled netCDF (by station)

    met_class: full path to the CLASS-CTEM met file to write.

    station_nr: station_number, as given in the stations .csv file to identify
                the station.

    The columns in CLASS-CTEM MET files are:
    1)  Hour
    2)  Minute
    3)  Day of year
    4)  Year YYYY
    5)  Shortwave Radiation (W/m2)
    6)  Longwave Radiation (W/m2)
    7)  Precip (mm/s)
    8)  Temp.(°C)
    9)  Specific Humidity (Kg/Kg)
    10) Wind Speed (m/s)
    11) Pressure (Pa)

    """

    # columns to export
    columns = ['time', 'SW_ERA_Wm2_sur', 'LW_ERA_Wm2_sur', 'PREC_ERA_mmsec_sur',
               'AIRT_ERA_C_sur', 'SH_ERA_kgkg_sur', 'WSPD_ERA_ms_sur',
               'AIRT_PRESS_Pa_pl']

    # output ASCII formatting
    formatters={"time": "  {:%H %M  %j  %Y}".format,
                "SW_ERA_Wm2_sur": "{:8.2f}".format,
                "LW_ERA_Wm2_sur": "{:8.2f}".format,
                "PREC_ERA_mmsec_sur": "{:13.4E}".format,
                "AIRT_ERA_C_sur": "{:8.2f}".format,
                "SH_ERA_kgkg_sur": "{:11.3E}".format,
                "WSPD_ERA_ms_sur": "{:7.2f}".format,
                "AIRT_PRESS_Pa_pl": "{:11.2f}".format}

    # get data
    df = globsimScaled2Pandas(ncdf_globsim, station_nr)

    # convery precipitation
    df["PREC_ERA_mmsec_sur"] = df["PREC_ERA_mm_sur"] / 1800.0

    # write FORTRAN formatted ASCII file
    with open(met_class, 'w') as f:
        f.write(' ')
        f.write(df.to_string(columns = columns,
                formatters=formatters,
                header=False, index=False))
    f.close()


def nc_to_clsmet(ncd, out_dir, src, start=None, end=None):
    """
    @args
    ncd: netcdf dataset
    src: data source ("ERAI","ERA5" "MERRA2", "JRA55")
    """
    # open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)

    # find number of stations
    nstn = len(n['station'][:])

    # get date / time columns
    time = nc.num2date(n['time'][:],
                       units=n['time'].units,
                       calendar=n['time'].calendar)
    time_step = (time[1] - time[0]).seconds # in second


    HH =   [x.timetuple().tm_hour for x in time]
    MM =   [x.timetuple().tm_min  for x in time]
    DDD =  [x.timetuple().tm_yday for x in time]
    YYYY = [x.timetuple().tm_year for x in time]

    TIME = np.stack((HH, MM, DDD, YYYY))

    # get shortwave
    SW = "SW_{}_Wm2_sur".format(src)
    SW = n[SW][:]

    # get longwave
    LW = "LW_{}_Wm2_sur".format(src)
    LW = n[LW][:]

    # get precip
    PREC = "PREC_{}_mm_sur".format(src)
    PREC = n[PREC][:]
    PREC /= time_step # [mm/0.5h] to [mm/s]

    # get temp
    AIRT = "AIRT_{}_C_sur".format(src)
    AIRT = n[AIRT][:]

    # get specific humidity
    SH = "SH_{}_kgkg_sur".format(src)
    SH = n[SH][:]

    # get wind speed
    WSPD = "WSPD_{}_ms_sur".format(src)
    WSPD = n[WSPD][:]

    # get pressure
    PRESS = "PRESS_{}_Pa_pl".format(src)
    PRESS = n[PRESS][:]

    # get site names
    NAMES = nc.chartostring(n['station_name'][:])

    data = np.stack((SW, LW, PREC, AIRT, SH, WSPD, PRESS))

    # write output files
    for i in range(nstn):

        # massage data into the right shape
        out_array = data[:, :, i]
        out_array = np.concatenate((TIME, out_array), 0)
        out_array = np.transpose(out_array)

        # get station name
        st_name = NAMES[i]

        filename = "{}_{}_{}.MET".format(i, st_name, src)
        savepath = path.join(out_dir, filename)

        # create file
        np.savetxt(savepath, out_array,
                   fmt = [" %02u", "%02u", " %03u",
                          " %2u", "%8.2f", "%8.2f",
                          "%13.4E", "%8.2f", "%11.3E",
                          "%7.2f", "%11.2f" ])


def nc_to_gtmet(ncd, out_dir, src, start=None, end=None):
    """
    @args
    ncd: netcdf dataset
    src: data source ("ERAI","ERA5" "MERRA2", "JRA55")
    """
    # open netcdf if string provided
    if type(ncd) is str:
        n = nc.Dataset(ncd)

    # find number of stations
    nstn = len(n['station'][:])

    # get date / time column
    time = nc.num2date(n['time'][:],
                       units=n['time'].units,
                       calendar=n['time'].calendar)

    time = [x.strftime('%d/%m/%Y %H:%M') for x in time]
    time = pd.DataFrame(time)

    # get precip
    PREC = "PREC_{}_mm_sur".format(src)
    PREC = n[PREC][:]

    # get wind velocity
    WSPD = "WSPD_{}_ms_sur".format(src)
    WSPD = n[WSPD][:]

    # get wind direction
    WDIR = "WDIR_{}_deg_sur".format(src)
    WDIR = n[WDIR][:]

    # get windx and windy

    # get RH
    RH = "RH_{}_per_sur".format(src)
    RH = n[RH][:]

    # get air temp
    AIRT = "AIRT_{}_C_sur".format(src)
    AIRT = n[AIRT][:]

    # get dew temp (missing?)

    # get air pressure
    PRESS = "PRESS_{}_Pa_pl".format(src)
    PRESS = n[PRESS][:]
    PRESS *= 1e-5      # convert to bar for geotop

    # get shortwave solar global (direct / diffuse missing?)
    SW = "SW_{}_Wm2_sur".format(src)
    SW = n[SW][:]

    # get longwave incoming
    LW = "LW_{}_Wm2_sur".format(src)
    LW = n[LW][:]

    # get site names
    NAMES = nc.chartostring(n['station_name'][:])

    # combine data variables into array
    data = np.stack((PREC, WSPD, WDIR, RH, AIRT, PRESS, SW, LW))

    # write output files
    for i in range(nstn):

        # massage data into the right shape
        out_df = pd.DataFrame(np.transpose(data[:, :, i]))
        out_df = pd.concat([time, out_df], axis=1)
        out_df.columns = ["Date", "IPrec", "WindVelocity", "WindDirection", "RH",
                            "AirTemp", "AirPress", "SWglobal", "LWin"]

        # get station name
        st_name = NAMES[i]

        # prepare paths
        filename = "{}-{}_{}_Forcing_0001.txt".format(i, st_name, src)
        savepath = path.join(out_dir, filename)

        # create file
        out_df.to_csv(savepath, index=False, float_format="%10.5f")