#TimeScaler
from datetime import datetime
import numpy as np
import netCDF4 as nc


"""
@param timestep in seconds
@param num_times how many
"""
def build_time_array(start_time: datetime, timestep: int, num_times: int, output_units: str, output_calendar: str):
    print(start_time)
    time_array = np.arange(num_times, dtype='int64') * timestep
    
    date_array = nc.num2date(time_array,
                             units=f"seconds since {start_time.strftime('%Y-%m-%d %H:%M:%S.%f')}",
                             calendar="gregorian")
    
    result = nc.date2num(date_array,
                         units=output_units,
                         calendar=output_calendar)

    return result


if __name__ == "__main__":
    f = r"/home/nbr512/glob-sg/interpolated/jra_pl_ldg_surface.nc"
    d = nc.Dataset(f, "r")
    units = d['time'].units
    cal = d['time'].calendar
    start = nc.num2date(d['time'][:][0], units, cal)
    res  = build_time_array(start, 3600, 400, units, cal)
    print(res)