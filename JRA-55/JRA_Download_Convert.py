from datetime import date, datetime
from dateutil.rrule import rrule, DAILY
from ftplib import FTP
from netCDF4 import Dataset
import netCDF4 as nc
import pygrib
import time
import sys


"""
Download ranges of data from the JRA-55 server
Data available from 1958 to present date (give 1 or 2 days for upload)
~Author: Christopher Molnar
~Date: July 19, 2017
"""
class JRA_Download:
  
  
  """
  Get the time sets that you want to download (start day to final day)
  returns:
  -start_date: A tuple containing the start year, month and day
  -end_date: A tuple containing the end year, month and day
  """  
  def TimeSet(self):
    x = True
    while(x == True):
      try: 
        x = False
        print "\nEnter the beginning date of the data you want"
        start_year = int(raw_input("\nStart Year: ")) # Get the start year
        start_month = int(raw_input("\nStart Month: ")) # Get the start month
        start_day = int(raw_input("\nStart Day: ")) # Get the start day
        start_data = date(start_year, start_month, start_day) # Put the start day in date format
       
        print "\nEnter the ending date of the data you want"
        end_year = int(raw_input("\nEnd Year: ")) # Get the end year
        end_month = int(raw_input("\nEnd Month: ")) # Get the end month
        end_day = int(raw_input("\nEnd Day: ")) # Get the end day
        end_data = date(end_year, end_month, end_day) # Put the end day in date format
                       
      except:
        x = True
        print"\nInvalid data entered"
        print "*Please enter valid digits for Year, Month, and Day"
        print "*Ex: Jan = 01"
       
      try:
        if (start_data > end_data): # Make sure the start day is before the end day
          print"\nInvalid data entered"
          print "*The end day you chose is before the start day you chose" 
          x = True
        elif (start_year < 1958 or end_year < 1958): # Make sure that the data selected is 1958 or later 
          print"\nInvalid data entered"
          print "*Data is only available from 1958 to the current day (data in the last couple of days may not be up instantly)" 
          x = True
        elif ((end_year > datetime.now().year) or (end_year >= datetime.now().year and end_month > datetime.now().month) or (end_year >= datetime.now().year and end_month >= datetime.now().month and end_day > datetime.now().day)): # Make sure none of the dates are in the future
          print"\nInvalid data entered"
          print "*You entered a date in the future"
          print "*Data is only available from 1958 to the current day (data in the last couple of days may not be up instantly)"
          x = True
        
      except:
        x = True
    
    return(start_data, end_data)
  
  
  """
  Choose your Latitude and Longitude restraints
  Return:
  -latLowerPosition: Lower position of the latitude
  -latUpperPosition: Upper position of the latitude
  -lonLowerPosition: Lower position of the longitude
  -lonUpperPosition: Upper position of the longitude
  """
  def ChooseLatLon(self):
      #Convert
      latRange = []
      lonRange = []
      
      x = 0.0
      while x < 360:
        lonRange.append(x)
        x += 1.25    
        
      y = 90.0
      while y >= -90:
        latRange.append(y)
        y -= 1.25   
      
      print "\nPlease specify the upper and lower bounds for Latitude and Longitude"
      print "\nLatitude is from -90 to 90"
      run1 = True
      while (run1 == True):
        try:
          latLow = float(raw_input("Enter the lower bound of Latitude: "))
          latUp = float(raw_input("Enter the upper bound of Latitude: "))
        
          latLowerPosition = abs(latLow - 90) // 1.25
          latUpperPosition = abs(latUp - 90) // 1.25
          
          if (-90 <= latLow < latUp and latLow < latUp <= 90):
            run1 = False
          else:
            print "\nThe lower bound must be less then the upper bound"
            print "And the bounds must be  between -90 and 90 (inclusively)"
        
        except:
          print "\nInvalid Entry"
          print "Please make sure your latitude position is between -90 and 90 (inclusively)"
          
      print "\nLongitude is from 0 to 358.75"
      run2 = True
      while (run2 == True):
        try:
          lonLow = float(raw_input("Enter the lower bound of Longitude: "))
          lonUp = float(raw_input("Enter the upper bound of Longitude: "))
        
          lonLowerPosition = lonLow // 1.25
          lonUpperPosition = lonUp // 1.25
          
          if (0 <= lonLow < lonUp and lonLow < lonUp <= 358.75):
            run2 = False
          else:
            print "\nThe lower bound must be less then the upper bound"
            print "And the bounds must be  between 0 and 358.75 (inclusively)"
       
        except:
          print "\nInvalid Entry"
          print "Please make sure your longitude position is between 0 and 358.75 (inclusively)"
  
      return (latLowerPosition, latUpperPosition, lonLowerPosition, lonUpperPosition)      
  
  
  """
  Download the files from the server
  Parameters:
  -start_day: The first day to start downloading from
  -end_day: The last day that you want to download
  -data_type: The file that is being downloaded
  -hourly_increment: The hour increments used to record the data
  -middle: Extra text in the middle of the filename needed for downloading certain files
  """
  def ftp_Download(self, start_day, end_day, data_type, hourly_increment, middle):
    
    tries = 1
    server = False
    
    while (tries <= 5 and server == False): # Try to connect to server (try 5 times)
      try: # Try to connect to JRA Site
        ftp = FTP("ds.data.jma.go.jp")
        server = True 
      except EOFError: # Catch server disconection error
        print "\nConnection to server terminated :("
        print "Try: " + str(tries) + " " + data_type + middle
        server = False
      tries += 1
    
    if (server == True): # If the connection works download files
      username = "jra04266"
      password = "xUhVnbjT"
      # print "\nEnter your JRA-55 Data"
      # username = raw_input("Enter your JRA-55 username: ")
      # password = raw_input("Enter your JRA-55 password: ")
      ftp.login(username, password) # Login (username, password)
      
      for dt in rrule(DAILY, dtstart = start_day, until = end_day): # Loop from start day to end day 
        path = "/JRA-55/Hist/Daily/" + data_type + "/" + dt.strftime("%Y") + dt.strftime("%m") # Generate the path  
        ftp.cwd(path) # Change Working Directory
      
        for x in range(0, 24, hourly_increment): # Loop through the hourly increments
      
          if (x <10): # Add a 0 infront any of the hour times that are less 10
            ending = "0" + str(x)
          else:
            ending = str(x)       

          filename = data_type + middle + "." +  dt.strftime("%Y") + dt.strftime("%m") + dt.strftime("%d") + ending # Generate the filename
          
          localfile = open(filename, 'wb')
     
          try: # try to download the file
            ftp.retrbinary("RETR %s" % filename , localfile.write) # Download file
          except:
            print "Error downloading file: " + filename
          
          localfile.close() # Close File 
        
      print data_type + middle + " Downloaded" # Download complete
        
      ftp.quit() # Close Connection   
    
    else:
      print "\nAttempted to download" + data_type + middle
      print "Tried to connect 5 times but failed"   
      print "Please retry later"
      sys.exit(0)  
    

"""
Convert the grib files downloaded from the JRA-55 into netCDF formats
"""
class Grib2CDF:


    """
    Convert lat and lon into the needed format for latitude and longitude
    Parameters:
    -lat: All of the latitude coordinates from the grib file
    -lon: All of the longitude coordinates from the grib file
    Returns:
    -lats: An array with all of the latitude coordinates
    -lons: An array with all of the longitude coordinates
    """   
    def ConvertLatLon(self, lat, lon):
      lats = []
      lons = []

      for x in range (0,145):
        for y in range (0,288):
          if (x == 1):
            lons.append(lon[x][y])
        lats.append(lat[x][0]) 
      
      return(lats, lons) 


    """
    Extracts the needed latitudes and longitudes from a Grib file
    Parameters:
    -filename: The name of the Grib file that the latitude and longitude are being extracted from
    Returns:
    -The result from the ConvertLatLon function
    """
    def GetLatLon(self, filename):
      try:
        grbs = pygrib.open(filename)
      except:
        print "file: " + filename +  " not found :("
        print "Quiting program3"
        sys.exit(0)
      
      for g in grbs:
        lat, lon = g.latlons() # Get Latitude and Longitutde
        break
        
      grbs.close()
      return (self.ConvertLatLon(lat,lon))
      
    
    """
    Extracts the needed values and levels from a Grib file
    Parameters:
    -filename: The name of the Grib file that the values and levels are being extracted from
    -date: An array with all of the dates in it
    -dataName: The name of the variable that you want the values and level for
    Returns:
    -value: An array with all of the values
    -levels: An int with the datas level
    """
    def ExtractData(self, filename, date, dataName):    
      value = []     
      for z in range(0, len(date)):
        try:
          grbs = pygrib.open(filename + date[z])
        except:
          print "file: " + filename +  " not found :("
          print "Quiting program"
          sys.exit(0)
                    
          # Loop through the Grib file and extract the data you need
        for g in grbs:
          if (g.shortName == dataName):
            value.append(g.values)
            level = g.level
            break
          
        grbs.close()
        
      return (value, level)


    """
    Creates the Latitude and Longitude dimensions, variables, standard names and units
    Parameter:
    -f: The netCDF file that is being created
    Returns:
    -latitude: The dimension latitude
    -longitude: The dimesion longitude
    -latitudes: The variable latitude
    -longitudes: The variable longitude
    """ 
    def CreateLatLon(self, f):
        # Create Dimensions for variables 
        latitude = f.createDimension("latitude", 145)
        longitude = f.createDimension("longitude", 288)
        # Create Variables
        latitudes = f.createVariable("latitude", "f4", "latitude")
        longitudes = f.createVariable("longitude", "f4", "longitude")
        # Standard Names
        latitudes.standard_name = "latitude"
        longitudes.standard_name = "longitude"
        # Units
        latitudes.units = "degree_north"
        longitudes.units = "degree_east"
        return(latitude, longitude, latitudes, longitudes)


    """
    Creates the time series for the netCDF file
    Adds the hour segments to the days
    Parameters:
    -date: An array with the dates
    -timeSize: The total amount of time segments in the netCDF file
    -hours: The 3 or 6 hours depending on the reanlysis data
    -numSegments: The amount of segments in each day
    Returns:
    -date: The updates array with all of the dates and hours
    """
    def TimeSeries(self, date, timeSize, hours, numSegments): 
      for size in range(0, timeSize, numSegments):
          YMD = date[size]
          date[size] = YMD + "00"

          for h in range(1, numSegments):
            end = h * hours
            if (end < 10):
              end = "0" + str(end)
            else:
              end = str(end)
            date.insert(size + h, YMD + end)
                  
      return(date)
    
    
    """
    Subsets the data in to the desired Latitude and Longitude coordinates
    Parameters:
    -dataValues: The array containing the reanlysis data
    -numDays: The total amount of time segments in the netCDF file
    -latMin: The minimum latitude coordinate
    -latMax: The maxmum latitude coordinate
    -lonMin: The minimum longitude coordinate
    -lonMax: The maximum longitude coordinate
    Returns:
    -ssDay: The subseted data set
    """
    def SubsetTheData(self, dataValues, numDays, latMin, latMax, lonMin, lonMax):
      ssDay = []
      for a in range(0,numDays):
        ssDay.append([])
        for b in range(0,145):
          ssDay[a].append([])
          for c in range(0,288):
            if (latMin <= b <= latMax and lonMin <= c <= lonMax):
              ssDay[a][b].append(dataValues[a][b][c])
            else:
              ssDay[a][b].append(None)
      
      return (ssDay)

   
"""
Class for Forcast data
""" 
class fcst_phy2m:
    
    
    """
    Creates the netCDF file with all of the forcasted variables
    Parameters:
    -startDay: Start day for the netCDF file
    -endDay: End day for the netCDF file
    -numDays: Total number of days being saved to the netCDF file
    -date: All of the dates in an array being saved to the netCDF file
    -lowLat: The minimum latitude
    -upLat The maximum latitude
    -lowLon: The minimum longitide
    -upLon: The maximum longitude
    """
    def Main(self, startDay, endDay, numDays, date, lowLat, upLat, lowLon, upLon):
             
        f = Dataset("subset2_fcst_" + startDay + "-" + endDay + ".nc", "w", format = "NETCDF4") # Name of the netCDF being created 
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon("fcst_phy2m125." + startDay + "00")
        latitude, longitude, latitudes, longitudes = Grib2CDF().CreateLatLon(f)
        latitudes[:] = lats
        longitudes[:] = lons

        
        # Create the Time
        timeSize = numDays * 8
        dateTimes = Grib2CDF().TimeSeries(date, timeSize, 3, 8)
        time = f.createDimension("time", timeSize)
        time = f.createVariable("time", "i4", "time")
        time.standard_name = "time"
        time.units = "Y/M/D/H"
        time.calendar = "gregorian"
        time[:] = dateTimes
        
        # Create the Clear Sky Downward Solar Radiation Flux
        csdsrf, csdsrfLevel = Grib2CDF().ExtractData("fcst_phy2m125.", date, "csdsf")
        csdsrfLevel = f.createDimension("csdsrfLevel", 1)
        clearSkySolarRF = f.createVariable("clear_sky_downward_solar_radiation_flux", "f4", ("time","csdsrfLevel", "latitude", "longitude"))
        clearSkySolarRF.standard_name = "net_downward_shortwave_flux_in_air_assuming_clear_sky"
        clearSkySolarRF.units = "W/(m^2)"
        clearSkySolarRF[:,:,:,:] = Grib2CDF().SubsetTheData(csdsrf, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the Clear Sky Downward Long Wave Radiation Flux
        csdlrf, csdlrfLevel = Grib2CDF().ExtractData("fcst_phy2m125.", date, "csdlf")
        csdlrfLevel = f.createDimension("csdlrfLevel", 1)
        clearSkyLongRF = f.createVariable("clear_sky_downward_long_wave_radiation_flux", "f4", ("time","csdlrfLevel",  "latitude", "longitude"))
        clearSkyLongRF.standard_name = "net_downward_longwave_flux_in_air_assuming_clear_sky"
        clearSkyLongRF.units = "W/(m^2)"
        clearSkyLongRF[:,:,:,:] = Grib2CDF().SubsetTheData(csdlrf, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the Downward Solar Radiation Flux
        dsrf, dsrfLevel = Grib2CDF().ExtractData("fcst_phy2m125.", date,"dswrf")
        dsrfLevel = f.createDimension("dsrfLevel", 1) 
        downwardSolarRF= f.createVariable("downward_solar_radiation_flux", "f4", ("time", "dsrfLevel", "latitude", "longitude"))
        downwardSolarRF.standard_name = "net_downward_shortwave_flux_in_air"
        downwardSolarRF.units = "W/(m^2)"
        downwardSolarRF[:,:,:,:] = Grib2CDF().SubsetTheData(dsrf, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the Downward Long Wave Radiation Flux
        dlrf, dlrfLevel = Grib2CDF().ExtractData("fcst_phy2m125.", date, "dlwrf")
        dlrfLevel = f.createDimension("dlrfLevel", 1) 
        downwardLongRF = f.createVariable("downward_long_wave_radiation_flux", "f4", ("time", "dlrfLevel", "latitude", "longitude"))
        downwardLongRF.standard_name = "net_downward_longwave_flux_in_air"
        downwardLongRF.units = "W/(m^2)"
        downwardLongRF[:,:,:,:] = Grib2CDF().SubsetTheData(dlrf, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the Total Precipitation
        tp, tpLevel = Grib2CDF().ExtractData("fcst_phy2m125.", date, "tpratsfc")
        tpLevel = f.createDimension("tpLevel", 1)
        totalPercipitation = f.createVariable("total_percipitation", "f4", ("time","tpLevel", "latitude", "longitude"))
        totalPercipitation.standard_name = "precipitation_amount"
        totalPercipitation.units =  "mm/day"
        totalPercipitation[:,:,:,:] = Grib2CDF().SubsetTheData(tp, timeSize, upLat, lowLat, lowLon, upLon)

        # Descriptions
        f.description = "fcst example"
        f.history = "Created today"
        f.source = "netCDF4 practice"
 
        f.close() 


"""
Class for Surface data
"""
class anl_surf:
    
    
    """
    Creates the netCDF file with all of the surface variables
    Parameters:
    -startDay: Start day for the netCDF file
    -endDay: End day for the netCDF file
    -numDays: Total number of days being saved to the netCDF file
    -date: All of the dates in an array being saved to the netCDF file
    -lowLat: The minimum latitude
    -upLat The maximum latitude
    -lowLon: The minimum longitide
    -upLon: The maximum longitude
    """
    def Main(self, startDay, endDay, numDays, date, lowLat, upLat, lowLon, upLon):
    
        f = Dataset("subset2_surf_" + startDay + "-" + endDay + ".nc", "w", format = "NETCDF4") # Name of the netCDF being created 
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon("anl_surf125." + startDay + "00")
        latitude, longitude, latitudes, longitudes = Grib2CDF().CreateLatLon(f)
        latitudes[:] = lats
        longitudes[:] = lons

        # Create the Time
        timeSize = numDays * 4  
        dateTimes = Grib2CDF().TimeSeries(date, timeSize, 6, 4)
        time = f.createDimension("time", timeSize)
        time = f.createVariable("time", "i4", "time")
        time.standard_name = "time"
        time.units = "Y/M/D/H"
        time.calendar = "gregorian"
        time[:] = dateTimes
        
        # Create the Dew-Point Depression
        d, dLevel = Grib2CDF().ExtractData("anl_surf125.", date, "depr")
        dLevel = f.createDimension("dLevel", 1) 
        dewPoint = f.createVariable("dew_point", "f4", ("time", "dLevel", "latitude", "longitude"))
        dewPoint.standard_name = "dew_point_depression"
        dewPoint.units = "K"
        dewPoint[:,:,:,:] = Grib2CDF().SubsetTheData(d, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the Relative Humidity
        r, rLevel = Grib2CDF().ExtractData("anl_surf125.", date, "r")
        rLevel = f.createDimension("rLevel", 1)
        rHumidity = f.createVariable("relative_humidity", "f4", ("time", "rLevel", "latitude", "longitude"))
        rHumidity.standard_name = "relative_humidity"
        rHumidity.units = "%"
        rHumidity[:,:,:,:] = Grib2CDF().SubsetTheData(r, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the Temperature
        t, tLevel = Grib2CDF().ExtractData("anl_surf125.", date, "t")
        tlevel = f.createDimension("tLevel", 1)
        temp = f.createVariable("air_temperature", "f4", ("time", "tLevel", "latitude", "longitude"))
        temp.standard_name = "air_temperature"
        temp.units = "K" 
        temp[:,:,:,:] = Grib2CDF().SubsetTheData(t, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the U-Wind Component
        u, uLevel = Grib2CDF().ExtractData("anl_surf125.", date, "u")
        uLevel = f.createDimension("uLevel", 1)
        uWind = f.createVariable("u_wind", "f4", ("time", "uLevel", "latitude", "longitude"))
        uWind.standard_name = "u_component_of_wind"
        uWind.units = "m/s"
        uWind[:,:,:,:] = Grib2CDF().SubsetTheData(u, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Create the V-Wind Component
        v, vLevel = Grib2CDF().ExtractData("anl_surf125.", date, "v")
        vLevel = f.createDimension("vLevel", 1)
        vWind = f.createVariable("v_wind", "f4", ("time", "vLevel", "latitude", "longitude"))
        vWind.standard_name = "v_component_of_wind"
        vWind.units = "m/s"        
        vWind[:,:,:,:] = Grib2CDF().SubsetTheData(v, timeSize, upLat, lowLat, lowLon, upLon)
        
        # Descriptions
        f.description = "anl_surf125 temp example"
        f.history = "Created today"
        f.source = "netCDF4 practice"
        
        f.close() 


"""
Class for Isobaric data
"""
class Isobaric:   


    """
    Extracts the needed values and levels from a Grib file
    Parameters:
    -filename: The name of the Grib file that the values and levels are being extracted from
    -date: An array with all of the dates in it
    -dataName: The name of the variable that you want the values and level for
    Returns:
    -allData: An array with all of the values
    -levels: An array with all of the datas level
    """
    def ExtractData(self, filename, date):    
      allData = []
      for z in range(0, len(date)):
        levels = []
        data = []
        try:
          grbs = pygrib.open(filename + date[z])
        except:
          print "file: " + filename +  " not found :("
          print "Quiting program"
          sys.exit(0)
                    
          # Loop through the Grib file and extract the data you need
        for g in grbs:
          levels.append(g.level)
          data.append(g.values)
        allData.append(data)
         
        grbs.close()
        
      return (allData, levels)
     
     
    """
    Subsets the data in to the desired Latitude and Longitude coordinates
    Parameters:
    -dataValues: The array containing the reanlysis data
    -numDays: The total amount of time segments in the netCDF file
    -latMin: The minimum latitude coordinate
    -latMax: The maxmum latitude coordinate
    -lonMin: The minimum longitude coordinate
    -lonMax: The maximum longitude coordinate
    Returns:
    -ssDay: The subseted data set
    """
    def SubsetTheData(self, dataValues, numDays, levels, latMin, latMax, lonMin, lonMax):
        ssDay = []
        for a in range(0,numDays):
          ssDay.append([])
          for b in range(0,levels):
            ssDay[a].append([])
            for c in range(0,145):
              ssDay[a][b].append([])
              for d in range(0,288):
                if (latMin <= c <= latMax and lonMin <= d <= lonMax):
                  ssDay[a][b][c].append(dataValues[a][b][c][d])
                else:
                  ssDay[a][b][c].append(None)
        
        return (ssDay)
       
    
    """
    Creates the netCDF file with all of the Isobaric variables
    Parameters:
    -startDay: Start day for the netCDF file
    -endDay: End day for the netCDF file
    -numDays: Total number of days being saved to the netCDF file
    -date: All of the dates in an array being saved to the netCDF file
    -lowLat: The minimum latitude
    -upLat The maximum latitude
    -lowLon: The minimum longitide
    -upLon: The maximum longitude
    """
    def Main(self, startDay, endDay, numDays, date, lowLat, upLat, lowLon, upLon):
         
        f = Dataset("subset2_Isobaric_" + startDay + "-" + endDay + ".nc", "w", format = "NETCDF4") # Name of the netCDF being created 
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon("anl_p125_hgt." + startDay + "00") 
        latitude, longitude, latitudes, longitudes = Grib2CDF().CreateLatLon(f)
        latitudes[:] = lats
        longitudes[:] = lons
                
        # Create the Time 
        timeSize = numDays * 4
        dateTimes = Grib2CDF().TimeSeries(date, timeSize, 6, 4)
        time = f.createDimension("time", timeSize)
        time = f.createVariable("time", "i4", "time")
        time.standard_name = "time"
        time.units = "Y/M/D/H"
        time.calendar = "gregorian"
        time[:] = dateTimes
        
        # Create the Geopotential Height
        hValues, hLevels = self.ExtractData("anl_p125_hgt.", date)
        hgtLevel = f.createDimension("hgtLevel", len(hLevels))
        hgtValues = f.createVariable("geopotential_height", "f4", ("time", "hgtLevel", "latitude", "longitude"))
        hgtValues.standard_name = "geopotential_height"
        hgtValues.units = "gpm"
        hgtValues[:,:,:,:] = self.SubsetTheData(hValues, timeSize, 37, upLat, lowLat, lowLon, upLon)
        
        # Create Relative Humidity
        rValues, rLevels = self.ExtractData("anl_p125_rh.", date)
        rhLevel = f.createDimension("rhLevel", len(rLevels))
        relativeHumidity = f.createVariable("relative_humidity", "f4", ("time", "rhLevel", "latitude", "longitude"))
        relativeHumidity.standard_name = "relative_humidity"
        relativeHumidity.units = "%"
        relativeHumidity[:,:,:,:] = self.SubsetTheData(rValues, timeSize, 27, upLat, lowLat, lowLon, upLon)
        
        # Create Temperature
        tValues, tLevels = self.ExtractData("anl_p125_tmp.", date)
        tmpLevel = f.createDimension("tmpLevel", len(tLevels))
        temperature = f.createVariable("temperature", "f4", ("time", "tmpLevel", "latitude", "longitude"))
        temperature.standard_name = "surface_temperature"
        temperature.units = "K"
        temperature[:,:,:,:] = self.SubsetTheData(tValues, timeSize, 37, upLat, lowLat, lowLon, upLon)
        
        # Create U-Wind Component
        uValues, uLevels = self.ExtractData("anl_p125_ugrd.", date)
        uLevel = f.createDimension("uLevel", len(uLevels))
        uWind = f.createVariable("u_component_of_wind", "f4", ("time", "uLevel", "latitude", "longitude"))
        uWind.standard_name = "u_component_of_wind"
        uWind.units = "m/s"
        uWind[:,:,:,:] = self.SubsetTheData(uValues, timeSize, 37, upLat, lowLat, lowLon, upLon)
        
        # Create V-Wind Component
        vValues, vLevels = self.ExtractData("anl_p125_vgrd.", date)
        vLevel = f.createDimension("vLevel", len(vLevels))
        vWind = f.createVariable("v_component_of_wind", "f4", ("time", "vLevel", "latitude", "longitude"))
        vWind.standard_name = "v_component_of_wind"
        vWind.units = "m/s"
        vWind[:,:,:,:] = self.SubsetTheData(vValues, timeSize, 37, upLat, lowLat, lowLon, upLon) 
  
        # Descriptions
        f.description = "Geography example"
        f.history = "Created today"
        f.source = "netCDF4 practice"
        
        f.close()


###################################
# Run the program
t0 = time.time()   
startDay, endDay = JRA_Download().TimeSet() # Get the time period for the data
latLowerPosition, latUpperPosition, lonLowerPosition, lonUpperPosition = JRA_Download().ChooseLatLon()
#startDay = date(2000, 01, 01) # (Y, M, D)
#endDay = date(2000, 01, 05) # (Y, M, D)
  

JRA_Download().ftp_Download(startDay, endDay, "anl_p125", 6, "_hgt") # Download anl_p125 for geopotential height (6 hour increments) 

JRA_Download().ftp_Download(startDay, endDay, "anl_p125", 6, "_tmp") # Download anl_p125 for temperatures (6 hour increments) 

JRA_Download().ftp_Download(startDay, endDay, "anl_p125", 6, "_ugrd") # Download anl_p125 for u-component of wind (6 hour increments)

JRA_Download().ftp_Download(startDay, endDay, "anl_p125", 6, "_vgrd") # Download anl_p125 for v-component of wind (6 hour increments)

JRA_Download().ftp_Download(startDay, endDay, "anl_p125", 6, "_rh") # Download anl_p125 for relative humidity (6 hour increments)

JRA_Download().ftp_Download(startDay, endDay, "anl_surf125", 6, "") # Download anl_surf125 (6 hour increments) 

JRA_Download().ftp_Download(startDay, endDay, "fcst_phy2m125", 3, "") # Download fcst_phy2m125 (3 hour increments) 

print "\nAll Downloads Finished :)"

timeSize = 3
Fdate = []
Adate = []
Idate = []
x = 0
finalDay = str(endDay)
for dt in rrule(DAILY, dtstart = startDay, until = endDay): 
  
  currentDay = str(dt.strftime("%Y") + "-" + dt.strftime("%m") + "-" + dt.strftime("%d"))
  if (x < timeSize):
    x += 1
    Fdate.append(str(dt.strftime("%Y")) + str(dt.strftime("%m")) + str(dt.strftime("%d")))
    Adate.append(str(dt.strftime("%Y")) + str(dt.strftime("%m")) + str(dt.strftime("%d")))
    Idate.append(str(dt.strftime("%Y")) + str(dt.strftime("%m")) + str(dt.strftime("%d")))
    
  if (x == timeSize or currentDay == finalDay):
    fcst_phy2m().Main(Fdate[0], Fdate[x-1], x, Fdate, latLowerPosition, latUpperPosition, lonLowerPosition, lonUpperPosition)
    anl_surf().Main(Adate[0], Adate[x-1], x, Adate, latLowerPosition, latUpperPosition, lonLowerPosition, lonUpperPosition)
    Isobaric().Main(Idate[0], Idate[x-1], x, Idate, latLowerPosition, latUpperPosition, lonLowerPosition, lonUpperPosition)
    x = 0
    Fdate = []
    Adate = []
    Idate = []   


print "\nAll Conversions Finished!"
print "Have a nice day! "  
t1 = time.time()
total = t1 - t0
print total  
