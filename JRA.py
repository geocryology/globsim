#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from datetime        import date, datetime, timedelta
from dateutil.rrule  import rrule, DAILY
from ftplib          import FTP
from netCDF4         import Dataset
from generic         import ParameterIO, StationListRead
from os              import path, listdir
from math            import exp
from fnmatch         import filter

import netCDF4       as nc
import numpy         as np
import pygrib
import time
import sys
import os.path
import shutil

try:
    import ESMF
except ImportError:
    print("*** ESMF not imported, interpolation not possible. ***")
    pass   


"""
Download ranges of data from the JRA-55 server
Data available from 1958 to present date (give 1 or 2 days for upload)
~Author: Christopher Molnar
~Date: October 25, 2017
"""
class JRA_Download:
      
    """
    Take the string containg the YYYY-MM_DD HH:MM:SS and convert it to a datetime.date tuple containg YYYY-MM-DD
    Get the time sets that you want to download (start day to final day)
    Parameters:
    -start: A string containing the start day in YYYY-MM-DD HH:MM:SS
    -end: A string containing the end day in YYYY-MM-DD HH:MM:SS
    Returns:
    -start_date: A tuple containing the start year, month and day
    -end_date: A tuple containing the end year, month and day
    """  
    def TimeSet(self, start, end):
        # Convert start and end date into datetime format 
        try:
            start_data = date(int(start.rsplit("-")[0]), int(start.rsplit("-")[1]), int(start.rsplit("-")[2].rsplit(" ")[0]))
            end_data = date(int(end.rsplit("-")[0]), int(end.rsplit("-")[1]), int(end.rsplit("-")[2].rsplit(" ")[0]))
            first_data = date(1979, 12, 31)
        except:
            print "Invalid start day or end day"
            sys.exit(0)
        # Make sure the day is suitable for our reanalysis       
        try:          
            if (start_data > end_data): # Make sure the start day is before the end day
                print"\nInvalid data entered"
                print "*The end day you chose is before the start day you chose" 
                sys.exit(0)
            elif (start_data < first_data or end_data < first_data):
                print"\nInvalid data entered"
                print "*There is no data available for this time frame"  
                sys.exit(0)
        except:
            print "***Invalid start day or end day"
            sys.exit(0)
    
        return(start_data, end_data)
      
    
    """
    Convert the latitude and longitude restraints into position in a list of JRA coordinates
    Parameter:
    -bottomLat: Bottom latitude coordinate
    -topLat: Top latitude coordinate
    -leftLon: Left longitude coordinate
    -rightLon: Right longitude coordinate
    Returns:
    -latTopPosition: The position of the top latitude position
    -latBottomPosition: The position of the bottom latitude position
    -lonLeftPosition: The position of the left longitude position
    -lonRightPosition: The position of the right longitude position
    """
    def ConvertLatLon(self, bottomLat, topLat, leftLon, rightLon):
        try:          
            latBottomPosition = int(abs(bottomLat - 90) / 1.25 + 1)
            latTopPosition = int (abs(topLat - 90) / 1.25)
            
            if (-90 <= bottomLat <= 90 and -90 <= topLat <= 90 and bottomLat != topLat):
                if (latBottomPosition > 144):
                    latBottomPosition = 144
                run1 = False
            else:
                print "\nThe bounds must be  between -90 and 90 (inclusively)"
                sys.exit(0)  
        except:
            print "\nInvalid Entry"
            print "Please make sure your latitude position is between -90 and 90 (inclusively)"
            sys.exit(0)
            
        try:
            lonLeftPosition = int(leftLon / 1.25)
            lonRightPosition = int(rightLon / 1.25 + 1)
          
            if (0 <= leftLon <= 358.75 and 0 <= rightLon <= 358.75 and leftLon != rightLon):
                if (lonRightPosition > 287):
                    lonRightPosition = 287
                run2 = False
            else:
                print "\nThe bounds must be  between 0 and 358.75 (inclusively)"
                sys.exit(0)
        except:
            print "\nInvalid Entry"
            print "Please make sure your longitude position is between 0 and 358.75 (inclusively)"
            sys.exit(0)
    
        return (latTopPosition, latBottomPosition, lonLeftPosition, lonRightPosition) # Return latTopPosition before latBottomPosition because JRA goes from 90 to -90   
    
    
    """
    Convert elevation into air pressure using barometric formula
    Paramters:
    -elevation: the height in m 
    Returns:
    -The elevation hight in air pressure
    """
    def getPressure(self, elevation):
        g  = 9.80665   #Gravitational acceleration [m/s2]
        R  = 8.31432   #Universal gas constant for air [N.m /(mol.K)]    
        M  = 0.0289644 #Molar mass of Earth's air [kg/mol]
        P0 = 101325    #Pressure at sea level [Pa]
        T0 = 288.15    #Temperature at sea level [K]
        #http://en.wikipedia.org/wiki/Barometric_formula
        return P0 * exp((-g * M * elevation) / (R * T0)) / 100 #[hPa] or [bar]
    
    """
    Takes the minimum and maximum elevation and finds their position in the total_elevation list
    Parameters:
    -elevationMin: An int of the minimum elevation
    -elevationMax: An int of the maximum elevation
    Returns:
    -elevationMinRange: The position of the minimum elevation from the total_elevation list
    -elevationMaxRange: The position of the maximum elevation from the total_elevation list
    """
    def ElevationCalculator(self, eMin, eMax):
        total_elevations = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
        
        # flip max and min because 1000 is the bottom and 0 is the top
        elevationMax = self.getPressure(eMin)
        elevationMin = self.getPressure(eMax)
        
        minNum = min(total_elevations, key=lambda x:abs(x-elevationMin))
        maxNum = min(total_elevations, key=lambda x:abs(x-elevationMax))
        
        if (minNum > elevationMin and total_elevations.index(minNum) > 0 ):
        	elevationMinRange = total_elevations.index(minNum) - 1
        else:
        	elevationMinRange = total_elevations.index(minNum)
        	
        if (maxNum < elevationMin and total_elevations.index(maxNum) < 36 ):
        	elevationMaxnRange = total_elevations.index(maxNum) - 1
        else:
        	elevationMaxRange = total_elevations.index(maxNum)
        
        return (elevationMinRange, elevationMaxRange)
    
    
    """
    Download the files from the server
    Parameters:
    -start_day: The first day to start downloading from
    -end_day: The last day that you want to download
    -download_list: A list of all the variables that need to be donwload, with data_type and hour increment
    -savePath: The directory to save the GRIB files from the JRA server
    -ftp: The ftp connection for the JRA website
    """
    def ftp_Download(self, start_day, end_day, download_list, savePath, ftp):   
          
        print "\nDownloading GRIB Files....."  
        for dt in rrule(DAILY, dtstart = start_day, until = end_day): # Loop from start day to end day
            for var in range(0,len(download_list)):
                path = "/JRA-55/Hist/Daily/" + download_list[var][0] + "/" + dt.strftime("%Y") + dt.strftime("%m") # Generate the path  
                ftp.cwd(path) # Change Working Directory
              
                for x in range(0, 24, download_list[var][2]): # Loop through the hourly increments
                    if (x <10): # Add a 0 infront any of the hour times that are less 10
                        ending = "0" + str(x)
                    else:
                        ending = str(x)  
                       
                    filename = download_list[var][0] + download_list[var][1] + "." +  dt.strftime("%Y") + dt.strftime("%m") + dt.strftime("%d") + ending
                    
                    try:
                        completeName = os.path.join(savePath, filename) # Generate the filename and save it to the proper folder
                    except:
                        print "Make sure you have a Grib and netCDF folder in your directory"
                        sys.exit(0)
               
                    localfile = open(completeName, 'wb')
               
                    try: # try to download the file
                        ftp.retrbinary("RETR %s" % filename , localfile.write) # Download file
                    except:
                        print "Error downloading file: " + filename
                    
                    localfile.close() # Close File
                     
            print "Downloaded all the data for:", str(dt.strftime("%Y")) + "-" + str(dt.strftime("%m")) + "-" + str(dt.strftime("%d"))  
        print "\nAll Downloads Finished :) \n"
            
                         
    """
    Connect to server and download the necessary Grib Files
    Parameters:
    -username: The username for the JRA server
    -password: The password for the JRA sever
    -startDay: The starting download day
    -endDay: The ending download day
    -save_path: The save path
    -variable_data: All of the data that we want to download
    -fcst_list: The forcast data that we need to download
    -surf_list: The surf data that we need to download
    """  
    def DownloadGribFile(self, username, password, startDay, endDay, save_path, isobaric_list, fcst_list, surf_list):
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
            
                ftp.login(username, password) # Login (username, password)
                
                download_list = []

                if ("geopotential_height" in isobaric_list): # Add hgt info
                    download_list.append(["anl_p125", "_hgt", 6])
                if ("air_temperature" in isobaric_list): # Add tmp info
                    download_list.append(["anl_p125", "_tmp", 6])
                if ("eastward_wind" in isobaric_list): # Add ugrd info
                    download_list.append(["anl_p125", "_ugrd", 6])
                if ("northward_wind" in isobaric_list): # Add vgrd info
                    download_list.append(["anl_p125", "_vgrd", 6]) 
                if ("relative_humidity" in isobaric_list):  # Add RH info
                    download_list.append(["anl_p125", "_rh", 6])
                if (len(surf_list) > 0): # Add anl_surf info
                    download_list.append(["anl_surf125", "", 6])
                if (len(fcst_list) > 0):  # Add fcst_phy2m info
                    download_list.append(["fcst_phy2m125", "", 3])
            
                if (len(download_list) > 0):
                    # Download all the files needed for convertion
                    self.ftp_Download(startDay, endDay, download_list, save_path + 'Grib', ftp) 
                    
                ftp.quit() # Close Connection   
              
            else:     
                print "\nAttempted to download" + data_type + middle
                print "Tried to connect 5 times but failed"   
                print "Please retry later"
                sys.exit(0)   


"""
Convert the grib files downloaded from the JRA-55 into netCDF format
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
    Set up the needed lattitude and longitude range
    Parameters:
    -latLow: Lower bound of the latitude
    -latHigh: Upper bound of the latitude
    -lonLow: Lower bound of the longitude
    -lonHigh: Upper bound of the longitude
    Returns:
    -lat: Latitude range
    -lon: Longitude range
    """
    def GetLatLon(self, latLow, latHigh, lonLow, lonHigh):
        
        lat = []
        lon = []
        for x in range(latLow, latHigh + 1):
            lat.append( 90 - (x * 1.25))
        for y in range(lonLow, lonHigh + 1):
            lon.append(y * 1.25)
        
        return (lat, lon)
      
    
    """
    Extracts the needed values from a Grib file
    Parameters:
    -filename: The name of the Grib file that the values are being extracted from
    -date: An array with all of the dates in it
    -dataName: The name of the variable that you want the values for
    -filePath: The direcory where you will find the folder with the GRIB files
    Returns:
    -value: An array with all of the values
    """
    def ExtractData(self, filename, date, dataName, filePath):    
        value = []    
        fileLocation = filePath + 'Grib/' 
        for z in range(0, len(date)):
            try:
                name =(str(date[z].strftime("%Y")) + str(date[z].strftime("%m")) + str(date[z].strftime("%d")) + str(date[z].strftime("%H")))
                grbs = pygrib.open(fileLocation + filename + name)
            except:
                print "file: " + filename + name +  " not found :("
                print "Quiting program"
                sys.exit(0)
                        
            # Loop through the Grib file and extract the data you need
            for g in grbs:
                if (g.shortName == dataName):
                    value.append(g.values)
                    break 
            grbs.close()
          
        return (value)


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
    def CreateLatLon(self, f, latSize, lonSize):
        # Create Dimensions for variables 
        latitude = f.createDimension("latitude", latSize + 1)
        longitude = f.createDimension("longitude", lonSize + 1)
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
        
        days = []
        
        for d in range(0, len(date)):
          for h in range(0, 24, hours):
              tempDate = date[d] + timedelta(hours = h)
              days.append(tempDate)
        
        return days

   
    """
    Creates the elevation series for the netCDF file
    Parameters:
    -elevationMin: The minimum elevation
    -elevationMax: The maximum elevation
    Returns:
    -elevation: A list of all the elvations
    """
    def ElevationSeries(self, elevationMin, elevationMax): 
        total_elevations = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
        elevation = []
        for e in range(elevationMin, elevationMax + 1):
            elevation.append(total_elevations[e])
                      
        return(elevation)
        
    
    """
    Subsets the data in to the desired Latitude and Longitude coordinates
    Parameters:
    -dataValues: The array containing the reanlysis data
    -numDays: The total amount of time segments in the netCDF file
    -latBottom: The bottom latitude coordinate
    -latTop: The top latitude coordinate
    -lonLeft: The left longitude coordinate
    -lonRight: The right longitude coordinate
    Returns:
    -ssDay: The subseted data set
    """
    def SubsetTheData(self, dataValues, numDays, latBottom, latTop, lonLeft, lonRight):
        ssDay = []
        for a in range(0,numDays):
          ssDay.append([])
          for b in range(latBottom, latTop + 1):
            ssDay[a].append([])
            for c in range(lonLeft, lonRight + 1):
              ssDay[a][b - latBottom].append(dataValues[a][b][c])

        return (ssDay)
        
    """
    Deletes the old Grib folder (removing all of the already used GRIB files) and
    replaces it with a new empty GRIB folder
    Parameters:
    -location: The directory where the GRIB folder is
    """
    def EmptyFolder(self, location):
        "Try to remove the GRIB folder"
        newlocation = location + "Grib"
        try:
            shutil.rmtree(newlocation) 
        except:
            print "Unable to delete the used Grib files"
            sys.exit(0)
        try:
            if not os.path.exists(newlocation):
                os.makedirs(newlocation)
        except:
            print "Unable to create the new Grib folder"
            sys.exit(0)

   
"""
Class for Forcast data
""" 
class fcst_phy2m:
    
    """
    The main function of the fcst_phy2m Class
    Parameters:
    -startDay: The first day for the data
    -endDay: The last day for the data
    -numDays: The total number of days
    -date: A list with all the days in it
    -bottomLat: The bottom latitude coordinate
    -topLat: The top latitude coordinate
    -leftLon: The left longitude coordinate
    -rightLon: The right longitude coordinate
    -savePath: The directory with the GRIB and netCDF folders 
    -JRA_Dictionary: A dictionary containing all of the JRA variables with there short_name, file, units, and elevations
    -fcst_data: A list of all the fcst variables that need to be downloaded
    """
    def Main(self, startDay, endDay, numDays, date, bottomLat, topLat, leftLon, rightLon, savePath, JRA_Dictionary, fcst_data):    
        
        #startName = (str(startDay.strftime("%Y")) + str(startDay.strftime("%m")) + str(startDay.strftime("%d")) + str(startDay.strftime("%H")))
        #endName = (str(endDay.strftime("%Y")) + str(endDay.strftime("%m")) + str(endDay.strftime("%d")) + str(endDay.strftime("%H")))
        
        startName = (str(startDay.strftime("%Y")) + str(startDay.strftime("%m")) + str(startDay.strftime("%d")))
        endName = (str(endDay.strftime("%Y")) + str(endDay.strftime("%m")) + str(endDay.strftime("%d")))
        
        try:
            completeName = os.path.join(savePath + 'jra55', "JRA_fcst_" + startName + "_" + endName + ".nc")
            f = Dataset(completeName, "w", format = "NETCDF4_CLASSIC") # Name of the netCDF being created  
        except:
            print "Make sure you have a netCDF folder in your directory"
            print "Once the netCDF files are created they will be stored in there"
            sys.exit(0)
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon(bottomLat, topLat, leftLon, rightLon)
        latitude, longitude, latitudes, longitudes = Grib2CDF().CreateLatLon(f, topLat - bottomLat, rightLon - leftLon)
        latitudes[:] = lats
        longitudes[:] = lons
        
        # Create the Time
        timeSize = numDays * 8
        dateTimes = Grib2CDF().TimeSeries(date, timeSize, 3, 8)
        time = f.createDimension("time", timeSize)
        time = f.createVariable("time", "i4", "time")
        time.standard_name = "time"
        time.units =  "hours since 1900-01-01 00:00:00"
        time.calendar = "gregorian"
        
        # Convert datetime object using the netCDF4 date2num function
        t = []
        for tt in range(0,len(dateTimes)):
          t.append(nc.date2num(dateTimes[tt], units = time.units, calendar = time.calendar))  
        time[:] = t

        for dataName in fcst_data: # Loop through all the needed varibales and make netCDF variables 
            data = Grib2CDF().ExtractData(JRA_Dictionary[dataName][1], dateTimes, JRA_Dictionary[dataName][0], savePath)

            dataVariable = f.createVariable(dataName, "f4", ("time", "latitude", "longitude"))
            dataVariable.standard_name = dataName
            dataVariable.units = JRA_Dictionary[dataName][3]
            dataVariable[:,:,:] = Grib2CDF().SubsetTheData(data, timeSize, bottomLat, topLat, leftLon, rightLon)
            print "Converted:", dataName

        # Description
        f.source = "JRA converted data"
 
        f.close() 


"""
Class for Surface data
"""
class anl_surf:
    
    """
    The main function of the anl_surf Class
    Parameters:
    -startDay: The first day for the data
    -endDay: The last day for the data
    -numDays: The total number of days
    -date: A list with all the days in it
    -bottomLat: The bottom latitude coordinate
    -topLat: The top latitude coordinate
    -leftLon: The left longitude coordinate
    -rightLon: The right longitude coordinate
    -savePath: The directory with the GRIB and netCDF folders 
    -JRA_Dictionary: A dictionary containing all of the JRA variables with there short_name, file, units, and elevations
    -surf_data: A list of all the surf variables that need to be downloaded
    """
    def Main(self, startDay, endDay, numDays, date, bottomLat, topLat, leftLon, rightLon, savePath, JRA_Dictionary, surf_data):   
        
        # startName = (str(startDay.strftime("%Y")) + str(startDay.strftime("%m")) + str(startDay.strftime("%d")) + str(startDay.strftime("%H")))
        # endName = (str(endDay.strftime("%Y")) + str(endDay.strftime("%m")) + str(endDay.strftime("%d")) + str(endDay.strftime("%H")))
        
        startName = (str(startDay.strftime("%Y")) + str(startDay.strftime("%m")) + str(startDay.strftime("%d")))
        endName = (str(endDay.strftime("%Y")) + str(endDay.strftime("%m")) + str(endDay.strftime("%d")))

        try:
            completeName = os.path.join(savePath + 'jra55', "JRA_surf_" + startName + "_" + endName + ".nc") 
            f = Dataset(completeName, "w", format = "NETCDF4_CLASSIC") # Name of the netCDF being created 
        except:
            print "Make sure you have a netCDF folder in your directory"
            print "Once the netCDF files are created they will be stored in there"
            sys.exit(0)
        
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon(bottomLat, topLat, leftLon, rightLon)
        latitude, longitude, latitudes, longitudes = Grib2CDF().CreateLatLon(f, topLat - bottomLat, rightLon - leftLon)
        latitudes[:] = lats
        longitudes[:] = lons

        # Create the Time
        timeSize = numDays * 4  
        dateTimes = Grib2CDF().TimeSeries(date, timeSize, 6, 4)
        time = f.createDimension("time", timeSize)
        time = f.createVariable("time", "i4", "time")
        time.standard_name = "time"
        time.units =  "hours since 1900-01-01 00:00:00"
        time.calendar = "gregorian"
        
        # Convert datetime object using the netCDF4 date2num function
        t = []
        for tt in range(0,len(dateTimes)):
          t.append(nc.date2num(dateTimes[tt], units = time.units, calendar = time.calendar))  
        time[:] = t
        
        for dataName in surf_data: # Loop through all the needed varibales and make netCDF variables 
            data = Grib2CDF().ExtractData(JRA_Dictionary[dataName][1], dateTimes, JRA_Dictionary[dataName][0], savePath)

            dataVariable = f.createVariable(dataName, "f4", ("time", "latitude", "longitude"))
            dataVariable.standard_name = dataName
            dataVariable.units = JRA_Dictionary[dataName][3]
            dataVariable[:,:,:] = Grib2CDF().SubsetTheData(data, timeSize, bottomLat, topLat, leftLon, rightLon)
            print "Converted:", dataName
        
        # Description
        f.source = "JRA converted data"
        
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
    -savePath: The directory where the GRIB files are stored
    Returns:
    -allData: An array with all of the values
    -levels: An array with all of the datas level
    """
    def ExtractData(self, filename, date, savePath):    
        allData = []
        fileLocation = savePath + 'Grib/'
        for z in range(0, len(date)):
            levels = []
            data = []
            try:
                name =(str(date[z].strftime("%Y")) + str(date[z].strftime("%m")) + str(date[z].strftime("%d")) + str(date[z].strftime("%H")))
                grbs = pygrib.open(fileLocation + filename + name)
            except:
                print "file: " + filename + name +  " not found :("
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
    -latBottom: The bottom latitude coordinate
    -latTop: The top latitude coordinate
    -lonLeft: The left longitude coordinate
    -lonRight: The right longitude coordinate
    -elevationMinRange: The minimum elevation needed to be downloaded
    -elevationMaxRange: The maximum elevation needed to be downloaded
    Returns:
    -ssDay: The subseted data set
    """
    def SubsetTheData(self, dataValues, numDays, latBottom, latTop, lonLeft, lonRight, elevationMinRange, elevationMaxRange):
        ssDay = []
        for a in range(0,numDays):
          ssDay.append([])
          for b in range(elevationMinRange, elevationMaxRange + 1): 
            ssDay[a].append([])
            for c in range(latBottom, latTop + 1):
              ssDay[a][b - elevationMinRange].append([])
              for d in range(lonLeft, lonRight + 1):
                ssDay[a][b - elevationMinRange][c - latBottom].append(dataValues[a][b][c][d])      
        return (ssDay)
       
       
    """
    The main function of the Isobaric Class
    Parameters:
    -startDay: The first day for the data
    -endDay: The last day for the data
    -numDays: The total number of days
    -date: A list with all the days in it
    -bottomLat: The bottom latitude coordinate
    -topLat: The top latitude coordinate
    -leftLon: The left longitude coordinate
    -rightLon: The right longitude coordinate
    -savePath: The directory with the GRIB and netCDF folders 
    -JRA_Dictionary: A dictionary containing all of the JRA variables with there short_name, file, units, and elevations
    -isobaric_data: A list of all the isobaric variables that need to be downloaded
    -elevationMinRange: The lower bound of the elevation
    -elevationMaxRange: The upper bound of the elevation
    """
    def Main(self, startDay, endDay, numDays, date, bottomLat, topLat, leftLon, rightLon, savePath, JRA_Dictionary, isobaric_data, elevationMinRange, elevationMaxRange):
        
        # startName = (str(startDay.strftime("%Y")) + str(startDay.strftime("%m")) + str(startDay.strftime("%d")) + str(startDay.strftime("%H")))
        # endName = (str(endDay.strftime("%Y")) + str(endDay.strftime("%m")) + str(endDay.strftime("%d")) + str(endDay.strftime("%H")))
        
        startName = (str(startDay.strftime("%Y")) + str(startDay.strftime("%m")) + str(startDay.strftime("%d")))
        endName = (str(endDay.strftime("%Y")) + str(endDay.strftime("%m")) + str(endDay.strftime("%d")))

        
        try:
            completeName = os.path.join(savePath + 'jra55', "JRA_Isobaric_" + startName + "_" + endName + ".nc") 
            f = Dataset(completeName, "w", format = "NETCDF4_CLASSIC") # Name of the netCDF being created
        except:
            print "Make sure you have a netCDF folder in your directory"
            print "Once the netCDF files are created they will be stored in there"
            sys.exit(0) 
        
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon(bottomLat, topLat, leftLon, rightLon)
        latitude, longitude, latitudes, longitudes = Grib2CDF().CreateLatLon(f, topLat - bottomLat, rightLon - leftLon)
        latitudes[:] = lats
        longitudes[:] = lons
                
        # Create the Time 
        timeSize = numDays * 4
        dateTimes = Grib2CDF().TimeSeries(date, timeSize, 6, 4)
        time = f.createDimension("time", timeSize)
        time = f.createVariable("time", "i4", "time")
        time.standard_name = "time"
        time.units =  "hours since 1900-01-01 00:00:00"
        time.calendar = "gregorian"
       
        # Convert datetime object using the netCDF4 date2num function
        t = []
        for tt in range(0,len(dateTimes)):
          t.append(nc.date2num(dateTimes[tt], units = time.units, calendar = time.calendar))
          
        time[:] = t
        
        x=0
        for dataName in isobaric_data: # Loop through all the needed varibales and make netCDF variables 
            data, level = self.ExtractData(JRA_Dictionary[dataName][1], dateTimes, savePath)
            
            # *Special Case: Relative Humidity only has 27 levels (if data from 1 - 99 hPa is needed make a seperate level for relative humidity 
            if (dataName == "relative_humidity" and elevationMinRange <= 9):
                if (elevationMinRange < 9):
                    tempMinRange = 0
                else:
                    tempMinRange = elevationMinRange - 9
                Levels = level
                Levels =f.createDimension("Levels", elevationMaxRange - 9 - tempMinRange)
                Levels = f.createVariable("Levels", "i4", "Levels")
                Levels.long_name = "pressure_level"
                Levels.units = "mbar"   
                Levels[:] = Grib2CDF().ElevationSeries(tempMinRange + 10, elevationMaxRange)
                dataVariable = f.createVariable(dataName, "f4", ("time", "Levels", "latitude", "longitude"))
                dataVariable.standard_name = dataName
                dataVariable.units = JRA_Dictionary[dataName][3]
                dataVariable[:,:,:,:] = self.SubsetTheData(data, timeSize, bottomLat, topLat, leftLon,  rightLon, tempMinRange, elevationMaxRange - 10)
            
            else:
                if (x==0):
                    levels = f.createDimension("level", elevationMaxRange - elevationMinRange + 1)
                    levels = f.createVariable("level", "i4", "level")
                    levels.long_name = "pressure_level"
                    levels.units = "mbar"
                    levels[:] = Grib2CDF().ElevationSeries(elevationMinRange, elevationMaxRange)
                    x =1
                if (dataName == "relative_humidity"):
                    tempMinRange = elevationMinRange - 10 
                    tempMaxRange = elevationMaxRange - 10 
                else:
                    tempMinRange = elevationMinRange
                    tempMaxRange = elevationMaxRange
                    
                dataVariable = f.createVariable(dataName, "f4", ("time","level", "latitude", "longitude"))
                dataVariable.standard_name = dataName
                dataVariable.units = JRA_Dictionary[dataName][3]
                dataVariable[:,:,:,:] = self.SubsetTheData(data, timeSize, bottomLat, topLat, leftLon,  rightLon, tempMinRange, tempMaxRange)
            print "Converted:", dataName
        
        # Description
        f.source = "JRA converted data"
        
        f.close()


"""
*From Dr. Grubers ERA file
Class for accessing the parameter file for downloading specified variables, latitude and longitude coordinates, start to end date and the chunk size
"""
class JRAdownload(object):

    """
    Initialize the JRAdownload class
    Parameters:
    -pfile: Full path to a Globsim Download Parameter file. 
    """
    def __init__(self, pfile):
        # read parameter file
        self.pfile = pfile
        par = ParameterIO(self.pfile)
        
        # assign bounding box
        self.area  = {'north':  par.bbN,
                      'south':  par.bbS,
                      'west' :  par.bbW,
                      'east' :  par.bbE}
                 
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}

        # elevation
        self.elevation = {'min' : par.ele_min, 
                          'max' : par.ele_max}
        
        # data directory for ERA-Interim  
        self.directory = par.project_directory
        #self.directory = path.join(par.project_directory, "eraint")  
        #if path.isdir(self.directory) == False:
            #raise ValueError("Directory does not exist: " + self.directory)
            
        self.credential = path.join(par.credentials_directory, ".jrarc")
        #print self.credential
        self.account = open(self.credential, "r")
        self.inf = self.account.readlines()
        self.username = ''.join(self.inf[0].split())
        self.password = ''.join(self.inf[1].split()) 
     
        # variables
        self.variables = par.variables
            
        # chunk size for downloading and storing data [days]        
        self.chunk_size = par.chunk_size           


    """
    Run the JRA program
    -Call the JRAdownload class to build a jra object with the parameters(area, date, elevation, directory, varianles, and chunk_size)
    -Convert the area data into latitude and longitude positions
    -Get the chunk size
    -Convert the start and end date into datetime.date tuples containing the YYYY-MM-DD
    -Convert the min and max elevation into there positions in the elevation_list
    -Get the directory
    -Extract the needed variables
    -Download the needed GRIB files
    -Use all of the retrieved information to convert the downloaded GRIB files into netCDF with the proper area resitrictions, dates, elevation, and chunk_size
    """
    """
    Start the reanalysis
    Parameter:
    -pfile: The parameter file
    """
    def retrieve(self):  #pfile
        # Run the program
        t0 = time.time()  
        
        #jra = JRAdownload("Parameter_Stuff.txt")
        #jra = JRAdownload(pfile)
        
        # Area data 
        try:
            latBottom = float(self.area["south"])
            latTop = float(self.area["north"])
            lonLeft = float(self.area["west"])
            lonRight = float(self.area["east"])
            
            if (lonLeft < 0):
                lonLeft = 360.0 + lonLeft
                
            if (lonRight < 0):
                lonRight = 360.0 + lonRight
                
            latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition = JRA_Download().ConvertLatLon(latBottom, latTop, lonLeft,lonRight)
        except:
            print "Invalid area format"
            sys.exit(0)
        
        # Chunk data
        try:
            chunk_size = int(self.chunk_size)
        except:
            print "Invalid chunk size"
            sys.exit(0)
        
        # Date data
        try:
            start = str(self.date["beg"])
            end = str(self.date["end"])
            startDay, endDay = JRA_Download().TimeSet(start, end) # Get the time period for the data
        except:
            print "Invalid date"
            print "Please make sure your date is in YYYY-MM-DD format"
            sys.exit(0)
        
        # Elevation data 
        try:
            elevationMin = self.elevation["min"]
            elevationMax = self.elevation["max"]
            elevationMinRange, elevationMaxRange = JRA_Download().ElevationCalculator(elevationMin, elevationMax) 
        except:
            print "Invalid elevation entered"
            sys.exit(0)
    
        # Directory Information
        directory = self.directory
        save_path = directory
        #save_path = '/home/cmolnar/FinishedCode/'
        
        # Get username and password
        username = self.username
        password = self.password
        
        # Create Grib and netCDF folders if necessary
        gribFolder = save_path + "Grib"
        netFolder = save_path + "jra55"
        try:
            if not os.path.exists(gribFolder):
                os.makedirs(gribFolder)
        except:
            print "Unable to create a Grib folder"
            sys.exit(0)
            
        try:
            if not os.path.exists(netFolder):
                os.makedirs(netFolder)
        except:
            print "Unable to create a netCDF folder"
            sys.exit(0)
        ######### TO FINISH
        
        # List of Variables
        variables = self.variables
        
        shared_data = {
                      "air_temperature"                                      : ["air_temperature", "surface_temperature", "geopotential_height"],
                      "relative_humidity"                                    : ["relative_humidity", "geopotential_height"],
                      "precipitation_amount"                                 : ["total_precipitation"],
                      "downwelling_longwave_flux_in_air"                     : ["downwelling_longwave_flux_in_air"],
                      "downwelling_longwave_flux_in_air_assuming_clear_sky"  : ["downwelling_longwave_flux_in_air_assuming_clear_sky"],
                      "downwelling_shortwave_flux_in_air"                    : ["downwelling_shortwave_flux_in_air" ],
                      "downwelling_shortwave_flux_in_air_assuming_clear_sky" : ["downwelling_shortwave_flux_in_air_assuming_clear_sky"],
                      "wind_from_direction"                                  : ["northward_wind", "eastward_wind","geopotential_height"],
                      "wind_speed"                                           : ["northward_wind", "eastward_wind","geopotential_height"],
                      "geopotential_height"                                  : ["geopotential_height"]
                      }
                      
        variable_data = []
        for x in variables:
            if (x in shared_data):
                for y in range(0, len(shared_data[x])):
                    variable_data.append(shared_data[x][y])
        
        # A dictionary for each file with all the variables available for donwloading with there standard name as the key and the values being a list of short-name, filename, number of levels and units                
        fcst_dictionary = {
                          "precipitation_amount"                                   : ["tpratsfc", "fcst_phy2m125.", 1, "mm/day"],
                          "downwelling_shortwave_flux_in_air_assuming_clear_sky"  : ["csdsf", "fcst_phy2m125.", 1, "W/(m^2)"],
                          "downwelling_longwave_flux_in_air_assuming_clear_sky"   : ["csdlf", "fcst_phy2m125.", 1, "W/(m^2)"],
                          "downwelling_shortwave_flux_in_air"                     : ["dswrf", "fcst_phy2m125.", 1, "W/(m^2)"],
                          "downwelling_longwave_flux_in_air"                      : ["dlwrf", "fcst_phy2m125.", 1, "W/(m^2)"]
                          }  
                          
        surf_dictionary = {
                          "surface_temperature"  : ["t", "anl_surf125.", 1, "K"],
                          "relative_humidity"    : ["r", "anl_surf125.", 1, "%"],
                          "eastward_wind"        : ["u", "anl_surf125.", 1, "m/s"],
                          "northward_wind"       : ["v", "anl_surf125.", 1, "m/s"]
                          }
      
        isobaric_dictionary = { 
                              "geopotential_height"  : ["gh", "anl_p125_hgt.", 37, "gpm"],
                              "air_temperature"     : ["t", "anl_p125_tmp.", 37, "K"],
                              "eastward_wind"        : ["u", "anl_p125_ugrd.", 37, "m/s"],
                              "northward_wind"       : ["v", "anl_p125_vgrd.", 37, "m/s"],
                              "relative_humidity"    : ["r", "anl_p125_rh.", 27, "%"]
                              }
        
        # Check to see which variables data need to be downloaded for fcst, surf, and isobaric list
        fcst_list = list(set(variable_data) & set(fcst_dictionary))
        surf_list = list(set(variable_data) & set(surf_dictionary))
        isobaric_list = list(set(variable_data) & set(isobaric_dictionary))
        
        timeSize = chunk_size # The number of days you want saved together
        Fdate = []
        Adate = []
        Idate = []
        x = 0
        finalDay = str(endDay)
        
        for dt in rrule(DAILY, dtstart = startDay, until = endDay): # Loop through the days to create netCDF files
            
            currentDay = str(dt.strftime("%Y") + "-" + dt.strftime("%m") + "-" + dt.strftime("%d"))
            if (x < timeSize): # If x is less than timeSize append the current day to the date lists
                x += 1
                Fdate.append(dt)
                Adate.append(dt)
                Idate.append(dt)
              
            if (x == timeSize or currentDay == finalDay): # If time size reached or last day reached build netCDF files
                JRA_Download().DownloadGribFile(username, password, Fdate[0], Fdate[x-1], save_path, isobaric_list, fcst_list, surf_list)
                fcst_phy2m().Main(Fdate[0], Fdate[x-1] + timedelta(hours = 21), x, Fdate, latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition, save_path, fcst_dictionary, fcst_list)
                anl_surf().Main(Adate[0], Adate[x-1] + timedelta(hours = 18), x, Adate, latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition, save_path, surf_dictionary, surf_list)
                Isobaric().Main(Idate[0], Idate[x-1] + timedelta(hours = 18), x, Idate, latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition, save_path, isobaric_dictionary, isobaric_list, elevationMinRange, elevationMaxRange)
                x = 0
                Fdate = []
                Adate = []
                Idate = []  
                
                Grib2CDF().EmptyFolder(save_path)                 
        
        # Empty out the GRIB files remaining in the folder
        Grib2CDF().EmptyFolder(save_path) 
        print "\nAll Conversions Finished!"
        print "Have a nice day! "  
        t1 = time.time()
        total = t1 - t0
        print "Total run time:", total  
        
        # Write time to textfile
        file = open("time.txt", "w")
        file.write("Time: " + str(total))
        file.close

class JRAinterpolate(object):
    """
    Algorithms to interpolate JRA55 netCDF files to station coordinates. 
    All variables retains their original units and time-steps. 
    
    Referenced from era_interim.py (Dr.Stephan Gruber): Class ERAinterpolate()     
    
    Args:
        ifile: Full path to a Globsim Interpolate Paramter file
        JRAinterpolate(ifile)


    Example:
        ifile = '/home/xquan/src/globsim/examples/par/examples.globsim_interpolate'
        JRAinterpolate(ifile)
      
    """

    def __init__(self, ifile):
        #read parameter file
        self.ifile = ifile
        par = ParameterIO(self.ifile)
        self.dir_inp = path.join(par.project_directory,'jra55') 
        self.dir_out = path.join(par.project_directory,'station')
        self.variables = par.variables
        self.list_name = par.list_name
        self.stations_csv = path.join(par.project_directory,
                                      'par', par.station_list)
        
        #read station points 
        self.stations = StationListRead(self.stations_csv)  
        
        # time bounds
        self.date  = {'beg' : par.beg,
                      'end' : par.end}

    def JRA2station(self, ncfile_in, ncfile_out, points,
                    variables=None, date=None):    
        """
        Biliner interpolation from fields on regular grid (latitude, longitude) 
        to individual point stations (latitude, longitude). This works for
        surface and for pressure level files (all jra55 files).
          
        Args:
            ncfile_in: Full path to an JRA-55 derived netCDF file. This can
                        contain wildcards to point to multiple files if temporal
                        chunking was used.
              
            ncfile_out: Full path to the output netCDF file to write.  
              
            points: A dictionary of locations. See method StationListRead in
                    generic.py for more details.
        
            variables:  List of variable(s) to interpolate such as 
                        [air_temperature, easteard_wind, northward_wind, relative_humidy, surface_temperature, 
                        downwelling_shortwave_flux_in_air, downwelling_longwave_flux_in_air,
                        downwelling_shortwave_flux_in_air_assuming_clear_sky, 
                        downwelling_longwave_flux_in_air_assuming_clear_sky].
                        Defaults to using all variables available.
        
            date: Directory to specify begin and end time for the derived time 
                  series. Defaluts to using all times available in ncfile_in.
              
        Example:
            from datetime import datetime
            date  = {'beg' : datetime(2008, 1, 1),
                      'end' : datetime(2008,12,31)}
            variables  = [air_temperature, easteard_wind, northward_wind]       
            stations = StationListRead("points.csv")      
            MERRA2station('jra_surf_*.nc', 'jra_sa_inter.nc', stations, 
                        variables=variables, date=date)        
        """   

        # open netcdf file handle, can be one file of several with wildcards
        ncf = nc.MFDataset(ncfile_in, 'r', aggdim ='time') 
        
        # is it a file with pressure levels?
        pl = 'level' in ncf.dimensions.keys()

        # get spatial dimensions
        lat  = ncf.variables['latitude'][:]
        lon  = ncf.variables['longitude'][:]
        if pl: # only for pressure level files
            lev  = ncf.variables['level'][:]
            nlev = len(lev)
    
        # get time and convert to datetime object
        nctime = ncf.variables['time'][:]
        t_unit = ncf.variables['time'].units #"hours since 1900-01-01 00:00:0.0"
        try :
            t_cal = ncf.variables['time'].calendar
        except AttributeError : # Attribute doesn't exist
            t_cal = u"gregorian" # or standard
        time = nc.num2date(nctime, units = t_unit, calendar = t_cal)
        
        # restrict to date/time range if given
        if date is None:
            tmask = time < datetime(3000, 1, 1)
        else:
            tmask = (time <= date['end']) * (time >= date['beg'])
          
        # test if time steps to interpolate remain
        nt = sum(tmask)
        if nt == 0:
            raise ValueError('No time steps from netCDF file selected.')
    
        # get variables
        varlist = [x.encode('UTF8') for x in ncf.variables.keys()]
        varlist.remove('time')
        varlist.remove('latitude')
        varlist.remove('longitude')
        if pl: #only for pressure level files
            varlist.remove('level')
    
        #list variables that should be interpolated
        if variables is None:
            variables = varlist
        #test is variables given are available in file
        if (set(variables) < set(varlist) == 0):
            raise ValueError('One or more variables not in netCDF file.')
        
        # Create source grid from a SCRIP formatted file. As ESMF needs one
        # file rather than an MFDataset, give first file in directory.
        ncsingle = filter(listdir(self.dir_inp), path.basename(ncfile_in))[0]
        ncsingle = path.join(self.dir_inp, ncsingle)
        sgrid = ESMF.Grid(filename=ncsingle, filetype=ESMF.FileFormat.GRIDSPEC)

        # create source field on source grid
        if pl: #only for pressure level files
            sfield = ESMF.Field(sgrid, name='sgrid',
                                staggerloc=ESMF.StaggerLoc.CENTER,
                                ndbounds=[len(variables), nt, nlev])
        else: # 2D files
            sfield = ESMF.Field(sgrid, name='sgrid',
                                staggerloc=ESMF.StaggerLoc.CENTER,
                                ndbounds=[len(variables), nt])
                                                 
        # assign data from ncdf: (variable, time, latitude, longitude) 
        for n, var in enumerate(variables):
            if pl: # only for pressure level files
                sfield.data[n,:,:,:,:] = ncf.variables[var][tmask,:,:,:].transpose((0,1,3,2)) 
            else:
                sfield.data[n,:,:,:] = ncf.variables[var][tmask,:,:].transpose((0,2,1))

        # create locstream, CANNOT have third dimension!!!
        locstream = ESMF.LocStream(len(self.stations), coord_sys=ESMF.CoordSys.SPH_DEG)
        locstream["ESMF:Lon"] = list(self.stations['longitude_dd'])
        locstream["ESMF:Lat"] = list(self.stations['latitude_dd'])

        # create destination field
        if pl: # only for pressure level files
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt, nlev])
        else:
            dfield = ESMF.Field(locstream, name='dfield', 
                                ndbounds=[len(variables), nt])    

        # regridding function, consider ESMF.UnmappedAction.ERROR
        regrid2D = ESMF.Regrid(sfield, dfield,
                                regrid_method=ESMF.RegridMethod.BILINEAR,
                                unmapped_action=ESMF.UnmappedAction.IGNORE,
                                dst_mask_values=None)
                          
        # regrid operation, create destination field (variables, times, points)
        dfield = regrid2D(sfield, dfield)        
        sfield.destroy() #free memory                  
		
        # === write output netCDF file =========================================
        # dimensions: station, time OR station, time, level
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, level, station) or (time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4_CLASSIC')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'JRA-55, interpolated bilinearly to stations'
        rootgrp.featureType = "timeSeries"

        # dimensions
        station = rootgrp.createDimension('station', len(self.stations))
        time    = rootgrp.createDimension('time', nt)
        if pl: # only for pressure level files
            level = rootgrp.createDimension('level', nlev)

        # base variables
        time           = rootgrp.createVariable('time',     'i4',('time'))
        time.long_name = 'time'
        time.units     = 'hours since 1900-01-01 00:00:0.0'
        time.calendar  = 'gregorian'
        station             = rootgrp.createVariable('station',  'i4',('station'))
        station.long_name   = 'station for time series data'
        station.units       = '1'
        latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
        latitude.long_name  = 'latitude'
        latitude.units      = 'degrees_north'    
        longitude           = rootgrp.createVariable('longitude','f4',('station'))
        longitude.long_name = 'longitude'
        longitude.units     = 'degrees_east'  
        height           = rootgrp.createVariable('height','f4',('station'))
        height.long_name = 'height_above_reference_ellipsoid'
        height.units     = 'm'  
        if pl: # only for pressure level files
            level           = rootgrp.createVariable('level','i4',('level'))
            level.long_name = 'pressure_level'
            level.units     = 'hPa'  
        
        # assign base variables
        time[:] = nctime[tmask]
        if pl: # only for pressure level files
            level[:] = lev
        station[:]   = list(self.stations['station_number'])
        latitude[:]  = list(self.stations['latitude_dd'])
        longitude[:] = list(self.stations['longitude_dd'])
        height[:]    = list(self.stations['elevation_m'])
      
        # create and assign variables from input file
        for n, var in enumerate(variables):
            vname = ncf.variables[var].standard_name.encode('UTF8')
            if pl: # only for pressure level files
                tmp   = rootgrp.createVariable(vname,
                                                'f4',('time', 'level', 'station'))
            else:
                tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))   
                  
            tmp.long_name = ncf.variables[var].standard_name.encode('UTF8')
            tmp.units     = ncf.variables[var].units.encode('UTF8')  
            # assign values
            if pl: # only for pressure level files
                tmp[:] = dfield.data[n,:,:,:]
            else:
                tmp[:] = dfield.data[n,:,:]    
    
        rootgrp.close()
        ncf.close()
        
        # closed file ==========================================================
    
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
        varlist = [x.encode('UTF8') for x in ncf.variables.keys()]
        varlist.remove('time')
        varlist.remove('station')
        varlist.remove('latitude')
        varlist.remove('longitude')
        varlist.remove('level')
        varlist.remove('height')
        varlist.remove('geopotential_height')

        # === open and prepare output netCDF file ==============================
        # dimensions: station, time
        # variables: latitude(station), longitude(station), elevation(station)
        #            others: ...(time, station)
        # stations are integer numbers
        # create a file (Dataset object, also the root group).
        rootgrp = nc.Dataset(ncfile_out, 'w', format='NETCDF4')
        rootgrp.Conventions = 'CF-1.6'
        rootgrp.source      = 'JRA55, interpolated (bi)linearly to stations'
        rootgrp.featureType = "timeSeries"

        # dimensions
        station = rootgrp.createDimension('station', len(height))
        time    = rootgrp.createDimension('time', nt)

        # base variables
        time           = rootgrp.createVariable('time',     'i4',('time'))
        time.long_name = 'time'
        time.units     = 'hours since 1900-01-01 00:00:0.0'
        time.calendar  = 'gregorian'
        station             = rootgrp.createVariable('station',  'i4',('station'))
        station.long_name   = 'station for time series data'
        station.units       = '1'
        latitude            = rootgrp.createVariable('latitude', 'f4',('station'))
        latitude.long_name  = 'latitude'
        latitude.units      = 'degrees_north'    
        longitude           = rootgrp.createVariable('longitude','f4',('station'))
        longitude.long_name = 'longitude'
        longitude.units     = 'degrees_east'  
        height           = rootgrp.createVariable('height','f4',('station'))
        height.long_name = 'height_above_reference_ellipsoid'
        height.units     = 'm'  
       
        # assign base variables
        time[:] = ncf.variables['time'][:]
        station[:]   = ncf.variables['station'][:]
        latitude[:]  = ncf.variables['latitude'][:]
        longitude[:] = ncf.variables['longitude'][:]
        height[:]    = ncf.variables['height'][:]
        
        # create and assign variables from input file
        for var in varlist:
            vname = ncf.variables[var].long_name.encode('UTF8')
            tmp   = rootgrp.createVariable(vname,'f4',('time', 'station'))    
            tmp.long_name = ncf.variables[var].long_name.encode('UTF8')
            tmp.units     = ncf.variables[var].units.encode('UTF8')  
        # end file prepation ===================================================
    
                                                                                                
        # loop over stations
        for n, h in enumerate(height): 
            # convert geopotential [millibar] to height [m]
            # shape: (time, level)
            ele = ncf.variables['geopotential_height'][:,:,n] / 9.80665
            # TODO: check if height of stations in data range (+50m at top, lapse r.)
            
            # difference in elevation. 
            # level directly above will be >= 0
            dele = ele - h
            # vector of level indices that fall directly above station. 
            # Apply after ravel() of data.
            va = np.argmin(dele + (dele < 0) * 100000, axis=1) 
            # mask for situations where station is below lowest level
            mask = va < (nl-1)
            va += np.arange(ele.shape[0]) * ele.shape[1]
            
            # Vector level indices that fall directly below station.
            # Apply after ravel() of data.
            vb = va + mask # +1 when OK, +0 when below lowest level
            
            # weights
            wa = np.absolute(dele.ravel()[vb]) 
            wb = np.absolute(dele.ravel()[va])
            
            wt = wa + wb
            
            wa /= wt # Apply after ravel() of data.
            wb /= wt # Apply after ravel() of data.
            
            #loop over variables and apply interpolation weights
            for v, var in enumerate(varlist):
                #read data from netCDF
                data = ncf.variables[var][:,:,n].ravel()
                ipol = data[va]*wa + data[vb]*wb   # interpolated value                    
                rootgrp.variables[var][:,n] = ipol # assign to file   
    
        rootgrp.close()
        # closed file ==========================================================    


    def TranslateCF2short(self, dpar):
        """
        Map CF Standard Names into short codes used in JRA-55 netCDF files.
        """
        varlist = [] 
        for var in self.variables:
            varlist.append(dpar.get(var))
        # drop none
        varlist = [item for item in varlist if item is not None]      
        # flatten
        varlist = [item for sublist in varlist for item in sublist]         
        return(varlist) 

    def process(self):
        """
        Interpolate point time series from downloaded data. Provides access to 
        the more generically JRA-like interpolation functions.
        """                       

        # === 2D Interpolation for Surface  Data ===    
        # dictionary to translate CF Standard Names into JRA55
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['surface_temperature'],  # [K] 2m values
                'relative_humidity' : ['relative_humidity'], # [%]                                                       
                'wind_speed' : ['eastward_wind', 'northward_wind']}   # [m s-1] 2m & 10m values   
        varlist = self.TranslateCF2short(dpar)                      
        self.JRA2station(path.join(self.dir_inp,'JRA_surf_*.nc'), 
                            path.join(self.dir_out,'jra_sa_' + 
                                      self.list_name + '.nc'), self.stations,
                                      varlist, date = self.date)          
        
        # 2D Interpolation for Radiation Data 
        # dictionary to translate CF Standard Names into JRA55
        # pressure level variable keys.       
        dpar = {'downwelling_shortwave_flux_in_air' : ['downwelling_shortwave_flux_in_air'], # [W/m2] short-wave downward
                'downwelling_longwave_flux_in_air'  : ['downwelling_longwave_flux_in_air'], # [W/m2] long-wave downward
                'downwelling_shortwave_flux_in_air_assuming_clear_sky': ['downwelling_shortwave_flux_in_air_assuming_clear_sky'], # [W/m2] short-wave downward assuming clear sky
                'downwelling_longwave_flux_in_air_assuming_clear_sky': ['downwelling_longwave_flux_in_air_assuming_clear_sky']} # [W/m2] long-wave downward assuming clear sky
        varlist = self.TranslateCF2short(dpar)                           
        self.JRA2station(path.join(self.dir_inp,'JRA_fcst_*.nc'), 
                          path.join(self.dir_out,'jra_sr_' + 
                                    self.list_name + '.nc'), self.stations,
                                    varlist, date = self.date)          

        # === 2D Interpolation for Pressure-Level, Analyzed Meteorological DATA ===
        # dictionary to translate CF Standard Names into MERRA2
        # pressure level variable keys. 
        dpar = {'air_temperature'   : ['air_temperature'],  # [K]
                'relative_humidity' : ['relative_humidity'], # [%]
                'wind_speed'        : ['eastward_wind', 'northward_wind']}  # [m s-1]
        varlist = self.TranslateCF2short(dpar).append('geopotential_height')
        self.JRA2station(path.join(self.dir_inp,'JRA_Isobaric_*.nc'), 
                          path.join(self.dir_out,'jra_pl_' + 
                                    self.list_name + '.nc'), self.stations,
                                    varlist, date = self.date)  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        # 1D Interpolation for Pressure Level Data 
        self.levels2elevation(path.join(self.dir_out,'jra_pl_' + 
                                        self.list_name + '.nc'), 
                              path.join(self.dir_out,'jra_pl_' + 
                                        self.list_name + '_surface.nc'))
        

