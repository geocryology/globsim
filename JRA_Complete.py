from datetime        import date, datetime
from dateutil.rrule  import rrule, DAILY
from ftplib          import FTP
from netCDF4         import Dataset
from generic         import ParameterIO, StationListRead
from os              import path, listdir
import netCDF4       as nc
import pygrib
import time
import sys
import os.path


"""
Download ranges of data from the JRA-55 server
Data available from 1958 to present date (give 1 or 2 days for upload)
~Author: Christopher Molnar
~Date: August 08, 2017
"""
class JRA_Download:
      
    """
    Get the time sets that you want to download (start day to final day)
    returns:
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
       
        except:
            print "\nInvalid Entry"
            print "Please make sure your longitude position is between 0 and 358.75 (inclusively)"
            sys.exit(0)
    
        return (latTopPosition, latBottomPosition, lonLeftPosition, lonRightPosition)      
    
    
    """
    """
    def ElevationCalculator(self, elevationMin, elevationMax):
        total_elevations = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]

        minNum = min(total_elevations, key=lambda x:abs(x-elevationMin))
        maxNum = min(total_elevations, key=lambda x:abs(x-elevationMax))
        
        print "YEAH"
        
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
    -data_type: The file that is being downloaded
    -hourly_increment: The hour increments used to record the data
    -middle: Extra text in the middle of the filename needed for downloading certain files
    """
    def ftp_Download(self, start_day, end_day, data_type, hourly_increment, middle, savePath):         
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
                           
                        filename = data_type + middle + "." +  dt.strftime("%Y") + dt.strftime("%m") + dt.strftime("%d") + ending
                        completeName = os.path.join(savePath, filename) # Generate the filename and save it to the proper folder
                   
                        localfile = open(completeName, 'wb')
                   
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
    Download the necessary Grib Files
    Parameters:
    -startDay: The starting download day
    -endDay: The ending download day
    -save_path: The save path
    -variable_data: All of the data that we want to download
    -fcst_list: The forcast data that we need to download
    -surf_list: The surf data that we need to download
    """  
    def DownloadGribFile(self, startDay, endDay, save_path, variable_data, fcst_list, surf_list):

        if ("Geopotential Height" in variable_data):
            # Download anl_p125 for geopotential height (6 hour increments) 
            self.ftp_Download(startDay, endDay, "anl_p125", 6, "_hgt", save_path + 'Grib') 
        
        if ("Temperature" in variable_data):
            # Download anl_p125 for temperatures (6 hour increments) 
            self.ftp_Download(startDay, endDay, "anl_p125", 6, "_tmp", save_path + 'Grib')
        
        if ("U-Component of Wind" in variable_data):
            # Download anl_p125 for u-component of wind (6 hour increments)
            self.ftp_Download(startDay, endDay, "anl_p125", 6, "_ugrd", save_path + 'Grib') 
        
        if ("V-Component of Wind" in variable_data):
            # Download anl_p125 for v-component of wind (6 hour increments)
            self.ftp_Download(startDay, endDay, "anl_p125", 6, "_vgrd", save_path + 'Grib') 
        
        if ("Relative Humidity" in variable_data):
            # Download anl_p125 for relative humidity (6 hour increments)
            self.ftp_Download(startDay, endDay, "anl_p125", 6, "_rh", save_path + 'Grib') 
        
        if (len(surf_list) > 0):
            # Download anl_surf125 (6 hour increments) 
            self.ftp_Download(startDay, endDay, "anl_surf125", 6, "", save_path + 'Grib') 
        
        if (len(fcst_list) > 0):
            # Download fcst_phy2m125 (3 hour increments) 
            self.ftp_Download(startDay, endDay, "fcst_phy2m125", 3, "", save_path + 'Grib') 
        
        print "\nAll Downloads Finished :)"


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
    -filePath: The direcory where you will find the folder with the GRIB files
    Returns:
    -The result from the ConvertLatLon function
    """
    def GetLatLon(self, filename, filePath):
        fileLocation = filePath+ 'Grib/'
        try:
            grbs = pygrib.open(fileLocation + filename)
        except:
            print "file: " + fileLocation + filename +  " not found :("
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
    -filePath: The direcory where you will find the folder with the GRIB files
    Returns:
    -value: An array with all of the values
    -levels: An int with the datas level
    """
    def ExtractData(self, filename, date, dataName, filePath):    
        value = []    
        fileLocation = filePath + 'Grib/' 
        for z in range(0, len(date)):
            try:
                grbs = pygrib.open(fileLocation + filename + date[z])
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
    -latBottom: The bottom latitude coordinate
    -latTop: The top latitude coordinate
    -lonLeft: The left longitude coordinate
    -lonRight: The right longitude coordinate
    Returns:
    -ssDay: The subseted data set
    """
    def SubsetTheData(self, dataValues, numDays, latBottom, latTop, lonLeft, lonRight):
        # Subset the data into boxes for easy transfer
        if (latBottom > latTop):
            lat1 = 0
            lat2 = latTop
            lat3 = latBottom
            lat4 = 145
        else:
            lat1 = latBottom
            lat2 = latTop
            lat3 = latBottom
            lat4 = latTop
        if (lonLeft > lonRight):
            lon1 = 0
            lon2 = lonRight
            lon3 = lonLeft
            lon4 = 288
        else:
            lon1 = lonLeft
            lon2 = lonRight
            lon3 = lonLeft
            lon4 = lonRight
        ssDay = []
        for a in range(0,numDays):
          ssDay.append([])
          for b in range(0,145):
            ssDay[a].append([])
            for c in range(0,288):
              if ((lat1 <= b <= lat2 or lat3 <= b <= lat4) and (lon1 <= c <= lon2 or lon3 <= c <= lon4)):
                ssDay[a][b].append(dataValues[a][b][c])
              else:
                ssDay[a][b].append(None)
        
        return (ssDay)

   
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
    """
    def Main(self, startDay, endDay, numDays, date, bottomLat, topLat, leftLon, rightLon, savePath, JRA_Dictionary, fcst_data):    
        
        completeName = os.path.join(savePath + 'netCDF', "Finalfcst_" + startDay + "-" + endDay + ".nc") 
    
        f = Dataset(completeName, "w", format = "NETCDF4") # Name of the netCDF being created 
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon("fcst_phy2m125." + startDay + "00", savePath)
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
        
        #########
        x=0
        for dataName in fcst_data:
        
            # Create netCDF Variable
            
            data, level = Grib2CDF().ExtractData(JRA_Dictionary[dataName][1], date, JRA_Dictionary[dataName][0], savePath)

            if (x == 0):
                levels = f.createDimension("level", JRA_Dictionary[dataName][2])
                x = 1

            dataVariable = f.createVariable(dataName, "f4", ("time","level", "latitude", "longitude"))
            dataVariable.standard_name = dataName
            dataVariable.units = JRA_Dictionary[dataName][3]
            dataVariable[:,:,:,:] = Grib2CDF().SubsetTheData(data, timeSize, bottomLat, topLat, leftLon, rightLon)

            print "Converted:", dataName

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
    """
    def Main(self, startDay, endDay, numDays, date, bottomLat, topLat, leftLon, rightLon, savePath, JRA_Dictionary, surf_data):   
        
        completeName = os.path.join(savePath + 'netCDF', "Finalsurf_" + startDay + "-" + endDay + ".nc") 
    
        f = Dataset(completeName, "w", format = "NETCDF4") # Name of the netCDF being created 
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon("anl_surf125." + startDay + "00", savePath)
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
        
        #########
        x=0
        for dataName in surf_data:
        
            # Create netCDF Variable
            data, level = Grib2CDF().ExtractData(JRA_Dictionary[dataName][1], date, JRA_Dictionary[dataName][0], savePath)
            if (x==0):
              levels = f.createDimension("level", JRA_Dictionary[dataName][2])
              x =1
            dataVariable = f.createVariable(dataName, "f4", ("time","level", "latitude", "longitude"))
            dataVariable.standard_name = dataName
            dataVariable.units = JRA_Dictionary[dataName][3]
            dataVariable[:,:,:,:] = Grib2CDF().SubsetTheData(data, timeSize, bottomLat, topLat, leftLon, rightLon)
        
            print "Converted:", dataName
        
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
                grbs = pygrib.open(fileLocation + filename + date[z])
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
    -latBottom: The bottom latitude coordinate
    -latTop: The top latitude coordinate
    -lonLeft: The left longitude coordinate
    -lonRight: The right longitude coordinate
    Returns:
    -ssDay: The subseted data set
    """
    def SubsetTheData(self, dataValues, numDays, latBottom, latTop, lonLeft, lonRight, elevationMinRange, elevationMaxRange):
        # Subset the data into boxes for easy transfer
        if (latBottom > latTop):
            lat1 = 0
            lat2 = latTop
            lat3 = latBottom
            lat4 = 145
        else:
            lat1 = latBottom
            lat2 = latTop
            lat3 = latBottom
            lat4 = latTop
        if (lonLeft > lonRight):
            lon1 = 0
            lon2 = lonRight
            lon3 = lonLeft
            lon4 = 288
        else:
            lon1 = lonLeft
            lon2 = lonRight
            lon3 = lonLeft
            lon4 = lonRight
        ssDay = []
        for a in range(0,numDays):
          ssDay.append([])
          for b in range(elevationMinRange,elevationMaxRange + 1):
            ssDay[a].append([])
            for c in range(0,145):
              ssDay[a][b].append([])
              for d in range(0,288):
                if ((lat1 <= c <= lat2 or lat3 <= c <= lat4) and (lon1 <= d <= lon2 or lon3 <= d <= lon4)):
                  ssDay[a][b][c].append(dataValues[a][b][c][d])
                else:
                  ssDay[a][b][c].append(None)
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
    """
    def Main(self, startDay, endDay, numDays, date, bottomLat, topLat, leftLon, rightLon, savePath, JRA_Dictionary, isobaric_data, elevationMinRange, elevationMaxRange):
        
        completeName = os.path.join(savePath + 'netCDF', "FinalIsobaric_" + startDay + "-" + endDay + ".nc") 
         
        f = Dataset(completeName, "w", format = "NETCDF4") # Name of the netCDF being created 
        
        # Create the Latitude and Longitude
        lats, lons = Grib2CDF().GetLatLon("anl_p125_hgt." + startDay + "00", savePath) 
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
        
        x=0
        for dataName in isobaric_data:
        
            # Create netCDF Variable
            data, level = self.ExtractData(JRA_Dictionary[dataName][1], date, savePath)
            
            if (dataName == "Relative Humidity" and elevationMaxRange >= 9):
                Levels =f.createDimension("Levels", elevationMaxRange - 9)
                dataVariable = f.createVariable(dataName, "f4", ("time", "Levels", "latitude", "longitude"))
                dataVariable.standard_name = dataName
                dataVariable.units = JRA_Dictionary[dataName][3]
                if (elevationMinRange < 9):
                    tempMinRange = 0
                else:
                    tempMinRange = elevationMinRange - 9
                  
                print tempMinRange
                print elevationMaxRange - 10
                dataVariable[:,:,:,:] = self.SubsetTheData(data, timeSize, bottomLat, topLat, leftLon,  rightLon, tempMinRange, elevationMaxRange - 10)
            
            elif (dataName != "Relative Humidity"):
                if (x==0):
                    levels = f.createDimension("level", elevationMaxRange + 1)
                    x =1
                dataVariable = f.createVariable(dataName, "f4", ("time","level", "latitude", "longitude"))
                dataVariable.standard_name = dataName
                dataVariable.units = JRA_Dictionary[dataName][3]
                dataVariable[:,:,:,:] = self.SubsetTheData(data, timeSize, bottomLat, topLat, leftLon,  rightLon, elevationMinRange, elevationMaxRange)
                
           print "Converted:", dataName
        # Descriptions
        f.description = "Geography example"
        f.history = "Created today"
        f.source = "netCDF4 practice"
        
        f.close()


"""
From Dr. Grubers ERA file
"""
class ERAdownload(object):
    """
    Class for ERA-Interim data that has methods for querying 
    the ECMWF server, returning all variables usually needed.
       
    Args:
        pfile: Full path to a Globsim Download Parameter file. 
              
    Example:          
        ERAd = ERAdownload(pfile) 
        ERAd.retrieve()
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
        self.directory = path.join(par.project_directory, "eraint")  
        if path.isdir(self.directory) == False:
            raise ValueError("Directory does not exist: " + self.directory)   
     
        # variables
        self.variables = par.variables
            
        # chunk size for downloading and storing data [days]        
        self.chunk_size = par.chunk_size           


###################################
# Run the program
t0 = time.time()  

jra = ERAdownload("Parameter_Stuff.txt")

# Area data 
try:
    latBottom = float(jra.area["south"])
    latTop = float(jra.area["north"])
    lonLeft = float(jra.area["west"])
    lonRight = float(jra.area["east"])
    
    if (lonLeft < 0):
        lonLeft = 180 + abs(lonLeft)
        
    if (lonRight < 0):
        lonRight = 180 + abs(lonRight)
        
    latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition = JRA_Download().ConvertLatLon(latBottom, latTop, lonLeft,lonRight)
except:
    print "Invalid area format"
    sys.exit(0)

# Chunk data
try:
    chunk_size = int(jra.chunk_size)
except:
    print "Invalid chunk size"
    sys.exit(0)

# Date data
try:
    start = str(jra.date["beg"])
    end = str(jra.date["end"])
    
    startDay, endDay = JRA_Download().TimeSet(start, end) # Get the time period for the data
    
except:
    print "Invalid date"
    print "Please make sure your date is in YYYY-MM-DD format"
    sys.exit(0)

# Elevation data 
try:
    elevationMin = jra.elevation["min"]
    elevationMax = jra.elevation["max"]

    elevationMinRange, elevationMaxRange = JRA_Download().ElevationCalculator(elevationMin, elevationMax)

except:
    print "Invalid elevation entered"
    sys.exit(0)

elevationMinRange = 0
elevationMaxRange = 8


# Directory Information
directory = jra.directory
print directory
# TO FINISH

# List of Variables
variables = jra.variables


shared_data = {
              "air_temperature"                                      : ["Temperature", "Temperature 2D"],
              "relative_humidity"                                    : ["Relative Humidity", "Relative Humidity 2D"],
              "precipitation_amount"                                 : ["Total Precipitation"],
              "downwelling_longwave_flux_in_air"                     : ["Downward Long Wave Radiation Flux"],
              "downwelling_longwave_flux_in_air_assuming_clear_sky"  : ["Clear Sky Downward Long Wave Radiation Flux"],
              "downwelling_shortwave_flux_in_air"                    : ["Downward Solar Radiation Flux"],
              "downwelling_shortwave_flux_in_air_assuming_clear_sky" : ["Clear Sky Downward Solar Radiation Flux"],
              "wind_from_direction"                                  : ["U-Component of Wind", "V-Component of Wind" ],
              "wind_speed"                                           : ["U-Component of Wind 2D", "V-Component of Wind 2D"],
              "geopotential_height"                                  : ["Geopotential Height"]
              }
              
variable_data = []
for x in variables:
    if (x in shared_data):
        for y in range(0, len(shared_data[x])):
            variable_data.append(shared_data[x][y])


JRA_Dictionary = {
                  "Geopotential Height"                         : ["gh", "anl_p125_hgt.", 37, "gpm"],
                  "Temperature"                                 : ["t", "anl_p125_tmp.", 37, "K"],
                  "U-Component of Wind"                         : ["u", "anl_p125_ugrd.", 37, "m/s"],
                  "V-Component of Wind"                         : ["v", "anl_p125_vgrd.", 37, "m/s"],
                  "Relative Humidity"                           : ["r", "anl_p125_rh.", 27, "%"],
                  "Temperature 2D"                              : ["t", "anl_surf125.", 1, "K"],
                  "Relative Humidity 2D"                        : ["r", "anl_surf125.", 1, "%"],
                  "U-Component of Wind 2D"                      : ["u", "anl_surf125.", 1, "m/s"],
                  "V-Component of Wind 2D"                      : ["v", "anl_surf125.", 1, "m/s"],
                  "Total Precipitation"                         : ["tpratsfc", "fcst_phy2m125.", 1, "mm/day"],
                  "Clear Sky Downward Solar Radiation Flux"     : ["csdsf", "fcst_phy2m125.", 1, "W/(m^2)"],
                  "Clear Sky Downward Long Wave Radiation Flux" : ["csdlf", "fcst_phy2m125.", 1, "W/(m^2)"],
                  "Downward Solar Radiation Flux"               : ["dswrf", "fcst_phy2m125.", 1, "W/(m^2)"],
                  "Downward Long Wave Radiation Flux"           : ["dlwrf", "fcst_phy2m125.", 1, "W/(m^2)"],
                  }


fcst_variables = ["Total Precipitation" , "Clear Sky Downward Solar Radiation Flux", "Clear Sky Downward Long Wave Radiation Flux", "Downward Solar Radiation Flux", "Downward Long Wave Radiation Flux"]

surf_variables = ["Temperature 2D", "Relative Humidity 2D", "U-Component of Wind 2D", "V-Component of Wind 2D"]

isobaric_variables = ["Geopotential Height", "Temperature", "U-Component of Wind", "V-Component of Wind", "Relative Humidity"]

fcst_list = list(set(variable_data) & set(fcst_variables))
surf_list = list(set(variable_data) & set(surf_variables))
isobaric_list = list(set(variable_data) & set(isobaric_variables))

save_path = '/home/cmolnar/Test/' # The current directory

#JRA_Download().DownloadGribFile(startDay, endDay, save_path, variable_data, fcst_list, surf_list)

timeSize = chunk_size # The number of days you want saved together
Fdate = []
Adate = []
Idate = []
x = 0
finalDay = str(endDay)
print finalDay

for dt in rrule(DAILY, dtstart = startDay, until = endDay): # Loop through the days to create netCDF files
    currentDay = str(dt.strftime("%Y") + "-" + dt.strftime("%m") + "-" + dt.strftime("%d"))
    if (x < timeSize): # If x is less than timeSize append the current day to the date lists
        x += 1
        Fdate.append(str(dt.strftime("%Y")) + str(dt.strftime("%m")) + str(dt.strftime("%d")))
        Adate.append(str(dt.strftime("%Y")) + str(dt.strftime("%m")) + str(dt.strftime("%d")))
        Idate.append(str(dt.strftime("%Y")) + str(dt.strftime("%m")) + str(dt.strftime("%d")))
      
    if (x == timeSize or currentDay == finalDay): # If time size reached or last day reached send data chuncks to there respective classes
        print "here"
       #fcst_phy2m().Main(Fdate[0], Fdate[x-1], x, Fdate, latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition, save_path, JRA_Dictionary, fcst_list)
        #anl_surf().Main(Adate[0], Adate[x-1], x, Adate, latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition, save_path, JRA_Dictionary, surf_list)
        Isobaric().Main(Idate[0], Idate[x-1], x, Idate, latBottomPosition, latTopPosition, lonLeftPosition, lonRightPosition, save_path, JRA_Dictionary, isobaric_list, elevationMinRange, elevationMaxRange)
        x = 0
        Fdate = []
        Adate = []
        Idate = []   


print "\nAll Conversions Finished!"
print "Have a nice day! "  
t1 = time.time()
total = t1 - t0
print "Total run time:", total  

