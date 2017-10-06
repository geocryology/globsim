# The original names and standard names of variables from MERRA-2 database are listed. 
# For obtaining the selected variables from MERRA-2 by conducting merra-2.py, please follow the steps as following: 
 
  - In project_merra.globsim_download:  
    
    1. To set up the project directories for keeping scripts and output results
    
    2. To set up the directory for keeping credential file (.merrarc) (in '/home/user/' directory generally, but depend on personal preference)
    
    3. To set up the wanted chunk size (Number of Days for each individual output files) 
    
    4. To set up the area bounding for wanted area:
       Range of bbN to bbS : -90.0 ~ 90.0
       Range of bbW to bbE : -180.0 ~ 180.0 
    
    5. To set up the wanted range of elevation above the ground (in Unit: meter):
       Available range of elevation : 0 ~ 11,000 metres (0.1 ~ 1000 hPa in atmospheric pressure)
       
    6. To set up the wanted time range (to make sure in YYYY/MM/DD format)
       Available range of time : 1980/01/01 ~ 2017/07/31 
       (Latest updated on 2017/09, need to keep tracking the database at : https://disc.sci.gsfc.nasa.gov/datasets?page=1&keywords=merra-2)
       
    7. To set up the wanted variables by giving a list of standard names in the Table as following (from MERRA-2 original dataset):
        
       CF Standard Names                         MERRA-2 Standard Names	                                   Names in MERRA-2 Database
    
       air_temperature                           air_temperature                                            T
       Wind Speed & Wind Direction               eastward_wind                                              U
       Wind Speed & Wind Direction               northward_wind                                             V
       geopotential_height                       geopotential_height                                        H
       relative_humidity                         relative_humidity                                          RH
       air_temperature                           2-meter_air_temperature                                    T2M
       Wind Speed & Wind Direction               2-meter_eastward_wind                                      U2M           
       Wind Speed & Wind Direction               2-meter_northward_wind                                     V2M
       Wind Speed & Wind Direction               10-meter_eastward_wind                                     U10M
       Wind Speed & Wind Direction               10-meter_northward_wind                                    V10M
       precipitation_flux                        precipitation_flux                                         PRECTOT
       precipitation_amount                      NONE                                                       NONE                                                   
       downwelling_shortwave_flux_in_air         NONE                                                       NONE
       downwelling_longwave_flux_in_air          NONE                                                       NONE
       downwelling_shortwave_flux_in_air
       _assuming_clear_sky                       NONE                                                       NONE
       downwelling_longwave_flux_in_air
       _assuming_clear_sky                       NONE                                                       NONE
       surface_net_downward_longwave_flux        surface_net_downward_longwave_flux                         LWGNT
       surface_net_downward_longward_flux        surface_net_downward_longwave_flux
       _assuming_clear_sky                       _assuming_clear_sky                                        LWGNTCLR
       NONE                                      surface_incoming_shortwave_flux                            SWGDN
       NONE                                      surface_incoming_shortewave_flux_assuming_clear_sky        SWGDNCLR
       NONE                                      longwave_flux_emitted_from_surface                         LWGEM


   
Important Notes:

1. Set up credential file (in .merrarc text file):
   Username (given in the first line without any symbols directly)
   password (given in the second line without any symbols directly)

2. For obtaining single variable, to set up the standard name of variable with ‘,’ in the end
   Example: variables = air_temperature,

3. To the variables in category of CF Variables Names, BUT which are NOT AVAILABLE in MERRA-2 Database, to get the similar variables instead as following table:
     
     Wanted Variables                                      Replaced Variables
     
     precipitation_amount                                  precipitation_flux
     downwelling_shortwave_flux_in_air                     surface_incoming_shortwave_flux
     downwelling_longwave_flux_in_air                      surface_net_downward_longwave_flux, longwave flux emitted from surface
     downwelling_shortwave_flux_in_air_assuming_clear_sky  surface_incoming_shortwave_flux_assuming_clear_sky
     downwelling_longwave_flux_in_air_assuming_clear_sky   surface_net_downward_longwave_flux_assuming_clear_sky, longwave_flux_emitted_from_surface
     

4. Update and extend the list of Standard Names if it is needed. 

5. Refer the CF Standard Names at the website: 
   http://cfconventions.org/Data/cf-standard-names/46/build/cf-standard-name-table.html



    
