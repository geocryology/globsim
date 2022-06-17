.. _stationlist:

The station_list file
=====================

The station list CSV file is used to describe the point locations at which globsim output is needed. 


Headers
^^^^^^^
The station  list CSV file uses the following headers for the data columns:

=========================  ==============      ===============
   **header**              **required?**       **Description** 
-------------------------  --------------      ---------------
``station_number``              yes            An integer. Must be unique within the file.
``station_name``                yes            The name of the location. Must be unique within the file.
``latitude_dd``                 yes            The WGS84 latitude of the location in decimal degrees.
``longitude_dd``                yes            The WGS84 latitude of the location in decimal degrees.
``elevation_m``                 yes            The elevation of the location in metres.
``sky_view``                    no             The sky-view factor (between 0 and 1) at the location. Used to correct for topographic effects (e.g. if the TOPOscale kernel is used).
``slope``                       no             The slope of the ground at the location in degrees. Used to correct for topographic effects (e.g. if the TOPOscale kernel is used).
``aspect``                      no             The slope of the ground at the location in degrees. Used to correct for topographic effects (e.g. if the TOPOscale kernel is used).
=========================  ==============      ===============

Example Siteslist file
^^^^^^^^^^^^^^^^^^^^^^
Here is an example of a station list file

::

   station_number,station_name,longitude_dd,latitude_dd,elevation_m,sky_view,slope,aspect
   1,'hilltop',-129.1,64.1,2500,0.9,70,90
   2,'swamp_2',-129.4,61.1,1000,0.5,0,0
   3,'research_station',-126.1,64.1,1500,1.0,10,180

