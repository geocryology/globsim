import unittest
from globsim import meteorology as met
from globsim.scale import toposcale as top
import numpy as np
import datetime


class SatVap(unittest.TestCase):

    def test_zero(self):
        self.assertAlmostEqual(met.satvapp_kPa_fT(0), 0.6113)


class PressureFromElevation(unittest.TestCase):
    
    def setUp(self) -> None:
        self.test_array = np.array([100,200,300])

    def test_zero(self):
        self.assertAlmostEqual(met.pressure_from_elevation(0), 1013.25)

    def test_array(self):
        self.assertEqual(4, len(met.pressure_from_elevation(np.array([10,150, 300,1000]))))


class WaterVaporPressure(unittest.TestCase):

    def test_zero(self):
        self.assertAlmostEqual(met.vapp_kPa_fTd(0), 6.112)
        
    def test_increasing(self):
        self.assertGreater(met.vapp_kPa_fTd(10), met.vapp_kPa_fTd(0))


class SpecificHumidity(unittest.TestCase):

    def test_increase_pressure(self):
        self.assertGreater(met.spec_hum_kgkg(10, 101325), met.spec_hum_kgkg(10, 105325))


class ZenithAngle(unittest.TestCase):

    def setUp(self) -> None:
        self.lat = np.array([50, 60, 70])
        self.lon = np.array([-110, -120, -130])
        self.time = np.array([datetime.datetime(year=2020, month=5, day=10, hour=hour, tzinfo=datetime.timezone.utc) for hour in [1, 9, 17]])

    def test_array(self):
        self.assertEqual(3, len(top.solar_zenith(lat=self.lat, lon=self.lon, time=self.time)))


class LocalAirmass(unittest.TestCase):

    pass

class ToposcaleShortwaveElevation(unittest.TestCase):

    def setUp(self) -> None:
        self.sw = np.array([500, 600, 700])
        self.lat = np.array([50, 50, 50])
        self.lon = np.array([-110, -110, -110])

        self.time = np.array([datetime.datetime(year=2020, month=5, day=10, hour=hour, tzinfo=datetime.timezone.utc) for hour in [1, 9, 17]])
        self.se = np.array([500, 500, 500])
        self.ge = np.array([100, 500, 1000])

    def test_array(self):
        self.assertEqual(3, len(top.elevation_corrected_sw(grid_sw=self.sw, lat=self.lat,
                                                           lon=self.lon, time=self.time,
                                                           grid_elevation=self.ge, sub_elevation=self.se)))

    def test_elevation_effect(self):
        self.assertGreater(top.elevation_corrected_sw(grid_sw=500, lat=50,
                                                      lon=-110,
                                                      time=datetime.datetime(year=2020, month=5, day=10, hour=18, tzinfo=datetime.timezone.utc),
                                                      grid_elevation=100, sub_elevation=2000)[0],
                           
                           top.elevation_corrected_sw(grid_sw=500, lat=50,
                                                      lon=-110,
                                                      time=datetime.datetime(year=2020, month=5, day=10, hour=18, tzinfo=datetime.timezone.utc),
                                                      grid_elevation=100, sub_elevation=1000)[0])
        

if __name__ == '__main__':
    unittest.main()
