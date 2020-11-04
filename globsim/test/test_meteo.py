import unittest
from globsim import meteorology as met


class SatVap(unittest.TestCase):

    def test_zero(self):
        self.assertAlmostEqual(met.satvapp_kPa_fT(0), 0.6113)


class PressureFromElevation(unittest.TestCase):
    
    def test_zero(self):
        self.assertAlmostEqual(met.pressure_from_elevation(0), 1013.25)


class WaterVaporPressure(unittest.TestCase):

    def test_zero(self):
        self.assertAlmostEqual(met.vapp_kPa_fTd(0), 6.112)
        
    def test_increasing(self):
        self.assertGreater(met.vapp_kPa_fTd(10), met.vapp_kPa_fTd(0))


class SpecificHumidity(unittest.TestCase):

    def test_increase_pressure(self):
        self.assertGreater(met.spec_hum_kgkg(10, 101325), met.spec_hum_kgkg(10, 105325))


if __name__ == '__main__':
    unittest.main()
