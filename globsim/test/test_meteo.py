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


class RainFractionJennings(unittest.TestCase):
    """Tests for rain_fraction_jennings (Jennings et al. 2018)."""

    def test_cold_dry_is_snow(self):
        """At very cold temperatures and low humidity, almost all precip should be snow."""
        rf = met.rain_fraction_jennings(-15.0, 50.0)
        self.assertLess(rf, 0.01)

    def test_warm_humid_is_rain(self):
        """At warm temperatures and high humidity, almost all precip should be rain."""
        rf = met.rain_fraction_jennings(10.0, 90.0)
        self.assertGreater(rf, 0.99)

    def test_monotonic_with_temperature(self):
        """Rain fraction should increase with temperature at constant RH."""
        T = np.array([-10.0, -5.0, 0.0, 5.0, 10.0])
        RH = np.full_like(T, 70.0)
        rf = met.rain_fraction_jennings(T, RH)
        self.assertTrue(np.all(np.diff(rf) > 0))

    def test_monotonic_with_rh(self):
        """Rain fraction should increase with RH at constant temperature."""
        RH = np.array([20.0, 40.0, 60.0, 80.0, 100.0])
        T = np.full_like(RH, 1.0)
        rf = met.rain_fraction_jennings(T, RH)
        self.assertTrue(np.all(np.diff(rf) > 0))

    def test_output_range(self):
        """Output should always be between 0 and 1."""
        T = np.linspace(-30, 30, 50)
        RH = np.linspace(10, 100, 50)
        rf = met.rain_fraction_jennings(T, RH)
        self.assertTrue(np.all(rf >= 0))
        self.assertTrue(np.all(rf <= 1))

    def test_scalar_input(self):
        """Should work with scalar inputs."""
        rf = met.rain_fraction_jennings(0.0, 50.0)
        self.assertEqual(rf.shape, (1,))

    def test_snow_fraction_complement(self):
        """snow_fraction_jennings should be 1 - rain_fraction_jennings."""
        T = np.array([-5.0, 0.0, 5.0])
        RH = np.array([50.0, 70.0, 90.0])
        rf = met.rain_fraction_jennings(T, RH)
        sf = met.snow_fraction_jennings(T, RH)
        np.testing.assert_allclose(rf + sf, 1.0)


class WetBulbTemperature(unittest.TestCase):
    """Tests for wet_bulb_temperature (Stull 2011 approximation)."""

    def test_returns_kelvin(self):
        """Output should be in Kelvin (above 273.15 K for 0 °C)."""
        Tw = met.wet_bulb_temperature(0.0, 50.0)
        self.assertGreater(Tw[0], 200.0)  # sanity check on Kelvin scale
        self.assertLess(Tw[0], 400.0)

    def test_wet_bulb_le_air_temperature(self):
        """Wet-bulb temperature should never exceed air temperature (for RH < 100%)."""
        T = np.array([-10.0, 0.0, 10.0, 20.0])
        RH = np.array([30.0, 50.0, 70.0, 90.0])
        Tw_K = met.wet_bulb_temperature(T, RH)
        Tw_C = Tw_K - 273.15
        self.assertTrue(np.all(Tw_C <= T + 0.01))  # allow tiny floating-point tolerance

    def test_saturated_air_equals_air_temperature(self):
        """At 100% RH, wet-bulb temperature should equal air temperature (within approximation error)."""
        T = np.array([0.0, 10.0, 20.0])
        RH = np.full_like(T, 100.0)
        Tw_K = met.wet_bulb_temperature(T, RH)
        Tw_C = Tw_K - 273.15
        np.testing.assert_allclose(Tw_C, T, atol=0.2)  # Stull (2011) accuracy is within ~1.5 °C

    def test_scalar_input(self):
        """Should work with scalar inputs."""
        Tw = met.wet_bulb_temperature(5.0, 60.0)
        self.assertEqual(Tw.shape, (1,))


class SnowFractionWang(unittest.TestCase):
    """Tests for snow_fraction_wang and rain_fraction_wang (Wang et al. 2019)."""

    def test_cold_is_snow(self):
        """Very cold wet-bulb temperature should yield near-100% snow."""
        Tw_cold = 258.15  # -15 °C in K
        sf = met.snow_fraction_wang(Tw_cold)
        self.assertGreater(sf[0], 0.99)

    def test_warm_is_rain(self):
        """Warm wet-bulb temperature should yield near-100% rain."""
        Tw_warm = 283.15  # 10 °C in K
        rf = met.rain_fraction_wang(Tw_warm)
        self.assertGreater(rf[0], 0.99)

    def test_monotonic_snow_fraction(self):
        """Snow fraction should decrease with increasing wet-bulb temperature."""
        Tw = np.linspace(263.15, 283.15, 20)  # -10 to 10 °C in K
        sf = met.snow_fraction_wang(Tw)
        self.assertTrue(np.all(np.diff(sf) < 0))

    def test_output_range(self):
        """Snow and rain fractions should always be between 0 and 1."""
        Tw = np.linspace(253.15, 293.15, 50)
        sf = met.snow_fraction_wang(Tw)
        rf = met.rain_fraction_wang(Tw)
        self.assertTrue(np.all(sf >= 0))
        self.assertTrue(np.all(sf <= 1))
        self.assertTrue(np.all(rf >= 0))
        self.assertTrue(np.all(rf <= 1))

    def test_complement(self):
        """rain_fraction_wang + snow_fraction_wang should equal 1."""
        Tw = np.array([263.15, 273.15, 278.15])
        sf = met.snow_fraction_wang(Tw)
        rf = met.rain_fraction_wang(Tw)
        np.testing.assert_allclose(sf + rf, 1.0)

    def test_scalar_input(self):
        """Should work with scalar inputs."""
        sf = met.snow_fraction_wang(273.15)
        self.assertEqual(sf.shape, (1,))

    def test_requires_kelvin(self):
        """Should raise an error if temperature is not in Kelvin."""
        with self.assertRaises(ValueError):
            met.snow_fraction_wang(-5.0)  # Celsius, not Kelvin

    def test_rain_fraction_requires_kelvin(self):
        """rain_fraction_wang should also raise an error for non-Kelvin input."""
        with self.assertRaises(ValueError):
            met.rain_fraction_wang(-5.0)  # Celsius, not Kelvin


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
        self.zenith = np.array([45, 45, 45])

    def test_array(self):
        self.assertEqual(3, len(top.elevation_corrected_sw(zenith=self.zenith, grid_sw=self.sw, lat=self.lat,
                                                           lon=self.lon, time=self.time,
                                                           grid_elevation=self.ge, sub_elevation=self.se)[0]))

    def test_elevation_effect(self):
        self.assertGreater(top.elevation_corrected_sw(zenith=45, grid_sw=500, lat=50,
                                                      lon=-110,
                                                      time=datetime.datetime(year=2020, month=5, day=10, hour=18, tzinfo=datetime.timezone.utc),
                                                      grid_elevation=100, sub_elevation=2000)[1],
                           
                           top.elevation_corrected_sw(zenith=45, grid_sw=500, lat=50,
                                                      lon=-110,
                                                      time=datetime.datetime(year=2020, month=5, day=10, hour=18, tzinfo=datetime.timezone.utc),
                                                      grid_elevation=100, sub_elevation=1000)[1])
        

if __name__ == '__main__':
    unittest.main()
