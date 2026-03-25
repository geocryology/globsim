import numpy as np
from cfunits import Units

from globsim.scale.JRAscale import JRAscale
from globsim.scale.scalenames import ScaleNames as SN


class JRA55(JRAscale):
      NAME = "JRA-55"
      REANALYSIS = "jra55"
      SCALING = {"sf": {"Total precipitation": (1 / (24 * 3600), 0)},
                  "sa": {},
                  "pl": {},
                  "to": {},
                  "pl_sur": {}}
      
      CONVERTERS = {
        ("sf", SN.precipitation_rate): "_daily_precip_to_rate",
            }
      
      def _daily_precip_to_rate(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        """mm day-1 → kg m-2 s-1"""
        input_units = Units(nc_var.units)
        converted_units = input_units / Units("s")
        converted_data = data / (24 * 3600)

        return converted_data, converted_units.units


class J3QSjma(JRAscale):
    NAME = "JRA-3Q-JMA"
    REANALYSIS = "jra3q-jma"

    # translator for JRA-55 names
    J55T = {"sa":{"Temperature": "2t",
                  "Specific humidity": "sp",
                  "Relative humidity": "2r",
                  "u-component of wind": "10u",
                  "v-component of wind": "10v"},
            "sf":{"Total precipitation": "mtpf",
                  'Downward solar radiation flux': "dswrf",
                  'Downward longwave radiation flux': "dswrf",},
            "pl":{"Temperature": "t",
                  },
            "pl_sur":{},
            "to":{}}
    CONVERTERS = {}

    SCALING = {"sa": {},
               "sf": {},
               "pl": {},
               "pl_sur": {},
               "to": {}}


class J3QS(JRAscale):
    NAME = "JRA-3Q"
    REANALYSIS = "jra3q"
    SCALING = {"sa": {},
            "sf": {},
            "pl": {},
            "pl_sur": {},
            "to": {}}
    CONVERTERS = {}



class J3QgS(JRAscale):
    NAME = "JRA-3QG"
    REANALYSIS = "jra3qg"
    SCALING = {"sa": {},
               "sf": {},
               "pl": {},
               "pl_sur": {},
               "to": {}}
    VARNAMES = {
        "sa":     {SN.time:        "time",
                   SN.longitude:     "longitude",
                   SN.latitude:      "latitude",
                   SN.temperature: "Temperature",
                   SN.rh:          "Relative humidity",
                   SN.u_wind:      "u-component of wind",
                   SN.v_wind:      "v-component of wind",
                   SN.pressure:     "air_pressure",
                   SN.specific_humidity: "Specific humidity"},
        "sf":     {SN.time:          "time",
                   SN.longitude:     "longitude",
                   SN.latitude:      "latitude",
                   SN.sw_down_flux:       "Downward solar radiation flux",
                   SN.lw_down_flux:       "Downward longwave radiation flux",
                   SN.precipitation_total: "Total precipitation",
                   SN.precipitation_rate: "Total precipitation"},
        "pl":     {SN.time:          "time",
                   SN.longitude:     "longitude",
                   SN.latitude:      "latitude",
                   SN.temperature:   "Temperature",
                   SN.elevation:  "Geopotential height"},  # note: pl already in m, but JRA-3QG has geopotential height instead of elevation
        "pl_sur": {SN.time:          "time",
                   SN.longitude:     "longitude",
                   SN.latitude:      "latitude",
                   SN.elevation:     "height",
                   SN.temperature:   "Temperature",
                   SN.rh:            "Relative humidity",
                   SN.pressure:      "air_pressure"},
        "to":     {SN.time:          "time",
                   SN.longitude:     "longitude",
                   SN.latitude:      "latitude",
                   SN.elevation:     "Geopotential",
                   SN.geopotential:  "Geopotential"},
    }

    CONVERTERS = {
                  ("to", SN.elevation): "_geopotential_to_m",
                  ("pl", SN.elevation): "_gpm_to_m_approx",
                  }  # note: pl already in m
    
    def _gpm_to_m_approx(self, data, nc_var, _slice) -> tuple[np.ndarray, str]:
        """Convert geopotential height in gpm to elevation in m, using an approximation that divides by g. This is not exact but is close enough for scaling purposes and avoids needing to read the surface pressure to get the actual gravity at each point."""
        converted_data = data  # geopotential height in gpm is approximately equal to elevation in m
        converted_units = "m"

        return converted_data, converted_units


