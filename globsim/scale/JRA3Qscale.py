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
    
    CONVERTERS = {
                  ("to", SN.elevation): "_geopotential_to_m"}  # note: pl already in m


