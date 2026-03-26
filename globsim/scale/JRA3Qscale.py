import numpy as np

from globsim.scale.JRAscale import JRAscale
from globsim.scale.scalenames import ScaleNames as SN


class J3QSjma(JRAscale):
    NAME = "JRA-3Q-JMA"
    REANALYSIS = "jra3q-jma"

    # translator for JRA-55 names
    VARNAMES = {"sa":{SN.temperature: "2t",
                      SN.specific_humidity: "sp",
                      SN.rh: "2r",
                      SN.u_wind: "10u",
                      SN.v_wind: "10v"},
            "sf":{SN.precipitation_total: "mtpf",
                  SN.sw_down_flux: "dswrf",
                  SN.lw_down_flux: "dswrf",},
            "pl":{SN.temperature: "t",
                  },
            "pl_sur":{},
            "to":{}}
    CONVERTERS = {}


class J3QS(JRAscale):
    NAME = "JRA-3Q"
    REANALYSIS = "jra3q"
    CONVERTERS = {}


class J3QgS(JRAscale):
    NAME = "JRA-3QG"
    REANALYSIS = "jra3qg"

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
        """Convert geopotential height in gpm to elevation in m, using an approximation that 
        they are equal. This is not exact but is close enough for scaling purposes and avoids 
        needing to read the surface pressure to get the actual gravity at each point."""
        converted_data = data  # geopotential height in gpm is approximately equal to elevation in m
        converted_units = "m"

        return converted_data, converted_units


