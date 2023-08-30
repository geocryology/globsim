from globsim.scale.JRAscale import JRAscale


class J3QS(JRAscale):
    NAME = "JRA-3Q"
    REANALYSIS = "jra3q"

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
                  }}

    SCALING = {"sa": {},
               "sf": {},
               "pl": {}}
    
    def get_name(self, file:str, jra55name:str) -> str:
        n = self.J55T[file].get(jra55name)
        if n is None:
            return jra55name
        else:
            return n
