from globsim.scale.JRAscale import JRAscale


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

    SCALING = {"sa": {},
               "sf": {},
               "pl": {},
               "pl_sur": {},
               "to": {}}


class J3QS(JRAscale):
    NAME = "JRA-3Q"
    REANALYSIS = "jra3q"



class J3QgS(JRAscale):
    NAME = "JRA-3QG"
    REANALYSIS = "jra3qg"

