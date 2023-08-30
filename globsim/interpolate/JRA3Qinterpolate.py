from globsim.interpolate.JRAinterpolate import JRAinterpolate
from globsim.download.JRA3Qdownload import J3QD

import logging

logger = logging.getLogger('globsim.interpolate')


class J3QI(JRAinterpolate):

    REANALYSIS = "jra3q"
    SA_INTERVAL = 6
    PL_INTERVAL = 6
    SF_INTERVAL = 1

    T_UNITS = J3QD._tunits

    GEOPOTENTIAL = "gh"

