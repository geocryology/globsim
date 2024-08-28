from globsim.interpolate.JRAinterpolate import JRAinterpolate
from globsim.download.JRA3Qdownload import J3QD

import logging

logger = logging.getLogger('globsim.interpolate')


class J3QI(JRAinterpolate):

    REANALYSIS = "jra3q"
    SA_INTERVAL = 6
    PL_INTERVAL = 6
    SF_INTERVAL = 1

    T_UNITS = "hours since 1900-01-01 00:00:00"


class J3QgI(J3QI):

    REANALYSIS = "jra3qg"



