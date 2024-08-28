import logging

from globsim.download.JRAdownload import JRAdownload
from globsim.download.RDA import Rdams 
from globsim.download.jra_dict_formatters import J3QDictFormatter, J3QGDictFormatter
from globsim.download.JraDownloadHandler import J3QDownloadHandler

logger = logging.getLogger('globsim.download')


class J3QD(JRAdownload):

    JRA_VERSION = 'jra3q'
    API = Rdams
    DICT_FORMATTER = J3QDictFormatter
    FILE_HANDLER = J3QDownloadHandler
    dsID = 'd640000'
    timeName = 'time'
    
class J3QgD(J3QD):

    JRA_VERSION = 'jra3qg'
    API = Rdams
    DICT_FORMATTER = J3QGDictFormatter
    FILE_HANDLER = J3QDownloadHandler
    dsID = 'd640000'
    timeName = 'time'
