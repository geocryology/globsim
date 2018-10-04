# -*- coding: utf-8 -*-

# ==== DOWNLOAD ===============================================================

from era5 import ERAdownload

pfile = 'C:/OneDrive/Bitbucket/era5/examples/par/examples.globsim_download'
ERAd = ERAdownload(pfile)
ERAd.retrieve()
ERAd.inventory()