# -*- coding: utf-8 -*-

# ==== DOWNLOAD ===============================================================

from era5 import ERAdownload, ERAscale, ERAinterpolate

pfile = "/home/nbrown/globsim/era5test/par/examples.globsim_download"
ifile = "/home/nbrown/globsim/era5test/par/examples.globsim_interpolate"
sfile =  "/home/nbrown/globsim/era5test/par/examples.globsim_scale"
#ERAd = ERAdownload(pfile)
#ERAd.retrieve()
#ERAd.inventory()

# ERAint = ERAinterpolate(ifile)
# ERAint.process()

ERAsc = ERAscale(sfile)
ERAsc.process()