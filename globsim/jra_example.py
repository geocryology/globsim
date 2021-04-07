from globsim.download import JRAdownload
from globsim.scale import JRAscale
from globsim.interpolate import JRAinterpolate

from os   import path

# ==== SETTING-UP =============================================================

dir_path = '/Users/bincao/OneDrive/GitHub/globsim/examples/par'
#dir_path = 'C:/OneDrive/GitHub/globsim/examples/par'

pfile = 'examples.globsim_download'
ifile = 'examples.globsim_interpolate'
sfile = 'examples.globsim_scale'

pfile = path.join(dir_path, pfile)
ifile = path.join(dir_path, ifile)
sfile = path.join(dir_path, sfile)

# ==== Download ===============================================================

JRAdownl = JRAdownload(pfile)
JRAdownl.retrieve()

'''
from JRA import RDA
rda = RDA(JRAdownl.username, JRAdownl.password)
#rda.getStatus()
#dsIndex = rda.getDSindex()

ds = '328839'
rda.download(JRAdownl.directory, ds)
JRAdownl.makeNCF(ds)
'''

# ==== Interpolation ==========================================================

JRAinterp = JRAinterpolate(ifile)
#JRAinterp.process()

# ==== Scale ==================================================================
JRAsc = JRAscale(sfile)
JRAsc.process()