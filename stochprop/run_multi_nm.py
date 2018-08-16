# prog_bar.py
#
# Progress bar methods to build, increment, and close
# a progress bar that looks like [>>>>>>>>  ]

import os

from .prog_bar import *

ncpa_prop_dir = "/Users/pblom/Research/Coding/Packages/ncpaprop-1.3.2/ncpaprop-1.3.2/"

def run_nm(profile, azimuth, freq, id, prog_step=0, grnd_elev = 0.0):
    os.system(ncpa_prop_dir + "bin/Modess --atmosfile " + profile + " --atmosfileorder ztuvdp --skiplines 0 --azimuth " + str(azimuth) + " --freq " + str(freq) + " --zground_km " + str(grnd_elev) + " > /dev/null")
    os.system("mv tloss_1d-%.3f.nm result." % freq + str(id) + ".dat")
    os.system("mv tloss_1d-%.3f.lossless.nm result.lossless." % freq + str(id) + ".dat")
    increment_progress_bar(n=prog_step)

def run_nm_wrapper(args): return run_nm(*args)





