# prop_analysis.py
#
# Run ray tracing and normal mode propagation simulations
# for the sampled atmospheric states and build propagation
# statistics for use in BISL and SpYE.
#
# Philip Blom (pblom@lanl.gov)
#
# Note: this analysis can be run via the stochprop CLI using the following steps:
# -------------------------------------------------------------------------------
#
# stochprop prop build-pgm --atmos-dir samples/winter --output-path prop/winter/winter --src-loc '[30.0, -120.0, 0.0]'  --cpu-cnt 6
# stochprop prop plot --model-file prop/winter/winter.pgm
# 
# stochprop prop build-tlm --atmos-dir samples/winter --output-path prop/winter/winter --freq 0.1  --cpu-cnt 6
# stochprop prop plot --model-file prop/winter/winter_0.100Hz.tlm

# ... (repeat tlm for other frequencies and all models for other seasons)


import os 

import subprocess

import numpy as np

from multiprocessing import cpu_count

from stochprop import propagation

if __name__ == '__main__':

    # ############################## #
    #      Define parameters for     #
    #  running infraga and building  #
    #     propagation statistics     #
    # ############################## #

    sample_dirs = "samples"
    results_dir = "prop"
    season_labels = ["winter", "spring", "summer"]

    freqs = [0.1, 0.2, 0.5, 1.0]
    inclinations = [2.5, 45.0, 2.5]
    azimuths = [0.0, 360.0, 6.0]
    src_loc = [30.0, -120, 0.0]

    cpu_cnt = int(cpu_count() / 2)

    # ########################################## #
    #   Build directories for results, run ray   #
    #      tracing analysis, and build path      #
    #            geometry statistics             #
    # ########################################## #
    if not os.path.isdir(results_dir):
        subprocess.call("mkdir " + results_dir, shell=True)

    for season in season_labels:
        if not os.path.isdir(results_dir + "/" + season):
            subprocess.call("mkdir " + results_dir + "/" + season, shell=True)

        results_id = results_dir + "/" + season + "/" + season

        propagation.run_infraga(sample_dirs + "/" + season, results_id + ".arrivals.dat", cpu_cnt=cpu_cnt, geom="sph", inclinations=inclinations, azimuths=azimuths, src_loc=src_loc)

        pgm = propagation.PathGeometryModel()
        pgm.build(results_id + ".arrivals.dat", results_id + ".pgm", show_fits=False, geom="sph", src_loc=src_loc)
        pgm.load(results_id + ".pgm", smooth=True)
        pgm.display(file_id=(results_id), subtitle=season)

        for fn in freqs:
            propagation.run_modess(sample_dirs + "/" + season, results_id + "_%.3f" % fn + "Hz", azimuths=azimuths, freq=fn, clean_up=True, cpu_cnt=cpu_cnt)

            tlm = propagation.TLossModel()
            tlm.build(results_id + "_%.3f" %fn + "Hz.nm", results_id + "_%.3f" %fn + "Hz.tlm", show_fits=False)
            tlm.load(results_id + "_%.3f" %fn + "Hz.tlm")
            tlm.display(file_id=(results_dir + "/" + season + "/" + season + "_%.3f" %fn), title=("Transmission Loss Statistics" + '\n' + season + ", %.3f" %fn + " Hz"))
