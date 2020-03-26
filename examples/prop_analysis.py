# run_stochprop_eofs.py
#
# Cosntruct the EOF coefficient statistics
# for a set of profiles and generate a set
# of seasonal atmosphere samples
#
# Philip Blom (pblom@lanl.gov)

import os

import numpy as np

from stochprop import propagation

# season_labels = ["winter", "spring", "summer", "fall"]
season_labels = ["winter"]


if __name__ == '__main__':

    # ############################## #
    #      Define parameters for     #
    #  running infraga and building  #
    #     propagation statistics     #
    # ############################## #

    sample_dirs = "samples"
    results_dir = "prop"
    freqs = [0.1, 1.0, 10]

    cpu_cnt = 6

    # ########################################## #
    #   Build directories for results, run ray   #
    #      tracing analysis, and build path      #
    #            geometry statistics             #
    # ########################################## #
    if not os.path.isdir(results_dir):
        os.system("mkdir " + results_dir)

    for season in season_labels:
        if not os.path.isdir(results_dir + "/" + season):
            os.system("mkdir " + results_dir + "/" + season)

        '''
        propagation.run_infraga(sample_dirs + "/" + season, results_dir + "/" + season + ".arrivals.dat", cpu_cnt=cpu_cnt, geom="3d", inclinations=[5.0, 45.0, 1.5], azimuths=[-180.0, 180.0, 10.0])

        pgm = propagation.PathGeometryModel()
        pgm.build(results_dir + "/" + season + ".arrivals.dat", results_dir + "/" + season + ".pgm", show_fits=False)
        pgm.load(results_dir + "/" + season + ".pgm", smooth=True)
        pgm.display(file_id=(results_dir + "/" + season))
        '''

        propagation.run_modess(sample_dirs + "/" + season, results_dir + "/" + season + "/" + season, azimuths=[-180.0, 180.0, 10.0], freqs=freqs)

        tlm = propagation.TLossModel()
        for fn in np.logspace(np.log10(freqs[0]), np.log10(freqs[1]), freqs[2]):
            tlm.build(results_dir + "/" + season + "/" + season + "_%.3f" %fn + ".nm", results_dir + "/" + season + "/" + season + "_%.3f" %fn + ".tlm", show_fits=False)
            tlm.load(results_dir + "/" + season + "/" + season + "_%.3f" %fn + ".tlm")
            tlm.display(file_id=(results_dir + "/" + season + "/" + season + "_%.3f" %fn), title=("Transmission Loss Statistics" + '\n' + season + ", %.3f" %fn + " Hz"))

