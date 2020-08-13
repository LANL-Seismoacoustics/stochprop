# run_stochprop_eofs.py
#
# Cosntruct the EOF coefficient statistics
# for a set of profiles and generate a set
# of seasonal atmosphere samples
#
# Philip Blom (pblom@lanl.gov)

import subprocess

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
    azimuths = [-180.0, 180.0, 10.0]
    src_loc = [30.0, -120, 0.0]

    cpu_cnt = 7

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

        propagation.run_infraga(sample_dirs + "/" + season, results_dir + "/" + season + "/" + season + ".arrivals.dat", cpu_cnt=cpu_cnt, geom="sph", inclinations=[5.0, 45.0, 1.5], azimuths=azimuths, src_loc=src_loc)

        pgm = propagation.PathGeometryModel()
        pgm.build(results_dir + "/" + season + "/" + season + ".arrivals.dat", results_dir + "/" + season + "/" + season + ".pgm", show_fits=False, geom="sph", src_loc=src_loc)
        pgm.load(results_dir + "/" + season + "/" + season + ".pgm", smooth=True)
        pgm.display(file_id=(results_dir + "/" + season + "/" + season), subtitle=season)

        '''
        propagation.run_modess(sample_dirs + "/" + season, results_dir + "/" + season + "/" + season, azimuths=azimuths, freqs=freqs)
        tlm = propagation.TLossModel()
        for fn in np.logspace(np.log10(freqs[0]), np.log10(freqs[1]), freqs[2]):
            tlm.build(results_dir + "/" + season + "/" + season + "_%.3f" %fn + ".nm", results_dir + "/" + season + "/" + season + "_%.3f" %fn + ".tlm", show_fits=False)
            tlm.load(results_dir + "/" + season + "/" + season + "_%.3f" %fn + ".tlm")
            tlm.display(file_id=(results_dir + "/" + season + "/" + season + "_%.3f" %fn), title=("Transmission Loss Statistics" + '\n' + season + ", %.3f" %fn + " Hz"))
        '''
