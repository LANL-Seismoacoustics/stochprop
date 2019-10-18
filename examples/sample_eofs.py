# run_stochprop_eofs.py
#
# Cosntruct the EOF coefficient statistics
# for a set of profiles and generate a set
# of seasonal atmosphere samples
#
# Philip Blom (pblom@lanl.gov)

import sys
import os
import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt

import pathos.multiprocessing as mp
from multiprocessing import cpu_count

from stochprop import eofs
from infrapy.association import hjl

if __name__ == '__main__':

    # ######################### #
    #   Define parameters for   #
    #    atmosphere sampling    #
    # ######################### #

    loc_id = "label"

    labels = ["winter",
              "spring",
              "summer",
              "fall"]
              
    season_months = [["10", "11", "12", "01", "02", "03"],
                     ["04", "05"],
                     ["06", "07"],
                     ["08", "09"]]
    
    # ################################ #
    #  Load coefficients and generate  #
    #      samples for each season     #
    # ################################ #
    for nS in range(4):
        print('\n' + "Building samples for " + labels[nS])
        coeffs = np.load(loc_id + "/coeffs/" + loc_id + "_" + season_months[nS][0] + "-coeffs.npy")
        for M in season_months[nS][1:]:
            coeffs = np.vstack((coeffs, np.load(loc_id + "/coeffs/" + loc_id + "_" + M + "-coeffs.npy")))

        eofs.sample_atmo(coeffs, loc_id + "/eofs/" + loc_id, loc_id + "/samples/" + labels[nS] + "/" + loc_id + "-" + labels[nS], eof_cnt=100, prof_cnt=250)
        eofs.maximum_likelihood_profile(coeffs, loc_id + "/eofs/" + loc_id, loc_id + "/samples/" + labels[nS] + "/" + loc_id + "-" + labels[nS], eof_cnt=100)




