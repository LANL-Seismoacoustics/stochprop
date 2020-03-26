# run_stochprop_eofs.py
#
# Cosntruct the EOF coefficient statistics
# for a set of profiles and generate a set
# of seasonal atmosphere samples
#
# Philip Blom (pblom@lanl.gov)

import os

from stochprop import eofs

if __name__ == '__main__':

    # ######################### #
    #   Define parameters for   #
    #   analysis of a g2s file  #
    # ######################### #

    prof_path = "profs/01/g2stxt_2010010100_39.7393_-104.9900.dat"
    eofs_path = "eofs/example"

    fit_eof_cnts = [1, 10, 25, 50]

    uncertainty = 15.0
    eof_cnt, eof_max = 50, 50
    smpl_cnt = 10

    # ####################################### #
    #  Fit the profile using limited numbers  #
    #    of EOFs and perturb the atmosphere   #
    #          with some uncertainty          #
    # ####################################### #
    if not os.path.isfile(eofs_path + "-singular_values.dat"):
        print("eofs/" + run_id + "-singular_values.dat doesn't exist  --->  run eof_analysis.py first to build EOFs")
    else:
        for N in fit_eof_cnts:
            eofs.fit_atmo(prof_path, eofs_path, "eof_fit-N=" + str(N) + ".met", eof_cnt=N)

        eofs.perturb_atmo(prof_path, eofs_path, "eof_perturb", uncertainty=uncertainty, eof_max=eof_max, eof_cnt=eof_cnt, sample_cnt=smpl_cnt)




