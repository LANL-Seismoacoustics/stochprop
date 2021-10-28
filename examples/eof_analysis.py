# run_stochprop_eofs.py
#
# Cosntruct the EOF coefficient statistics
# for a set of profiles and generate a set
# of seasonal atmosphere samples
#
# Philip Blom (pblom@lanl.gov)

import os
import subprocess
import numpy as np

import pathos.multiprocessing as mp
from multiprocessing import cpu_count

from stochprop import eofs

eof_dirs = ["coeffs", "eofs", "samples"]
season_labels = ["winter", "spring", "summer", "fall"]

if __name__ == '__main__':

    # ######################### #
    #   Define parameters for   #
    #      EOF construction     #
    # ######################### #
    prof_dir = "dir/of/g2s/"
    prof_prefix = "g2stxt_"
    year_lims = [2010, 2016]

    run_id = "example"

    eof_cnt = 50
    smpl_cnt = 25

    # update season labels with clustering results for actual analysis
    season_months = [["10", "11", "12", "01", "02", "03"],
                     ["04", "05"],
                     ["06", "07", "08"],
                     ["09"]]

    # ######################### #
    #     Remove results and    #
    #  directories of examples  #
    # ######################### #
    subprocess.call("rm -rf eofs/ coeffs/ samples/", shell=True)

    # ######################### #
    #    Build directories to   #
    #      organize results     #
    # ######################### #
    print("Creating directories...")
    for dir in eof_dirs:
        if not os.path.isdir(dir):
            subprocess.call("mkdir " + dir, shell=True)

    for season in season_labels:
        if not os.path.isdir("samples/" + season):
            subprocess.call("mkdir samples/" + season, shell=True)

    # UNCOMMENT TO POPULATE DIRECTORIES and QC
    '''
    print('\n' + "Moving profiles into profs directory...")
    for year in range(year_lims[0], year_lims[1] + 1):
        print('\t' + "Moving " + "{:04d}".format(year) + " profiles...")
        for M in range(1, 13):
            command = "mv " + prof_dir + prof_prefix + "{:04d}".format(year) + "{:02d}".format(M) + "* profs/" + "{:02d}".format(M)
            subprocess.call(command, shell=True)

    for M in range(1, 13):
        eofs.profs_check("profs/{:02d}/".format(M), prof_prefix + "*")
    '''

    # ######################### #
    #  Build atmosphere matrix  #
    #   for EOF construction    #
    # ######################### #
    if os.path.isfile("eofs/" + run_id + "-singular_values.dat"):
        print("eofs/" + run_id + "-singular_values.dat exists  --->  Skipping construction of EOFs")
    else:
        print("Building atmosphere matrix for all profiles...")
        A, z0 = eofs.build_atmo_matrix("profs/{:02d}/".format(1), prof_prefix + "*")
        for n in range(1, 12):
            A_temp, _ = eofs.build_atmo_matrix("profs/{:02d}/".format(n + 1), prof_prefix + "*", ref_alts=z0)
            A = np.vstack((A, A_temp))
        eofs.compute_eofs(A, z0, "eofs/" + run_id, eof_cnt=eof_cnt)

    # #################################### #
    #  Compute EOFs coefficients for each  #
    #  month and build season definitions  #
    # #################################### #
    if os.path.isfile("coeffs/" + run_id + "_12-coeffs.npy"):
        print("coeffs/" + run_id + "_12-coeffs.npy exists  --->  Skipping calcuation of EOF coefficients")
    else:
        print('\n' + "Computing EOF coefficients for each month...")
        pl = mp.ProcessingPool(int(cpu_count() - 1))
        for m in range(12):
            Am, zm = eofs.build_atmo_matrix("profs/{:02d}/".format(m + 1), prof_prefix + "*")
            coeffs = eofs.compute_coeffs(Am, zm, "eofs/" + run_id, "coeffs/" + run_id + "_{:02d}".format(m + 1), eof_cnt=eof_cnt, pool=pl)
        pl.close()
        pl.terminate()

    if os.path.isfile("coeffs/" + run_id + "-overlap.npy"):
        print("coeffs/" + run_id + "-overlap.npy exists  --->  Skipping calculation of seasonality")
    else:
        coeffs = [0] * 12
        for m in range(12):
            coeffs[m] = np.load("coeffs/" + run_id + "_{:02d}-coeffs.npy".format(m + 1))

        overlap = eofs.compute_overlap(coeffs, "eofs/" + run_id, eof_cnt=eof_cnt, method="kde") 
        np.save("coeffs/" + run_id + "-overlap_kde", overlap)

    eofs.compute_seasonality("coeffs/" + run_id + "-overlap_kde.npy", "coeffs/" + run_id + "_kde")

    # ################################ #
    #  Load coefficients and generate  #
    #      samples for each season     #
    # ################################ #
    for nS in range(4):
        print('\n' + "Building samples for " + season_labels[nS])
        coeffs = np.load("coeffs/" + run_id + "_" + season_months[nS][0] + "-coeffs.npy")
        for M in season_months[nS][1:]:
            coeffs = np.vstack((coeffs, np.load("coeffs/" + run_id + "_" + M + "-coeffs.npy")))

        eofs.sample_atmo(coeffs, "eofs/" + run_id, "samples/" + season_labels[nS] + "/" + run_id + "-" + season_labels[nS], eof_cnt=eof_cnt, prof_cnt=smpl_cnt, coeff_label=season_labels[nS], output_mean=True)
        eofs.maximum_likelihood_profile(coeffs, "eofs/" + run_id, "samples/" + run_id + "-" + season_labels[nS], eof_cnt=eof_cnt, coeff_label=season_labels[nS])
