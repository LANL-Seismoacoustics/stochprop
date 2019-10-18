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

eof_dirs = ["coeffs", "eofs", "profs", "samples"]
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
    
    clustering_thresh = 0.3
    eof_cnt = 50
    
    # update season labels with clustering results for actual analysis
    season_months = [["10", "11", "12", "01", "02", "03"],
                     ["04", "05"],
                     ["06", "07"],
                     ["08", "09"]]
              
    # ######################### #
    #    Build directories to   #
    #      organize results     #
    # ######################### #
    print("Creating directories...")
    for dir in eof_dirs:
        if not os.path.isdir(dir):
            os.system("mkdir " + dir)
    
    for M in range(1, 13):
        if not os.path.isdir("profs/" + "{:02d}".format(M)):
            os.system("mkdir profs/" + "{:02d}".format(M))

    for season in season_labels:
        if not os.path.isdir("samples/" + season):
            os.system("mkdir samples/" + season)

    # UNCOMMENT TO POPULATE DIRECTORIES and QC
    '''
    print('\n' + "Moving profiles into profs directory...")
    for year in range(year_lims[0], year_lims[1] + 1):
        print('\t' + "Moving " + "{:04d}".format(year) + " profiles...")
        for M in range(1, 13):
            os.system("mv " + prof_dir + prof_prefix + "{:04d}".format(year) + "{:02d}".format(M) + "* profs/" + "{:02d}".format(M))

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
        A0, z0 = eofs.build_atmo_matrix("profs/{:02d}/".format(1), prof_prefix + "*")
        A_list = [A0]
        for n in range(1, 12):
            A_temp, z_temp = eofs.build_atmo_matrix("profs/{:02d}/".format(n + 1), prof_prefix + "*", ref_alts=z0)
            if np.any(A_temp):
                A_list += [A_temp]
    
        A = A_list[0]
        for n in range(1, len(A_list)):
            A = np.vstack((A, A_list[n]))
    
        print("Computing SVD...")
        eofs.compute_svd(A, z0, "eofs/" + run_id, eof_cnt=eof_cnt)

    # ################################### #
    #  Compute EOF coefficients for each  #
    #     month and build season stats    #
    # ################################### #
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

        overlap = eofs.compute_overlap(coeffs, eof_cnt=eof_cnt)
        np.save("coeffs/" + run_id + "-overlap", overlap)

    overlap = np.load("coeffs/" + run_id + "-overlap.npy")
    eof_weight = np.loadtxt("eofs/" + run_id + "-singular_values.dat")[:, 1]
    overlap_proj = np.average(overlap, weights=eof_weight[:overlap.shape[0]], axis=0)

    dist_mat = -np.log(overlap_proj)
    np.fill_diagonal(dist_mat, 0.0)
    hjl.cluster(dist_mat, clustering_thresh, file_id="coeffs/" + run_id + "-clustering")

    M1, M2 = np.meshgrid(range(12), range(12))
    plt.figure(1, figsize=(4, 4))
    plt.scatter(M1.flatten() + 1, M2.flatten() + 1, c=overlap_proj.flatten(), cmap=cm.Blues, marker='s', s=200.0)
    plt.savefig("coeffs/" + run_id + "-seasonality.png", dpi=300)

    # ################################ #
    #  Load coefficients and generate  #
    #      samples for each season     #
    # ################################ #
    for nS in range(4):
        print('\n' + "Building samples for " + season_labels[nS])
        coeffs = np.load("coeffs/" + run_id + "_" + season_months[nS][0] + "-coeffs.npy")
        for M in season_months[nS][1:]:
            coeffs = np.vstack((coeffs, np.load("coeffs/" + run_id + "_" + M + "-coeffs.npy")))
        
        eofs.sample_atmo(coeffs, "eofs/" + run_id, "samples/" + season_labels[nS] + "/" + run_id + "-" + season_labels[nS], eof_cnt=eof_cnt, prof_cnt=25)
        eofs.maximum_likelihood_profile(coeffs, "eofs/" + run_id, "samples/" + season_labels[nS] + "/" + run_id + "-" + season_labels[nS], eof_cnt=eof_cnt)


