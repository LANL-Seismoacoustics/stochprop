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
    #      EOF construction     #
    # ######################### #
    loc_id = "label"
    year_lims = [2010, 2016]
    clustering_thresh = 0.2

    prof_prefix = "g2stxt_"

    # ######################### #
    #  Build atmosphere matrix  #
    #   for EOF construction    #
    # ######################### #
    if os.path.isfile(loc_id + "/eofs/" + loc_id + "-singular_values.dat"):
        print(loc_id + "/eofs/" + loc_id + "-singular_values.dat exists" + '\n\t' + "Skipping construction of EOFs")
    else:
        print("Building atmosphere matrix for all profiles...")
        A0, z0 = eofs.build_atmo_matrix(loc_id + "/profs/{:02d}/".format(1), prof_prefix + "*")
        A_list = [A0]
        for n in range(1, 12):
            A_temp, z_temp = eofs.build_atmo_matrix(loc_id + "/profs/{:02d}/".format(n + 1), prof_prefix + "*", ref_alts=z0)
            if np.any(A_temp):
                A_list += [A_temp]
    
        A = A_list[0]
        for n in range(1, len(A_list)):
            A = np.vstack((A, A_list[n]))
    
        print("Computing SVD...")
        eofs.compute_svd(A, z0, loc_id + "/eofs/" + loc_id)

    # ################################### #
    #  Compute EOF coefficients for each  #
    #     month and build season stats    #
    # ################################### #
    if os.path.isfile(loc_id + "/coeffs/" + loc_id + "_12-coeffs.npy"):
        print(loc_id + "/coeffs/" + loc_id + "_12-coeffs.npy exists" + '\n\t' + "Skipping calcuation of EOF coefficients")
    else:
        print('\n' + "Computing EOF coefficients for each month...")
        pl = mp.ProcessingPool(int(cpu_count() - 1))
        for m in range(12):
            Am, zm = eofs.build_atmo_matrix(loc_id + "/profs/{:02d}/".format(m + 1), prof_prefix + "*")
            coeffs = eofs.compute_coeffs(Am, zm, loc_id + "/eofs/" + loc_id, loc_id + "/coeffs/" + loc_id + "_{:02d}".format(m + 1), eof_cnt=100, pool=pl)
        pl.close()
        pl.terminate()

    if os.path.isfile(loc_id + "/coeffs/" + loc_id + "-overlap.npy"):
        print(loc_id + "/coeffs/" + loc_id + "-overlap.npy exists" + '\n\t' + "Skipping calculation of seasonality")
    else:
        coeffs = [0] * 12
        for m in range(12):
            coeffs[m] = np.load(loc_id + "/coeffs/" + loc_id + "_{:02d}-coeffs.npy".format(m + 1))

        overlap = eofs.compute_overlap(coeffs)
        np.save(loc_id + "/coeffs/" + loc_id + "-overlap", overlap)

    overlap = np.load(loc_id + "/coeffs/" + loc_id + "-overlap.npy")
    eof_weight = np.loadtxt(loc_id + "/eofs/" + loc_id + "-singular_values.dat")[:, 1]
    overlap_proj = np.average(overlap, weights=eof_weight[:overlap.shape[0]], axis=0)

    dist_mat = -np.log(overlap_proj)
    np.fill_diagonal(dist_mat, 0.0)
    hjl.cluster(dist_mat, clustering_thresh, show_result=True, file_id=loc_id + "/coeffs/" + loc_id + "-clustering")

    M1, M2 = np.meshgrid(range(12), range(12))
    plt.figure(1, figsize=(4, 4))
    plt.scatter(M1.flatten() + 1, M2.flatten() + 1, c=overlap_proj.flatten(), cmap=cm.Blues, marker='s', s=200.0)
    plt.savefig(loc_id + "/coeffs/" + loc_id + "-seasonality.png", dpi=300)



