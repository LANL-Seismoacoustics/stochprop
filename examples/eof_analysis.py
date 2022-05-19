# run_stochprop_eofs.py
#
# Cosntruct the EOF coefficient statistics
# for a set of profiles and generate a set
# of seasonal atmosphere samples
#
# Philip Blom (pblom@lanl.gov)
#
#
# Note: these steps can all be run via the CLI
# --------------------------------------------
# stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_low_alt --max-alt 80.0 --eof-cnt 50
# stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_low_alt --coeff-path coeffs/example_low_alt_01 --month-selection '01' --eof-cnt 50
# stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_low_alt --coeff-path coeffs/example_low_alt_02 --month-selection '02' --eof-cnt 50
# ... (repeat for each month)
#
# stochprop eof seasonality --eofs-path eofs/example_low_alt --coeff-path coeffs/example_low_alt --eof-cnt 50
#
# stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_winter --month-selection '[10, 11, 12, 01, 02, 03]' --eof-cnt 50
# stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --month-selection '[10, 11, 12, 01, 02, 03]' --eof-cnt 50
# stochprop eof sample --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --sample-path samples/winter/example_winter --sample-cnt 50 --eof-cnt 50
# 
# ... (repeat for spring and summer with modified --month-selection definitions)

import os
import subprocess
import numpy as np

from stochprop import eofs

eof_dirs = ["coeffs", "eofs", "samples"]
season_labels = ["winter", "spring", "summer"]

if __name__ == '__main__':

    # ######################### #
    #   Define parameters for   #
    #      EOF construction     #
    # ######################### #
    prof_dir = "dir/of/g2s/"
    prof_prefix = "g2stxt_"

    run_id = "example"

    eof_cnt = 50
    smpl_cnt = 25

    # update season labels with clustering results for actual analysis
    season_months = [["10", "11", "12", "01", "02", "03"],
                     ["04", "09"],
                     ["05", "06", "07", "08"]]

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

    # ######################### #
    #  Build atmosphere matrix  #
    #   for EOF construction    #
    # ######################### #
    print("Building atmosphere matrix for all profiles...")
    A, z0 = eofs.build_atmo_matrix("profs/", "*.dat", max_alt=80.0)
    eofs.compute_eofs(A, z0, "eofs/" + run_id + "_low_alt", eof_cnt=eof_cnt)

    # #################################### #
    #  Compute EOFs coefficients for each  #
    #  month and build season definitions  #
    # #################################### #
    print('\n' + "Computing EOF coefficients for each month...")
    for m in range(12):
        Am, zm = eofs.build_atmo_matrix("profs/", "*.dat", months = ['%02d' % (m + 1)])
        coeffs = eofs.compute_coeffs(Am, zm, "eofs/" + run_id + "_low_alt", "coeffs/" + run_id + "_low_alt" + "_{:02d}".format(m + 1), eof_cnt=eof_cnt)

    coeffs = [0] * 12
    for m in range(12):
        coeffs[m] = np.load("coeffs/" + run_id + "_low_alt" + "_{:02d}-coeffs.npy".format(m + 1))

    overlap = eofs.compute_overlap(coeffs, "eofs/" + run_id + "_low_alt", eof_cnt=eof_cnt, method="kde") 
    np.save("coeffs/" + run_id + "_low_alt" + "-overlap", overlap)
    eofs.compute_seasonality("coeffs/" + run_id + "_low_alt" + "-overlap_kde.npy", "coeffs/" + run_id)

    # ################################ #
    #  Load coefficients and generate  #
    #      samples for each season     #
    # ################################ #
    for nS in range(len(season_labels)):
        print('\n' + "Building EOF for " + season_labels[nS])
        A, z0 = eofs.build_atmo_matrix("profs/", "*.dat", months=season_months[nS])
        eofs.compute_eofs(A, z0, "eofs/" + run_id + "_" + season_labels[nS], eof_cnt=eof_cnt)
        coeffs = eofs.compute_coeffs(A, z0, "eofs/" + run_id+ "_" + season_labels[nS], "coeffs/" + run_id + "_" + season_labels[nS], eof_cnt=eof_cnt)

        print('\n' + "Sampling atmospheric state for " + season_labels[nS])
        eofs.sample_atmo(coeffs, "eofs/" + run_id + "_" + season_labels[nS], "samples/" + season_labels[nS] + "/" + run_id + "-" + season_labels[nS], eof_cnt=eof_cnt, prof_cnt=smpl_cnt, coeff_label=season_labels[nS], output_mean=True)
        eofs.maximum_likelihood_profile(coeffs, "eofs/" + run_id+ "_" + season_labels[nS], "samples/" + run_id + "-" + season_labels[nS], eof_cnt=eof_cnt, coeff_label=season_labels[nS])

