# gen-eofs.py
#
# Cycle through .met files and extract wind profiles to
# construct empirical orthogonal function basis
#
# Philip Blom (pblom@lanl.gov)

import sys
import os
import numpy as np

from code import eofs

# ######################### #
#   Define parameters for   #
#      EOF construction     #
# ######################### #
profs_path = "../Locations/DPRK_Priors/Profiles-old/"
eofs_path = "DPRK/eofs/DPRK"
coeffs_path = "DPRK/coeffs/DPRK"
seasonality_path = "DPRK/results/DPRK"

overlap_file = None
# overlap_file = seasonality_path + "-seasonality.npy"

# eofs.compute_svd(profs_path, eofs_path)
# eofs.compute_coeffs(profs_path, eofs_path, coeffs_path)
eofs.cluster_seasonality(coeffs_path, eofs_path, seasonality_path, overlap_file=overlap_file)




