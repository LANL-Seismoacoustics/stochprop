# gen-eofs.py
#
# Cycle through .met files and extract wind profiles to
# construct empirical orthogonal function basis
#
# Philip Blom (pblom@lanl.gov)

import sys
import os
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
# months = ["08"]

# load profile to be fit
for month in months:
    for day_index in range(31):
        day = "%02d" % (day_index + 1)
        for hour_index in range(4):
            hour = "%02d" % (hour_index * 6)

            datetime_id = "2014" + month + day + hour
            filename = "Profiles/2014/" + datetime_id + ".met"

            if os.path.isfile(filename):
                profile = np.loadtxt(filename)
                dTdz = np.gradient(profile[:,1]) / np.gradient(profile[:,0])
                if max(abs(dTdz)) > 50.0:
                    print(filename)

                    plt.figure(1)
                    plt.clf()
                    plt.plot(profile[:, 1], profile[:, 0])
                    plt.plot(dTdz, profile[:, 0])
                    plt.title(filename)
                    plt.show()



