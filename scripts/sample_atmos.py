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

from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import bisect
from scipy.stats import gaussian_kde

# #################### #
#    Set parameters    #
# #################### #
eof_path, eof_id = "UTTR/eofs/", "UTTR-all"
coeff_path = "UTTR/coeffs/"
output_path = "UTTR/results/profs/hybrid/"

eof_cnt = 50
prof_cnt = 100

'''
# winter (October - February)
months = ["10", "11", "12", "01", "02"]
output_id = "winter"

# early-spring (March)
months = ["03"]
output_id = "early-spring"

# mid-spring (April)
months = ["04"]
output_id = "mid-spring"

# late-spring (May)
months = ["05"]
output_id = "late-spring"

# summer (June - August)
months = ["06", "07", "08"]
output_id = "summer"
'''
# fall (September)
months = ["09"]
output_id = "fall"

##################
# Define methods #
##################
def define_limits(coeff_vals):
    coeff_min, coeff_max = min(coeff_vals), max(coeff_vals)
    return np.array([(coeff_max + coeff_min) / 2.0 - 1.5 * (coeff_max - coeff_min),
                     (coeff_max + coeff_min) / 2.0 + 1.5 * (coeff_max - coeff_min)])

def build_cdf(pdf, lims, pnts=100):
    norm = quad(pdf, lims[0], lims[1])[0]

    x_vals = np.linspace(lims[0], lims[1], pnts)
    cdf_vals = np.empty_like(x_vals)
    for n in range(pnts):
        cdf_vals[n] = quad(pdf, lims[0], x_vals[n])[0] / norm

    return interp1d(x_vals, cdf_vals)

def draw_from_pdf(pdf, lims, cdf=None, size=1):
    if not cdf:
        cdf = build_cdf(pdf, lims)

    def cdf_w_shift(x, d):
        return cdf(x) - d

    rnd_vals = np.random.random(size)
    drawn_vals = np.empty_like(rnd_vals)
    for n in range(size):
        drawn_vals[n] = bisect(cdf_w_shift, lims[0], lims[1], args=(rnd_vals[n],))

    return drawn_vals

def eff_gas_constant(z):
    return 287.058 + 0.75 * (z - 80.0) / (1.0 + np.exp(-(z - 102.5) / 6.6))


# ############################ #
#  Read in coefficient values  #
#    for months of interest    #
# ############################ #

# load mean profile and eofs
print("Loading mean profile info and eofs...")
means = np.loadtxt(eof_path + eof_id + "_means.dat")
T_eofs = np.loadtxt(eof_path + eof_id + "_temperature.eofs")
u_eofs = np.loadtxt(eof_path + eof_id + "_zonal_winds.eofs")
v_eofs = np.loadtxt(eof_path + eof_id + "_merid_winds.eofs")
d_eofs = np.loadtxt(eof_path + eof_id + "_density.eofs")

# load and stack coefficient values for all months of interest
print("Loading coefficient values and using KDE to build pdf's...")
T_coeffs = np.empty((0, T_eofs.shape[1] - 1))
u_coeffs = np.empty((0, u_eofs.shape[1] - 1))
v_coeffs = np.empty((0, v_eofs.shape[1] - 1))
d_coeffs = np.empty((0, d_eofs.shape[1] - 1))

for month in months:
    T_coeffs = np.vstack((T_coeffs, np.load(coeff_path + month + "-T_coeffs.npy")))
    u_coeffs = np.vstack((u_coeffs, np.load(coeff_path + month + "-u_coeffs.npy")))
    v_coeffs = np.vstack((v_coeffs, np.load(coeff_path + month + "-v_coeffs.npy")))
    d_coeffs = np.vstack((d_coeffs, np.load(coeff_path + month + "-d_coeffs.npy")))

# use kernel density estimates to define coefficient pdf's
# and use the mean and span to define the limits for sampling
T_kernels, T_lims = [0] * eof_cnt, np.empty((eof_cnt, 2))
u_kernels, u_lims = [0] * eof_cnt, np.empty((eof_cnt, 2))
v_kernels, v_lims = [0] * eof_cnt, np.empty((eof_cnt, 2))
d_kernels, d_lims = [0] * eof_cnt, np.empty((eof_cnt, 2))

for eof_id in range(eof_cnt):
    T_kernels[eof_id] = gaussian_kde(T_coeffs[:, eof_id])
    u_kernels[eof_id] = gaussian_kde(u_coeffs[:, eof_id])
    v_kernels[eof_id] = gaussian_kde(v_coeffs[:, eof_id])
    d_kernels[eof_id] = gaussian_kde(d_coeffs[:, eof_id])

    T_lims[eof_id] = define_limits(T_coeffs[:, eof_id])
    u_lims[eof_id] = define_limits(u_coeffs[:, eof_id])
    v_lims[eof_id] = define_limits(v_coeffs[:, eof_id])
    d_lims[eof_id] = define_limits(d_coeffs[:, eof_id])

# compute the cumulative density functions for all pdf's
T_cdfs = [0] * eof_cnt
u_cdfs = [0] * eof_cnt
v_cdfs = [0] * eof_cnt
d_cdfs = [0] * eof_cnt
print("Pre-computing cdf's for random sampling...")
for eof_id in range(eof_cnt):
    print('\t' + "eof_id: " + str(eof_id) + "...")
    T_cdfs[eof_id] = build_cdf(T_kernels[eof_id].pdf, T_lims[eof_id])
    u_cdfs[eof_id] = build_cdf(u_kernels[eof_id].pdf, u_lims[eof_id])
    v_cdfs[eof_id] = build_cdf(v_kernels[eof_id].pdf, v_lims[eof_id])
    d_cdfs[eof_id] = build_cdf(d_kernels[eof_id].pdf, d_lims[eof_id])

# Generate a random atmosphere sample
print("Generating profiles...")
profs = np.array([means] * prof_cnt)
for eof_id in range(eof_cnt):
    T_coeff_vals = draw_from_pdf(T_kernels[eof_id].pdf, T_lims[eof_id], cdf=T_cdfs[eof_id], size=prof_cnt)
    u_coeff_vals = draw_from_pdf(u_kernels[eof_id].pdf, u_lims[eof_id], cdf=u_cdfs[eof_id], size=prof_cnt)
    v_coeff_vals = draw_from_pdf(v_kernels[eof_id].pdf, v_lims[eof_id], cdf=v_cdfs[eof_id], size=prof_cnt)
    d_coeff_vals = draw_from_pdf(d_kernels[eof_id].pdf, d_lims[eof_id], cdf=d_cdfs[eof_id], size=prof_cnt)

    for pn in range(prof_cnt):
        profs[pn][:, 1] += T_coeff_vals[pn] * T_eofs[:, eof_id + 1]
        profs[pn][:, 2] += u_coeff_vals[pn] * u_eofs[:, eof_id + 1]
        profs[pn][:, 3] += v_coeff_vals[pn] * v_eofs[:, eof_id + 1]
        profs[pn][:, 4] += d_coeff_vals[pn] * d_eofs[:, eof_id + 1]

# convert density back to physical units (eof analysis is log10 space) and use
# the effective gas constant to compute physical pressure
for pn in range(prof_cnt):
    profs[pn][:, 4] = 10.0**profs[pn][:, 4]
    profs[pn][:, 5] = profs[pn][:, 4] * eff_gas_constant(profs[pn][:, 0]) * profs[pn][:, 1] * 10.0

'''
# Plot the results
print("Plotting profiles...")
fig, ax = plt.subplots(1, 5, sharey=True)
for n in range(3):
    for pn in range(prof_cnt):
        ax[n].plot(profs[pn][:, n + 1], profs[pn][:, 0], linewidth=0.5, color='0.5')
    ax[n].plot(means[:, n + 1], means[:, 0], linewidth=1.0, color='k')
    ax[n].plot(np.average(profs, axis=0)[:, n + 1], np.average(profs, axis=0)[:, 0], linewidth=1.5, color='b')
for n in range(3, 5):
    for pn in range(prof_cnt):
        ax[n].semilogx(profs[pn][:, n + 1], profs[pn][:, 0], linewidth=0.5, color='0.5')
    ax[n].semilogx(10.0**means[:, n + 1], means[:, 0], linewidth=1.0, color='k')
    ax[n].semilogx(np.average(profs, axis=0)[:, n + 1], np.average(profs, axis=0)[:, 0], linewidth=1.5, color='b')
plt.savefig(output_path + output_id + "/" + output_id + ".png", dpi=300)
# plt.show()
'''

# Write the profiles to file
print("Writing profiles to file...")
for pn in range(prof_cnt):
    np.savetxt(output_path + output_id + "/" + output_id + "-" + "%02d" % pn + ".met", profs[pn])

for pn in range(prof_cnt):
    np.savetxt(output_path + output_id + "/" + output_id + "-mean.met", np.average(profs, axis=0))


