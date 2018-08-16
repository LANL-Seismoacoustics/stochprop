# eofs.py
#
# Methods to construct a singular value decomposition (SVD) from a set of
# G2S profiles and define empirical orthogonal functions.  Also included
# are methods to compute coefficient values for the profiles and define
# statistics of those coefficients to randomly sample the atmosphere state
#
# Philip Blom (pblom@lanl.gov)

import sys
import os
import numpy as np
import calendar

import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import bisect
from scipy.stats import gaussian_kde

from infrapy.association import hjl

def profile_qc(g2s_file):
    profile = np.loadtxt(g2s_file)
    dTdz = np.gradient(profile[:,1]) / np.gradient(profile[:,0])
    if max(abs(dTdz)) > 50.0:
        return False
    else:
        return True

def profs_check(profs_path, months=range(1, 13), years = range(2007, 2015), plot_bad_profs=False):
    for month_index in months:
        month = "%02d" % month_index
        print("Loading profiles for (month):", calendar.month_name[month_index])
        for year_index in years:
            year = "%04d" % year_index
            for day_index in range(32):
                day = "%02d" % (day_index + 1)
                for hour_index in range(4):
                    hour = "%02d" % (hour_index * 6)

                    datetime_id = year + month + day + hour
                    filename = profs_path + year + "/" + datetime_id + ".met"

                    if os.path.isfile(filename):
                        profile = np.loadtxt(filename)
                        dTdz = np.gradient(profile[:,1]) / np.gradient(profile[:,0])
                        if max(abs(dTdz)) > 50.0:
                            print(filename)
                            os.system("mv " + filename + " " + profs_path + "bad_profs")

                            if plot_bad_profs:
                                plt.figure(1)
                                plt.clf()
                                plt.plot(profile[:, 1], profile[:, 0])
                                plt.plot(dTdz, profile[:, 0])
                                plt.title(filename)
                                plt.show()

def compute_svd(profs_path, output_path, eof_cnt=50, g2s_file_len=901, months=range(1, 13), years = range(2007, 2015)):
    print("-" * 50 + '\n' + "Computing SVD to build EOFs")
    # Prep the arrays to hold the profile information
    T_vals = np.empty((0, g2s_file_len))
    u_vals = np.empty((0, g2s_file_len))
    v_vals = np.empty((0, g2s_file_len))
    d_vals = np.empty((0, g2s_file_len))
    p_vals = np.empty((0, g2s_file_len))

    # Read in data from G2S profiles
    for month_index in months:
        month = "%02d" % month_index
        print('\t' + "Loading profiles for " + calendar.month_name[month_index])
        for year_index in years:
            year = "%04d" % year_index
            for day_index in range(32):
                day = "%02d" % (day_index + 1)
                for hour_index in range(4):
                    hour = "%02d" % (hour_index * 6)

                    datetime_id = year + month + day + hour
                    filename = profs_path + year + "/" + datetime_id + ".met"

                    if os.path.isfile(filename):
                        profile = np.loadtxt(filename)
                        alt_vals = profile[:, 0]
                        T_vals = np.vstack((T_vals, profile[:, 1]))
                        u_vals = np.vstack((u_vals, profile[:, 2]))
                        v_vals = np.vstack((v_vals, profile[:, 3]))
                        d_vals = np.vstack((d_vals, np.log10(profile[:, 4])))
                        p_vals = np.vstack((p_vals, np.log10(profile[:, 5])))

    print('\t' + "Building SVD expansion...")
    T_mean = np.mean(T_vals, axis=0)
    u_mean = np.mean(u_vals, axis=0)
    v_mean = np.mean(v_vals, axis=0)
    d_mean = np.mean(d_vals, axis=0)
    p_mean = np.mean(p_vals, axis=0)

    _, T_singular_vals, T_eofs = np.linalg.svd(T_vals - np.array([T_mean] * T_vals.shape[0]), full_matrices=False)
    _, u_singular_vals, u_eofs = np.linalg.svd(u_vals - np.array([u_mean] * u_vals.shape[0]), full_matrices=False)
    _, v_singular_vals, v_eofs = np.linalg.svd(v_vals - np.array([v_mean] * v_vals.shape[0]), full_matrices=False)
    _, d_singular_vals, d_eofs = np.linalg.svd(d_vals - np.array([d_mean] * d_vals.shape[0]), full_matrices=False)
    _, p_singular_vals, p_eofs = np.linalg.svd(p_vals - np.array([p_mean] * p_vals.shape[0]), full_matrices=False)

    print('\t' + "Normalizing EOFs...")
    for n in range(eof_cnt):
        T_eofs[n] = T_eofs[n] / np.sqrt(quad(interp1d(alt_vals, T_eofs[n]**2), alt_vals[0], alt_vals[-1])[0])
        u_eofs[n] = u_eofs[n] / np.sqrt(quad(interp1d(alt_vals, u_eofs[n]**2), alt_vals[0], alt_vals[-1])[0])
        v_eofs[n] = v_eofs[n] / np.sqrt(quad(interp1d(alt_vals, v_eofs[n]**2), alt_vals[0], alt_vals[-1])[0])
        d_eofs[n] = d_eofs[n] / np.sqrt(quad(interp1d(alt_vals, d_eofs[n]**2), alt_vals[0], alt_vals[-1])[0])
        p_eofs[n] = p_eofs[n] / np.sqrt(quad(interp1d(alt_vals, p_eofs[n]**2), alt_vals[0], alt_vals[-1])[0])

    # Write the mean profiles and desired number of EOFs to file
    print('\t' + "Writing to file...")
    file_out = open(output_path + "-singular_values.dat", 'w')
    for n in range(len(T_singular_vals)):
        print(n, '\t', T_singular_vals[n],
                 '\t', u_singular_vals[n],
                 '\t', v_singular_vals[n],
                 '\t', d_singular_vals[n],
                 '\t', p_singular_vals[n], file=file_out)
    file_out.close()

    file_out = open(output_path + "-means.dat", 'w')
    for n in range(len(alt_vals)):
        print(alt_vals[n], '\t', T_mean[n], '\t', u_mean[n], '\t', v_mean[n], '\t', d_mean[n], '\t', p_mean[n], file=file_out)
    file_out.close()

    file_out = open(output_path + "-temperature.eofs", 'w')
    for n in range(len(alt_vals)):
        print(alt_vals[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', T_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()

    file_out = open(output_path + "-zonal_winds.eofs", 'w')
    for n in range(len(alt_vals)):
        print(alt_vals[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', u_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()

    file_out = open(output_path + "-merid_winds.eofs", 'w')
    for n in range(len(alt_vals)):
        print(alt_vals[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', v_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()

    file_out = open(output_path + "-density.eofs", 'w')
    for n in range(len(alt_vals)):
        print(alt_vals[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', d_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()

    file_out = open(output_path + "-pressure.eofs", 'w')
    for n in range(len(alt_vals)):
        print(alt_vals[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', p_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()


def compute_coeffs(profs_path, eofs_path, output_path, months=range(1, 13), years = range(2007, 2015)):
    # load means and eofs
    print("-" * 50 + '\n' + "Loading eofs...")
    means = np.loadtxt(eofs_path + "-means.dat")
    T_eofs = np.loadtxt(eofs_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eofs_path + "-density.eofs")

    # load profile to be fit
    for month_index in months:
        month = "%02d" % month_index
        for hour_index in range(4):
            hour = "%02d" % (hour_index * 6)

            T_coeffs = np.empty((0, T_eofs.shape[1] - 1))
            u_coeffs = np.empty((0, u_eofs.shape[1] - 1))
            v_coeffs = np.empty((0, v_eofs.shape[1] - 1))
            d_coeffs = np.empty((0, d_eofs.shape[1] - 1))

            for year_index in years:
                year = "%04d" % year_index
                for day_index in range(32):
                    day = "%02d" % (day_index + 1)

                    datetime_id = year + month + day + hour
                    filename = profs_path + year + "/" + datetime_id + ".met"

                    if os.path.isfile(filename):
                        print('\t' + "Computing EOF coefficients for " + filename)
                        profile = np.loadtxt(filename)

                        # apply the shifts from means (note that density is log-scaled first)
                        T_diff = profile[:, 1] - means[:, 1]
                        u_diff = profile[:, 2] - means[:, 2]
                        v_diff = profile[:, 3] - means[:, 3]
                        d_diff = np.log10(profile[:, 4]) - means[:, 4]

                        # interpolate and compute the projection onto each of the eofs
                        T_coeffs_temp = np.empty(T_eofs.shape[1] - 1)
                        u_coeffs_temp = np.empty(u_eofs.shape[1] - 1)
                        v_coeffs_temp = np.empty(v_eofs.shape[1] - 1)
                        d_coeffs_temp = np.empty(d_eofs.shape[1] - 1)

                        for n in range(T_eofs.shape[1] - 1):
                            T_coeffs_temp[n] = quad(interp1d(profile[:, 0], T_diff * T_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            u_coeffs_temp[n] = quad(interp1d(profile[:, 0], u_diff * u_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            v_coeffs_temp[n] = quad(interp1d(profile[:, 0], v_diff * v_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            d_coeffs_temp[n] = quad(interp1d(profile[:, 0], d_diff * d_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]

                        T_coeffs = np.vstack((T_coeffs, T_coeffs_temp))
                        u_coeffs = np.vstack((u_coeffs, u_coeffs_temp))
                        v_coeffs = np.vstack((v_coeffs, v_coeffs_temp))
                        d_coeffs = np.vstack((d_coeffs, d_coeffs_temp))

            # store coefficient values as numpy array file
            np.save(output_path + "/" + hour + "00/" + month + "-T_coeffs", T_coeffs)
            np.save(output_path + "/" + hour + "00/" + month + "-u_coeffs", u_coeffs)
            np.save(output_path + "/" + hour + "00/" + month + "-v_coeffs", v_coeffs)
            np.save(output_path + "/" + hour + "00/" + month + "-d_coeffs", d_coeffs)

def cluster_seasonality(coeff_path, eofs_path, output_path, cluser_thresh=0.15, eof_cnt=50, eof_cluster_cnt=15, overlap_file=None):
    print("-" * 50 + '\n' + "Computing coefficient PDF clustering...")
    if overlap_file:
        print('\t' + "Loading overlap values from file...")
        overlap_values = np.load(overlap_file)
    else:
        print('\t' + "Computing overlap values from coefficient distributions...")
        overlap_values = np.empty((eof_cnt, 4, 12, 12))
        for eof_id in range(eof_cnt):
            print('\t\t' + "Analyzing seasonality of EOF id: " + str(eof_id))

            # define coefficient pdfs for each month
            T_kernels = [0] * 12
            u_kernels = [0] * 12
            v_kernels = [0] * 12
            d_kernels = [0] * 12

            T_lims = np.array([np.inf, -np.inf])
            u_lims = np.array([np.inf, -np.inf])
            v_lims = np.array([np.inf, -np.inf])
            d_lims = np.array([np.inf, -np.inf])

            for month_index in range(1, 13):
                month = "%02d" % month_index

                T_coeffs = np.load(coeff_path + month + "-T_coeffs.npy")
                u_coeffs = np.load(coeff_path + month + "-u_coeffs.npy")
                v_coeffs = np.load(coeff_path + month + "-v_coeffs.npy")
                d_coeffs = np.load(coeff_path + month + "-d_coeffs.npy")

                T_kernels[int(month) - 1] = gaussian_kde(T_coeffs[:, eof_id])
                u_kernels[int(month) - 1] = gaussian_kde(u_coeffs[:, eof_id])
                v_kernels[int(month) - 1] = gaussian_kde(v_coeffs[:, eof_id])
                d_kernels[int(month) - 1] = gaussian_kde(d_coeffs[:, eof_id])

                T_lims = np.array([min(T_lims[0], np.min(T_coeffs[:, eof_id])), max(T_lims[1], np.max(T_coeffs[:, eof_id]))])
                u_lims = np.array([min(u_lims[0], np.min(u_coeffs[:, eof_id])), max(u_lims[1], np.max(u_coeffs[:, eof_id]))])
                v_lims = np.array([min(v_lims[0], np.min(v_coeffs[:, eof_id])), max(v_lims[1], np.max(v_coeffs[:, eof_id]))])
                d_lims = np.array([min(d_lims[0], np.min(d_coeffs[:, eof_id])), max(d_lims[1], np.max(d_coeffs[:, eof_id]))])

            T_lims = np.array([(T_lims[1] + T_lims[0]) / 2.0 - 1.5 * (T_lims[1] - T_lims[0]), (T_lims[1] + T_lims[0]) / 2.0 + 1.5 * (T_lims[1] - T_lims[0])])
            u_lims = np.array([(u_lims[1] + u_lims[0]) / 2.0 - 1.5 * (u_lims[1] - u_lims[0]), (u_lims[1] + u_lims[0]) / 2.0 + 1.5 * (u_lims[1] - u_lims[0])])
            v_lims = np.array([(v_lims[1] + v_lims[0]) / 2.0 - 1.5 * (v_lims[1] - v_lims[0]), (v_lims[1] + v_lims[0]) / 2.0 + 1.5 * (v_lims[1] - v_lims[0])])
            d_lims = np.array([(d_lims[1] + d_lims[0]) / 2.0 - 1.5 * (d_lims[1] - d_lims[0]), (d_lims[1] + d_lims[0]) / 2.0 + 1.5 * (d_lims[1] - d_lims[0])])

            for m1 in range(12):
                overlap_values[eof_id][0][m1][m1] = 1.0
                overlap_values[eof_id][1][m1][m1] = 1.0
                overlap_values[eof_id][2][m1][m1] = 1.0
                overlap_values[eof_id][3][m1][m1] = 1.0

                for m2 in range(m1 + 1, 12):
                    print('\t\t\t' + "Computing overlap for " + calendar.month_name[m1 + 1] + " vs. " + calendar.month_name[m2 + 1])

                    # compute temperature coefficient overlap
                    norm1 = quad(lambda x: T_kernels[m1].pdf(x)**2, T_lims[0], T_lims[1])[0]
                    norm2 = quad(lambda x: T_kernels[m2].pdf(x)**2, T_lims[0], T_lims[1])[0]
                    overlap = quad(lambda x: T_kernels[m1].pdf(x) * T_kernels[m2].pdf(x), T_lims[0], T_lims[1])[0] / np.sqrt(norm1 * norm2)
                    overlap_values[eof_id][0][m1][m2] = overlap
                    overlap_values[eof_id][0][m2][m1] = overlap

                    # compute zonal wind coefficient overlap
                    norm1 = quad(lambda x: u_kernels[m1].pdf(x)**2, u_lims[0], u_lims[1])[0]
                    norm2 = quad(lambda x: u_kernels[m2].pdf(x)**2, u_lims[0], u_lims[1])[0]
                    overlap = quad(lambda x: u_kernels[m1].pdf(x) * u_kernels[m2].pdf(x), u_lims[0], u_lims[1])[0] / np.sqrt(norm1 * norm2)
                    overlap_values[eof_id][1][m1][m2] = overlap
                    overlap_values[eof_id][1][m2][m1] = overlap

                    # compute meridional wind coefficient overlap
                    norm1 = quad(lambda x: v_kernels[m1].pdf(x)**2, v_lims[0], v_lims[1])[0]
                    norm2 = quad(lambda x: v_kernels[m2].pdf(x)**2, v_lims[0], v_lims[1])[0]
                    overlap = quad(lambda x: v_kernels[m1].pdf(x) * v_kernels[m2].pdf(x), v_lims[0], v_lims[1])[0] / np.sqrt(norm1 * norm2)
                    overlap_values[eof_id][2][m1][m2] = overlap
                    overlap_values[eof_id][2][m2][m1] = overlap

                    # compute density coefficient overlap
                    norm1 = quad(lambda x: d_kernels[m1].pdf(x)**2, d_lims[0], d_lims[1])[0]
                    norm2 = quad(lambda x: d_kernels[m2].pdf(x)**2, d_lims[0], d_lims[1])[0]
                    overlap = quad(lambda x: d_kernels[m1].pdf(x) * d_kernels[m2].pdf(x), d_lims[0], d_lims[1])[0] / np.sqrt(norm1 * norm2)
                    overlap_values[eof_id][3][m1][m2] = overlap
                    overlap_values[eof_id][3][m2][m1] = overlap

        overlap_values = np.nan_to_num(overlap_values)
        np.save(output_path + "-seasonality", overlap_values)

    singular_values = np.loadtxt(eofs_path + "-singular_values.dat")
    T_weights = singular_values[:, 1][:eof_cluster_cnt] / np.sum(singular_values[:, 1][:eof_cluster_cnt])
    u_weights = singular_values[:, 2][:eof_cluster_cnt] / np.sum(singular_values[:, 2][:eof_cluster_cnt])
    v_weights = singular_values[:, 3][:eof_cluster_cnt] / np.sum(singular_values[:, 3][:eof_cluster_cnt])
    d_weights = singular_values[:, 4][:eof_cluster_cnt] / np.sum(singular_values[:, 4][:eof_cluster_cnt])

    overlap_values[overlap_values > 1.0] = 1.0
    T_dist_mat = -np.log10(np.average(overlap_values[:eof_cluster_cnt, 0, :, :], axis=0, weights=T_weights))
    u_dist_mat = -np.log10(np.average(overlap_values[:eof_cluster_cnt, 1, :, :], axis=0, weights=u_weights))
    v_dist_mat = -np.log10(np.average(overlap_values[:eof_cluster_cnt, 2, :, :], axis=0, weights=v_weights))
    d_dist_mat = -np.log10(np.average(overlap_values[:eof_cluster_cnt, 3, :, :], axis=0, weights=d_weights))

    for n in range(12):
        T_dist_mat[n][n] = 0.0
        u_dist_mat[n][n] = 0.0
        v_dist_mat[n][n] = 0.0
        d_dist_mat[n][n] = 0.0

    T_results = hjl.cluster(T_dist_mat, cluser_thresh, linkage_method='weighted', show_result=True, file_id=output_path + "-T_overlap")
    u_results = hjl.cluster(u_dist_mat, cluser_thresh, linkage_method='weighted', show_result=True, file_id=output_path + "-u_overlap")
    v_results = hjl.cluster(v_dist_mat, cluser_thresh, linkage_method='weighted', show_result=True, file_id=output_path + "-v_overlap")
    d_results = hjl.cluster(d_dist_mat, cluser_thresh, linkage_method='weighted', show_result=True, file_id=output_path + "-d_overlap")


################################
#   Define methods to sample   #
#     the coefficient PDFs     #
################################
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


def sample_atmo(months, eof_path, coeff_path, output_path, eof_cnt=50, prof_cnt=50):
    print("-" * 50 + '\n' + "Generating atmosphere states from coefficient PDFs...")
    # load mean profile and eofs
    print('\t' + "Loading mean profile info and eofs...")
    means = np.loadtxt(eof_path + "-means.dat")
    T_eofs = np.loadtxt(eof_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eof_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eof_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eof_path + "-density.eofs")

    # load and stack coefficient values for all months of interest
    print('\t' + "Loading coefficient values and using KDE to build pdf's...")
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

    print('\t' + "Pre-computing cdf's for random sampling...")
    for eof_id in range(eof_cnt):
        print('\t' + "eof_id: " + str(eof_id) + "...")
        T_cdfs[eof_id] = build_cdf(T_kernels[eof_id].pdf, T_lims[eof_id])
        u_cdfs[eof_id] = build_cdf(u_kernels[eof_id].pdf, u_lims[eof_id])
        v_cdfs[eof_id] = build_cdf(v_kernels[eof_id].pdf, v_lims[eof_id])
        d_cdfs[eof_id] = build_cdf(d_kernels[eof_id].pdf, d_lims[eof_id])

    # Generate a random atmosphere sample
    print('\t' + "Generating profiles...")
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

    # Write the profiles to file
    for pn in range(prof_cnt):
        np.savetxt(output_path + "-" + "%02d" % pn + ".met", profs[pn])

    for pn in range(prof_cnt):
        np.savetxt(output_path + "-mean.met", np.average(profs, axis=0))

