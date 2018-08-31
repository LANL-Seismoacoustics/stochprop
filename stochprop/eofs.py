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

def compute_svd(profs_path, output_path, eof_cnt=100, g2s_file_len=901, months=range(1, 13), years = range(2007, 2015)):
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
                        d_vals = np.vstack((d_vals, profile[:, 4]))
                        p_vals = np.vstack((p_vals, profile[:, 5]))

    print('\t' + "Building SVD expansion...")

    stacked_vals = np.hstack((T_vals, u_vals))
    stacked_vals = np.hstack((stacked_vals, v_vals))
    stacked_vals = np.hstack((stacked_vals, d_vals))
    stacked_vals = np.hstack((stacked_vals, p_vals))

    stacked_means = np.mean(stacked_vals, axis=0)
    T_mean = stacked_means[g2s_file_len * 0:g2s_file_len * 1]
    u_mean = stacked_means[g2s_file_len * 1:g2s_file_len * 2]
    v_mean = stacked_means[g2s_file_len * 2:g2s_file_len * 3]
    d_mean = stacked_means[g2s_file_len * 3:g2s_file_len * 4]
    p_mean = stacked_means[g2s_file_len * 4:g2s_file_len * 5]

    stacked_diffs = stacked_vals[:, :g2s_file_len * 3] - np.array([stacked_means] * stacked_vals.shape[0])[:, :g2s_file_len * 3]
    stacked_diffs = np.hstack((stacked_diffs, np.log(stacked_vals[:, g2s_file_len * 3:]) - np.log(np.array([stacked_means] * stacked_vals.shape[0])[:, g2s_file_len * 3:])))

    _, singular_vals, eofs = np.linalg.svd(stacked_diffs, full_matrices=False)

    eofs = eofs / np.sqrt(0.2)
    '''
    print('\t' + "Normalizing EOFs...")
    for n in range(len(eofs)):
        norm = 0.0
        norm = norm + quad(interp1d(profile[:, 0], eofs[n][g2s_file_len * 0:g2s_file_len * 1]**2), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        norm = norm + quad(interp1d(profile[:, 0], eofs[n][g2s_file_len * 1:g2s_file_len * 2]**2), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        norm = norm + quad(interp1d(profile[:, 0], eofs[n][g2s_file_len * 2:g2s_file_len * 3]**2), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        norm = norm + quad(interp1d(profile[:, 0], eofs[n][g2s_file_len * 3:g2s_file_len * 4]**2), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        norm = norm + quad(interp1d(profile[:, 0], eofs[n][g2s_file_len * 4:g2s_file_len * 5]**2), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        eofs[n] = eofs[n] / np.sqrt(norm)
    '''

    T_eofs = eofs[:, g2s_file_len * 0:g2s_file_len * 1]
    u_eofs = eofs[:, g2s_file_len * 1:g2s_file_len * 2]
    v_eofs = eofs[:, g2s_file_len * 2:g2s_file_len * 3]
    d_eofs = eofs[:, g2s_file_len * 3:g2s_file_len * 4]
    p_eofs = eofs[:, g2s_file_len * 4:g2s_file_len * 5]

    # Write the mean profiles and desired number of EOFs to file
    print('\t' + "Writing to file...")
    file_out = open(output_path + "-singular_values.dat", 'w')
    for n in range(len(singular_vals)):
        print(n, '\t', singular_vals[n], file=file_out)
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

def test_coeffs(prof_path, eofs_path, coeff_cnt):
    # load means and eofs
    means = np.loadtxt(eofs_path + "-means.dat")
    T_eofs = np.loadtxt(eofs_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eofs_path + "-density.eofs")
    p_eofs = np.loadtxt(eofs_path + "-pressure.eofs")

    # load profile to fit and compute differences from mean
    profile = np.loadtxt(prof_path)

    T_diff = profile[:, 1] - means[:, 1]
    u_diff = profile[:, 2] - means[:, 2]
    v_diff = profile[:, 3] - means[:, 3]
    d_diff = np.log(profile[:, 4]) - np.log(means[:, 4])
    p_diff = np.log(profile[:, 5]) - np.log(means[:, 5])

    T_coeffs = np.empty(coeff_cnt)
    u_coeffs = np.empty(coeff_cnt)
    v_coeffs = np.empty(coeff_cnt)
    d_coeffs = np.empty(coeff_cnt)
    p_coeffs = np.empty(coeff_cnt)

    # integrate to compute EOF coefficients
    for n in range(coeff_cnt):
        T_coeffs[n] = quad(interp1d(profile[:, 0], T_diff * T_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        u_coeffs[n] = quad(interp1d(profile[:, 0], u_diff * u_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        v_coeffs[n] = quad(interp1d(profile[:, 0], v_diff * v_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        d_coeffs[n] = quad(interp1d(profile[:, 0], d_diff * d_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]
        p_coeffs[n] = quad(interp1d(profile[:, 0], p_diff * p_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=len(profile[:, 0]), epsrel=1.0e-3)[0]

    coeffs = T_coeffs + u_coeffs + v_coeffs + d_coeffs + p_coeffs

    # Generate profile
    prof_fit = np.copy(means)
    prof_fit[:, 4:] = np.log(prof_fit[:, 4:])

    for n in range(coeff_cnt):
        prof_fit[:, 1] = prof_fit[:, 1] + coeffs[n] * T_eofs[:, n + 1]
        prof_fit[:, 2] = prof_fit[:, 2] + coeffs[n] * u_eofs[:, n + 1]
        prof_fit[:, 3] = prof_fit[:, 3] + coeffs[n] * v_eofs[:, n + 1]
        prof_fit[:, 4] = prof_fit[:, 4] + coeffs[n] * d_eofs[:, n + 1]
        prof_fit[:, 5] = prof_fit[:, 5] + coeffs[n] * p_eofs[:, n + 1]

    # convert density and pressure back to physical units
    prof_fit[:, 4] = np.exp(prof_fit[:, 4])
    prof_fit[:, 5] = np.exp(prof_fit[:, 5])

    return prof_fit


def compute_coeffs(profs_path, eofs_path, output_path, months=range(1, 13), years = range(2007, 2015)):
    # load means and eofs
    print("-" * 50 + '\n' + "Loading eofs...")
    means = np.loadtxt(eofs_path + "-means.dat")
    T_eofs = np.loadtxt(eofs_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eofs_path + "-density.eofs")
    p_eofs = np.loadtxt(eofs_path + "-pressure.eofs")

    # load profile to be fit
    for month_index in months:
        month = "%02d" % month_index
        for hour_index in range(4):
            hour = "%02d" % (hour_index * 6)

            coeffs = np.empty((0, T_eofs.shape[1] - 1))
            coeff_cnt = T_eofs.shape[1] - 1

            for year_index in years:
                year = "%04d" % year_index
                for day_index in range(32):
                    day = "%02d" % (day_index + 1)

                    datetime_id = year + month + day + hour
                    filename = profs_path + year + "/" + datetime_id + ".met"

                    if os.path.isfile(filename):
                        print('\t' + "Computing EOF coefficients for " + filename)
                        profile = np.loadtxt(filename)

                        # apply the shifts from means
                        # note that density and pressure are log-scaled first
                        T_diff = profile[:, 1] - means[:, 1]
                        u_diff = profile[:, 2] - means[:, 2]
                        v_diff = profile[:, 3] - means[:, 3]
                        d_diff = np.log(profile[:, 4]) - np.log(means[:, 4])
                        p_diff = np.log(profile[:, 5]) - np.log(means[:, 5])

                        # interpolate and compute the projection onto each of the eofs
                        # stacking the entire atmosphere specification requires summing the projections onto each field (c_n^{(T(} + c_n^{(u)} + ...)
                        T_coeffs = np.empty(coeff_cnt)
                        u_coeffs = np.empty(coeff_cnt)
                        v_coeffs = np.empty(coeff_cnt)
                        d_coeffs = np.empty(coeff_cnt)
                        p_coeffs = np.empty(coeff_cnt)

                        for n in range(coeff_cnt):
                            T_coeffs[n] = quad(interp1d(profile[:, 0], T_diff * T_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            u_coeffs[n] = quad(interp1d(profile[:, 0], u_diff * u_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            v_coeffs[n] = quad(interp1d(profile[:, 0], v_diff * v_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            d_coeffs[n] = quad(interp1d(profile[:, 0], d_diff * d_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]
                            p_coeffs[n] = quad(interp1d(profile[:, 0], p_diff * p_eofs[:, n + 1]), profile[0][0], profile[-1][0], limit=int(len(profile[:, 0]) / 2), epsrel=1.0e-3)[0]

                        coeffs = np.vstack((coeffs, T_coeffs + u_coeffs + v_coeffs + d_coeffs + p_coeffs))

            # store coefficient values as numpy array file
            np.save(output_path + "/" + hour + "00/" + month + "-coeffs", coeffs)

def cluster_seasonality(coeff_path, eofs_path, output_path, cluser_thresh=0.15, eof_cnt=50, eof_cluster_cnt=15, overlap_file=None):
    print("-" * 50 + '\n' + "Computing coefficient PDF clustering...")
    if overlap_file:
        print('\t' + "Loading overlap values from file...")
        overlap_values = np.load(overlap_file)
    else:
        print('\t' + "Computing overlap values from coefficient distributions...")
        overlap_values = np.empty((eof_cnt, 12, 12))
        for eof_id in range(eof_cnt):
            print('\t\t' + "Analyzing seasonality of EOF id: " + str(eof_id))

            # define coefficient pdfs for each month
            kernels = [0] * 12
            lims = np.array([np.inf, -np.inf])

            for month_index in range(1, 13):
                month = "%02d" % month_index

                coeffs = np.load(coeff_path + month + "-coeffs.npy")
                kernels[int(month) - 1] = gaussian_kde(coeffs[:, eof_id])
                lims = np.array([min(lims[0], np.min(coeffs[:, eof_id])), max(lims[1], np.max(coeffs[:, eof_id]))])

            lims = np.array([(lims[1] + lims[0]) / 2.0 - 1.5 * (lims[1] - lims[0]), (lims[1] + lims[0]) / 2.0 + 1.5 * (lims[1] - lims[0])])

            for m1 in range(12):
                overlap_values[eof_id][m1][m1] = 1.0

                for m2 in range(m1 + 1, 12):
                    print('\t\t\t' + "Computing overlap for " + calendar.month_name[m1 + 1] + " vs. " + calendar.month_name[m2 + 1])

                    # compute coefficient overlap
                    norm1 = quad(lambda x: kernels[m1].pdf(x)**2, lims[0], lims[1])[0]
                    norm2 = quad(lambda x: kernels[m2].pdf(x)**2, lims[0], lims[1])[0]
                    overlap = quad(lambda x: kernels[m1].pdf(x) * kernels[m2].pdf(x), lims[0], lims[1])[0] / np.sqrt(norm1 * norm2)

                    overlap_values[eof_id][m1][m2] = overlap
                    overlap_values[eof_id][m2][m1] = overlap

        overlap_values = np.nan_to_num(overlap_values)
        np.save(output_path + "-seasonality", overlap_values)

    singular_values = np.loadtxt(eofs_path + "-singular_values.dat")
    weights = singular_values[:eof_cluster_cnt] / np.sum(singular_values[:eof_cluster_cnt])

    overlap_values[overlap_values > 1.0] = 1.0
    dist_mat = -np.log10(np.average(overlap_values[:eof_cluster_cnt, :, :], axis=0, weights=weights))

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

def sample_atmo(months, hour, eof_path, coeff_path, output_path, eof_cnt=50, prof_cnt=50):
    print("-" * 50 + '\n' + "Generating atmosphere states from coefficient PDFs...")
    # load mean profile and eofs
    print('\t' + "Loading mean profile info and eofs...")
    means = np.loadtxt(eof_path + "-means.dat")
    T_eofs = np.loadtxt(eof_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eof_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eof_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eof_path + "-density.eofs")
    p_eofs = np.loadtxt(eof_path + "-pressure.eofs")

    # load and stack coefficient values for all months of interest
    print('\t' + "Loading coefficient values and using KDE to build pdf's...")
    coeffs = np.empty((0, T_eofs.shape[1] - 1))

    for month in months:
        coeffs = np.vstack((coeffs, np.load(coeff_path + "/" + hour + "00/" + month + "-coeffs.npy")))

    # use kernel density estimates to define coefficient pdf's
    # and use the mean and span to define the limits for sampling
    kernels, lims = [0] * eof_cnt, np.empty((eof_cnt, 2))

    for eof_id in range(eof_cnt):
        kernels[eof_id] = gaussian_kde(coeffs[:, eof_id])
        lims[eof_id] = define_limits(coeffs[:, eof_id])

    # compute the cumulative density functions for all pdf's
    cdfs = [0] * eof_cnt

    print('\t' + "Pre-computing cdf's for random sampling...")
    for eof_id in range(eof_cnt):
        print('\t' + "eof_id: " + str(eof_id) + "...")
        cdfs[eof_id] = build_cdf(kernels[eof_id].pdf, lims[eof_id])

    # Generate a random atmosphere sample
    print('\t' + "Generating profiles...")
    means[:, 4:] = np.log(means[:, 4:])
    profs = np.array([means] * prof_cnt)

    for eof_id in range(eof_cnt):
        coeff_vals = draw_from_pdf(kernels[eof_id].pdf, lims[eof_id], cdf=cdfs[eof_id], size=prof_cnt)

        for pn in range(prof_cnt):
            profs[pn][:, 1] += coeff_vals[pn] * T_eofs[:, eof_id + 1]
            profs[pn][:, 2] += coeff_vals[pn] * u_eofs[:, eof_id + 1]
            profs[pn][:, 3] += coeff_vals[pn] * v_eofs[:, eof_id + 1]
            profs[pn][:, 4] += coeff_vals[pn] * d_eofs[:, eof_id + 1]
            profs[pn][:, 5] += coeff_vals[pn] * p_eofs[:, eof_id + 1]

    # convert density and pressure back to physical units
    for pn in range(prof_cnt):
        profs[pn][:, 4] = np.exp(profs[pn][:, 4])
        profs[pn][:, 5] = np.exp(profs[pn][:, 5])

    # Write the profiles to file
    for pn in range(prof_cnt):
        np.savetxt(output_path + "-" + "%02d" % pn + ".met", profs[pn])

    for pn in range(prof_cnt):
        np.savetxt(output_path + "-mean.met", np.average(profs, axis=0))

