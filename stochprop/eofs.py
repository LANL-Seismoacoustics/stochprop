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
import fnmatch
import warnings

import matplotlib.pyplot as plt

from scipy.integrate import quad, simps
from scipy.interpolate import interp1d
from scipy.optimize import bisect, minimize
from scipy.stats import gaussian_kde

from infrapy.association import hjl
from infrapy.utils import prog_bar

gamma = 1.4
gamR = gamma * 287.06

def profile_qc(g2s_file):
    profile = np.loadtxt(g2s_file)
    dTdz = np.gradient(profile[:,1]) / np.gradient(profile[:,0])
    if max(abs(dTdz)) > 50.0:
        return False
    else:
        return True

def profs_check(path, pattern="*.met", skiprows=0, plot_bad_profs=False):
    print("Running QC on profiles in " + path + " matching pattern " + pattern + "...")
    file_list = []
    dir_files = os.listdir(path)
    for file in dir_files:
        if fnmatch.fnmatch(file, pattern):
            file_list += [file]

    if len(file_list) > 0:
        for file in file_list:
            atmo = np.loadtxt(path + file, skiprows=skiprows)
            dTdz = np.gradient(atmo[:,1]) / np.gradient(atmo[:,0])
            dudz = np.gradient(atmo[:,2]) / np.gradient(atmo[:,0])
            dvdz = np.gradient(atmo[:,3]) / np.gradient(atmo[:,0])

            if max(abs(dTdz)) > 50.0 or max(abs(dudz)) > 50.0 or max(abs(dvdz)) > 50.0:
                os.system("mv " + path + file + " " + path + "bad_profs")
                if plot_bad_profs:
                    plt.figure(1)
                    plt.clf()
                    plt.plot(profile[:, 1], profile[:, 0])
                    plt.plot(dTdz, profile[:, 0])
                    plt.title(filename)
                    plt.pause(0.1)
    else:
        print("WARNING!! No profiles matching specified pattern in path.")

def build_atmo_matrix(path, pattern="*.met", skiprows=0, ref_alts=None):
    print("Loading profiles from " + path + " with pattern: " + pattern)

    file_list = []
    dir_files = os.listdir(path)
    for file in dir_files:
        if fnmatch.fnmatch(file, pattern):
            file_list += [file]

    if len(file_list) > 0:
        atmo = np.loadtxt(path + file_list[0], skiprows=skiprows)
        if np.any(ref_alts) == None:
            z0 = atmo[:, 0]
        else:
            z0 = ref_alts

        if np.allclose(z0, atmo[:, 0]):
            T = atmo[:, 1]
            u = atmo[:, 2]
            v = atmo[:, 3]
            d = atmo[:, 4]
            p = atmo[:, 5]
        else:
            print("WARNING!!  Altitudes in " + path + file_list[0] + " don't match expected values.  Interpolating to resolve...")
            T = interp1d(atmo[:, 0], atmo[:, 1])(z0)
            u = interp1d(atmo[:, 0], atmo[:, 2])(z0)
            v = interp1d(atmo[:, 0], atmo[:, 3])(z0)
            d = interp1d(atmo[:, 0], atmo[:, 4])(z0)
            p = interp1d(atmo[:, 0], atmo[:, 5])(z0)

    
        for file in file_list[1:]:
            atmo = np.loadtxt(path + file, skiprows=skiprows)
                
            if np.allclose(z0, atmo[:, 0]):
                T = np.vstack((T, atmo[:, 1]))
                u = np.vstack((u, atmo[:, 2]))
                v = np.vstack((v, atmo[:, 3]))
                d = np.vstack((d, atmo[:, 4]))
                p = np.vstack((p, atmo[:, 5]))
            else:
                print("WARNING!!  Altitudes in " + path + file + " don't match expected values.  Interpolating to resolve...")
                T = np.vstack((T, interp1d(atmo[:, 0], atmo[:, 1])(z0)))
                u = np.vstack((u, interp1d(atmo[:, 0], atmo[:, 2])(z0)))
                v = np.vstack((v, interp1d(atmo[:, 0], atmo[:, 3])(z0)))
                d = np.vstack((d, interp1d(atmo[:, 0], atmo[:, 4])(z0)))
                p = np.vstack((p, interp1d(atmo[:, 0], atmo[:, 5])(z0)))
        

        A = np.hstack((T, u))
        A = np.hstack((A, v))
        A = np.hstack((A, d))
        A = np.hstack((A, p))

        return A, z0
    else:
        return None, None

def compute_svd(A, alts, output_path, eof_cnt=100):
    file_len = int(A.shape[1] / 5)
    alts_len = len(alts)
    if file_len != alts_len:
        print("Warning!  Altitudes don't match dimension of A matrix.")
        return 0

    if eof_cnt > A.shape[0]:
        eof_cnt = A.shape[0]
        msg = "Warning!  Atmosphere count less than eof_cnt.  Only returning " + str(eof_cnt) + " EOFs."
        raise ValueError(msg)
      
    # compute the differences to fit with EOFs
    T_vals = A[:, file_len * 0:file_len * 1]
    u_vals = A[:, file_len * 1:file_len * 2]
    v_vals = A[:, file_len * 2:file_len * 3]
    d_vals = A[:, file_len * 3:file_len * 4]
    p_vals = A[:, file_len * 4:file_len * 5]

    T_mean = np.mean(T_vals, axis=0)
    u_mean = np.mean(u_vals, axis=0)
    v_mean = np.mean(v_vals, axis=0)
    d_mean = np.mean(d_vals, axis=0)
    p_mean = np.mean(p_vals, axis=0)

    cT_vals =  np.sqrt(gamR * T_vals)
    cp_vals =  np.sqrt(gamma / 10.0 * p_vals / d_vals)

    u_diff  = u_vals - np.array([u_mean] * u_vals.shape[0])
    v_diff  = v_vals - np.array([v_mean] * v_vals.shape[0])
    cT_diff = cT_vals - np.array([np.mean(cT_vals, axis=0)] * cT_vals.shape[0])
    cp_diff = cp_vals - np.array([np.mean(cp_vals, axis=0)] * cp_vals.shape[0])
    
    # stack diffs and evaluate SVD
    diffs = np.hstack((cT_diff, cp_diff))
    diffs = np.hstack((diffs, u_diff))
    diffs = np.hstack((diffs, v_diff))

    _, singular_vals, eofs = np.linalg.svd(diffs, full_matrices=False)

    # normalize the eofs by the altitude resolution and separate into cT, cp, u, v
    eofs = eofs / np.sqrt(abs(alts[1] - alts[0]))
    
    cT_eofs = eofs[:, file_len * 0:file_len * 1]
    cp_eofs = eofs[:, file_len * 1:file_len * 2]
    u_eofs  = eofs[:, file_len * 2:file_len * 3]
    v_eofs  = eofs[:, file_len * 3:file_len * 4]

    # Write the mean profiles and desired number of EOFs to file
    file_out = open(output_path + "-singular_values.dat", 'w')
    for n in range(len(singular_vals)):
        print(n, '\t', singular_vals[n], file=file_out)
    file_out.close()

    file_out = open(output_path + "-mean_atmo.dat", 'w')
    for n in range(len(alts)):
        print(alts[n], '\t', T_mean[n], '\t', u_mean[n], '\t', v_mean[n], '\t', d_mean[n], '\t', p_mean[n], file=file_out)
    file_out.close()

    file_out = open(output_path + "-ideal_gas_snd_spd.eofs", 'w')
    for n in range(len(alts)):
        print(alts[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', cT_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()


    file_out = open(output_path + "-adiabatic_snd_spd.eofs", 'w')
    for n in range(len(alts)):
        print(alts[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', cp_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()

    file_out = open(output_path + "-zonal_winds.eofs", 'w')
    for n in range(len(alts)):
        print(alts[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', u_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()

    file_out = open(output_path + "-merid_winds.eofs", 'w')
    for n in range(len(alts)):
        print(alts[n], '\t', end=' ', file=file_out)
        for m in range(eof_cnt):
            print('\t', v_eofs[m][n], end=' ', file=file_out)
        print(' ', file=file_out)
    file_out.close()


def fit_profile(prof_path, eofs_path, eof_cnt, skiprows=0, pool=None):
    # load means and eofs
    means = np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")

    # load profile and identify common altitudes
    # Note: means and eofs are assumed to have identical sampling (means[:, 0] = {..}_eofs[:, 0])
    #       but the profile being fit may not have matching sampling
    profile = np.loadtxt(prof_path, skiprows=skiprows)
    
    alt_min = max(means[0, 0], profile[0, 0])
    alt_max = min(means[-1, 0], profile[-1, 0])

    eofs_mask = np.logical_and(alt_min <= means[:, 0], means[:, 0] <= alt_max)
    prof_mask = np.logical_and(alt_min <= profile[:, 0], profile[:, 0] <= alt_max)
    
    # interpolate the profile and evaluate differences at the sampled points of the mean/eofs
    T_interp = interp1d(profile[:, 0][prof_mask], profile[:, 1][prof_mask], kind='cubic')
    u_interp = interp1d(profile[:, 0][prof_mask], profile[:, 2][prof_mask], kind='cubic')
    v_interp = interp1d(profile[:, 0][prof_mask], profile[:, 3][prof_mask], kind='cubic')
    d_interp = interp1d(profile[:, 0][prof_mask], profile[:, 4][prof_mask], kind='cubic')
    p_interp = interp1d(profile[:, 0][prof_mask], profile[:, 5][prof_mask], kind='cubic')

    cT_diff =  np.sqrt(gamR * T_interp(means[:, 0][eofs_mask]))
    cT_diff -= np.sqrt(gamR * means[:, 1][eofs_mask])

    cp_diff =  np.sqrt(gamma / 10.0 * p_interp(means[:, 0][eofs_mask]) / d_interp(means[:, 0][eofs_mask]))
    cp_diff -= np.sqrt(gamma / 10.0 * means[:, 5][eofs_mask] / means[:, 4][eofs_mask])
    
    u_diff = u_interp(means[:, 0][eofs_mask]) - means[:, 2][eofs_mask]
    v_diff = v_interp(means[:, 0][eofs_mask]) - means[:, 3][eofs_mask]
   
    # define the integration to compute the EOF coefficients
    def calc_coeff(n):
        result  = simps(cT_diff * cT_eofs[:, n + 1][eofs_mask], cT_eofs[:, 0][eofs_mask])
        result += simps(cp_diff * cp_eofs[:, n + 1][eofs_mask], cp_eofs[:, 0][eofs_mask])
        result += simps(u_diff * u_eofs[:, n + 1][eofs_mask], u_eofs[:, 0][eofs_mask])
        result += simps(v_diff * v_eofs[:, n + 1][eofs_mask], v_eofs[:, 0][eofs_mask])

        return result

    # run the integration to define coefficients
    if pool:
        args = list(range(eof_cnt))
        coeffs = pool.map(calc_coeff, args)
    else:
        coeffs = np.empty(eof_cnt)
        for n in range(eof_cnt):
            coeffs[n] = calc_coeff(n)

    # Generate the profile fit
    cT_fit = np.sqrt(gamR * means[:, 1])
    cp_fit = np.sqrt(gamma / 10.0 * means[:, 5] / means[:, 4])
    for n in range(eof_cnt):
        cT_fit = cT_fit + coeffs[n] * cT_eofs[:, n + 1]
        cp_fit = cp_fit + coeffs[n] * cp_eofs[:, n + 1]

    fit = np.copy(means)
    for n in range(eof_cnt):
        fit[:, 2] = fit[:, 2] + coeffs[n] * u_eofs[:, n + 1]
        fit[:, 3] = fit[:, 3] + coeffs[n] * v_eofs[:, n + 1]
        fit[:, 4] = fit[:, 4] + 2.0 / (gamma - 1.0) * (means[:, 4] / cp_fit**2) * coeffs[n] * cp_eofs[:, n + 1]

    fit[:, 1] = cT_fit**2 / gamR
    fit[:, 5] = fit[:, 4] * cp_fit**2 / (gamma / 10.0)

    return profile, fit


def compute_coeffs(A, alts, eofs_path, output_path, eof_cnt=100, pool=None):
    print("Computing EOF coefficients for set of profiles...", '\t', end=' ')
    # load means and eofs
    means = np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs =  np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs =  np.loadtxt(eofs_path + "-merid_winds.eofs")
    
    # check A and specified alts are consistent
    file_len = int(A.shape[1] / 5)
    alts_len = len(alts)
    if file_len != alts_len:
        print("Warning!  Altitudes don't match dimension of A matrix.")
        return 0

    alt_min = max(means[0, 0], alts[0])
    alt_max = min(means[-1, 0], alts[-1])

    eofs_mask = np.logical_and(alt_min <= means[:, 0], means[:, 0] <= alt_max)
    A_mask = np.logical_and(alt_min <= alts, alts <= alt_max)

    A_alts = alts[A_mask]
    eof_alts = means[:, 0][eofs_mask]
                        
    # cycle through profiles in A to define coefficients
    coeffs = np.empty((A.shape[0], eof_cnt))
    prog_bar.prep(50)
    for n, An in enumerate(A):
        # interpolate the profile
        T_interp = interp1d(A_alts, An[file_len * 0:file_len * 1][A_mask], kind='cubic')
        u_interp = interp1d(A_alts, An[file_len * 1:file_len * 2][A_mask], kind='cubic')
        v_interp = interp1d(A_alts, An[file_len * 2:file_len * 3][A_mask], kind='cubic')
        d_interp = interp1d(A_alts, An[file_len * 3:file_len * 4][A_mask], kind='cubic')
        p_interp = interp1d(A_alts, An[file_len * 4:file_len * 5][A_mask], kind='cubic')
    
        # evaluate differences at the sampled points of the mean/eofs
        cT_diff =  np.sqrt(gamR * T_interp(eof_alts))
        cT_diff -= np.sqrt(gamR * means[:, 1][eofs_mask])
    
        cp_diff =  np.sqrt(gamma / 10.0 * p_interp(eof_alts) / d_interp(eof_alts))
        cp_diff -= np.sqrt(gamma / 10.0 * means[:, 5][eofs_mask] / means[:, 4][eofs_mask])
    
        u_diff = u_interp(eof_alts) - means[:, 2][eofs_mask]
        v_diff = v_interp(eof_alts) - means[:, 3][eofs_mask]

        # define the integration to compute the EOF coefficients
        def calc_coeff(n):
            result  = simps(cT_diff * cT_eofs[:, n + 1][eofs_mask], eof_alts)
            result += simps(cp_diff * cp_eofs[:, n + 1][eofs_mask], eof_alts)
            result += simps(u_diff * u_eofs[:, n + 1][eofs_mask], eof_alts)
            result += simps(v_diff * v_eofs[:, n + 1][eofs_mask], eof_alts)
    
            return result
    
        # run the integration to define coefficients
        if pool:
            args = list(range(eof_cnt))
            coeffs[n] = pool.map(calc_coeff, args)
        else:
            for m in range(eof_cnt):
                coeffs[n][m] = calc_coeff(m)

        prog_bar.increment(int(np.floor((50.0 * (n + 1)) / A.shape[0]) - np.floor((50.0 * n) / A.shape[0])))
    prog_bar.close()

    # store coefficient values as numpy array file
    np.save(output_path + "-coeffs", coeffs)

    return coeffs


def cluster_seasonality(coeffs, eofs_path, output_path, cluser_thresh=0.15, eof_cnt=100, eof_cluster_cnt=15, overlap_file=None):
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

def build_cdf(pdf, lims, pnts=250):
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

def sample_atmo(coeffs, eof_path, output_path, eof_cnt=100, prof_cnt=50):
    print("-" * 50 + '\n' + "Generating atmosphere states from coefficient PDFs...")
    # load mean profile and eofs
    print('\t' + "Loading mean profile info and eofs...")
    means = np.loadtxt(eof_path + "-means.dat")
    T_eofs = np.loadtxt(eof_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eof_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eof_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eof_path + "-density.eofs")
    p_eofs = np.loadtxt(eof_path + "-pressure.eofs")

    # use kernel density estimates to define coefficient pdf's
    # and use the mean and span to define the limits for sampling
    print('\t' + "Mapping coefficidnts onto PDF's using KDE...")
    kernels, lims = [0] * eof_cnt, np.empty((eof_cnt, 2))
    for eof_id in range(eof_cnt):
        kernels[eof_id] = gaussian_kde(coeffs[:, eof_id])
        lims[eof_id] = define_limits(coeffs[:, eof_id])

    # Generate prof_cnt random atmosphere samples
    print('\t' + "Generating sample atmosphere profiles...")
    sampled_profs = np.array([np.copy(means)] * prof_cnt)
    for pn in range(prof_cnt):
        sampled_profs[pn][:, 4] = np.log(sampled_profs[pn][:, 4])
        sampled_profs[pn][:, 5] = np.log(sampled_profs[pn][:, 5])

    for eof_id in range(eof_cnt):
        print('\t\t' + "Sampling EOF coefficient for eof_id = " + str(eof_id) + "...")
        sampled_coeffs = draw_from_pdf(kernels[eof_id].pdf, lims[eof_id], size=prof_cnt)

        for pn in range(prof_cnt):
            sampled_profs[pn][:, 1] = sampled_profs[pn][:, 1] + sampled_coeffs[pn] * T_eofs[:, eof_id + 1]
            sampled_profs[pn][:, 2] = sampled_profs[pn][:, 2] + sampled_coeffs[pn] * u_eofs[:, eof_id + 1]
            sampled_profs[pn][:, 3] = sampled_profs[pn][:, 3] + sampled_coeffs[pn] * v_eofs[:, eof_id + 1]
            sampled_profs[pn][:, 4] = sampled_profs[pn][:, 4] + sampled_coeffs[pn] * d_eofs[:, eof_id + 1]
            sampled_profs[pn][:, 5] = sampled_profs[pn][:, 5] + sampled_coeffs[pn] * p_eofs[:, eof_id + 1]

    # convert density and pressure back to physical units and save
    for pn in range(prof_cnt):
        sampled_profs[pn][:, 4] = np.exp(sampled_profs[pn][:, 4])
        sampled_profs[pn][:, 5] = np.exp(sampled_profs[pn][:, 5])
        np.savetxt(output_path + "-" + "%02d" % pn + ".met", sampled_profs[pn])

    np.savetxt(output_path + "-mean.met", np.average(sampled_profs, axis=0))

def maximum_likelihood_profile(coeffs, eof_path, output_path, eof_cnt=100):
    print("-" * 50 + '\n' + "Generating atmosphere states from coefficient PDFs...")
    # load mean profile and eofs
    print('\t' + "Loading mean profile info and eofs...")
    means = np.loadtxt(eof_path + "-means.dat")
    T_eofs = np.loadtxt(eof_path + "-temperature.eofs")
    u_eofs = np.loadtxt(eof_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eof_path + "-merid_winds.eofs")
    d_eofs = np.loadtxt(eof_path + "-density.eofs")
    p_eofs = np.loadtxt(eof_path + "-pressure.eofs")
    
    plt.show(block=False)

    # use kernel density estimates to define coefficient pdf's
    # and use the mean and span to define the limits for sampling
    print('\t' + "Determining maximum likelihood values of coefficidnts...")

    ml_prof = np.copy(means)
    ml_prof[:, 4] = np.log(ml_prof[:, 4])
    ml_prof[:, 5] = np.log(ml_prof[:, 5])
    
    for eof_id in range(eof_cnt):
        kernel = gaussian_kde(coeffs[:, eof_id])
        lims = define_limits(coeffs[:, eof_id])
    
        c_vals = np.linspace(lims[0], lims[1], 101)
        c_ml_est = c_vals[np.argmax(kernel.pdf(c_vals))]
        
        def temp(c):
            return -kernel.pdf(c)
        
        c_ml = minimize(temp, c_ml_est).x[0]
        print("EOF-" + str(eof_id) + '\t' + "ML Value: " + str(np.round(c_ml, 2)) + '\t' + "Range of values: [" + str(np.round(np.min(coeffs[:, eof_id]), 2)) + ", " + str(np.round(np.max(coeffs[:, eof_id]), 2)) + "]", end=' ')
        print('\t' + "Range of values: [" + str(np.round(np.min(coeffs[:, eof_id]) - c_ml, 2)) + ", " + str(np.round(np.max(coeffs[:, eof_id]) - c_ml, 2)) + "]")
        
        ml_prof[:, 1] = ml_prof[:, 1] + c_ml * T_eofs[:, eof_id + 1]
        ml_prof[:, 2] = ml_prof[:, 2] + c_ml * u_eofs[:, eof_id + 1]
        ml_prof[:, 3] = ml_prof[:, 3] + c_ml * v_eofs[:, eof_id + 1]
        ml_prof[:, 4] = ml_prof[:, 4] + c_ml * d_eofs[:, eof_id + 1]
        ml_prof[:, 5] = ml_prof[:, 5] + c_ml * p_eofs[:, eof_id + 1]

    ml_prof[:, 4] = np.exp(ml_prof[:, 4])
    ml_prof[:, 5] = np.exp(ml_prof[:, 5])

    print("Writing maximum likelihood atmosphere to file...")
    np.savetxt(output_path, ml_prof)

