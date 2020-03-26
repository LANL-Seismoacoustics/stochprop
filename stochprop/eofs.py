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

from scipy.cluster import hierarchy
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d
from scipy.optimize import bisect, minimize
from scipy.stats import gaussian_kde
from scipy.spatial.distance import squareform

from infrapy.association import hjl
from infrapy.utils import prog_bar

gam = 1.4
gamR = gam * 287.06


def profiles_qc(path, pattern="*.met", skiprows=0, plot_bad_profs=False):
    """
        Runs a quality control (QC) check on profiles in the path
        matching the pattern.  It can optionally plot the bad
        profiles.  If it finds any, it makes a new direcotry
        in the path location called "bad_profs" and moves those
        profiles into the directory for you to check

        Parameters
        ----------
        path : string
            Path to the profiles to be QC'd
        pattern : string
            Pattern defining the list of profiles in the path
        skiprows : int
            Number of header rows in the profiles
        plot_bad_profs : boolean
            Flag to plot bad profiles identified in the analysis
    """

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
                if not os.path.isdir(path + "bad_profs/"):
                    os.system("mkdir " + path + "bad_profs/")

                os.system("mv " + path + file + " " + path + "bad_profs/")
                print("Bad profile:", file)
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
    """
        Read in a list of atmosphere files from the path location
            matching a specified pattern for continued analysis.

        Parameters
        ----------
        path : string
            Path to the profiles to be loaded
        pattern : string
            Pattern defining the list of profiles in the path
        skiprows : int
            Number of header rows in the profiles
        ref_alts : 1darray
            Reference altitudes if comparison is needed

        Returns:
        A : 2darray
            Atmosphere array of size M x (5 * N) for M atmospheres
                where each atmosphere samples N altitudes
    """

    print('\t' + "Loading profiles from " + path + " with pattern: " + pattern)

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
    """
        Computes the singular value decomposition (SVD)
            of an atmosphere set read into an array by
            stochprop.eofs.build_atmo_matrix() and saves
            the basis functions (empirical orthogonal
            functions) and singular values to file

        Parameters
        ----------
        A : 2darray
            Suite of atmosphere specifications from build_atmo_matrix
        alts : 1darray
            Altitudes at which the atmosphere is sampled from build_atmo_matrix
        output_path : string
            Path to output the SVD results
        eof_cnt : int
            Number of basic functions to save
    """

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
    cp_vals =  np.sqrt((gam / 10.0) * p_vals / d_vals)

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
    cT_eofs = eofs[:, file_len * 0:file_len * 1] / np.sqrt(abs(alts[1] - alts[0]))
    cp_eofs = eofs[:, file_len * 1:file_len * 2] / np.sqrt(abs(alts[1] - alts[0]))
    u_eofs  = eofs[:, file_len * 2:file_len * 3] / np.sqrt(abs(alts[1] - alts[0]))
    v_eofs  = eofs[:, file_len * 3:file_len * 4] / np.sqrt(abs(alts[1] - alts[0]))

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


def compute_coeffs(A, alts, eofs_path, output_path, eof_cnt=100, pool=None):
    """
        Compute the EOF coefficients for a suite of atmospheres
            and store the coefficient values.

        Parameters
        ----------
        A : 2darray
            Suite of atmosphere specifications from build_atmo_matrix
        alts : 1darray
            Altitudes at which the atmosphere is sampled from build_atmo_matrix
        eofs_path : string
            Path to the .eof results from compute_svd
        output_path : string
            Path where output will be stored
        eof_cnt : int
            Number of EOFs to consider in computing coefficients
        pool : pathos.multiprocessing.ProcessingPool
            Multiprocessing pool for accelerating calculations            

        Returns:
        coeffs : 2darray
            Array containing coefficient values of size prof_cnt
                by eof_cnt.  Result is also written to file.
    """

    print('\t' + "Computing EOF coefficients for profiles...", '\t', end=' ')
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
        print("Warning!  Altitudes don't match dimension of atmosphere matrix.")
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
    
        cp_diff =  np.sqrt(gam / 10.0 * p_interp(eof_alts) / d_interp(eof_alts))
        cp_diff -= np.sqrt(gam / 10.0 * means[:, 5][eofs_mask] / means[:, 4][eofs_mask])
    
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
            coeffs[n] = pool.map(calc_coeff, list(range(eof_cnt)))
        else:
            for m in range(eof_cnt):
                coeffs[n][m] = calc_coeff(m)

        prog_bar.increment(int(np.floor((50.0 * (n + 1)) / A.shape[0]) - np.floor((50.0 * n) / A.shape[0])))
    prog_bar.close()

    # store coefficient values as numpy array file
    np.save(output_path + "-coeffs", coeffs)

    return coeffs


def compute_overlap(coeffs, eof_cnt=100):
    """
        Compute the overlap of EOF coefficient distributions

        Parameters
        ----------
        coeffs : list of 2darrays
            List of 2darrays containing coefficients to consider 
                overlap in PDF of values
        eof_cnt : int
            Number of EOFs to compute

        Returns:
        overlap : 3darray
            Array containing overlap values of size coeff_cnt
                by coeff_cnt by eof_cnt
    """

    print("-" * 50 + '\n' + "Computing coefficient PDF overlap...")
    overlap = np.empty((eof_cnt, 12, 12))
       
    for eof_id in range(eof_cnt):
        print('\t' + "Computing overlap values for EOF " + str(eof_id) + "...")

        for m1 in range(12):
            overlap[eof_id][m1][m1] = 1.0
            kernel1 = gaussian_kde(coeffs[m1][:, eof_id])
        
            for m2 in range(m1 + 1, 12):
                kernel2 = gaussian_kde(coeffs[m2][:, eof_id])

                # define coefficient value limits
                lims = np.array([min(min(coeffs[m1][:, eof_id]), min(coeffs[m2][:, eof_id])),
                                 max(max(coeffs[m1][:, eof_id]), max(coeffs[m2][:, eof_id]))])
                lims = np.array([(lims[1] + lims[0]) / 2.0 - 1.25 * (lims[1] - lims[0]), (lims[1] + lims[0]) / 2.0 + 1.25 * (lims[1] - lims[0])])
                c_vals = np.linspace(lims[0], lims[1], 1001)

                # compute coefficient overlap via integration
                norm1 = simps(kernel1.pdf(c_vals)**2, c_vals)
                norm2 = simps(kernel2.pdf(c_vals)**2, c_vals)
                overlap[eof_id][m1][m2] = simps(kernel1.pdf(c_vals) * kernel2.pdf(c_vals), c_vals) / np.sqrt(norm1 * norm2)
                overlap[eof_id][m2][m1] = overlap[eof_id][m1][m2]

    return overlap


def compute_seasonality(overlap_file, eofs_path, file_id=None):
    """
        Compute the overlap of EOF coefficients for

        Parameters
        ----------
        overlap_file : string
            Path and name of file containing results of stochprop.eofs.compute_overlap
        eofs_path : string
            Path to the .eof results from compute_svd
        file_id : string
            Path and ID to save the dendrogram result of the overlap analysis
    """

    overlap = np.load(overlap_file)
    eof_weights = np.loadtxt(eofs_path + "-singular_values.dat")[:, 1]

    dist_mat = -np.log(np.average(overlap, weights=eof_weights[:overlap.shape[0]], axis=0))
    np.fill_diagonal(dist_mat, 0.0)

    links = hierarchy.linkage(squareform(dist_mat), 'weighted')

    f, (ax1) = plt.subplots(1, 1)
    den = hierarchy.dendrogram(links, orientation="right", ax=ax1, labels=[calendar.month_abbr[n] for n in range(1, 13)])
    if file_id:
        plt.savefig(file_id + "-seasonality.png", dpi=300)
    plt.show(block=False)
    plt.pause(5.0)
    plt.close('all')


################################
#   Define methods to sample   #
#     the coefficient PDFs     #
################################
def define_coeff_limits(coeff_vals):
    """
        Compute upper and lower bounds for coefficient values

        Parameters
        ----------
        coeff_vals : 2darrays
            Coefficients computed with stochprop.eofs.compute_coeffs

        Returns:
        lims : 1darray
            Lower and upper bounds of coefficient value distribution
    """

    coeff_min, coeff_max = min(coeff_vals), max(coeff_vals)
    return np.array([(coeff_max + coeff_min) / 2.0 - 1.5 * (coeff_max - coeff_min),
                     (coeff_max + coeff_min) / 2.0 + 1.5 * (coeff_max - coeff_min)])

def build_cdf(pdf, lims, pnts=250):
    """
        Compute the cumulative distribution of a pdf within specified limits

        Parameters
        ----------
        pdf : function
            Probability distribution function (PDF) for a single variable
        lims : 1darray
            Iterable containing lower and upper bound for integration
        pnts : int
            Number of points to consider in defining the cumulative distribution

        Returns:
        cfd : interp1d
            Interpolated results for the cdf
    """

    norm = quad(pdf, lims[0], lims[1])[0]

    x_vals = np.linspace(lims[0], lims[1], pnts)
    cdf_vals = np.empty_like(x_vals)
    for n in range(pnts):
        cdf_vals[n] = quad(pdf, lims[0], x_vals[n])[0] / norm

    return interp1d(x_vals, cdf_vals)

def draw_from_pdf(pdf, lims, cdf=None, size=1):
    """
        Sample a number of values from a probability distribution 
            function (pdf) with specified limits

        Parameters
        ----------
        pdf : function
            Probability distribution function (PDF) for a single variable
        lims : 1darray
            Iterable containing lower and upper bound for integration
        cdf : function
            Cumulative distribution function (CDF) from stochprop.eofs.build_cfd
        size : int
            Number of samples to generate

        Returns:
        samples : 1darray
            Sampled values from the PDF
    """

    if not cdf:
        cdf = build_cdf(pdf, lims)

    def cdf_w_shift(x, d):
        return cdf(x) - d

    rnd_vals = np.random.random(size)
    samples = np.empty_like(rnd_vals)
    for n in range(size):
        samples[n] = bisect(cdf_w_shift, lims[0], lims[1], args=(rnd_vals[n],))

    return samples

def sample_atmo(coeffs, eofs_path, output_path, eof_cnt=100, prof_cnt=250, output_mean=False):
    """
        Generate atmosphere states using coefficient distributions for 
            a set of empirical orthogonal basis functions

        Parameters
        ----------
        coeffs : 2darrays
            Coefficients computed with stochprop.eofs.compute_coeffs
        eofs_path : string
            Path to the .eof results from compute_svd
        output_path : string
            Path where output will be stored
        eof_cnt : int
            Number of EOFs to use in building sampled specifications
        prof_cnt : int
            Number of atmospheric specification samples to generate
        output_mean : bool
            Flag to output the mean profile from the samples generated
    """
    
    print("-" * 50 + '\n' + "Generating atmosphere state samples from coefficient PDFs...")
    print('\t' + "Loading mean profile info and eofs...")
    means  =  np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs =  np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs =  np.loadtxt(eofs_path + "-merid_winds.eofs")

    # use kernel density estimates to define coefficient pdf's
    # and use the mean and span to define the limits for sampling
    print('\t' + "Mapping coefficidnts onto PDF's using KDE...")
    kernels, lims = [0] * eof_cnt, np.empty((eof_cnt, 2))
    for eof_id in range(eof_cnt):
        kernels[eof_id] = gaussian_kde(coeffs[:, eof_id])
        lims[eof_id] = define_coeff_limits(coeffs[:, eof_id])

    # Generate prof_cnt random atmosphere samples
    print('\t' + "Generating sample atmosphere profiles...")
    cT0 = np.sqrt(gamR * means[:, 1])
    cp0 = np.sqrt(gam / 10.0 * means[:, 5]  / means[:, 4])
    
    sampled_profs = np.array([np.copy(means)] * prof_cnt)
    sampled_cTs = np.array([cT0] * prof_cnt)
    sampled_cps = np.array([cp0] * prof_cnt)

    for eof_id in range(eof_cnt):
        print('\t\t' + "Sampling EOF coefficient for eof_id = " + str(eof_id) + "...")
        sampled_coeffs = draw_from_pdf(kernels[eof_id].pdf, lims[eof_id], size=prof_cnt)

        for pn in range(prof_cnt):
            sampled_cTs[pn] = sampled_cTs[pn] + sampled_coeffs[pn] * cT_eofs[:, eof_id + 1]
            sampled_cTs[pn] = sampled_cTs[pn] + sampled_coeffs[pn] * cp_eofs[:, eof_id + 1]

            sampled_profs[pn][:, 2] = sampled_profs[pn][:, 2] + sampled_coeffs[pn] * u_eofs[:, eof_id + 1]
            sampled_profs[pn][:, 3] = sampled_profs[pn][:, 3] + sampled_coeffs[pn] * v_eofs[:, eof_id + 1]

    sampled_profs[:, :, 1] = sampled_cTs**2 / gamR
    sampled_profs[:, :, 5] = sampled_profs[:, :, 4] * sampled_cps**2 / (gam / 10.0)

    # save the individual profiles and the mean profile
    print('\t' + "Writing sampled atmospheres to file...", '\n')
    for pn in range(prof_cnt):
        np.savetxt(output_path + "-" + "%02d" % pn + ".met", sampled_profs[pn])
    if output_mean:
        np.savetxt(output_path + "-mean.met", np.average(sampled_profs, axis=0))
        np.savetxt(output_path + "-stdev.met", np.std(sampled_profs, axis=0))


def maximum_likelihood_profile(coeffs, eofs_path, output_path, eof_cnt=100):
    """
        Use coefficient distributions for a set of empirical orthogonal
            basis functions to compute the maximum likelihood specification

        Parameters
        ----------
        coeffs : 2darrays
            Coefficients computed with stochprop.eofs.compute_coeffs
        eofs_path : string
            Path to the .eof results from compute_svd
        output_path : string
            Path where output will be stored
        eof_cnt : int
            Number of EOFs to use in building sampled specifications
    """

    print("-" * 50 + '\n' + "Generating maximum likelihood atmosphere states from coefficient PDFs...")
    # load mean profile and eofs
    print('\t' + "Loading mean profile info and eofs...")
    means  =  np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs =  np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs =  np.loadtxt(eofs_path + "-merid_winds.eofs")
    
    plt.show(block=False)

    # use kernel density estimates to define coefficient pdf's
    # and use the mean and span to define the limits for sampling
    print('\t' + "Determining maximum likelihood coefficient values...")

    # Generate the profile fit
    cT0 = np.sqrt(gamR * means[:, 1])
    cp0 = np.sqrt(gam / 10.0 * means[:, 5] / means[:, 4])
    
    ml_prof = np.copy(means)
    cT_ml = np.copy(cT0)
    cp_ml = np.copy(cp0)
    for n in range(eof_cnt):
        kernel = gaussian_kde(coeffs[:, n])
        lims = define_coeff_limits(coeffs[:, n])
        
        c_vals = np.linspace(lims[0], lims[1], 1001)
        coeff_ml = c_vals[np.argmax(kernel.pdf(c_vals))]
        
        cT_ml = cT_ml + coeff_ml * cT_eofs[:, n + 1]
        cp_ml = cp_ml + coeff_ml * cp_eofs[:, n + 1]
    
        ml_prof[:, 2] = ml_prof[:, 2] + coeff_ml * u_eofs[:, n + 1]
        ml_prof[:, 3] = ml_prof[:, 3] + coeff_ml * v_eofs[:, n + 1]
    
    ml_prof[:, 1] = cT_ml**2 / gamR
    ml_prof[:, 5] = ml_prof[:, 4] * cp_ml**2 / (gam / 10.0)
    
    print('\t' + "Writing maximum likelihood atmosphere to file...", '\n')
    np.savetxt(output_path+ "-maximum_likelihood.met", ml_prof)



################################
#   Define methods to fit or   #
#    perturb specific atmos    #
################################
def fit_atmo(prof_path, eofs_path, output_path, eof_cnt=100):
    """
        Compute a given number of EOF coefficients to fit a given
            atmophere specification using the basic functions.  Write
            the resulting approximated atmospheric specification to
            file.

        Parameters
        ----------
        prof_path : string
            Path and name of the specification to be fit
        eofs_path : string
            Path to the .eof results from compute_svd
        output_path : string
            Path where output will be stored
        eof_cnt : int
            Number of EOFs to use in building approximate specification
    """

    print("Generating EOF fit to " + prof_path + "...")

    # load means and eofs
    means = np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")
    
    # load profile and identify common altitudes
    # Note: means and eofs are assumed to have identical sampling (means[:, 0] = {..}_eofs[:, 0])
    #       but the profile being fit may not have matching sampling
    profile = np.loadtxt(prof_path)
    
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
    
    cp_diff =  np.sqrt(gam / 10.0 * p_interp(means[:, 0][eofs_mask]) / d_interp(means[:, 0][eofs_mask]))
    cp_diff -= np.sqrt(gam / 10.0 * means[:, 5][eofs_mask] / means[:, 4][eofs_mask])
    
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
    coeffs = np.empty(eof_cnt)
    for n in range(eof_cnt):
        coeffs[n] = calc_coeff(n)

    # Generate the profile fit
    cT0 = np.sqrt(gamR * means[:, 1])
    cp0 = np.sqrt(gam / 10.0 * means[:, 5] / means[:, 4])

    fit = np.copy(means)
    cT_fit = np.copy(cT0)
    cp_fit = np.copy(cp0)

    for n in range(eof_cnt):
        cT_fit = cT_fit + coeffs[n] * cT_eofs[:, n + 1]
        cp_fit = cp_fit + coeffs[n] * cp_eofs[:, n + 1]
    
        fit[:, 2] = fit[:, 2] + coeffs[n] * u_eofs[:, n + 1]
        fit[:, 3] = fit[:, 3] + coeffs[n] * v_eofs[:, n + 1]
    
    fit[:, 1] = cT_fit**2 / gamR
    fit[:, 5] = fit[:, 4] * cp_fit**2 / (gam / 10.0)

    np.savetxt(output_path, fit)

def perturb_atmo(prof_path, eofs_path, output_path, uncertainty=10.0, eof_max=100, eof_cnt=50, sample_cnt=1, alt_wt_pow=2.0, sing_val_wt_pow=0.25):
    """
        Use EOFs to perturb a specified profile using a given scale

        Parameters
        ----------
        prof_path : string
            Path and name of the specification to be fit
        eofs_path : string
            Path to the .eof results from compute_svd
        output_path : string
            Path where output will be stored
        uncertainty : float
            Estimate of uncertainty in wind speeds; 95% confidence is set to this value
        eof_max : int
            Higher numbered EOF to sample
        eof_cnt : int
            Number of EOFs to sample in the perturbation (can be less than eof_max)
        sample_cnt : int
            Number of perturbed atmospheric samples to generate
        alt_wt_pow : float
            Power raising relative mean altitude value in weighting
        sing_val_wt_pow : float
            Power raising relative singular value in weighting
    """

    print("Generating EOF perturbations to " + prof_path + "...")
    ref_atmo = np.loadtxt(prof_path)
    
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs  = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs  = np.loadtxt(eofs_path + "-merid_winds.eofs")
    sing_vals = np.loadtxt(eofs_path + "-singular_values.dat")
    
    # Define altitude limits
    alt_min = max(u_eofs[ 0, 0], ref_atmo[ 0, 0])
    alt_max = min(u_eofs[-1, 0], ref_atmo[-1, 0])
    
    eof_mask = np.logical_and(alt_min <= u_eofs[:, 0], u_eofs[:, 0] <= alt_max)
    ref_mask = np.logical_and(alt_min <= ref_atmo[:, 0], ref_atmo[:, 0] <= alt_max)
    z_vals = np.copy(ref_atmo[:, 0][ref_mask])

    # interpolate the eofs and use random coefficients to generate perturbations
    cT_perturb_all = np.empty((sample_cnt, len(u_eofs[:, 0][eof_mask])))
    cp_perturb_all = np.empty((sample_cnt, len(u_eofs[:, 0][eof_mask])))
    u_perturb_all = np.empty((sample_cnt, len(u_eofs[:, 0][eof_mask])))
    v_perturb_all = np.empty((sample_cnt, len(u_eofs[:, 0][eof_mask])))

    for m in range(sample_cnt):        
        wts = np.empty(eof_cnt)

        cT_perturb = np.zeros((eof_cnt, len(u_eofs[:, 0][eof_mask])))
        cp_perturb = np.zeros((eof_cnt, len(u_eofs[:, 0][eof_mask])))
        u_perturb = np.zeros((eof_cnt, len(u_eofs[:, 0][eof_mask])))
        v_perturb = np.zeros((eof_cnt, len(u_eofs[:, 0][eof_mask])))

        for j, n in enumerate(np.random.choice(range(eof_max), eof_cnt, replace=False)):
            coeff_val = np.random.randn();    
            
            cT_perturb[j] = coeff_val * interp1d(cT_eofs[:, 0][eof_mask], cT_eofs[:, n + 1][eof_mask])(z_vals)
            cp_perturb[j] = coeff_val * interp1d(cp_eofs[:, 0][eof_mask], cp_eofs[:, n + 1][eof_mask])(z_vals)
            u_perturb[j]  = coeff_val * interp1d(u_eofs[:, 0][eof_mask],  u_eofs[:, n + 1][eof_mask])(z_vals)
            v_perturb[j]  = coeff_val * interp1d(v_eofs[:, 0][eof_mask],  v_eofs[:, n + 1][eof_mask])(z_vals)

            wts[j]  = simps(u_eofs[:, 0][eof_mask] * abs(u_eofs[:, n + 1][eof_mask]), u_eofs[:, 0][eof_mask]) / simps(abs(u_eofs[:, n + 1][eof_mask]), u_eofs[:, 0][eof_mask])
            wts[j] += simps(u_eofs[:, 0][eof_mask] * abs(v_eofs[:, n + 1][eof_mask]), v_eofs[:, 0][eof_mask]) / simps(abs(v_eofs[:, n + 1][eof_mask]), u_eofs[:, 0][eof_mask])
            wts[j] /= (2.0 * max(u_eofs[:, 0][eof_mask]))
            wts[j] = wts[j]**alt_wt_pow

            wts[j] *= (sing_vals[n, 1] / sing_vals[0, 1])**sing_val_wt_pow

        cT_perturb_all[m] = np.average(cT_perturb, axis=0, weights=wts)
        cp_perturb_all[m] = np.average(cp_perturb, axis=0, weights=wts)
        u_perturb_all[m]  = np.average(u_perturb, axis=0, weights=wts)
        v_perturb_all[m]  = np.average(v_perturb, axis=0, weights=wts)

    wind_perturbation = np.sqrt(u_perturb_all**2 + v_perturb_all**2)
    mid_alt_mask = np.logical_and(30.0 <= u_eofs[:, 0][eof_mask], u_eofs[:, 0][eof_mask] <= 90.0)

    for m in range(sample_cnt):
        T_vals = np.copy(ref_atmo[:, 1][ref_mask])
        u_vals = np.copy(ref_atmo[:, 2][ref_mask])
        v_vals = np.copy(ref_atmo[:, 3][ref_mask])
        d_vals = np.copy(ref_atmo[:, 4][ref_mask])
        p_vals = np.copy(ref_atmo[:, 5][ref_mask])
        
        cT_vals = np.sqrt(gamR * T_vals)
        cp_vals = np.sqrt(gam / 10.0 * p_vals / d_vals)

        cT_vals = cT_vals + cT_perturb_all[m] * (uncertainty / (np.average(wind_perturbation[:, mid_alt_mask]) * 2.0))
        cp_vals = cp_vals + cp_perturb_all[m] * (uncertainty / (np.average(wind_perturbation[:, mid_alt_mask]) * 2.0))
        u_vals = u_vals + u_perturb_all[m] * (uncertainty / (np.average(wind_perturbation[:, mid_alt_mask]) * 2.0))
        v_vals = v_vals + v_perturb_all[m] * (uncertainty / (np.average(wind_perturbation[:, mid_alt_mask]) * 2.0))

        T_vals = cT_vals**2 / gamR
        p_vals = d_vals * cp_vals**2 / (gam / 10.0)
        
        np.savetxt(output_path + "-" + str(m) + ".met", np.vstack((z_vals, T_vals, u_vals, v_vals, d_vals, p_vals)).T)





