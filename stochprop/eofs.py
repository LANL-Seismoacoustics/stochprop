# eofs.py
#
# Methods to construct a singular value decomposition (SVD) of a set of
# G2S profiles and define empirical orthogonal functions (EOFs).  Also 
# included are methods to compute coefficient values for the profiles and
# estimate statistics of those coefficients to generate atmospheric 
# samples for different seasons.  Lastly, tools are included to perturb
# an initial atmospheric state given some specified level of uncertainty.
#
# Philip Blom (pblom@lanl.gov)

import os
import calendar
import fnmatch
import datetime
import imp
import subprocess
import pkg_resources
import re 

from netCDF4 import Dataset

import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster import hierarchy
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import bisect
from scipy.stats import gaussian_kde
from scipy.spatial.distance import squareform

gam = 1.4
gamR = gam * 287.06


gasR = 287.0
den0 = 0.001225
coeffs_A = np.array([-3.9082017e-2, -1.1526465e-3, 3.2891937e-5, -2.0494958e-7,
                        -4.7087295e-2, 1.2506387e-3, -1.5194498e-5, 6.581877e-8])
coeffs_B = np.array([-4.9244637e-3,  -1.2984142e-6, -1.5701595e-6, 1.5535974e-8,
                        -2.7221769e-2, 4.247473e-4, -3.9583181e-6, 1.7295795e-8])

def density(z):
    """
        Computes the atmospheric density according to the US standard atmosphere model using a polynomial fit

        Parameters
        ----------
        z: float
            Altitude above sea level [km]

        Returns
        -------
        density: float
            Density of the atmosphere at altitude z [g/cm^3] 
    """

    poly_A, poly_B = 0.0, 1.0
    for n in range(4):
        poly_A += coeffs_A[n] * z**(n + 1)
        poly_B += coeffs_B[n] * z**(n + 1)

    return den0 * 10.0**(poly_A / poly_B)

def pressure(z, T):
    """
        Computes the atmospheric pressure according to the US standard atmosphere model using a polynomial fit assuming an ideal gas

        Parameters
        ----------
        z: float
            Altitude above sea level [km]

        Returns
        -------
        pressure: float
            Pressure of the atmosphere at altitude :math:`z` [mbar] and temperature :math:`T` [K]
    """
     
    return density(z) * gasR * T * 10.0
    

def profiles_qc(path, pattern="*.dat", skiprows=0):
    """
        Runs a quality control (QC) check on profiles in the path
        matching the pattern.  It can optionally plot the bad
        profiles.  If it finds any, it makes a new direcotry
        in the path location called "bad_profs" and moves those
        profiles into the directory for you to check

        Parameters
        ----------
        path: string
            Path to the profiles to be QC'd
        pattern: string
            Pattern defining the list of profiles in the path
        skiprows: int
            Number of header rows in the profiles

        Returns
        -------

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
            dTdz = np.gradient(atmo[:, 1]) / np.gradient(atmo[:, 0])
            dudz = np.gradient(atmo[:, 2]) / np.gradient(atmo[:, 0])
            dvdz = np.gradient(atmo[:, 3]) / np.gradient(atmo[:, 0])

            if max(abs(dTdz)) > 50.0 or max(abs(dudz)) > 50.0 or max(abs(dvdz)) > 50.0:
                if not os.path.isdir(path + "bad_profs/"):
                    subprocess.call("mkdir " + path + "bad_profs/", shell=True)

                subprocess.call("mv " + path + file + " " + path + "bad_profs/", shell=True)
                print("Bad profile:", file)
    else:
        print("WARNING!! No profiles matching specified pattern in path.")




def build_atmo_matrix(path, pattern="*.dat", years=None, months=None, weeks=None, hours=None, max_alt=None, skiprows=0, ref_alts=None, prof_format="zTuvdp", latlon0=None, return_datetime=False):
    """
        Read in a list of atmosphere files from the path location
        matching a specified pattern for continued analysis.

        Parameters
        ----------
        path: string
            Path to the profiles to be loaded
        pattern: string
            Pattern defining the list of profiles in the path
        years: iterable
            Iterable of years to include in analysis (e.g., ['2010', '2011', '2012'])
        months: iterable
            Iterable of months to include in analysis (e.g., ['11', '12', '01', '02'])
        weeks: iterable
            Iterable of weeks to include in analysis (e.g., ['01', '02', '03'])
        hours: iterable
            Iterable of hours to include in analysis (e.g., ['00', '06', '18'])
        max_alt: float
            Altitude maximum to trim specifications
        skiprows: int
            Number of header rows in the profiles
        ref_alts: 1darray
            Reference altitudes if comparison is needed
        prof_format: string
            Profile format is either 'ECMWF' or column specifications (e.g., 'zTuvdp')
        return_datetime: bool
            Option to return the datetime info of ingested atmosphere files for future reference

        Returns
        -------
        A: 2darray
            Atmosphere array of size M x (5 * N) for M atmospheres where each atmosphere samples N altitudes
        z: 1darray
            Altitude reference values [km]
        datetime: 1darray
            List of dates and times for each specification in the matrix (optional output, see Parameters)
    """

    print('\t' + "Loading profiles from " + path + " with pattern: " + pattern)
    if years is not None:
        print('\t\t' + "Years filter:", years)
    if months is not None:
        print('\t\t' + "Months filter:", months)
    if weeks is not None:
        print('\t\t' + "Weeks filter:", weeks)
    if hours is not None:
        print('\t\t' + "Hours filter:", hours)

    file_list = []
    dir_files = os.listdir(path)
    for file in np.sort(dir_files):
        if fnmatch.fnmatch(file, pattern):
            date_check = False

            date_parse = re.search(r'\d{10}', file)[0]

            if len(date_parse) > 0:
                year = date_parse[:4]
                month = date_parse[4:6]
                day = date_parse[6:8]
                hour = date_parse[8:10]
                date_check = True

            if not date_check:
                file_test = open(path + file)
                for line in file_test:
                    if "# Model Time" in line:
                        year = line[15:19]
                        month = line[20:22]
                        day = line[23:25]
                        hour = line[26:28]

                        date_check = True
                        file_test.close()
                        break

            week = "%02d" % datetime.date(int(year), int(month), int(day)).isocalendar()[1]

            include_check = True 
            if years is not None:
                if year not in years:
                    include_check = False 

            if months is not None:
                if month not in months:
                    include_check = False
            
            if weeks is not None:
                if week not in weeks:
                    include_check = False

            if hours is not None:
                if hour not in hours:
                    include_check = False
            
            if include_check:
                file_list += [file]

    if len(file_list) > 0:
        file_list = np.sort(file_list)
        if prof_format == "ecmwf" or prof_format == "ECMWF":
            # Add a check and download here for the ETOPO file
            etopo_file = imp.find_module('stochprop')[1] + '/resources/ETOPO1_Ice_g_gmt4.grd'

            etopo1 = Dataset(etopo_file)

            grid_lons = etopo1.variables['x'][:]
            grid_lats = etopo1.variables['y'][:]
            grid_elev = etopo1.variables['z'][:]

            # Change underwater values to sea level
            grid_elev[grid_elev < 0.0] = 0.0

            # Interpolate and evaluate the ground level elevation at the specified latitude and longitude
            lat_mask = np.logical_and(latlon0[0] - 1.0 <= grid_lats, grid_lats <= latlon0[0] + 1.0).nonzero()[0]
            lon_mask = np.logical_and(latlon0[1] - 1.0 <= grid_lons, grid_lons <= latlon0[1] + 1.0).nonzero()[0]

            etopo_interp = interp2d(grid_lons[lon_mask], grid_lats[lat_mask], grid_elev[lat_mask,:][:,lon_mask] / 1000.0, kind='linear')

            # Load ECMWF file and identify indices of nearest node for specified loccation
            if return_datetime:
                dt_parse = re.search(r'\d{10}', file_list[0])[0]
                dt_list = [np.datetime64(dt_parse[0:4] + "-" + dt_parse[4:6] + "-" + dt_parse[6:8] + "T" + dt_parse[8:10] + ":00:00")]

            ecmwf = Dataset(path + file_list[0])

            lat0, dlat = float(ecmwf.variables['ylat0'][:].data), float(ecmwf.variables['dy'][:].data)
            lon0, dlon = float(ecmwf.variables['xlon0'][:].data), float(ecmwf.variables['dx'][:].data)
    
            lat_vals = np.arange(lat0, lat0 + dlat * ecmwf.variables['T'].shape[1], dlat)
            lon_vals = np.arange(lon0, lon0 + dlon * ecmwf.variables['T'].shape[2], dlon)

            n_lat = np.argmin(abs(lat_vals - latlon0[0]))
            n_lon = np.argmin(abs(lon_vals - latlon0[1]))
            z_gl = etopo_interp(lon_vals[n_lon], lat_vals[n_lat])[0]

            z0 = ecmwf.variables['height'][:].data / 1000.0 + z_gl    
            T = ecmwf.variables['T'][:, n_lat, n_lon].data   
            u = ecmwf.variables['U'][:, n_lat, n_lon].data   
            v = ecmwf.variables['V'][:, n_lat, n_lon].data   
            d = density(z0)
            p = pressure(z0, T)

            for file in file_list[1:]:
                if return_datetime:
                    dt_parse = re.search(r'\d{10}', file)[0]
                    dt_list = dt_list + [[np.datetime64(dt_parse[0:4] + "-" + dt_parse[4:6] + "-" + dt_parse[6:8] + "T" + dt_parse[8:10] + ":00:00")]]
                ecmwf = Dataset(path + file)
                
                lat0, dlat = float(ecmwf.variables['ylat0'][:].data), float(ecmwf.variables['dy'][:].data)
                lon0, dlon = float(ecmwf.variables['xlon0'][:].data), float(ecmwf.variables['dx'][:].data)
    
                lat_vals = np.arange(lat0, lat0 + dlat * ecmwf.variables['T'].shape[1], dlat)
                lon_vals = np.arange(lon0, lon0 + dlon * ecmwf.variables['T'].shape[2], dlon)

                n_lat = np.argmin(abs(lat_vals - latlon0[0]))
                n_lon = np.argmin(abs(lon_vals - latlon0[1]))

                T = np.vstack((T, ecmwf.variables['T'][:, n_lat, n_lon].data))
                u = np.vstack((u, ecmwf.variables['U'][:, n_lat, n_lon].data))
                v = np.vstack((v, ecmwf.variables['V'][:, n_lat, n_lon].data))
                d = np.vstack((d, density(z0)))
                p = np.vstack((p, pressure(z0, ecmwf.variables['T'][:, n_lat, n_lon].data)))

            A = np.hstack((T, u))
            A = np.hstack((A, v))
            A = np.hstack((A, d))
            A = np.hstack((A, p))

            A = np.atleast_2d(A)
            return A, z0
        else:
            if return_datetime:
                dt_parse = re.search(r'\d{10}', file_list[0])[0]
                dt_list = [np.datetime64(dt_parse[0:4] + "-" + dt_parse[4:6] + "-" + dt_parse[6:8] + "T" + dt_parse[8:10] + ":00:00")]

            # add parser to determine indices of fields of interest (T or p, u, v, d)
            atmo = np.loadtxt(path + file_list[0], skiprows=skiprows)
            if np.any(ref_alts) is None:
                z0 = atmo[:, 0]
            else:
                z0 = ref_alts

            if max_alt is not None:
                alt_mask = tuple([z0 <= max_alt])
                print('\t\t' + "Trimming above maximum altitude:", max_alt)
            else:
                alt_mask = np.ones_like(z0, dtype=bool)
            
            if np.allclose(z0, atmo[:, 0]):
                T = atmo[:, 1][alt_mask]
                u = atmo[:, 2][alt_mask]
                v = atmo[:, 3][alt_mask]
                d = atmo[:, 4][alt_mask]
                p = atmo[:, 5][alt_mask]
            else:
                print("WARNING!!  Altitudes in " + path + file_list[0] + " don't match expected values.  Interpolating to resolve...")
                T = interp1d(atmo[:, 0], atmo[:, 1])(z0)[alt_mask]
                u = interp1d(atmo[:, 0], atmo[:, 2])(z0)[alt_mask]
                v = interp1d(atmo[:, 0], atmo[:, 3])(z0)[alt_mask]
                d = interp1d(atmo[:, 0], atmo[:, 4])(z0)[alt_mask]
                p = interp1d(atmo[:, 0], atmo[:, 5])(z0)[alt_mask]

            for file in file_list[1:]:
                if return_datetime:
                    dt_parse = re.search(r'\d{10}', file)[0]
                    dt_list = dt_list + [np.datetime64(dt_parse[0:4] + "-" + dt_parse[4:6] + "-" + dt_parse[6:8] + "T" + dt_parse[8:10] + ":00:00")]

                atmo = np.loadtxt(path + file, skiprows=skiprows)
                if np.allclose(z0, atmo[:, 0]):
                    T = np.vstack((T, atmo[:, 1][alt_mask]))
                    u = np.vstack((u, atmo[:, 2][alt_mask]))
                    v = np.vstack((v, atmo[:, 3][alt_mask]))
                    d = np.vstack((d, atmo[:, 4][alt_mask]))
                    p = np.vstack((p, atmo[:, 5][alt_mask]))
                else:
                    print("WARNING!!  Altitudes in " + path + file + " don't match expected values.  Interpolating to resolve...")
                    T = np.vstack((T, interp1d(atmo[:, 0], atmo[:, 1])(z0)[alt_mask]))
                    u = np.vstack((u, interp1d(atmo[:, 0], atmo[:, 2])(z0)[alt_mask]))
                    v = np.vstack((v, interp1d(atmo[:, 0], atmo[:, 3])(z0)[alt_mask]))
                    d = np.vstack((d, interp1d(atmo[:, 0], atmo[:, 4])(z0)[alt_mask]))
                    p = np.vstack((p, interp1d(atmo[:, 0], atmo[:, 5])(z0)[alt_mask]))

            A = np.hstack((T, u))
            A = np.hstack((A, v))
            A = np.hstack((A, d))
            A = np.hstack((A, p))

        A = np.atleast_2d(A)
        if return_datetime:
            return A, z0[alt_mask], np.array(dt_list)
        else:
            return A, z0[alt_mask]
    else:
        return None, None


def compute_eofs(A, alts, output_path, eof_cnt=100):
    """
        Computes the singular value decomposition (SVD)
        of an atmosphere set read into an array by
        stochprop.eofs.build_atmo_matrix() and saves
        the basis functions (empirical orthogonal
        functions) and singular values to file

        Parameters
        ----------
        A: 2darray
            Suite of atmosphere specifications from build_atmo_matrix
        alts: 1darray
            Altitudes at which the atmosphere is sampled from build_atmo_matrix
        output_path: string
            Path to output the SVD results
        eof_cnt: int
            Number of basic functions to save
    """

    print('\t' + "Building EOFs using SVD...")

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

    cT_vals = np.sqrt(gamR * T_vals)
    cp_vals = np.sqrt((gam / 10.0) * p_vals / d_vals)

    u_diff = u_vals - np.array([u_mean] * u_vals.shape[0])
    v_diff = v_vals - np.array([v_mean] * v_vals.shape[0])
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
    u_eofs = eofs[:, file_len * 2:file_len * 3] / np.sqrt(abs(alts[1] - alts[0]))
    v_eofs = eofs[:, file_len * 3:file_len * 4] / np.sqrt(abs(alts[1] - alts[0]))

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
        A: 2darray
            Suite of atmosphere specifications from build_atmo_matrix
        alts: 1darray
            Altitudes at which the atmosphere is sampled from build_atmo_matrix
        eofs_path: string
            Path to the .eof results from compute_eofs
        output_path: string
            Path where output will be stored
        eof_cnt: int
            Number of EOFs to consider in computing coefficients
        pool: pathos.multiprocessing.ProcessingPool
            Multiprocessing pool for accelerating calculations

        Returns
        -------
        coeffs: 2darray
            Array containing coefficient values of size prof_cnt by eof_cnt.  Result is also written to file.
    """

    print('\t' + "Computing EOF coefficients for profiles...")
    # load means and eofs
    means = np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")

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
    for n, An in enumerate(A):
        if A.shape[0] <= 20 :
            print('\t\t' + "Currently on profile " + str(n + 1) + " of " + str(len(A)) + "...")
        elif (n + 1) % 10 == 0:
            print('\t\t' + "Currently on profile " + str(n + 1) + " of " + str(len(A)) + "...")

        # interpolate the profile
        T_interp = interp1d(A_alts, An[file_len * 0:file_len * 1][A_mask], kind='cubic')
        u_interp = interp1d(A_alts, An[file_len * 1:file_len * 2][A_mask], kind='cubic')
        v_interp = interp1d(A_alts, An[file_len * 2:file_len * 3][A_mask], kind='cubic')
        d_interp = interp1d(A_alts, An[file_len * 3:file_len * 4][A_mask], kind='cubic')
        p_interp = interp1d(A_alts, An[file_len * 4:file_len * 5][A_mask], kind='cubic')

        # evaluate differences at the sampled points of the mean/eofs
        cT_diff = np.sqrt(gamR * T_interp(eof_alts))
        cT_diff -= np.sqrt(gamR * means[:, 1][eofs_mask])

        cp_diff = np.sqrt(gam / 10.0 * p_interp(eof_alts) / d_interp(eof_alts))
        cp_diff -= np.sqrt(gam / 10.0 * means[:, 5][eofs_mask] / means[:, 4][eofs_mask])

        u_diff = u_interp(eof_alts) - means[:, 2][eofs_mask]
        v_diff = v_interp(eof_alts) - means[:, 3][eofs_mask]

        # define the integration to compute the EOF coefficients
        def calc_coeff(n):
            result = simps(cT_diff * cT_eofs[:, n + 1][eofs_mask], eof_alts)
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

    # store coefficient values as numpy array file
    np.save(output_path + "-coeffs", coeffs)

    return coeffs


def compute_overlap(coeffs, eofs_path, eof_cnt=100, method="mean"):
    """
        Compute the overlap of EOF coefficient distributions

        Parameters
        ----------
        coeffs: list of 2darrays
            List of 2darrays containing coefficients to consider
                overlap in PDF of values
        eofs_path: string
            Path to the .eof results from compute_eofs
        eof_cnt: int
            Number of EOFs to compute
        method : string
            Option to decide which overlap to use ("kde" or "mean")
            

        Returns
        -------
        overlap: 3darray
            Array containing overlap values of size coeff_cnt by coeff_cnt by eof_cnt
    """

    print("-" * 50 + '\n' + "Computing coefficient overlap...")

    overlap = np.ones((12, 12))
    
    eof_weights = np.loadtxt(eofs_path + "-singular_values.dat")[:, 1][:eof_cnt]**2
    eof_weights /= np.sum(eof_weights)  

    if method == "mean":
        for m1 in range(12):
            m1_C_means = np.mean(coeffs[m1], axis=0)          
            norm1 = np.dot(m1_C_means, m1_C_means)
            for m2 in range(m1 + 1, 12):
                m2_C_means = np.mean(coeffs[m2], axis=0)
                norm2 = np.dot(m2_C_means, m2_C_means)

                overlap[m1][m2] = (np.dot(eof_weights * m1_C_means, m2_C_means) / np.sqrt(norm1 * norm2) + 1) / 2.0
                overlap[m2][m1] = overlap[m1][m2]

    else: 
        overlap_temp = np.ones((eof_cnt, 12, 12))
        for eof_id in range(eof_cnt):
            print('\t' + "Computing overlap values for EOF " + str(eof_id) + " using KDE...")
            for m1 in range(12):
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
                    overlap_temp[eof_id][m1][m2] = simps(kernel1.pdf(c_vals) * kernel2.pdf(c_vals), c_vals) / np.sqrt(norm1 * norm2)
                    overlap_temp[eof_id][m2][m1] = overlap_temp[eof_id][m1][m2]

        overlap = np.average(overlap_temp, weights=eof_weights, axis=0)

    return overlap


def compute_seasonality(overlap_file, file_id=None):
    """
        Compute the overlap of EOF coefficients to identify seasonality

        Parameters
        ----------
        overlap_file: string
            Path and name of file containing results of stochprop.eofs.compute_overlap
        file_id: string
            Path and ID to save the dendrogram result of the overlap analysis
    """

    print("Generating seasonality plot...")

    dist_mat = -np.log(np.load(overlap_file))
    np.fill_diagonal(dist_mat, 0.0)   
    links = hierarchy.linkage(squareform(dist_mat), 'weighted')

    f, (ax1) = plt.subplots(1, 1)
    den = hierarchy.dendrogram(links, orientation="right", ax=ax1, labels=[calendar.month_abbr[n] for n in range(1, 13)])
    plt.title(file_id.rpartition('/')[-1])

    if file_id:
        plt.savefig(file_id + "-seasonality.png", dpi=300)

    plt.show(block=False)
    plt.pause(5.0)
    plt.close('all')


################################
#   Define methods to sample   #
#     the coefficient PDFs     #
################################
def _coeff_smpl_header_txt(coeff_label, eofs_path, eof_cnt, n, prof_cnt):
    result = "# Data Source: stochprop v" + pkg_resources.get_distribution("stochprop").version
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Coefficient KDE Sampling"
    result = result + '\n' + "# Coeff Label = " + coeff_label
    result = result + '\n' + "# EOF Set = " + eofs_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# EOF Cnt = " + str(eof_cnt)
    result = result + '\n' + "# Sample: " + str(n) + "/" + str(prof_cnt)
    result = result + '\n' + "# Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm3), P(mbar) ]"
    result = result + '\n' + "# The following lines are formatted input for ncpaprop"
    result = result + '\n' + "#% 0, Z0, km, 0.0"
    result = result + '\n' + "#% 1, Z, km"
    result = result + '\n' + "#% 2, T, degK"
    result = result + '\n' + "#% 3, U, m/s"
    result = result + '\n' + "#% 4, V, m/s"
    result = result + '\n' + "#% 5, RHO, g/cm3"
    result = result + '\n' + "#% 6, P, mbar"

    return result

def _coeff_smpl_mean_header_txt(coeff_label, eofs_path, eof_cnt):
    result = "# Data Source: stochprop v" + pkg_resources.get_distribution("stochprop").version
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Coefficient KDE Sampling"
    result = result + '\n' + "# Coeff Label = " + coeff_label
    result = result + '\n' + "# EOF Set = " + eofs_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# EOF Cnt = " + str(eof_cnt)
    result = result + '\n' + "# Sample: mean"
    result = result + '\n' + "# Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm3), P(mbar) ]"
    result = result + '\n' + "# The following lines are formatted input for ncpaprop"
    result = result + '\n' + "#% 0, Z0, km, 0.0"
    result = result + '\n' + "#% 1, Z, km"
    result = result + '\n' + "#% 2, T, degK"
    result = result + '\n' + "#% 3, U, m/s"
    result = result + '\n' + "#% 4, V, m/s"
    result = result + '\n' + "#% 5, RHO, g/cm3"
    result = result + '\n' + "#% 6, P, mbar"

    return result
    

def _coeff_smpl_max_header_txt(coeff_label, eofs_path, eof_cnt):
    result = "# Data Source: stochprop v" + pkg_resources.get_distribution("stochprop").version
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Coefficient KDE Sampling"
    result = result + '\n' + "# Coeff Label = " + coeff_label
    result = result + '\n' + "# EOF Set = " + eofs_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# EOF Cnt = " + str(eof_cnt)
    result = result + '\n' + "# Sample: max likelihood"
    result = result + '\n' + "# Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm3), P(mbar) ]"
    result = result + '\n' + "# The following lines are formatted input for ncpaprop"
    result = result + '\n' + "#% 0, Z0, km, 0.0"
    result = result + '\n' + "#% 1, Z, km"
    result = result + '\n' + "#% 2, T, degK"
    result = result + '\n' + "#% 3, U, m/s"
    result = result + '\n' + "#% 4, V, m/s"
    result = result + '\n' + "#% 5, RHO, g/cm3"
    result = result + '\n' + "#% 6, P, mbar"

    return result

def _fit_header_txt(prof_path, eofs_path, eof_cnt):
    result = "# Data Source: stochprop v" + pkg_resources.get_distribution("stochprop").version
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Fitting"
    result = result + '\n' + "# Reference Specification = " + prof_path
    result = result + '\n' + "# EOF Set = " + eofs_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# EOF Cnt = " + str(eof_cnt)
    result = result + '\n' + "# Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm3), P(mbar) ]"
    result = result + '\n' + "# The following lines are formatted input for ncpaprop"
    result = result + '\n' + "#% 0, Z0, km, 0.0"
    result = result + '\n' + "#% 1, Z, km"
    result = result + '\n' + "#% 2, T, degK"
    result = result + '\n' + "#% 3, U, m/s"
    result = result + '\n' + "#% 4, V, m/s"
    result = result + '\n' + "#% 5, RHO, g/cm3"
    result = result + '\n' + "#% 6, P, mbar"

    return result

def _perturb_header_txt(prof_path, eofs_path, eof_cnt, stdev, n, prof_cnt):
    result = "# Data Source: stochprop v" + pkg_resources.get_distribution("stochprop").version
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: EOF Perturbation"
    result = result + '\n' + "# Reference Specification = " + prof_path
    result = result + '\n' + "# EOF Set = " + eofs_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# EOF Cnt = " + str(eof_cnt)
    result = result + '\n' + "# Perturbation St Dev (winds) = " + str(stdev) + " m/s"
    result = result + '\n' + "# Sample: " + str(n) + "/" + str(prof_cnt)
    result = result + '\n' + "# Fields = [ Z(km), T(K), U(m/s), V(m/s), R(g/cm3), P(mbar) ]"
    result = result + '\n' + "# The following lines are formatted input for ncpaprop"
    result = result + '\n' + "#% 0, Z0, km, 0.0"
    result = result + '\n' + "#% 1, Z, km"
    result = result + '\n' + "#% 2, T, degK"
    result = result + '\n' + "#% 3, U, m/s"
    result = result + '\n' + "#% 4, V, m/s"
    result = result + '\n' + "#% 5, RHO, g/cm3"
    result = result + '\n' + "#% 6, P, mbar"

    return result

def define_coeff_limits(coeff_vals):
    """
        Compute upper and lower bounds for coefficient values

        Parameters
        ----------
        coeff_vals: 2darrays
            Coefficients computed with stochprop.eofs.compute_coeffs

        Returns
        -------
        lims: 1darray
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
        pdf: function
            Probability distribution function (PDF) for a single variable
        lims: 1darray
            Iterable containing lower and upper bound for integration
        pnts: int
            Number of points to consider in defining the cumulative distribution

        Returns
        -------
        cfd: interp1d
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
        pdf: function
            Probability distribution function (PDF) for a single variable
        lims: 1darray
            Iterable containing lower and upper bound for integration
        cdf: function
            Cumulative distribution function (CDF) from stochprop.eofs.build_cfd
        size: int
            Number of samples to generate

        Returns
        -------
        samples: 1darray
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


def sample_atmo(coeffs, eofs_path, output_path, eof_cnt=100, prof_cnt=250, output_mean=False, coeff_label="None"):
    """
        Generate atmosphere states using coefficient distributions for
        a set of empirical orthogonal basis functions

        Parameters
        ----------
        coeffs: 2darrays
            Coefficients computed with stochprop.eofs.compute_coeffs
        eofs_path: string
            Path to the .eof results from compute_eofs
        output_path: string
            Path where output will be stored
        eof_cnt: int
            Number of EOFs to use in building sampled specifications
        prof_cnt: int
            Number of atmospheric specification samples to generate
        output_mean: bool
            Flag to output the mean profile from the samples generated
    """

    print("-" * 50 + '\n' + "Generating atmosphere state samples from coefficient PDFs...")
    print('\t' + "Loading mean profile info and eofs...")
    means = np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")

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
    cp0 = np.sqrt(gam / 10.0 * means[:, 5] / means[:, 4])

    sampled_profs = np.array([np.copy(means)] * prof_cnt)
    sampled_cTs = np.array([cT0] * prof_cnt)
    sampled_cps = np.array([cp0] * prof_cnt)

    for eof_id in range(eof_cnt):
        print('\t\t' + "Sampling EOF coefficient for eof_id = " + str(eof_id) + "...")
        sampled_coeffs = draw_from_pdf(kernels[eof_id].pdf, lims[eof_id], size=prof_cnt)

        for pn in range(prof_cnt):
            sampled_cTs[pn] = sampled_cTs[pn] + sampled_coeffs[pn] * cT_eofs[:, eof_id + 1]
            sampled_cps[pn] = sampled_cps[pn] + sampled_coeffs[pn] * cp_eofs[:, eof_id + 1]

            sampled_profs[pn][:, 2] = sampled_profs[pn][:, 2] + sampled_coeffs[pn] * u_eofs[:, eof_id + 1]
            sampled_profs[pn][:, 3] = sampled_profs[pn][:, 3] + sampled_coeffs[pn] * v_eofs[:, eof_id + 1]

    sampled_profs[:, :, 1] = sampled_cTs**2 / gamR
    sampled_profs[:, :, 5] = sampled_profs[:, :, 4] * sampled_cps**2 / (gam / 10.0)

    # save the individual profiles and the mean profile
    print('\t' + "Writing sampled atmospheres to file...", '\n')
    for pn in range(prof_cnt):
        np.savetxt(output_path + "-" + "%02d" % pn + ".met", sampled_profs[pn], header=_coeff_smpl_header_txt(coeff_label, eofs_path, eof_cnt, pn, prof_cnt), comments='')

    if output_mean:
        np.savetxt(output_path + "-mean.met", np.average(sampled_profs, axis=0), header=_coeff_smpl_mean_header_txt(coeff_label, eofs_path, eof_cnt), comments='')


def maximum_likelihood_profile(coeffs, eofs_path, output_path, eof_cnt=100, coeff_label="None"):
    """
        Use coefficient distributions for a set of empirical orthogonal
        basis functions to compute the maximum likelihood specification

        Parameters
        ----------
        coeffs: 2darrays
            Coefficients computed with stochprop.eofs.compute_coeffs
        eofs_path: string
            Path to the .eof results from compute_eofs
        output_path: string
            Path where output will be stored
        eof_cnt: int
            Number of EOFs to use in building sampled specifications
    """

    print("-" * 50 + '\n' + "Generating maximum likelihood atmosphere states from coefficient PDFs...")
    # load mean profile and eofs
    print('\t' + "Loading mean profile info and eofs...")
    means = np.loadtxt(eofs_path + "-mean_atmo.dat")
    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")

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

        c_vals = np.linspace(lims[0], lims[1], 5001)
        coeff_ml = c_vals[np.argmax(kernel.pdf(c_vals))]

        cT_ml = cT_ml + coeff_ml * cT_eofs[:, n + 1]
        cp_ml = cp_ml + coeff_ml * cp_eofs[:, n + 1]

        ml_prof[:, 2] = ml_prof[:, 2] + coeff_ml * u_eofs[:, n + 1]
        ml_prof[:, 3] = ml_prof[:, 3] + coeff_ml * v_eofs[:, n + 1]

    ml_prof[:, 1] = cT_ml**2 / gamR
    ml_prof[:, 5] = ml_prof[:, 4] * cp_ml**2 / (gam / 10.0)

    print('\t' + "Writing maximum likelihood atmosphere to file...", '\n')
    np.savetxt(output_path + "-maximum_likelihood.met", ml_prof, header=_coeff_smpl_max_header_txt(coeff_label, eofs_path, eof_cnt), comments='')


# ############################ #
#   Define methods to fit or   #
#    perturb specific atmos    #
# ############################ #
def fit_atmo(prof_path, eofs_path, output_path, eof_cnt=100):
    """
        Compute a given number of EOF coefficients to fit a given
        atmophere specification using the basic functions.  Write
        the resulting approximated atmospheric specification to
        file.

        Parameters
        ----------
        prof_path: string
            Path and name of the specification to be fit
        eofs_path: string
            Path to the .eof results from compute_eofs
        output_path: string
            Path where output will be stored
        eof_cnt: int
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

    cT_diff = np.sqrt(gamR * T_interp(means[:, 0][eofs_mask]))
    cT_diff -= np.sqrt(gamR * means[:, 1][eofs_mask])

    cp_diff = np.sqrt(gam / 10.0 * p_interp(means[:, 0][eofs_mask]) / d_interp(means[:, 0][eofs_mask]))
    cp_diff -= np.sqrt(gam / 10.0 * means[:, 5][eofs_mask] / means[:, 4][eofs_mask])

    u_diff = u_interp(means[:, 0][eofs_mask]) - means[:, 2][eofs_mask]
    v_diff = v_interp(means[:, 0][eofs_mask]) - means[:, 3][eofs_mask]

    # define the integration to compute the EOF coefficients
    def calc_coeff(n):
        result = simps(cT_diff * cT_eofs[:, n + 1][eofs_mask], cT_eofs[:, 0][eofs_mask])
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

    np.savetxt(output_path, fit, header=_fit_header_txt(prof_path, eofs_path, eof_cnt), comments='')


def perturb_atmo(prof_path, eofs_path, output_path, stdev=10.0, eof_max=100, eof_cnt=50, sample_cnt=1, alt_wt_pow=2.0, sing_val_wt_pow=0.25):
    """
        Use EOFs to perturb a specified profile using a given scale

        Parameters
        ----------
        prof_path: string
            Path and name of the specification to be fit
        eofs_path: string
            Path to the .eof results from compute_eofs
        output_path: string
            Path where output will be stored
        stdev: float
            Standard deviation of wind speed used to scale perturbation
        eof_max: int
            Higher numbered EOF to sample
        eof_cnt: int
            Number of EOFs to sample in the perturbation (can be less than eof_max)
        sample_cnt: int
            Number of perturbed atmospheric samples to generate
        alt_wt_pow: float
            Power raising relative mean altitude value in weighting
        sing_val_wt_pow: float
            Power raising relative singular value in weighting
    """

    print("Generating EOF perturbations to " + prof_path + "...")
    ref_atmo = np.loadtxt(prof_path)

    cT_eofs = np.loadtxt(eofs_path + "-ideal_gas_snd_spd.eofs")
    cp_eofs = np.loadtxt(eofs_path + "-adiabatic_snd_spd.eofs")
    u_eofs = np.loadtxt(eofs_path + "-zonal_winds.eofs")
    v_eofs = np.loadtxt(eofs_path + "-merid_winds.eofs")
    sing_vals = np.loadtxt(eofs_path + "-singular_values.dat")

    # Define altitude limits
    alt_min = max(u_eofs[0, 0], ref_atmo[0, 0])
    alt_max = min(u_eofs[-1, 0], ref_atmo[-1, 0])

    eof_mask = np.logical_and(alt_min <= u_eofs[:, 0], u_eofs[:, 0] <= alt_max)
    ref_mask = np.logical_and(alt_min <= ref_atmo[:, 0], ref_atmo[:, 0] <= alt_max)

    z_vals = np.copy(ref_atmo[:, 0][ref_mask])

    # interpolate the eofs and use random coefficients to generate perturbations
    cT_perturb_all = np.empty((sample_cnt, len(z_vals)))
    cp_perturb_all = np.empty((sample_cnt, len(z_vals)))
    u_perturb_all = np.empty((sample_cnt, len(z_vals)))
    v_perturb_all = np.empty((sample_cnt, len(z_vals)))

    for m in range(sample_cnt):
        wts = np.empty(eof_cnt)

        cT_perturb = np.zeros((eof_cnt, len(z_vals)))
        cp_perturb = np.zeros((eof_cnt, len(z_vals)))
        u_perturb = np.zeros((eof_cnt, len(z_vals)))
        v_perturb = np.zeros((eof_cnt, len(z_vals)))

        for j, n in enumerate(np.random.choice(range(eof_max), eof_cnt, replace=False)):
            coeff_val = np.random.randn()

            cT_perturb[j] = coeff_val * interp1d(cT_eofs[:, 0][eof_mask], cT_eofs[:, n + 1][eof_mask])(z_vals)
            cp_perturb[j] = coeff_val * interp1d(cp_eofs[:, 0][eof_mask], cp_eofs[:, n + 1][eof_mask])(z_vals)
            u_perturb[j] = coeff_val * interp1d(u_eofs[:, 0][eof_mask], u_eofs[:, n + 1][eof_mask])(z_vals)
            v_perturb[j] = coeff_val * interp1d(v_eofs[:, 0][eof_mask], v_eofs[:, n + 1][eof_mask])(z_vals)

            wts[j] = simps(u_eofs[:, 0][eof_mask] * abs(u_eofs[:, n + 1][eof_mask]), u_eofs[:, 0][eof_mask]) / simps(abs(u_eofs[:, n + 1][eof_mask]), u_eofs[:, 0][eof_mask])
            wts[j] += simps(u_eofs[:, 0][eof_mask] * abs(v_eofs[:, n + 1][eof_mask]), v_eofs[:, 0][eof_mask]) / simps(abs(v_eofs[:, n + 1][eof_mask]), u_eofs[:, 0][eof_mask])
            wts[j] /= (2.0 * max(u_eofs[:, 0][eof_mask]))
            wts[j] = wts[j]**alt_wt_pow

            wts[j] *= (sing_vals[n, 1] / sing_vals[0, 1])**sing_val_wt_pow

        cT_perturb_all[m] = np.average(cT_perturb, axis=0, weights=wts)
        cp_perturb_all[m] = np.average(cp_perturb, axis=0, weights=wts)
        u_perturb_all[m] = np.average(u_perturb, axis=0, weights=wts)
        v_perturb_all[m] = np.average(v_perturb, axis=0, weights=wts)

    scaling = stdev / (np.average(np.sqrt(u_perturb_all**2 + v_perturb_all**2)))

    for m in range(sample_cnt):
        T_vals = np.copy(ref_atmo[:, 1][ref_mask])
        u_vals = np.copy(ref_atmo[:, 2][ref_mask])
        v_vals = np.copy(ref_atmo[:, 3][ref_mask])
        d_vals = np.copy(ref_atmo[:, 4][ref_mask])
        p_vals = np.copy(ref_atmo[:, 5][ref_mask])

        cT_vals = np.sqrt(gamR * T_vals)
        cp_vals = np.sqrt(gam / 10.0 * p_vals / d_vals)

        cT_vals = cT_vals + cT_perturb_all[m] * scaling
        cp_vals = cp_vals + cp_perturb_all[m] * scaling
        u_vals = u_vals + u_perturb_all[m] * scaling
        v_vals = v_vals + v_perturb_all[m] * scaling

        T_vals = cT_vals**2 / gamR
        p_vals = d_vals * cp_vals**2 / (gam / 10.0)

        np.savetxt(output_path + "-" + str(m) + ".met", np.vstack((z_vals, T_vals, u_vals, v_vals, d_vals, p_vals)).T, header=_perturb_header_txt(prof_path, eofs_path, eof_cnt, stdev, m, sample_cnt), comments='')
