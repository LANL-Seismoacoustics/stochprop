#!/usr/bin/env python
#
# Methods to use G2S profiles for inversion framework
#
# Philip Blom (pblom@lanl.gov)

import sys
import os
import warnings
import copy
import time
import itertools

import numpy as np

from scipy.interpolate import interp1d
from scipy.optimize import newton

from pyproj import Geod

from infrapy.utils import prog_bar

wgs84_proj = Geod(ellps='WGS84')

#############################################
### Read in a G2S profile and set of EOFs ###
##### to build a perturbable atmosphere #####
#############################################
class AtmosphereState(object):
    def __init__(self, g2s_file, zonal_eof_file, merid_eof_file, g2s_format="zTuvdp"):
        print("Loading G2S atmosphere from " + g2s_file + " using format '" + g2s_format + "' and interpolating...")
        if g2s_format == "zTuvdp":
            self.z_vals, self.T_vals, self.u_vals, self.v_vals, self.d_vals, self.p_vals = np.loadtxt(g2s_file, unpack=True)
        elif g2s_format == "zuvwTdp":
            self.z_vals, self.u_vals, self.v_vals, self.w_vals, self.T_vals, self.d_vals, self.p_vals = np.loadtxt(g2s_file, unpack=True)
        else:
            print("Error!  Bad G2S format provided.")

        self.T_spline = interp1d(self.z_vals, self.T_vals, kind='cubic')
        self.u_spline = interp1d(self.z_vals, self.u_vals, kind='cubic')
        self.v_spline = interp1d(self.z_vals, self.v_vals, kind='cubic')
        self.d_spline = interp1d(self.z_vals, self.d_vals, kind='cubic')
        self.p_spline = interp1d(self.z_vals, self.p_vals, kind='cubic')

        print("Loading zonal EOFs from " + zonal_eof_file + " and interpolating...", '\t\t', end = ' ')
        input = np.loadtxt(zonal_eof_file)
        self.zonal_z_vals = input[:,0]
        self.zonal_eof_vals = input[:,1:]

        self.zonal_eof_cnt = len(self.zonal_eof_vals[0, :])
        self.zonal_eof_splines = [0] * self.zonal_eof_cnt

        prog_bar_len, prog = min(50, self.zonal_eof_cnt), 0
        prog_bar.prep(prog_bar_len)
        
        for n in range(self.zonal_eof_cnt):
            self.zonal_eof_splines[n] = interp1d(self.zonal_z_vals, self.zonal_eof_vals[:, n], kind='cubic')
            while(float(n) / self.zonal_eof_cnt * prog_bar_len >= prog):
                prog_bar.increment()
                prog += 1
        prog_bar.close()

        print("Loading meridional EOFs from " + merid_eof_file + " and interpolating...", '\t', end = ' ')
        input = np.loadtxt(merid_eof_file)
        self.merid_z_vals = input[:, 0]
        self.merid_eof_vals = input[:, 1:]
        
        self.merid_eof_cnt = len(self.merid_eof_vals[0, :])
        self.merid_eof_splines = [0] * self.merid_eof_cnt

        prog_bar_len, prog = min(50, self.merid_eof_cnt), 0
        prog_bar.prep(prog_bar_len)
        
        for n in range(self.merid_eof_cnt):
            self.merid_eof_splines[n] = interp1d(self.merid_z_vals, self.merid_eof_vals[:, n], kind='cubic')
            while(float(n) / self.merid_eof_cnt * prog_bar_len >= prog):
                prog_bar.increment()
                prog += 1
        prog_bar.close()
        print(' ')

    def c(self, z): return self.c_spline(z)
    def u(self, z): return self.u_spline(z)
    def v(self, z): return self.v_spline(z)
    def d(self, z): return self.d_spline(z)
    def p(self, z): return self.p_spline(z)

    def write_met_file(self, zonal_coeffs, merid_coeffs, file_name):
        file_out = open(file_name, 'w')
        
        for z_val in self.z_vals:
            if z_val > max(self.zonal_z_vals): break
            if z_val > max(self.merid_z_vals): break
            
            perturbed_zonal = self.u(z_val)
            perturbed_merid = self.v(z_val)
            for coeff in zonal_coeffs: perturbed_zonal += self.zonal_eof_splines[coeff[0]](z_val) * coeff[1]
            for coeff in merid_coeffs: perturbed_merid += self.merid_eof_splines[coeff[0]](z_val) * coeff[1]

            print(z_val, self.T_spline(z_val), perturbed_zonal, perturbed_merid, self.d_spline(z_val), self.p_spline(z_val), file=file_out)

        file_out.close()


########################################
####### Define a g2s profile and #######
### identify arrival characteristics ###
########################################
def compute_rng_az_dt(detection_list, source_loc, unique_latlon=False):
    # convert detection lat/lon and source lat/lon to propagation range and azimuth
    latlons = np.array([[det.lat, det.lon] for det in detection_list])
    temp = wgs84_proj.inv([source_loc[1]] * len(latlons), [source_loc[0]] * len(latlons), latlons[:, 1], latlons[:, 0])

    azimuths = np.array(temp[0])
    ranges = np.array(temp[2]) / 1000.0

    # compute relative arrival times from earliest detection in the list
    times = np.array([det.t for det in detection_list])
    rel_times = (times - min(times)).astype('m8[s]').astype(float)

    if unique_latlon:
        _, ids = np.unique(latlons, axis=0, return_index=True)
        return ranges[np.sort(ids)], azimuths[np.sort(ids)], rel_times[np.sort(ids)]
    else:
        return ranges, azimuths, rel_times

def find_arrivals(ranges, azimuths, atmosphere_state, zonal_coeffs, merid_coeffs, incl_min=0.5, incl_max=45.0, incl_step=0.1, cpu_cnt=1, show_prog=False):
    if show_prog:
        print("Computing arrivals..."+ '\t', end=' ')
    atmosphere_state.write_met_file(zonal_coeffs, merid_coeffs, "temp.met")

    array_cnt = len(ranges)
    arrivals = np.empty((0, 4))
    if show_prog:
        prog_ref = 0
        prog_bar.prep(40)

    for az in np.unique(azimuths):
        if cpu_cnt > 1 :
            os.system("mpirun -np " + str(cpu_cnt) + " infraga-accel-2d -prop temp.met azimuth=" + str(az) + " bounces=0 incl_min=" + str(incl_min) + " incl_max=" + str(incl_max) + " incl_step=" + str(incl_step) + " calc_amp=false > /dev/null")
        else:
            os.system("infraga-2d -prop temp.met azimuth=" + str(az) + " bounces=0 incl_min=" + str(incl_min) + " incl_max=" + str(incl_max) + " incl_step=" + str(incl_step) + " calc_amp=false > /dev/null")
        if show_prog:
            prog_bar.increment(int(np.floor((40.0 * (prog_ref + 1)) / len(np.unique(azimuths))) - np.floor((40.0 * prog_ref) / len(np.unique(azimuths)))))
            prog_ref += 1

        if sum(1 for line in open("temp.results.dat")) > 10:
            # If ray paths returned to the ground, load ray tracing results and
            # interpolate the range and propagation time values vs. inclination
            results = np.loadtxt("temp.results.dat")
            arrival_range_interp = interp1d(results[:, 0], results[:, 3], kind='linear')
            prop_time_interp  = interp1d(results[:, 0], results[:, 4],  kind='linear')

            # For each detection range at this azimuth, define r_n - r(incl) and find
            # roots, identify all arrivals and stack onto arrivals list
            for rng in ranges[azimuths==az]:
                def dr(incl, bnc_cnt):
                    return arrival_range_interp(incl) * bnc_cnt - rng

                incl_cnt = 100
                incl_vals = np.linspace(results[0, 0], results[-1, 0], incl_cnt)
                incl_step = incl_vals[1] - incl_vals[0]
                for m in range(incl_cnt - 1):
                    for bncs in range(4):
                        if dr(incl_vals[m], bncs) * dr(incl_vals[m + 1], bncs) <= 0.0:
                            dr_deriv = (dr(incl_vals[m] + 1.0e-2, bncs) - dr(incl_vals[m], bncs)) / 1.0e-2
                            incl_correction = max(- 2.0 * dr(incl_vals[m], bncs) / dr_deriv, 0.0)
                            incl_correction = min(incl_correction, incl_step)

                            incl_intercept = incl_vals[m] + incl_correction

                            arrivals = np.vstack((arrivals, np.array([rng, az, prop_time_interp(incl_intercept) * bncs, incl_intercept])))

        os.system("rm atmo.dat temp.results.dat")
    os.system("rm temp.met")

    if show_prog:
        prog_bar.close()

    return arrivals

def compute_src_time_bnds(ranges, times):
    src_time_min = min(times) - np.timedelta64(1, 'D')
    src_time_max = max(times)

    for n in range(len(ranges)):
        src_time_min = max(src_time_min, times[n] - np.timedelta64(int(ranges[n] / 0.2 * 1e3), 'ms'))
        src_time_max = min(src_time_max, times[n] - np.timedelta64(int(ranges[n] / 0.4 * 1e3), 'ms'))

    return src_time_min, src_time_max


def comparator(observations, src_loc, arrivals, default_err=1.0e7):
    # Use comparator to compute residuals
    uniq_ranges, uniq_azimuths, _ = compute_rng_az_dt(observations, src_loc, unique_latlon=True)
    ranges, azimuths, _ = compute_rng_az_dt(observations, src_loc, unique_latlon=False)
    times = np.array([det.t for det in observations])

    src_time_min, src_time_max = compute_src_time_bnds(ranges, times)
    dt_list = np.arange(0.0, (src_time_max - src_time_min) / np.timedelta64(1, 's'))
    rss_list = np.array([])

    resid_min = np.array([default_err] * len(observations))
    for m in range(len(observations)):
        src_time = src_time_min + np.timedelta64(int(dt_list[m] * 1e3), 'ms')

        residuals = np.array([])
        for rng in uniq_ranges:
            # Extract the propagation times for arrivals predicted at a given array
            prop_tms = arrivals[:, 2][arrivals[:, 0] == rng]
            arr_tms = np.array([src_time] * len(prop_tms))
            for k in range(len(prop_tms)):
                arr_tms[k] += np.timedelta64(int(prop_tms[k] * 1e3), 'ms')
            arr_tms = np.sort(arr_tms)

            # Extract the observed arrival times
            obs_tms = np.sort(times[ranges == rng])

            ln, kn = len(arr_tms), len(obs_tms)
            if ln == kn:
                array_residuals = (obs_tms - arr_tms).astype('m8[s]').astype(float)
            elif kn < ln:
                array_residuals = np.array([default_err] * len(obs_tms))
                for vals in itertools.product([0, 1], repeat=ln):
                    if sum(vals) == kn:
                        mask = np.array(vals).astype(bool)
                        residual_test = (obs_tms - arr_tms[mask]).astype('m8[s]').astype(float)

                        if np.linalg.norm(residual_test) < np.linalg.norm(array_residuals):
                            array_residuals = residual_test

            elif kn > ln:
                array_residuals = np.array([default_err] * len(obs_tms))
                if ln > 0:
                    for vals in itertools.product([0, 1], repeat=kn):
                        if sum(vals) == ln:
                            mask = np.array(vals).astype(bool)
                            residual_test = np.array([default_err] * kn)
                            residual_test[mask] = (obs_tms[mask] - arr_tms).astype('m8[s]').astype(float)

                            if np.linalg.norm(residual_test) < np.linalg.norm(array_residuals):
                                array_residuals = residual_test

            residuals = np.concatenate((residuals, array_residuals))

        if np.linalg.norm(residuals) < np.linalg.norm(resid_min):
            resid_min = residuals

    return resid_min

def jacobian(observations, src_loc, atmosphere_state, zonal_coeffs, merid_coeffs, incl_min=0.5, incl_max=45.0, incl_step=0.1, coeff_eps=0.5, cpu_cnt=1, default_err=1.0e7, show_prog=False):
    if show_prog:
        print("Computing Jacobian..." + '\t', end=' ')

    coeff_cnt = len(zonal_coeffs) + len(merid_coeffs)
    uniq_ranges, uniq_azimuths, _ = compute_rng_az_dt(observations, src_loc, unique_latlon=True)

    D = np.empty([len(observations), coeff_cnt])
    if show_prog:
        prog_bar.prep(40)
    for n_coeff in range(coeff_cnt):
        zonal_coeffs_up, zonal_coeffs_dn = copy.deepcopy(zonal_coeffs), copy.deepcopy(zonal_coeffs)
        merid_coeffs_up, merid_coeffs_dn = copy.deepcopy(merid_coeffs), copy.deepcopy(merid_coeffs)

        if n_coeff < len(zonal_coeffs):
            zonal_coeffs_up[n_coeff][1] += coeff_eps
            zonal_coeffs_dn[n_coeff][1] -= coeff_eps
        else:
            merid_coeffs_up[n_coeff - len(zonal_coeffs)][1] += coeff_eps
            merid_coeffs_dn[n_coeff - len(zonal_coeffs)][1] -= coeff_eps

        arrivals_up = find_arrivals(uniq_ranges, uniq_azimuths, atmosphere_state, zonal_coeffs_up, merid_coeffs_up, incl_min=incl_min, incl_max=incl_max, incl_step=incl_step, cpu_cnt=cpu_cnt, show_prog=False)
        if show_prog:
            prog_bar.increment(int(np.floor((20.0 * (n_coeff + 1)) / coeff_cnt) - np.floor((20.0 * n_coeff) / coeff_cnt)))
        arrivals_dn = find_arrivals(uniq_ranges, uniq_azimuths, atmosphere_state, zonal_coeffs_dn, merid_coeffs_dn, incl_min=incl_min, incl_max=incl_max, incl_step=incl_step, cpu_cnt=cpu_cnt, show_prog=False)
        if show_prog:
            prog_bar.increment(int(np.floor((20.0 * (n_coeff + 1)) / coeff_cnt) - np.floor((20.0 * n_coeff) / coeff_cnt)))

        resid_up = comparator(observations, src_loc, arrivals_up, default_err=default_err)
        resid_dn = comparator(observations, src_loc, arrivals_dn, default_err=default_err)

        for n_obs in range(len(observations)):
            D[n_obs][n_coeff] = (resid_up[n_obs] - resid_dn[n_obs]) / (2.0 * coeff_eps)

    prog_bar.close()
    return D

def compute_corrections(jacobian, residual, damping=1.0e-2):
    DtD = np.dot(jacobian.T, jacobian)
    DtD_diag = np.mean(np.diag(DtD)) * np.eye(DtD.shape[0])
    corrections = np.dot(np.linalg.inv(DtD + damping * DtD_diag), np.dot(jacobian.T, residual))

    return corrections

