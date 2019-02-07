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

####################################################
#### The US standard atmosphere for temperature ####
#### and density with simplistic models for the ####
##### tropospheric and stratospheric wind jets #####
####################################################
class US_Standard_Atmosphere(object):
    gamma = 1.4
    R = 287.0
    
    T0, rho0 = 288.2, 0.001225
    A_coeffs = [-3.9082017e-2, -1.1526465e-3,  3.2891937e-5, -2.0494958e-7, -4.7087295e-2, 1.250639e-3, -1.5194498e-5, 6.581877e-8];
    B_coeffs = [-4.9244637e-3, -1.2984142e-6, -1.5701595e-6,  1.5535974e-8, -2.7221769e-2, 4.247473e-4, -3.9583181e-6, 1.729580e-8];

    tropo_jet_speed = 0.0
    tropo_jet_alt = 10.0
    tropo_jet_width = 2.5
    tropo_jet_az = 45.0
        
    strat_jet_speed = 60.0
    strat_jet_alt = 60.0
    strat_jet_width = 17.5
    strat_jet_az = -90.0
    
    def __init__(self, tropo_jet = [0.0, 10.0, 2.5, 45.0], strat_jet = [60.0, 60.0, 17.5, -90.0],):
        tropo_jet_speed = tropo_jet[0]
        tropo_jet_alt = tropo_jet[1]
        tropo_jet_width = tropo_jet[2]
        tropo_jet_az = tropo_jet[3]

        strat_jet_speed = strat_jet[0]
        strat_jet_alt = strat_jet[1]
        strat_jet_width = strat_jet[2]
        strat_jet_az = strat_jet[3]

    def rho(self, z):
        A, B = 0.0, 1.0

        for n in range(4):
            A = A + self.A_coeffs[n] * pow(z, n + 1);
            B = B + self.B_coeffs[n] * pow(z, n + 1);
    
        return self.rho0 * pow(10.0, A / B);


    def T(self, z):
        A, B = 1.0, 1.0
        for n in range(4, 8):
            A = A + self.A_coeffs[n] * pow(z, n - 3);
            B = B + self.B_coeffs[n] * pow(z, n - 3);

        return self.T0 * (A / B);


    def zonal_winds(self, z):
        tropo_jet_arg = (z - self.tropo_jet_alt) / self.tropo_jet_width
        strat_jet_arg = (z - self.strat_jet_alt) / self.strat_jet_width
    
        result = self.tropo_jet_speed * np.exp(-tropo_jet_arg**2) * np.sin(np.radians(self.tropo_jet_az))
        result = result + self.strat_jet_speed * np.exp(-strat_jet_arg**2) * np.sin(np.radians(self.strat_jet_az))
        return result

    def merid_winds(self, z):
        tropo_jet_arg = (z - self.tropo_jet_alt) / self.tropo_jet_width
        strat_jet_arg = (z - self.strat_jet_alt) / self.strat_jet_width

        result = self.tropo_jet_speed * np.exp(-tropo_jet_arg**2) * np.cos(np.radians(self.tropo_jet_az))
        result = result + self.strat_jet_speed * np.exp(-strat_jet_arg**2) * np.cos(np.radians(self.strat_jet_az))
        return result
    
    def write_profile(self, file_name, z_rng = [0, 160], dz=0.1):
        file_out = open(file_name, 'w')
        for z_val in np.arange(z_rng[0], z_rng[1], dz):
            print(z_val, self.T(z_val), self.zonal_winds(z_val), self.merid_winds(z_val), self.rho(z_val), self.rho(z_val) * self.R * self.T(z_val) * 10.0, file=file_out)
        file_out.close()


class AtmosphereState(object):
    
    z_min = 0.0
    z_max = 120.0
    
    z_eofs = []
    T_eofs = []
    u_eofs = []
    v_eofs = []
    d_eofs = []
    p_eofs = []
    
    def __init__(self, estimate_profile, eof_path, profile_format="zTuvdp"):
        print("Loading estimated atmospheric state from " + estimate_profile + " using format '" + profile_format + "' and interpolating...")
        if profile_format == "zTuvdp":
            z_vals, T_vals, u_vals, v_vals, d_vals, p_vals = np.loadtxt(estimate_profile, unpack=True)
        elif profile_format == "zuvwTdp":
            z_vals, u_vals, v_vals, w_vals, T_vals, d_vals, p_vals = np.loadtxt(estimate_profile, unpack=True)
        else:
            print("Error!  Bad profile format provided.")

        self.z_est = z_vals
        self.T_est = interp1d(z_vals, T_vals, kind='cubic')
        self.u_est = interp1d(z_vals, u_vals, kind='cubic')
        self.v_est = interp1d(z_vals, v_vals, kind='cubic')
        self.d_est = interp1d(z_vals, d_vals, kind='cubic')
        self.p_est = interp1d(z_vals, p_vals, kind='cubic')
        
        print('\t' + "Loading EOFs from " + eof_path + " and interpolating...")
        T_temp = np.loadtxt(eof_path + "-temperature.eofs")
        u_temp = np.loadtxt(eof_path + "-zonal_winds.eofs")
        v_temp = np.loadtxt(eof_path + "-merid_winds.eofs")
        d_temp = np.loadtxt(eof_path + "-density.eofs")
        p_temp = np.loadtxt(eof_path + "-pressure.eofs")

        self.z_eofs = T_temp[:, 0]
        self.T_eofs = []
        self.u_eofs = []
        self.v_eofs = []
        self.d_eofs = []
        self.p_eofs = []

        for m in range(T_temp.shape[1] - 1):
            self.T_eofs += [interp1d(T_temp[:, 0], T_temp[:, m + 1], kind='cubic')]
            self.u_eofs += [interp1d(u_temp[:, 0], u_temp[:, m + 1], kind='cubic')]
            self.v_eofs += [interp1d(v_temp[:, 0], v_temp[:, m + 1], kind='cubic')]
            self.d_eofs += [interp1d(d_temp[:, 0], d_temp[:, m + 1], kind='cubic')]
            self.p_eofs += [interp1d(p_temp[:, 0], p_temp[:, m + 1], kind='cubic')]

        self.z_min = max(z_vals[0],  T_temp[0][0])
        self.z_max = min(z_vals[-1], T_temp[-1][0])
        print('\t' + "Setting alititude limits to " + str(self.z_min) + " - " + str(self.z_max) + "km...")


    def perturb_atmosphere(self, coeffs, file_name, dz=0.2):
        print("Writing perturbed atmosphere state into " + file_name + "...")
        file_out = open(file_name, 'w')
        
        for z in np.arange(self.z_min, self.z_max, dz):
            T_val = self.T_est(z)
            u_val = self.u_est(z)
            v_val = self.v_est(z)
            d_val = np.log(self.d_est(z))
            p_val = np.log(self.p_est(z))
            
            for coeff in coeffs:
                T_val += coeff[1] * self.T_eofs[coeff[0]](z)
                u_val += coeff[1] * self.u_eofs[coeff[0]](z)
                v_val += coeff[1] * self.v_eofs[coeff[0]](z)
                d_val += coeff[1] * self.d_eofs[coeff[0]](z)
                p_val += coeff[1] * self.p_eofs[coeff[0]](z)
            
            print(z, T_val, u_val, v_val, np.exp(d_val), np.exp(p_val), file=file_out)

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

