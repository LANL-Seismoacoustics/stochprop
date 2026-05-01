# gravity_waves.py
#
# Methods to compute the gravity wave spectra that perturb a G2S specification
# based on the publication by Drob et al. (2013).  Additional insights from the
# Lalande & Waxler (2016) analysis and Fritts & Alexander (2003) manuscripts on
# atmospheric gravity waves.  The source and saturation spectra are from the 
# Warner & McIntyre (1996) work.
#  
# Drob, D. P., Broutman, D., Hedlin, M. A., Winslow, N. W., & Gibson, R. G. (2013). 
# A method for specifying atmospheric gravity wavefields for long‐range infrasound 
# propagation calculations. Journal of Geophysical Research: Atmospheres, 118(10), 
# 3933-3943.
#
# Lalande, J. M., & Waxler, R. (2016). The interaction between infrasonic waves
# and gravity wave perturbations: Application to observations using UTTR rocket 
# motor fuel elimination events. Journal of Geophysical Research: Atmospheres, 
# 121(10), 5585-5600.
#
# Fritts, D. C., & Alexander, M. J. (2003). Gravity wave dynamics and effects in 
# the middle atmosphere. Reviews of geophysics, 41(1).
#
# Warner, C. D., & McIntyre, M. E. (1996). On the propagation and dissipation of
# gravity wave spectra through a realistic middle atmosphere. Journal of 
# Atmospheric Sciences, 53(22), 3213-3235.
#
# Philip Blom (pblom@lanl.gov)

import os 
import sys
import time
import datetime

import numpy as np

from multiprocess import Pool
from importlib.metadata import version

import matplotlib.pyplot as plt 

try:
    from scipy.integrate import simpson
except:
    from scipy.integrate import simps as simpson

from scipy.interpolate import interp1d
from scipy.special import airy

from . import utils as sp_utils

# Progress bar methods
def prog_prep(bar_length):
    sys.stdout.write("[%s]" % (" " * bar_length))
    sys.stdout.flush()

    sys.stdout.write("\b" * (bar_length + 1))
    sys.stdout.flush()

def prog_increment(n=1):
    for j in range(n):
        sys.stdout.write(">")
        sys.stdout.flush()
        time.sleep(0.01)

def prog_close():
    sys.stdout.write("\n")
    sys.stdout.flush()

def prog_set_step(n, N, bar_length):
    return int(np.floor((float(bar_length) * (n + 1)) / N) - np.floor((float(bar_length) * n) / N))



def single_fourier_component(k, l, om_intr, atmo_info, src_index, om_min, t0_evol, figure_out=None, prog_step=0):
    """
        Compute the vertical structure of a specific Fourier component, :math:`\hat{w} \left( k, l, \omega, z \right), 
        by first identifying critical layers and turning heights then using the appropriate solution form (free or 
        trapped solution) to evaluate the component.

        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om: float
            Absolute frequency (relative to the ground) [Hz]
        atmo_specification: string
            Atmospheric specification file path
        t0: float
            Reference time for gravity wave propagation (typically 4 - 6 hours)
        src_index: int
            Index of the source height within the atmo_info z values
        m_star: float
            Source parameter m_* (default value, :math:`\frac{2 \pi}{2.5} \text{ km}^{-1}` is for 20 km altitude source)
        om_min: float
            Minimum absolute frequency used in analysis
        k_max: float
            Maximum horzintal wavenumber value used in 1 grid dimension
        random_phase: bool
            Flag to randomize initial phase of freely propagating solution
        figure_out: string
            Path to output figure showing component charcterisitcs 
        prop_step: int
            Progress bar increment

        Returns
        -------
        u_spec: 1darray
            Zonal wind perturbation spectrum, \hat{u}(k, l, z, omega)
        v_spec: 1darray
            Meridional wind perturbation spectrum, \hat{v}(k, l, z, omega)
        w_spec: 1darray
            Vertical wind perturbation spectrum, \hat{w}(k, l, z, omega)
        eta_spec: 1darray
            Displacement spectrum used to compute temperature and pressure perturbations 

    """
    # extract atmospheric information
    z0 = atmo_info[0]
    T0 = atmo_info[1]
    d0 = atmo_info[2]
    H = atmo_info[3]
    N = atmo_info[4]
    u0 = atmo_info[5]
    v0 = atmo_info[6]

    # Compute the horizontal wavenumber 
    kh = np.sqrt(k**2 + l**2)

    m_sqr_vals = (kh / om_intr)**2 * (N**2 - om_intr**2) + 1.0 / (4.0 * H**2)
    m_vals = np.sqrt(abs(m_sqr_vals))

    src_cg = (m_vals[src_index] * kh * N[src_index]) / (kh**2 + m_sqr_vals[src_index] + 1.0 / (4.0 * H[src_index]**2))**(3.0 / 2.0)

    # original cg and N / m limits
    # cg_check = abs(src_cg) > 2.0e-4
    # Nm_check = abs(np.max(N / m_vals)) < 0.4

    cg_lim = 1.0e-4
    Nm_lim = 0.6

    cg_check = abs(src_cg) > cg_lim
    Nm_check = abs(np.max(N / m_vals)) < Nm_lim

    if cg_check and Nm_check:
        # set source strength
        m_star = (2.0 * np.pi) / 2.5
        src_Omega = om_min**(2.0 / 3.0) / (1.0 - (om_min / N[src_index])**(2.0 / 3.0))
        w0 = np.sqrt(abs(2.7e-2 * (m_sqr_vals[src_index] / (m_star**4 + m_sqr_vals[src_index]**2)) * (om_intr**(1.0 / 3.0) / kh**2) * src_Omega))

        # free propagating phase integration
        w_phase_vals = np.zeros_like(z0)
        w_phase_vals[:src_index + 1] = np.array([m_vals[src_index] * (z0[zj] - z0[src_index]) for zj in range(src_index + 1)])
        w_phase_vals[src_index + 1:] = np.array([simpson(m_vals[src_index:zj], z0[src_index:zj]) for zj in range(src_index + 1, len(z0))])

        # losses above 100 km
        # original viscosity scaling is 1.0e-10
        w_losses = np.zeros_like(z0)

        env = 1.0 / (1.0 + np.exp(-(z0 - 100.0) / 7.5))
        visc = 3.58e-7 * (T0**0.69 / d0) * 1.0e-10
        m_imag = (visc * abs(m_sqr_vals)**(3.0 / 2.0) / om_intr) * env

        damping_index = np.argmin(abs(z0 - 75.0))
        w_losses[damping_index + 1:] = np.array([simpson(m_imag[damping_index:zj], z0[damping_index:zj]) for zj in range(damping_index + 1, len(z0))])

        # combine into the vertical spectra
        w_spec = np.zeros_like(z0, dtype=complex)
        d0_m_ratios = np.sqrt((d0[src_index] / d0) * (m_vals[src_index] / m_vals))
        w_spec = w0 * d0_m_ratios * np.exp(-w_losses) * (np.cos(w_phase_vals) - 1.0j * np.sin(w_phase_vals))

        # apply saturation threshold
        # Original coefficient is 2.7e-2
        w_sat = np.sqrt(2.7e-2 * om_intr**(1.0 / 3.0) / (m_sqr_vals * kh**2)) 

        sat_mask = abs(w_spec) > w_sat 
        w_spec[sat_mask] = w_spec[sat_mask] * (w_sat[sat_mask] / abs(w_spec[sat_mask]))

        # apply phasing for field evolution time (convert u0, v0 to km/s)
        t0_phase = (k * u0 + l * v0) * 1.0e-3 * t0_evol
        w_spec = w_spec * (np.cos(t0_phase) - 1.0j * np.sin(t0_phase))
    else:
        w_spec = np.zeros_like(z0, dtype=complex)

    # project onto zonal/meridional winds and temperature
    u_spec = - w_spec * (k / kh**2) * m_vals
    v_spec = - w_spec * (l / kh**2) * m_vals
    eta_spec = (-1.0j / om_intr) * w_spec * (np.gradient(T0) / np.gradient(z0))
    sp_utils.prog_increment(prog_step)

    if figure_out:
        cg_check = abs(src_cg) > cg_lim
        Nm_check = abs(np.max(N / m_vals)) < Nm_lim

        f, a = plt.subplots(1, 2, sharey=True)
        a[0].set_ylabel("Altitude [km]")
        a[0].set_xlabel("m^2")
        a[1].set_xlabel("spectral component")

        a[0].plot(m_sqr_vals, z0, '-k')
        a[1].plot(np.real(u_spec) * kh, z0, '-b')
        a[1].plot(np.real(v_spec) * kh, z0, '-r')

        a[0].set_title("(k_perp, om_intr) = (" + str(np.round(kh ,4)) + ", " + str(np.round(om_intr, 4)) + '\n' + "N/m max:" + str(np.round(abs(np.max(N / m_vals)), 4))
                        + "(" + str(Nm_check) + ")" + "), src_cg: "+ str(np.round(src_cg, 4)) + "(" + str(cg_check) + ")")
    
        plt.savefig(figure_out + "-" + str(k) + "-" + str(l) + "-" + str(om_intr) + ".png", dpi=250)
        # plt.show()
        plt.close()

        '''
        # Plot cleaner looking panel (for ASA meeting)
        plt.figure(2, figsize=(5, 8))
        plt.xlabel("spectral component")
        plt.ylabel("Altitude [km]")

        plt.plot(np.real(u_spec), z0, '-b')
        plt.plot(np.real(v_spec), z0, '-r')

        plt.title("k_perp = " + str(np.round(kh, 3)) + " km^{-1}, om_intr = " + str(np.round(om_intr, 3)) + " Hz")

        plt.savefig(figure_out + "_2-" + str(k) + "-" + str(l) + "-" + str(om_intr) + ".png", dpi=250)
        plt.close()
        '''

    return [u_spec, v_spec, w_spec, eta_spec]



def single_fourier_component_wrapper(args):
    return single_fourier_component(*args)


def _perturb_header_txt(k_max, fourier_cnt, n, smpl_cnt, src_lat, t0_evol, ref_header):
    result = "# Data Source: stochprop v" + version("stochprop")
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Gravity Wave Perturbation"
    result = result + '\n' + "# Maximum Wavenumber [km^{-1]}] = " + str(k_max)
    result = result + '\n' + "# Fourier Component Count = " + str(fourier_cnt)
    result = result + '\n' + "# WG-Wind Evolution Time [s] = " + str(t0_evol)
    result = result + '\n' + "# Sample: " + str(n) + "/" + str(smpl_cnt)

    if src_lat is not None:
        result = result + '\n' + "# Overwritten Reference Latitude = " + str(src_lat)

    for line in ref_header:
        result = result + '\n' + line

    return result


def perturb_atmo(atmo_file, output_path, k_max, fourier_cnt, smpl_cnt, t0_evol, gw_scaling=1.0, src_lat=None, debug_fig_out=None, pool=None):

    """
        Use gravity waves to perturb a specified profile using some of the methods in Drob et al. (2013)

        Parameters
        ----------
        atmo_file: string
            Path and name of the atmosphere file to be used as the reference
        output_path: string
            Path where output will be stored
        k_max: float 
            Maximum horizontal wavenumber for the gravity wave spectrum (0.25 - 0.5 km^{-1})
        fourier_cnt: int
            Number of Fourier coefficient sets to consider in approximating the integral            
        smpl_cnt: int
            Number of perturbed atmospheric samples to generate
        t0_evol: float
            Reference time for wind-gravity wave interaction (0 - 1200 seconds)
        gw_scaling: float
            Scaling factor to increase or decrease perturbations
        src_lat: float
            Custom source latitude for Coriolis force calculation
        debug_fig_out: string
            Path for debugging figures to be written
        pool: multiprocessing Pool
            Pool for multi-threading calculation

    """

    ref_header = ["# Reference Specification = " + atmo_file + " (cwd: " + os.getcwd() + ")"]
    with open(atmo_file, 'r') as file:
        for line in file:
            if "#" in line:
                if any(s in line for s in ["Data Source", "Model Time", "Location"]):
                    ref_header = ref_header + ["# Reference " + line[2:].replace('\n','')]
                else:
                    ref_header = ref_header + [line.replace('\n','')]
            else:
                break

    if src_lat is not None:
        ref_lat = src_lat
    else:
        try:
            for line in ref_header:
                if "Location" in line:
                    ref_lat = float(line.split(' ')[-4][:-1])
        except:
            ref_lat = 40.0

    # re-sample fields to higher resolution
    z0, T0, u0, v0, d0, p0 = np.loadtxt(atmo_file, unpack=True)

    # define dz needed to resolve phasing from wind interactions
    # Want 8 points/cycle at highest wind gradient
    dudz = np.gradient(u0) / np.gradient(z0)
    dvdz = np.gradient(v0) / np.gradient(z0)
    dPhi_dz = (k_max * t0_evol) * np.max(np.sqrt(dudz**2 + dvdz**2)) * 1.0e-3
    dz = min(abs(z0[1] - z0[0]) / 2.0, (np.pi / 4.0) / dPhi_dz)

    z = np.arange(z0[0], z0[-1], dz)
    T = interp1d(z0, T0, kind='cubic')(z)
    u = interp1d(z0, u0, kind='cubic')(z)
    v = interp1d(z0, v0, kind='cubic')(z)
    d = interp1d(z0, d0, kind='cubic')(z)
    p = interp1d(z0, p0, kind='cubic')(z)

    src_index = np.argmin(abs(z - 20.0))

    H = - d / (np.gradient(d) / np.gradient(z))
    N = np.sqrt(9.8e-3 / H)

    om_min = 2.0 * 7.29e-5 * np.sin(np.radians(ref_lat))
    om_max = np.max(N) / np.sqrt(5)

    # randomize Fourier components
    # sampling on circle requires square root scaling to get distribution correct
    rn_gen = np.random.default_rng()

    k_perp = k_max * np.sqrt(rn_gen.uniform(low=0.0, high=1.0, size=fourier_cnt))
    k_angle = rn_gen.uniform(low=0.0, high=2.0 * np.pi, size=fourier_cnt)   
    k, l = k_perp * np.cos(k_angle), k_perp * np.sin(k_angle)

    om = rn_gen.uniform(low=om_min, high=om_max, size=fourier_cnt)

    print('Computing Fourier components...\n\t', end='')
    sp_utils.prog_prep(50)
    if pool is not None:
        args = [[k[n], l[n], om[n], [z, T, d, H, N, u, v], src_index, om_min, t0_evol, debug_fig_out, sp_utils.prog_set_step(n, fourier_cnt, 50)] for n in range(fourier_cnt)]
        results = np.array(pool.map(single_fourier_component_wrapper, args))
    else:
        results = np.array([single_fourier_component(k[n], l[n], om[n], [z, T, d, H, N, u, v], src_index, om_min, t0_evol, figure_out=debug_fig_out, prog_step=sp_utils.prog_set_step(n, fourier_cnt, 50)) for n in range(fourier_cnt)])
    sp_utils.prog_close()

    print('\nApplying perturbations to reference atmosphere...\n\t', end='')
    u_spec = results[:, 0, :]
    v_spec = results[:, 1, :]
    
    # scale for cylindrical volume, convert km to meters, account for nulled components
    scaling = ((2.0 * np.pi) * k_max * (om_max - om_min)) * 1.0e3 * (fourier_cnt / np.count_nonzero(np.count_nonzero(u_spec, axis=1)))
    u_spec = u_spec * np.outer(k_perp, np.ones_like(z)) * scaling
    v_spec = v_spec * np.outer(k_perp, np.ones_like(z)) * scaling
    
    sp_utils.prog_prep(50)
    for m in range(smpl_cnt):
        # randomly phase and sum together
        rndm_phasing = np.outer(rn_gen.uniform(low=0.0, high=2.0 * np.pi, size=fourier_cnt), np.ones_like(z))
        u_perturb = u_spec * (np.cos(rndm_phasing) - 1.0j * np.sin(rndm_phasing))
        v_perturb = v_spec * (np.cos(rndm_phasing) - 1.0j * np.sin(rndm_phasing))

        u_out = u + gw_scaling * np.sqrt(fourier_cnt) * np.real(np.mean(u_perturb, axis=0))
        v_out = v + gw_scaling * np.sqrt(fourier_cnt) * np.real(np.mean(v_perturb, axis=0))

        np.savetxt(output_path + "-" + str(m) + ".met", np.vstack((z, T, u_out, v_out, d, p)).T, 
            header=_perturb_header_txt(k_max, fourier_cnt, m, smpl_cnt, src_lat, t0_evol, ref_header), comments='')
        sp_utils.prog_increment(sp_utils.prog_set_step(m, smpl_cnt, 50))
    sp_utils.prog_close()








