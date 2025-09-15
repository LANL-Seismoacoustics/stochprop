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



def single_fourier_component(k, l, om_intr, atmo_info, src_index, om_min, figure_out=None, prog_step=0):
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

    # Compute the horizontal wavenumber 
    kh = np.sqrt(k**2 + l**2)

    m_sqr_vals = (kh / om_intr)**2 * (N**2 - om_intr**2) + 1.0 / (4.0 * H**2)
    m_vals = np.sqrt(abs(m_sqr_vals))

    src_cg = (m_vals[src_index] * kh * N[src_index]) / (kh**2 + m_sqr_vals[src_index] + 1.0 / (4.0 * H[src_index]**2))**(3.0 / 2.0)
    cg_check = abs(src_cg) > 2.0e-4
    Nm_check = abs(np.max(N / m_vals)) < 0.4

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
        w_losses = np.zeros_like(z0)

        env = 1.0 / (1.0 + np.exp(-(z0 - 100.0) / 2.5))
        visc = 3.58e-7 * (T0**0.69 / d0) * 1.0e-10
        m_imag = (visc * abs(m_sqr_vals)**(3.0 / 2.0) / om_intr) * env

        damping_index = np.argmin(abs(z0 - 80.0))
        w_losses[damping_index + 1:] = np.array([simpson(m_imag[damping_index:zj], z0[damping_index:zj]) for zj in range(damping_index + 1, len(z0))])

        # combine into the vertical spectra
        w_spec = np.zeros_like(z0, dtype=complex)
        d0_m_ratios = np.sqrt((d0[src_index] / d0) * (m_vals[src_index] / m_vals))
        w_spec = w0 * d0_m_ratios * np.exp(-w_losses) * (np.cos(w_phase_vals) - 1.0j * np.sin(w_phase_vals))

        # apply saturation threshold
        w_sat = np.sqrt(2.7e-2 * om_intr**(1.0 / 3.0) / (m_sqr_vals * kh**2))
        sat_mask = abs(w_spec) > w_sat 
        w_spec[sat_mask] = w_spec[sat_mask] * (w_sat[sat_mask] / abs(w_spec[sat_mask]))

        if max(abs(w_spec)) > 25.0:
            w_spec = np.zeros_like(z0, dtype=complex)

    else:
        w_spec = np.zeros_like(z0, dtype=complex)

    # project onto zonal/meridional winds and temperature
    u_spec = - w_spec * (k / kh**2) * m_vals
    v_spec = - w_spec * (l / kh**2) * m_vals
    eta_spec = (-1.0j / om_intr) * w_spec * (np.gradient(T0) / np.gradient(z0))

    prog_increment(prog_step)

    if figure_out:
        cg_check = abs(src_cg) > 1.0e-3
        Nm_check = abs(np.max(N / m_vals)) < 0.4

        f, a = plt.subplots(1, 2, sharey=True)
        a[0].set_ylabel("Altitue [km]")
        a[0].set_xlabel("m^2")
        a[1].set_xlabel("spectral component")

        a[0].plot(m_sqr_vals, z0, '-k')
        a[1].plot(np.real(w_spec), z0, '-b')
        a[1].plot(np.imag(w_spec), z0, '-r')

        a[0].set_title("(k_perp, om_intr) = (" + str(np.round(kh ,4)) + ", " + str(np.round(om_intr, 4)) + '\n' + "N/m max:" + str(np.round(abs(np.max(N / m_vals)), 4)) + "(" + str(Nm_check) + ")" + "), src_cg: "+ str(np.round(src_cg, 4)) + "(" + str(cg_check) + ")")
    
        plt.show()
        plt.savefig(figure_out + "-" + str(k) + "-" + str(l) + "-" + str(om_intr) + ".png", dpi=250)
        plt.close()
    
    return [u_spec, v_spec, w_spec, eta_spec]



def single_fourier_component_wrapper(args):
    return single_fourier_component(*args)


def _perturb_header_txt(prof_path, src_lat, k_max, fourier_cnt, n, smpl_cnt):
    result = "# Data Source: stochprop v" + version("stochprop")
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Gravity Wave Perturbation"
    result = result + '\n' + "# Reference Specification = " + prof_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# Source Latitude = " + str(src_lat)
    result = result + '\n' + "# Maximum Wavenumber [km^{-1]}] = " + str(k_max)
    result = result + '\n' + "# Fourier Component Count = " + str(fourier_cnt)
    result = result + '\n' + "# Sample: " + str(n) + "/" + str(smpl_cnt)
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


def perturb_atmo(atmo_file, sample_path, k_max, fourier_cnt, smpl_cnt, src_lat=None, debug_fig_out=None, pool=None):

    """
        Use gravity waves to perturb a specified profile using some of the methods in Drob et al. (2013)

        Parameters
        ----------
        atmo_spec: string
            Path and name of the specification to be used as the reference
        output_path: string
            Path where output will be stored
        sample_cnt: int
            Number of perturbed atmospheric samples to generate
        t0: float
            Reference time for gravity wave propagation (typically 4 - 6 hours)
        dx: float
            Horizontal wavenumber resolution [km]
        dz: float
            Vertical resolution for integration steps [km]
        Nk: int
            Horizontal wavenumber grid dimensions (Nk x Nk)
        N_om: int
            Frequency resolution (typically 5)
        ref_lat: float
            Reference latitude used to define the Coriolis frequency used as the minimum frequency
        random_phase: boolean
            Controls inclusion of random initial phase shifts
        env_below: boolean
            Controls whether perturbations below the source height are included
        cpu_cnt: int
            Number of CPUs to use for parallel computation of Fourier components (defaults to None)

    """

    ref_lat = 40.0
    if src_lat is not None:
        ref_lat = src_lat
    else:
        try:
            temp = open(atmo_file, 'r')
            for line in temp:
                if "Location" in line:
                    ref_lat = float(line.split(' ')[-4][:-1])
        except:
            ref_lat = 40.0

    z0, T0, u0, v0, d0, p0 = np.loadtxt(atmo_file, unpack=True)

    z_src = 20.0
    src_index = np.argmin(abs(z0 - z_src))

    H = - d0 / (np.gradient(d0) / np.gradient(z0))
    N = np.sqrt(9.8e-3 / H)

    om_min = 2.0 * 7.29e-5 * np.sin(np.radians(ref_lat))
    om_max = np.max(N) / np.sqrt(5)

    # randomize Fourier components
    k_perp = k_max * np.random.random(fourier_cnt)
    k_angle = (2.0 * np.pi) * np.random.random(fourier_cnt)
    k, l = k_perp * np.cos(k_angle), k_perp * np.sin(k_angle)

    om = om_min + (om_max - om_min) * np.random.random(fourier_cnt)

    print('Computing Fourier components...\n\t', end='')
    prog_prep(50)
    if pool is not None:
        args = [[k[n], l[n], om[n], [z0, T0, d0, H, N], src_index, om_min, debug_fig_out, prog_set_step(n, fourier_cnt, 50)] for n in range(fourier_cnt)]
        results = np.array(pool.map(single_fourier_component_wrapper, args))
    else:
        results = np.array([single_fourier_component(k[n], l[n], om[n], [z0, T0, d0, H, N], src_index, om_min, figure_out=debug_fig_out, prog_step=prog_set_step(n, fourier_cnt, 50)) for n in range(fourier_cnt)])
    prog_close()

    print('\nApplying perturbations to reference atmosphere...\n\t', end='')
    u_spec = results[:, 0, :]
    v_spec = results[:, 1, :]

    scaling = (2.0 * np.pi) * k_max**2 * (om_max - om_min)
    scaling = scaling * (fourier_cnt / np.count_nonzero(np.count_nonzero(u_spec, axis=1)))

    prog_prep(50)
    for m in range(smpl_cnt):
        rndm_phasing = (2.0 * np.pi) * np.outer(np.random.random(fourier_cnt), np.ones_like(z0))

        u_perturb = np.mean(u_spec * np.exp(-1.0j * rndm_phasing), axis=0) * scaling
        v_perturb = np.mean(v_spec * np.exp(-1.0j * rndm_phasing), axis=0) * scaling

        u_vals = u0 + 2.0 * np.real(u_perturb) * 1.0e3
        v_vals = v0 + 2.0 * np.real(v_perturb) * 1.0e3

        np.savetxt(sample_path + "-" + str(m) + ".met", np.vstack((z0, T0, u_vals, v_vals, d0, p0)).T, 
            header=_perturb_header_txt(atmo_file, src_lat, k_max, fourier_cnt, m, smpl_cnt), comments='')
        prog_increment(prog_set_step(m, smpl_cnt, 50))
    prog_close()








