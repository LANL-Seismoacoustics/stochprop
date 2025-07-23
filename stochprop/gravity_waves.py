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

from scipy.integrate import simps
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


# Gravity wave methods
def BV_freq(H):
    """
        Compute the Brunt-Vaisala frequency defined as :math:`N = \sqrt{\frac{g}{H}}` where 
        :math:`H = \frac{rho0}{\frac{\partial \rho_0}{\partial z}` is the density scale height
    
        Parameters
        ----------
        H: float
            Scale height, :math:`H = \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}`

        Returns:
        f_BV: float
            Brunt-Vaisalla (bouyancy) frequency, :math:`f_BV = sqrt{\frac{g}{H}}`

    """
    return np.sqrt(9.8e-3 / H)


def m_sqr(k, l, om_intr, H):
    """
        Compute the vertical wavenumber dispersion relation for gravity wave propagation
        defined as :math:`m^2 = \frac{k_h^2}{\hat{\omega}^2} \left( N^2 - \hat{\omega}^2 \right) + \frac{1}{4 H^2}`

        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om_intr: float
            Intrinsic frequency (relative to winds), defined as :math:`\hat{\omega} = \omega - k u_0 - l v_0`
        H: float
            Scale height, :math:`H = \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}`

        Returns:
        m_sqr: float
            Vertical wave number squared, :math:`m^2 = \frac{k_h^2}{\hat{\omega}^2 \left( N^2 - \hat{\omega}^2 \right) + \frac{1}{4 H^2}}`


    """
    return ((k**2 + l**2) / om_intr**2) * (BV_freq(H)**2 - om_intr**2) + 1.0 / (4.0 * H**2)


def cg(k, l, om_intr, H):
    """
        Compute the vertical group velocity for gravity wave propagation 
        as :math:`cg = \frac{\partial \hat{omega}}{\partial m} = \frac{m k_h N}{ \left(k_h^2 + m^2 + \frac{1}{4H^2 \right)^{\frac{3}{2}}}`

        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om_intr: float
            Intrinsic frequency (relative to winds), defined as :math:`\hat{\omega} = \omega - k u_0 - l v_0`
        H: float
            Scale height, :math:`H = \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}`

        Returns:
        c_g: float
            Vertical group velocity of gravity waves


    """
    m_sqr_val = abs(m_sqr(k, l, om_intr, H))
    kh = np.sqrt(k**2 + l**2)

    return (np.sqrt(m_sqr_val) * kh * BV_freq(H)) / (kh**2 + m_sqr_val + 1.0 / (4.0 * H**2))**(3.0 / 2.0)


def m_imag(k, l, om_intr, z, H, T0, d0):
    """
        Compute the imaginary wave number component to add attenuation effects
        The imaginary component is defined as :math:`m_\text{im} = -\nu \frac{m^3}{\hat{\omega}}`
        where the viscosity is :math:`\nu = 3.563 \times 10^{-7} \frac{T_0 \left( z \right)}{\rho_0 \left( z \right)}`

        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om_intr: float
            Intrinsic frequency (relative to winds), defined as :math:`\hat{\omega} = \omega - k u_0 - l v_0`
        z: float
            Absolute height (used for turning attenuation "off" below 100 km)
        H: float
            Scale height, :math:`H = \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}`
        T0: float
            Ambient temperature in the atmosphere
        d0: float
            Ambient density in the atmosphere

        Returns:
        m_i: float
            Imaginary component of the wavenumber used for damping above 100 km (note: 100 km limit is applied elsewhere)
    """
    env = 1.0 / (1.0 + np.exp(-(z - 100.0) / 2.5))
    visc = 3.563e-7 * (T0**0.69 / d0)
    return visc * (abs(m_sqr(k, l, om_intr, H))**(3.0 / 2.0) / om_intr) * env * 1.0e-10



def single_fourier_component(k, l, om_intr, atmo_info, t0, src_index, m_star, om_min, k_max, random_phase=False, figure_out=None, prog_step=0):
    """
        Compute the vertical structure of a specific Fourier component, :math:`\hat{w} \left( k, l, \omega, z \right), 
        by first identifying critical layers and turning heights then using the appropriate solution form (free or 
        trapped solution) to evalute the component.

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
    z_vals = atmo_info[0]
    H_vals = atmo_info[1]
    u0_vals = atmo_info[2]
    v0_vals = atmo_info[3]
    T0_vals = atmo_info[4]
    d0_vals = atmo_info[5]

    dz = (z_vals[1] - z_vals[0])

    # Define combined horizontal wavenumber, intrinsic frequency 
    # values, and vertical wavenumber values (squared) 
    kh = np.sqrt(k**2 + l**2)
    m_sqr_vals = m_sqr(k, l, om_intr, H_vals)

    if random_phase:
        ph0 = np.random.random(1)[0] * (2.0 * np.pi)
    else:
        ph0 = 0.0

    prog_increment(prog_step)
    if kh < 1.0e-6 or kh > k_max:
        return [np.zeros_like(z_vals, dtype=complex)] * 4
    else:
        # check for turning height
        if np.all(m_sqr_vals > 0.0):
            turn_ht_index = len(z_vals) + 1
        else:
            jz = np.where(m_sqr_vals < 0.0)[0][0]
            if jz > src_index + 2:
                turn_ht_index = jz - 1
                refl_phase = simps(np.sqrt(m_sqr_vals[:turn_ht_index]), z_vals[:turn_ht_index])
                refl_time = simps(1.0 / abs(cg(k, l, om_intr, H_vals[:turn_ht_index])), z_vals[:turn_ht_index])
                refl_loss = simps(m_imag(k, l, om_intr, z_vals[:turn_ht_index], H_vals[:turn_ht_index], T0_vals[:turn_ht_index], d0_vals[:turn_ht_index]), z_vals[:turn_ht_index])
            else:
                return [np.zeros_like(z_vals, dtype=complex)] * 4
        
        # check velocity condition (N/m < 90 m/s)
        if turn_ht_index > 0:
            velocity_check = np.max(BV_freq(H_vals[:turn_ht_index]) / np.sqrt(m_sqr(k, l, om_intr, H_vals[:turn_ht_index])))
            src_cg = cg(k, l, om_intr, H_vals[src_index])

            if abs(velocity_check) < 0.09 and src_cg > 5.0e-4:
                m_sqr_src = m_sqr(k, l, om_intr, H_vals[src_index])
                w0 = 2.7e-2 * (m_sqr_src / (m_star**4 + m_sqr_src**2)) * (om_intr**(1.0 / 3.0) / kh**2)
                w0 = w0 * om_min**(2.0 / 3.0) / (1.0 - (om_min / BV_freq(H_vals[0]))**(2.0 / 3.0))
                w0 = np.sqrt(abs(w0))        
            else:
                return [np.zeros_like(z_vals, dtype=complex)] * 4

        # prep numpy arrays for spectral info
        u_spec = np.zeros_like(z_vals, dtype=complex)
        v_spec = np.zeros_like(z_vals, dtype=complex)
        w_spec = np.zeros_like(z_vals, dtype=complex)
        eta_spec = np.zeros_like(z_vals, dtype=complex)

        # prep information for spectral calculation
        m_vals = np.sqrt(abs(m_sqr_vals))
        d0_m_ratios = np.array([np.sqrt((d0_vals[src_index] / d0_vals[zj]) * (m_vals[src_index] / m_vals[zj])) for zj in range(len(z_vals))])

        w_sat_vals = np.sqrt(2.7e-2 * om_intr**(1.0 / 3.0) / (m_sqr_vals * kh**2))
        
        prop_times = np.zeros_like(z_vals)
        prop_times[src_index + 1:] = np.array([simps(1.0 / abs(cg(k, l, om_intr, H_vals[src_index:zj])), z_vals[src_index:zj]) for zj in range(src_index + 1, len(z_vals))])
        if np.all(prop_times < t0):
            prop_ht_index = len(z_vals) + 1
        else:
            prop_ht_index = np.where(prop_times > t0)[0][0]

        w_losses = np.zeros_like(z_vals)
        w_losses[src_index + 1:] = np.array([simps(m_imag(k, l, om_intr, z_vals[src_index:zj], H_vals[src_index:zj], T0_vals[src_index:zj], d0_vals[src_index:zj]), z_vals[src_index:zj]) for zj in range(src_index + 1, len(z_vals))])

        if prop_ht_index <= turn_ht_index:
            # free propagating solution
            w_phase_vals = np.zeros_like(z_vals)
            w_phase_vals[:src_index + 1] = np.array([m_vals[src_index] * (z_vals[zj] - z_vals[src_index]) for zj in range(src_index + 1)])
            w_phase_vals[src_index + 1: prop_ht_index] = np.array([simps(m_vals[src_index:zj], z_vals[src_index:zj]) for zj in range(src_index + 1, min(len(z_vals), prop_ht_index))])

            w_spec[:prop_ht_index] = w0 * d0_m_ratios[:prop_ht_index] * (np.cos(ph0 + w_phase_vals[:prop_ht_index]) - 1.0j * np.sin(ph0 + w_phase_vals[:prop_ht_index]))
            w_spec[:prop_ht_index] = w_spec[:prop_ht_index] * np.exp(-w_losses[:prop_ht_index])

            u_spec[:prop_ht_index] = - w_spec[:prop_ht_index] * (k / kh**2) * m_vals[:prop_ht_index]
            v_spec[:prop_ht_index] = - w_spec[:prop_ht_index] * (l / kh**2) * m_vals[:prop_ht_index]
            eta_spec[:prop_ht_index] = -1.0j * w_spec[:prop_ht_index] / om_intr
        else:
            # trapped (Airy function) solution 
            # set source region below source altitude
            refl_cnt = int(np.floor(t0 / (2.0 * refl_time)))
            refl_ph_shift = np.exp(-refl_cnt * (2.0 * refl_loss))
            refl_ph_shift = refl_ph_shift * (1.0 + np.sum(np.array([np.exp(-1.0j * (N - 1) * (2.0 * refl_phase - np.pi / 2.0)) for N in range(2, refl_cnt)])))
            
            if np.all(m_sqr_vals[turn_ht_index:turn_ht_index + int(5.0 / dz)] < 0.0):
                airy_lim_index = turn_ht_index + int(5.0 / dz)
            else:
                airy_lim_index = turn_ht_index + (np.where(m_sqr_vals[turn_ht_index:turn_ht_index + int(5.0 / dz)] > 0.0)[0][0] - 1)

            airy_args = np.zeros_like(z_vals)
            airy_scaling0 = np.zeros_like(z_vals)
            airy_scaling1 = np.zeros_like(z_vals)

            m_integral_above = simps(m_vals[src_index:turn_ht_index], z_vals[src_index:turn_ht_index])
            airy_args[:src_index + 1] = np.array([-((3.0 / 2.0) * (m_integral_above + m_vals[src_index] * abs(z_vals[src_index] - z_vals[zj])))**(2.0 / 3.0) for zj in range(src_index + 1)])
            airy_args[src_index + 1:turn_ht_index - 1] = np.array([-((3.0 / 2.0) * simps(m_vals[zj:turn_ht_index], z_vals[zj:turn_ht_index]))**(2.0 / 3.0) for zj in range(src_index, turn_ht_index - 1)])
            airy_args[turn_ht_index + 1:airy_lim_index] = np.array([((3.0 / 2.0) * simps(m_vals[turn_ht_index:zj], z_vals[turn_ht_index:zj]))**(2.0 / 3.0) for zj in range(turn_ht_index + 1, airy_lim_index)])

            temp = (-airy_args[:airy_lim_index])**0.25 * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift
            airy_scaling0[:airy_lim_index] = temp * airy(airy_args[:airy_lim_index])[0]
            airy_scaling1[:airy_lim_index] = temp * airy(airy_args[:airy_lim_index])[1]

            u_spec[:airy_lim_index] = -2.0j * np.sqrt(np.pi) * w0 * (k / kh**2) * d0_m_ratios[:airy_lim_index] * airy_scaling1[:airy_lim_index]
            v_spec[:airy_lim_index] = -2.0j * np.sqrt(np.pi) * w0 * (l / kh**2) * d0_m_ratios[:airy_lim_index] * airy_scaling1[:airy_lim_index]
            w_spec[:airy_lim_index] = -2.0j * np.sqrt(np.pi) * w0 * d0_m_ratios[:airy_lim_index] * airy_scaling0[:airy_lim_index]
            eta_spec[:airy_lim_index] = -1.0j * w_spec[:airy_lim_index] / om_intr

        sat_mask = abs(w_spec) > w_sat_vals
        u_spec[sat_mask] = u_spec[sat_mask] * (w_sat_vals[sat_mask] / abs(w_spec[sat_mask]))
        v_spec[sat_mask] = v_spec[sat_mask] * (w_sat_vals[sat_mask] / abs(w_spec[sat_mask]))
        w_spec[sat_mask] = w_spec[sat_mask] * (w_sat_vals[sat_mask] / abs(w_spec[sat_mask]))

        wind_phasor = np.cos(-(k * u0_vals + l * v0_vals)) - 1.0j * np.sin(-(k * u0_vals + l * v0_vals))
        u_spec = u_spec * wind_phasor
        v_spec = v_spec * wind_phasor
        w_spec = w_spec * wind_phasor
        eta_spec = eta_spec * wind_phasor
        
        if figure_out:
            prop_times = np.array([simps(1.0 / abs(cg(k, l, om_intr, H_vals[src_index:zk])), z_vals[src_index:zk]) for zk in range(src_index + 1, len(z_vals))])

            f, a = plt.subplots(1, 3, sharey=True)
            a[0].set_ylabel("Altitue [km]")
            a[0].set_xlabel("prop. time [hrs]")
            a[1].set_xlabel("m^2")
            a[2].set_xlabel("spectral component")

            a[0].plot(prop_times / 3600.0, z_vals[src_index + 1:], '-k')
            a[0].axvline(t0 / 3600.0, color="k")
            a[1].semilogx(m_sqr_vals, z_vals, '-k')
            a[2].plot(np.real(w_spec), z_vals, '-b')
            a[2].plot(np.imag(w_spec), z_vals, '-r')
            for n in range(3):
                a[n].axhline(z_vals[min(turn_ht_index, len(z_vals) - 1)], color="g")
            a[1].set_title("(k, l, om_intr) = (" + str(k) + ", " + str(l) + ", " + str(om_intr) + ")" + '\n' + "Source cg: " + str(cg(k, l, om_intr, H_vals[src_index])))
        
            plt.savefig(figure_out + "-" + str(k) + "-" + str(l) + "-" + str(om_intr) + ".png", dpi=250)
            plt.close()
        
        return [u_spec, v_spec, w_spec, eta_spec]


def single_fourier_component_wrapper(args):
    return single_fourier_component(*args)


def perturbations(atmo_specification, t0=4.0 * 3600.0, dx=2.0, dz=0.2, Nk=128, N_om=5, ref_lat=40.0, random_phase=False, z_src=20.0, m_star=(2.0*np.pi)/2.5, figure_out=None, pool=None):
    """
        Loop over Fourier components :math:`\left(k, l, \omega \right)` and compute the spectral components for :math:`\hat{u} \left(k, l, \omega, z \right)`, 
        :math:`\hat{v}\left(k, l, \omega, z \right)`, and :math:`\hat{w} \left(k, l, \omega, z \right)`.  Once computed, apply inverse Fourier transforms to 
        obtain the space and time domain forms.

        Parameters
        ----------
        atmo_specification: string
            Atmospheric specification file path
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
        figure_out: string
            Option to output a figure with each component's structure (slows down calculations notably, useful for debugging)
        pool: multiprocessing.Pool
            Multprocessing option for parallel computation of Fourier components
                       
        Returns
        -------
        z_vals: 1darray
            Altitudes of output 
        du_vals: 3darray
            Zonal (E/W) wind perturbations, du(x, y, z, t0)
        dv_vals: 3darray
            Meridional (N/S) wind perturbations, dv(x, y, z, t0)
        dw_vals: 3darray
            Vertical wind perturbations, dw(x, y, z, t0)
            
    """

    print("Computing gravity wave perturbations for " + atmo_specification)

    # Set up atmosphere info
    z, Temp0, u0, v0, d0, _ = np.loadtxt(atmo_specification, unpack=True)

    T0_interp = interp1d(z, Temp0)
    u0_interp = interp1d(z, u0 * 1.0e-3)
    v0_interp = interp1d(z, v0 * 1.0e-3)
    d0_interp = interp1d(z, d0)

    d0_finite_diffs = np.array([(d0[min(n + 1, len(z) - 1)] - d0[max(n - 1, 0)]) / (z[min(n + 1, len(z) - 1)] - z[max(n - 1, 0)]) for n in range(len(z))])
    sc_ht = interp1d(z, -d0 / d0_finite_diffs)

    z_vals = np.arange(0.0, z[-1], dz)
    atmo_info = [z_vals, sc_ht(z_vals), u0_interp(z_vals), v0_interp(z_vals), T0_interp(z_vals), d0_interp(z_vals)]
    src_index = np.argmin(abs(z_vals - z_src))

    # Define spectral component grid
    om_min = 2.0 * 7.292e-5 * np.sin(np.radians(ref_lat))
    om_max = np.max(BV_freq(sc_ht(z))) / np.sqrt(5)

    k_vals = np.fft.fftfreq(Nk, d=1.0 / dx)
    l_vals = np.fft.fftfreq(Nk, d=1.0 / dx)
    intr_om_vals = np.linspace(om_min, om_max, N_om)

    # Define spectral information for each Fourier component
    u_spec = np.zeros((Nk, Nk, N_om, len(z_vals)), dtype=complex)
    v_spec = np.zeros_like(u_spec, dtype=complex)
    w_spec = np.zeros_like(u_spec, dtype=complex)
    eta_spec = np.zeros_like(u_spec, dtype=complex)
               
    # Compute each Fourier component
    print('\t' + "Integrating Fourier components to compute spectra...")
    print('\t' + "Progress: ", end='')
    prog_len = 50
    prog_prep(prog_len)
    if pool:
        args = []
        M = 0
        for nk, k in enumerate(k_vals):
            for nl, l in enumerate(l_vals):
                for n_om, om in enumerate(intr_om_vals):
                    step = prog_set_step(M, Nk * Nk * N_om, float(prog_len))
                    args = args + [[k, l, om, atmo_info, t0, src_index, m_star, intr_om_vals[0], max(abs(k_vals)), random_phase, figure_out, step]]
                    M = M + 1

        results = pool.map(single_fourier_component_wrapper, args)

        M = 0
        for nk in range(Nk):
            for nl in range(Nk):
                for n_om  in range(N_om):
                    u_spec[nk][nl][n_om] = results[M][0]
                    v_spec[nk][nl][n_om] = results[M][1]
                    w_spec[nk][nl][n_om] = results[M][2]
                    eta_spec[nk][nl][n_om] = results[M][3]
                    M = M + 1
    else:
        M = 0
        for nk, k in enumerate(k_vals):
            for nl, l in enumerate(l_vals):
                for n_om, om in enumerate(intr_om_vals):
                    step = prog_set_step(M, Nk* Nk * N_om, float(prog_len))
                    results = single_fourier_component(k, l, om, atmo_info, t0, src_index, m_star, intr_om_vals[0], max(abs(k_vals)), random_phase, figure_out, step)

                    u_spec[nk][nl][n_om] = results[0]
                    v_spec[nk][nl][n_om] = results[1]
                    w_spec[nk][nl][n_om] = results[2]
                    eta_spec[nk][nl][n_om] = results[3]

    prog_close()

    print('\t' + "Running inverse Fourier transforms to obtain space and time domain solutions.")

    dk = Nk / dx
    # dk = 1.0 / dx 

    du_vals = np.fft.ifft(u_spec, axis=0) * dk 
    du_vals = np.fft.ifft(du_vals, axis=1) * dk
    du_vals = simps(du_vals, intr_om_vals, axis=2)

    dv_vals = np.fft.ifft(v_spec, axis=0) * dk 
    dv_vals = np.fft.ifft(dv_vals, axis=1) * dk 
    dv_vals = simps(dv_vals, intr_om_vals, axis=2)

    dw_vals = np.fft.ifft(w_spec, axis=0) * dk 
    dw_vals = np.fft.ifft(dw_vals, axis=1) * dk 
    dw_vals = simps(dw_vals, intr_om_vals, axis=2)

    eta_vals = np.fft.ifft(eta_spec, axis=0) * dk 
    eta_vals = np.fft.ifft(eta_vals, axis=1) * dk 
    eta_vals = simps(eta_vals, intr_om_vals, axis=2)

    return z_vals, np.real(du_vals), np.real(dv_vals), np.real(dw_vals), np.real(eta_vals)


def _perturb_header_txt(prof_path, t0, dx, Nk, N_om, random_phase, z_src, m_star, n, prof_cnt, indices):
    result = "# Data Source: stochprop v" + version("stochprop")
    result = result + '\n' + "# Calculated: " + str(datetime.datetime.now())
    result = result + '\n' + "# Method: Gravity Wave Perturbation"
    result = result + '\n' + "# Reference Specification = " + prof_path + " (cwd: " + os.getcwd() + ")"
    result = result + '\n' + "# Gravity Wave Propagation Time [hr] = " + str(t0 / 3600.0)
    result = result + '\n' + "# Gravity Wave Spatial Scale Factor (dx) [km] = " + str(dx)
    result = result + '\n' + "# Gravity Wave Spatial Resolution (N_k) = " + str(Nk)
    result = result + '\n' + "# Gravity Wave Frequency Resolution (N_om) = " + str(N_om)
    result = result + '\n' + "# Gravity Wave Phase Randomization = " + str(random_phase)
    result = result + '\n' + "# Gravity Wave Source Elevation [km] = " + str(z_src)
    result = result + '\n' + "# Gravity Wave Source Wave Number (m_star) [km^{-1}] = " + str(m_star)
    result = result + '\n' + "# Gravity Wave Sample indices = [" + str(indices[0]) + ", " + str(indices[1]) + "]"
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

def perturb_atmo(atmo_spec, output_path, sample_cnt=50, t0=8.0 * 3600.0, dx=4.0, dz=0.2, Nk=128, N_om=5, random_phase=False, z_src=20.0, m_star=(2.0*np.pi)/2.5, env_below=True, cpu_cnt=None, fig_out=None):
    """
        Use gravity waves to perturb a specified profile using the methods in Drob et al. (2013)

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

    ref_atmo = np.loadtxt(atmo_spec)
    z0_vals = np.copy(ref_atmo[:, 0])
    T0_vals = np.copy(ref_atmo[:, 1])
    u0_vals = np.copy(ref_atmo[:, 2])
    v0_vals = np.copy(ref_atmo[:, 3])
    d0_vals = np.copy(ref_atmo[:, 4])
    p0_vals = np.copy(ref_atmo[:, 5])

    ref_lat = 40.0
    try:
        temp = open(atmo_spec, 'r')
        for line in temp:
            if "Location" in line:
                ref_lat = float(line.split(' ')[-4][:-1])
    except:
        ref_lat = 40.0

    if cpu_cnt:
        pl = Pool(cpu_cnt)
    else:
        pl = None

    z, du, dv, _, eta = perturbations(atmo_spec, t0=t0, dx=dx, dz=dz, Nk=Nk, N_om=N_om, ref_lat=ref_lat, random_phase=random_phase, z_src=z_src, m_star=m_star, figure_out=fig_out, pool=pl)

    np.save(output_path + ".z_vals", z)
    np.save(output_path + ".du_vals", du)
    np.save(output_path + ".dv_vals", dv)
    np.save(output_path + ".eta_vals", eta)

    if cpu_cnt:
        pl.terminate()
        pl.close()

    # Envelope perturbations off below the source height
    if env_below:
        temp = (1.0 + np.exp(-(z - z_src) / 2.5))
        du = du / temp
        dv = dv / temp
        eta = eta / temp

    # Sample atmosphere spatially
    print("Generating perturbed atmospheric samples...")
    n1_vals = np.random.default_rng().choice(Nk, size=sample_cnt, replace=False)
    n2_vals = np.random.default_rng().choice(Nk, size=sample_cnt, replace=False)
    for m in range(sample_cnt):
        du_interp = interp1d(z, du[n1_vals[m]][n2_vals[m]], fill_value="extrapolate", bounds_error=False)
        dv_interp = interp1d(z, dv[n1_vals[m]][n2_vals[m]], fill_value="extrapolate", bounds_error=False)
        eta_interp = interp1d(z, eta[n1_vals[m]][n2_vals[m]], fill_value="extrapolate", bounds_error=False)

        dT_vals = np.array([(T0_vals[min(n + 1, len(z0_vals) - 1)] - T0_vals[max(n - 1, 0)]) / (z0_vals[min(n + 1, len(z0_vals) - 1)] - z0_vals[max(n - 1, 0)]) for n in range(len(z0_vals))])
        dp_vals = np.array([(p0_vals[min(n + 1, len(z0_vals) - 1)] - p0_vals[max(n - 1, 0)]) / (z0_vals[min(n + 1, len(z0_vals) - 1)] - z0_vals[max(n - 1, 0)]) for n in range(len(z0_vals))])

        u_vals = u0_vals + du_interp(z0_vals) * 1000.0
        v_vals = v0_vals + dv_interp(z0_vals) * 1000.0
        T_vals = T0_vals # + dT_vals * eta_interp(z0_vals)
        p_vals = p0_vals # + dp_vals * eta_interp(z0_vals)

        np.savetxt(output_path + "-" + str(m) + ".met", np.vstack((z0_vals, T_vals, u_vals, v_vals, d0_vals, p_vals)).T, 
            header=_perturb_header_txt(atmo_spec, t0, dx, Nk, N_om, random_phase, z_src, m_star, m, sample_cnt, [n1_vals[m], n2_vals[m]]), comments='')

