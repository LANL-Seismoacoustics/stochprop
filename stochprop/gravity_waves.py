# gravity_saves.py
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

import itertools

import numpy as np

import matplotlib.pyplot as plt 

from scipy.integrate import simps
from scipy.interpolate import interp1d
from scipy.special import airy


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
        as :math:`cg = \frac{\partial \hat{omega}}{\partial m} = \frac{m k_h N}{k_h^2 + m^2 + \frac{1}{4H^2}`

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
    return visc * (abs(m_sqr(k, l, om_intr, H))**(3.0 / 2.0) / om_intr) * env 


def prep(k, l, om, atmo_info, om_min, k_max, src_index, m_star):
    """
        Identify critical layers ($ \hat{\omega}(z) = 0$) and turning points ($m(z) = 0$) and compute relevant quantities 
        for the trapped GW solution.  Also checks for valid Fourier component status, group velocity minimum, and other 
        conditions that would produce a bad Fourier component

        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om: float
            Absolute frequency (relative to the ground)
        atmo_info: 2darray
            Array with columns containing altitude , scale height, zonal wind, meridional wind, temperature, and density
        om_min: float
            Minimum absolute frequency used in analysis
        k_max: float
            Maximum horzintal wavenumber value used in 1 grid dimension
        src_index: int
            Index of the source height within the atmo_info z values
        m_star: float
            Source parameter m_* (default value, :math:`\frac{2 \pi}{2.5} \text{ km}^{-1}` is for 20 km altitude source)
        
        Returns:
        crit_lyr_index: int
            Index of the lowest critical layer altitude (set to len(z_vals) from atmo_info if no critical layer)
        turn_ht_index: int
            Index of the lowest turning height (set to len(z_vals) from atmo_info if no critical layer)
        w0: float
            Source spectral amplitude for this Fourier component
        refl_phase: float
            Reflection phase from ground to turning height (relevant if a turning point is identified)
        refl_time: float
            Reflection propagation time from ground to turning height (relevant if a turning point is identified)
        refl_loss: float
            Reflection attenuation losses from 100 km to turning height (relevant if a turning point is identified, and only computed if the turning height is above 100 km)

    """

    # extract atmospheric information
    z_vals = atmo_info[0]
    H_vals = atmo_info[1]
    u0_vals = atmo_info[2]
    v0_vals = atmo_info[3]
    T0_vals = atmo_info[4]
    d0_vals = atmo_info[5]

    # define parameters to return
    crit_lyr_index, turn_ht_index = 0, 0
    w0, refl_phase, refl_time, refl_loss = 0.0, 0.0, 0.0, 0.0

    # compute combined horizontal wavenumber, intrinsic frequency, and vertical wavenumber values
    kh = np.sqrt(k**2 + l**2)
    om_intr_vals = om - k * u0_vals - l * v0_vals
    m_sqr_vals = m_sqr(k, l, om_intr_vals, H_vals)

    if 0.0 < kh and kh <= k_max:
        # check for critical layer
        if om_intr_vals[0] <= 0.0:
            crit_lyr_index = 0
        elif np.all(om_intr_vals > 0.0):
            crit_lyr_index = len(z_vals) + 1
        else:
            jz = np.where(om_intr_vals < 0.0)[0][0]
            if jz > src_index + 2:
                crit_lyr_index = jz - 1
            else:
                crit_lyr_index = 0

        # check for turning height
        if m_sqr_vals[0] <= 0.0:
            turn_ht_index = 0
        elif np.all(m_sqr_vals > 0.0):
            turn_ht_index = len(z_vals) + 1
        else:
            jz = np.where(m_sqr_vals < 0.0)[0][0]
            if jz > src_index + 2:
                turn_ht_index = jz - 1
                refl_phase = simps(np.sqrt(m_sqr_vals[:turn_ht_index]), z_vals[:turn_ht_index])
                refl_time = simps(1.0 / abs(cg(k, l, om_intr_vals[:turn_ht_index], H_vals[:turn_ht_index])), z_vals[:turn_ht_index])
                refl_loss = simps(m_imag(k, l, om_intr_vals[:turn_ht_index], z_vals[:turn_ht_index], H_vals[:turn_ht_index], T0_vals[:turn_ht_index], d0_vals[:turn_ht_index]), z_vals[:turn_ht_index])
            else:
                turn_ht_index = 0
        
        # check velocity (N/m < 90 m/s)
        min_index = min(crit_lyr_index, turn_ht_index)
        if min_index > 0:
            velocity_check = np.max(BV_freq(H_vals[:min_index]) / np.sqrt(m_sqr(k, l, om_intr_vals[:min_index], H_vals[:min_index])))
            if abs(velocity_check) > 0.09:
                turn_ht_index = 0
                crit_lyr_index = 0
            else:
                m_sqr_src = m_sqr(k, l, om_intr_vals[src_index], H_vals[src_index])
                w0 = 2.7e-2 * (m_sqr_src / (m_star**4 + m_sqr_src**2)) * (om_intr_vals[0]**(1.0 / 3.0) / kh**2)
                w0 = w0 * om_min**(2.0 / 3.0) / (1.0 - (om_min / BV_freq(H_vals[0]))**(2.0 / 3.0))
                w0 = np.sqrt(abs(w0))
        
    return [crit_lyr_index, turn_ht_index, w0, refl_phase, refl_time, refl_loss]


def prop_upward(k, l, om, zj, t0, crit_lyr_index, turning_ht_index, w0, atmo_info, trapped_info, src_index):
    """
        Integrate the gravity wave spectrum upward using the freely propagating solution for Fourier
        combinations without turning points and the Airy function form for trapped components.


        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om: float
            Absolute frequency (relative to the ground)
        zj: int
            Index of the current altitude point
        t0: float
            Reference propagation time used for free solution maximum altitude and trapped solution reflection count
        crit_lyr_index: int
            Index of the shallowest critical layer
        turning_ht_index: int
            Index of the shallowest turning height
        w0: Float
            Spectral amplitude at the source
        atmo_info: 2darray
            Array with columns containing altitude , scale height, zonal wind, meridional wind, temperature, and density
        trapped_info: float
            Reflection phase, propagation time, and losses for trapped solution
        src_index: int
            Index of the source altitude in atmo_info z_vals

        Returns
        u_spec: Complex float
            Fourier component for the zonal wind perturbation
        v_spec: Complex float
            Fourier component for the meridional wind perturbation
        w_spec: Complex float
            Fourier component of the vertical wind perturbation
        w_phase: Float
            Cumulative vertical Fourier component phase
        w_loss: Float
            Cumulative vertical Fourier component attenuation

    """

    # extract atmospheric and trapped solution information
    z_vals = atmo_info[0]
    H_vals = atmo_info[1]
    u0_vals = atmo_info[2]
    v0_vals = atmo_info[3]
    T0_vals = atmo_info[4]
    d0_vals = atmo_info[5]

    refl_phase = trapped_info[0]
    refl_time = trapped_info[1]
    refl_loss = trapped_info[2]

    dz = (z_vals[zj] - z_vals[zj - 1])

    # define returned values
    u_spec, v_spec, w_spec = 0.0, 0.0, 0.0

    # Check that m was valid at the source, altitude is below critical layer, and no more than 2.0 km above a turning height
    if turning_ht_index > 0 and zj < crit_lyr_index and zj < min(len(z_vals) + 1, turning_ht_index + int(2.0 / dz)):
        kh = np.sqrt(k**2 + l**2)
        om_intr = om - k * u0_vals - l * v0_vals

        # Check reference time vs. propagation time to this altitude
        if zj <= crit_lyr_index:
            prop_time = simps(1.0 / abs(cg(k, l, om_intr[src_index:zj], H_vals[src_index:zj])), z_vals[src_index:zj])
        else:
            prop_time = np.inf 
                            
        if prop_time < t0:
            if z_vals[zj] > 80.0:
                m_imag_vals = m_imag(k, l, om_intr[src_index:zj], z_vals[src_index:zj], H_vals[src_index:zj], T0_vals[src_index:zj], d0_vals[src_index:zj])
                w_loss = simps(m_imag_vals, z_vals[src_index:zj])
            else:
                w_loss = 0.0

            # Define saturation spectra for this Fourier component
            w_sat_sqr = 2.7e-2 * om_intr[zj]**(1.0 / 3.0) / (m_sqr(k, l, om_intr[zj], H_vals[zj]) * kh**2)

            # integrate using free solution if there is no turning point or if the 
            # evaluation time is less than the trapped propagation time
            m_src = np.sqrt(abs(m_sqr(k, l, om_intr[src_index], H_vals[src_index])))
            m_curr = np.sqrt(abs(m_sqr(k, l, om_intr[zj], H_vals[zj])))
            d0_m_ratio = np.sqrt((d0_vals[src_index] / d0_vals[zj]) * (m_src / m_curr))

            if turning_ht_index > len(z_vals) or t0 < refl_time:
                m_vals = np.sqrt(abs(m_sqr(k, l, om_intr[src_index:zj], H_vals[src_index:zj])))
                w_phase = simps(m_vals, z_vals[src_index:zj])

                w_spec = w0 * d0_m_ratio * (np.cos(w_phase) - 1.0j * np.sin(w_phase)) * np.exp(-w_loss)

                if abs(w_spec) > np.sqrt(abs(w_sat_sqr)):
                    w_spec = w_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))

                u_spec = - w_spec * (k / kh**2) * m_curr
                v_spec = - w_spec * (l / kh**2) * m_curr

            else:
                # Alternately, use the Airy function form in the case of trapped waves
                if zj < turning_ht_index:
                    m_vals = np.sqrt(abs(m_sqr(k, l, om_intr[zj:turning_ht_index], H_vals[zj:turning_ht_index])))
                    airy_arg = -((3.0 / 2.0) * simps(m_vals, z_vals[zj:turning_ht_index]))**(2.0 / 3.0)
                elif zj > turning_ht_index:
                    # check that m(z) doesn't become real again at some altitude 
                    # above the first reflection point and stop integration if it does
                    m_sqr_vals = m_sqr(k, l, om_intr[turning_ht_index:zj], H_vals[turning_ht_index:zj])
                    if np.all(m_sqr_vals < 0.0):
                        j_temp = zj
                    else:
                        j_temp = turning_ht_index + (np.where(m_sqr_vals > 0.0)[0][0] - 1)

                    m_vals = np.sqrt(-m_sqr(k, l, om_intr[turning_ht_index:j_temp], H_vals[turning_ht_index:j_temp]))
                    airy_arg = ((3.0 / 2.0) * simps(m_vals, z_vals[turning_ht_index:j_temp]))**(2.0 / 3.0)
                else:
                    airy_arg = 0.0

                refl_cnt = int(np.floor(t0 / (2.0 * refl_time)))
                refl_ph_shift = np.exp(-refl_cnt * (2.0 * refl_loss) - w_loss)
                refl_ph_shift = refl_ph_shift * (1.0 + np.sum(np.array([np.exp(-1.0j * (N - 1) * (2.0 * refl_phase - np.pi / 2.0)) for N in range(2, refl_cnt)])))

                airy_scaling = np.array([(-airy_arg)**0.25 * airy(airy_arg)[n] * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift for n in range(2)])

                u_spec = -2.0j * np.sqrt(np.pi) * w0 * (k / kh**2) * d0_m_ratio * airy_scaling[1]
                v_spec = -2.0j * np.sqrt(np.pi) * w0 * (l / kh**2) * d0_m_ratio * airy_scaling[1]
                w_spec = -2.0j * np.sqrt(np.pi) * w0 * d0_m_ratio * airy_scaling[0]

                if abs(w_spec) > np.sqrt(abs(w_sat_sqr)):
                    u_spec = u_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))
                    v_spec = v_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))
                    w_spec = w_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))

    return [u_spec, v_spec, w_spec]


def prop_downward(k, l, om, zj, t0, crit_lyr_index, turning_ht_index, w0, atmo_info, trapped_info, src_index):
    """
        Compute the gravity wave spectum below the source height ignoring refraction effects

        Parameters
        ----------
        k: float
            Zonal wave number [km^{-1}]
        l: float
            Meridional wave number [km^{-1}]
        om: float
            Absolute frequency (relative to the ground)
        zj: int
            Index of the current altitude point
        t0: float
            Reference propagation time used for free solution maximum altitude and trapped solution reflection count
        crit_lyr_index: int
            Index of the shallowest critical layer
        turning_ht_index: int
            Index of the shallowest turning height
        w0: Float
            Spectral amplitude at the source
        atmo_info: 2darray
            Array with columns containing altitude , scale height, zonal wind, meridional wind, temperature, and density
        trapped_info: float
            Reflection phase, propagation time, and losses for trapped solution
        src_index: int
            Index of the source altitude in atmo_info z_vals

        Returns
        u_spec: Complex float
            Fourier component for the zonal wind perturbation
        v_spec: Complex float
            Fourier component for the meridional wind perturbation
        w_spec: Complex float
            Fourier component of the vertical wind perturbation
        w_phase: Float
            Cumulative vertical Fourier component phase
        w_loss: Float
            Cumulative vertical Fourier component attenuation
    """

    # extract atmospheric and trapped solution information
    z_vals = atmo_info[0]
    H_vals = atmo_info[1]
    u0_vals = atmo_info[2]
    v0_vals = atmo_info[3]
    d0_vals = atmo_info[5]

    refl_phase = trapped_info[0]
    refl_time = trapped_info[1]
    refl_loss = trapped_info[2]

    # define returned values
    u_spec, v_spec, w_spec = 0.0, 0.0, 0.0

    # Check that m was valid at the source, altitude is below critical layer, and no more than 2.0 km above a turning height
    if turning_ht_index > 0 and zj < crit_lyr_index:
        
        kh = np.sqrt(k**2 + l**2)
        om_intr = om - k * u0_vals - l * v0_vals
    
        # Define saturation spectra for this Fourier component
        w_sat_sqr = 2.7e-2 * om_intr[zj]**(1.0 / 3.0) / (m_sqr(k, l, om_intr[zj], H_vals[zj]) * kh**2)

        # integrate using free solution if there is no turning point or if the 
        # evaluation time is less than the trapped propagation time
        m_val = np.sqrt(abs(m_sqr(k, l, om_intr[src_index], H_vals[src_index])))
        w_phase = m_val * (z_vals[src_index] - z_vals[zj])

        d0_ratio = np.sqrt(d0_vals[src_index] / d0_vals[zj])

        if turning_ht_index > len(z_vals) or t0 < refl_time:
            w_spec = w0 * d0_ratio * (np.cos(w_phase) - 1.0j * np.sin(w_phase))

            if abs(w_spec) > np.sqrt(abs(w_sat_sqr)):
                w_spec = w_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))

            u_spec = - w_spec * (k / kh**2) * m_val
            v_spec = - w_spec * (l / kh**2) * m_val

        else:
            # Alternately, extend the Airy function form in the case of trapped waves
            m_above = np.sqrt(abs(m_sqr(k, l, om_intr[src_index:turning_ht_index], H_vals[src_index:turning_ht_index])))
            airy_arg = -((3.0 / 2.0) * simps(m_above, z_vals[src_index:turning_ht_index]) + w_phase)**(2.0 / 3.0)

            refl_cnt = int(np.floor(t0 / (2.0 * refl_time)))
            refl_ph_shift = np.exp(-refl_cnt * (2.0 * refl_loss))
            refl_ph_shift = refl_ph_shift * (1.0 + np.sum(np.array([np.exp(-1.0j * (N - 1) * (2.0 * refl_phase - np.pi / 2.0)) for N in range(2, refl_cnt)])))

            airy_scaling = np.array([(-airy_arg)**0.25 * airy(airy_arg)[n] * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift for n in range(2)])

            u_spec = -2.0j * np.sqrt(np.pi) * w0 * (k / kh**2) * d0_ratio * airy_scaling[1]
            v_spec = -2.0j * np.sqrt(np.pi) * w0 * (l / kh**2) * d0_ratio * airy_scaling[1]
            w_spec = -2.0j * np.sqrt(np.pi) * w0 * d0_ratio * airy_scaling[0]

            if abs(w_spec) > np.sqrt(abs(w_sat_sqr)):
                u_spec = u_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))
                v_spec = v_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))
                w_spec = w_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))

    return [u_spec, v_spec, w_spec]


def build_spec(atmo_specification, t0=4.0 * 3600.0, dx=4.0, dz=0.2, Nk=128, N_om=5, ref_lat=30.0, random_phase=False, z_src=20.0, m_star=2.0*np.pi/2.5):
    """
        Function definition...

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

    # Set up atmosphere info
    z, Temp0, u0, v0, d0, p0 = np.loadtxt(atmo_specification, unpack=True)

    T0_interp = interp1d(z, Temp0)
    u0_interp = interp1d(z, u0 * 1.0e-3)
    v0_interp = interp1d(z, v0 * 1.0e-3)
    d0_interp = interp1d(z, d0 * 1.0e3)

    d0_finite_diffs = np.array([(d0[min(n + 1, len(z) - 1)] - d0[max(n - 1, 0)]) / (z[min(n + 1, len(z) - 1)] - z[max(n - 1, 0)]) for n in range(len(z))])
    sc_ht = interp1d(z, -d0 / d0_finite_diffs)

    z_vals = np.arange(0.0, z[-1], dz)
    atmo_info = [z_vals, sc_ht(z_vals), u0_interp(z_vals), v0_interp(z_vals), T0_interp(z_vals), d0_interp(z_vals)]
    src_index = np.argmin(abs(z_vals - z_src))

    # Define spectral component grid
    om_min = 2.0 * 7.292e-5 * np.sin(np.radians(ref_lat))
    om_max = np.max(BV_freq(sc_ht(z))) / np.sqrt(5)

    k_vals = np.fft.fftfreq(Nk, d=(2.0 * np.pi) / dx)
    l_vals = np.fft.fftfreq(Nk, d=(2.0 * np.pi) / dx)
    om_vals = np.linspace(om_min, om_max, N_om)

    # Define propagation metadata (turning points, critical layers, etc.)
    turn_ht_indices = np.zeros((Nk, Nk, N_om), dtype=int)
    crit_lyr_indices = np.zeros_like(turn_ht_indices, dtype=int)
    w0_vals = np.zeros_like(turn_ht_indices, dtype=float)
    refl_info = np.zeros((Nk, Nk, N_om, 3), dtype=float)

    # Define spectral information for each Fourier component
    u_spec = np.zeros((Nk, Nk, N_om, len(z_vals)), dtype=complex)
    v_spec = np.zeros_like(u_spec, dtype=complex)
    w_spec = np.zeros_like(u_spec, dtype=complex)

    w_phase = np.zeros_like(w_spec, dtype=float)

    print("Running pre-analysis to identify turning heights and critical layers...")
    results = np.array([[[prep(k, l, om, atmo_info, om_vals[0], max(abs(k_vals)), src_index, m_star) for om in om_vals] for l in l_vals] for k in k_vals])
    for kn, ln in itertools.product(range(Nk), repeat=2):
        for om_n in range(N_om):
            crit_lyr_indices[kn][ln][om_n] = results[kn][ln][om_n][0]
            turn_ht_indices[kn][ln][om_n] = results[kn][ln][om_n][1]
            w0_vals[kn][ln][om_n] = results[kn][ln][om_n][2]
            refl_info[kn][ln][om_n] = results[kn][ln][om_n][3:]

    # Randomize initial phase and set initial solution
    if random_phase:
        w_phase[:,:,:,src_index] = np.random.rand((Nk, Nk, N_om)) * (2.0 * np.pi)
    w_spec[:, :, :, src_index] = w0_vals * (np.cos(w_phase[:,:,:,0]) - 1.0j * np.sin(w_phase[:,:,:,0]))   
    
    # set solution below source
    for j in range(0, src_index):
        if j % 5 == 0:
            print("Defining solution below source height.  Currently at z = " + str(np.round(z_vals[j], 2)) + " km...")

        results = np.array([[[prop_downward(k_vals[kn], l_vals[ln], om_vals[om_n], j, t0, crit_lyr_indices[kn][ln][om_n], turn_ht_indices[kn][ln][om_n], 
        w0_vals[kn][ln][om_n], atmo_info, refl_info[kn][ln][om_n], src_index) for om_n in range(N_om)] for ln in range(Nk)] for kn in range(Nk)])

        for kn, ln in itertools.product(range(Nk), repeat=2):
            for om_n in range(N_om):
                u_spec[kn][ln][om_n][j] = results[kn][ln][om_n][0]
                v_spec[kn][ln][om_n][j] = results[kn][ln][om_n][1]
                w_spec[kn][ln][om_n][j] = results[kn][ln][om_n][2]
                
    # integrate upward
    for j in range(src_index + 1, len(z_vals)):
        if j % 5 == 0:
            print("Integrating upward.  Currently at z = " + str(np.round(z_vals[j], 2)) + " km...")
        results = np.array([[[prop_upward(k_vals[kn], l_vals[ln], om_vals[om_n], j, t0, crit_lyr_indices[kn][ln][om_n], turn_ht_indices[kn][ln][om_n], 
        w0_vals[kn][ln][om_n], atmo_info, refl_info[kn][ln][om_n], src_index) for om_n in range(N_om)] for ln in range(Nk)] for kn in range(Nk)])
        
        for kn, ln in itertools.product(range(Nk), repeat=2):
            for om_n in range(N_om):
                u_spec[kn][ln][om_n][j] = results[kn][ln][om_n][0]
                v_spec[kn][ln][om_n][j] = results[kn][ln][om_n][1]
                w_spec[kn][ln][om_n][j] = results[kn][ln][om_n][2]

    du_vals = np.fft.ifft(u_spec, axis=0) * Nk / dx
    du_vals = np.fft.ifft(du_vals, axis=1) * Nk / dx
    du_vals = simps(du_vals, om_vals, axis=2)

    dv_vals = np.fft.ifft(v_spec, axis=0) * Nk / dx
    dv_vals = np.fft.ifft(dv_vals, axis=1) * Nk / dx
    dv_vals = simps(dv_vals, om_vals, axis=2)

    dw_vals = np.fft.ifft(w_spec, axis=0) * Nk / dx
    dw_vals = np.fft.ifft(dw_vals, axis=1) * Nk / dx
    dw_vals = simps(dw_vals, om_vals, axis=2)

    return z_vals, np.real(du_vals), np.real(dv_vals), np.real(dw_vals)



def single_fourier_component(k, l, om, atmo_info, t0, src_index, m_star, om_min, k_max):
    """
        Function definition...

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

        Returns
        -------
        u_spec: 1darray
            Zonal wind perturbation spectrum, \hat{u}(k, l, z, omega)
        v_spec: 1darray
            Meridional wind perturbation spectrum, \hat{v}(k, l, z, omega)
        w_spec: 1darray
            Vertical wind perturbation spectrum, \hat{w}(k, l, z, omega)

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
    om_intr_vals = om - k * u0_vals - l * v0_vals
    m_sqr_vals = m_sqr(k, l, om_intr_vals, H_vals)

    if kh < 1.0e-6 or kh > k_max:
        return [np.zeros_like(z_vals, dtype=complex)] * 3
    else:
        # check for critical layer and turning height, then 
        # set the source spectral component
        if np.all(om_intr_vals > 0.0):
            crit_lyr_index = len(z_vals) + 1
        else:
            jz = np.where(om_intr_vals < 0.0)[0][0]
            if jz > src_index + 2:
                crit_lyr_index = jz - 1
            else:
                return [np.zeros_like(z_vals, dtype=complex)] * 3

        # check for turning height
        if np.all(m_sqr_vals > 0.0):
            turn_ht_index = len(z_vals) + 1
        else:
            jz = np.where(m_sqr_vals < 0.0)[0][0]
            if jz > src_index + 2:
                turn_ht_index = jz - 1
                refl_phase = simps(np.sqrt(m_sqr_vals[:turn_ht_index]), z_vals[:turn_ht_index])
                refl_time = simps(1.0 / abs(cg(k, l, om_intr_vals[:turn_ht_index], H_vals[:turn_ht_index])), z_vals[:turn_ht_index])
                refl_loss = simps(m_imag(k, l, om_intr_vals[:turn_ht_index], z_vals[:turn_ht_index], H_vals[:turn_ht_index], T0_vals[:turn_ht_index], d0_vals[:turn_ht_index]), z_vals[:turn_ht_index])
            else:
                return [np.zeros_like(z_vals, dtype=complex)] * 3
        
        # check velocity condition (N/m < 90 m/s)
        min_index = min(crit_lyr_index, turn_ht_index)
        if min_index > 0:
            velocity_check = np.max(BV_freq(H_vals[:min_index]) / np.sqrt(m_sqr(k, l, om_intr_vals[:min_index], H_vals[:min_index])))
            if abs(velocity_check) < 0.09:
                m_sqr_src = m_sqr(k, l, om_intr_vals[src_index], H_vals[src_index])
                w0 = 2.7e-2 * (m_sqr_src / (m_star**4 + m_sqr_src**2)) * (om_intr_vals[0]**(1.0 / 3.0) / kh**2)
                w0 = w0 * om_min**(2.0 / 3.0) / (1.0 - (om_min / BV_freq(H_vals[0]))**(2.0 / 3.0))
                w0 = np.sqrt(abs(w0))        
            else:
                return [np.zeros_like(z_vals, dtype=complex)] * 3

        # prep numpy arrays for spectral info
        u_spec = np.zeros_like(z_vals, dtype=complex)
        v_spec = np.zeros_like(z_vals, dtype=complex)
        w_spec = np.zeros_like(z_vals, dtype=complex)

        # set source elevation values
        m_src = np.sqrt(m_sqr_vals[src_index])
        w_spec[src_index] = w0
        u_spec[src_index] = - w0 * (k / kh**2) * m_src
        v_spec[src_index] = - w0 * (l / kh**2) * m_src

        # Evaluate values below the source (ignore refraction)
        if turn_ht_index > len(z_vals):
            # free propagating solution
            w_phase = np.array([m_src * (z_vals[zj] - z_vals[src_index]) for zj in range(src_index)])
            w_spec[:src_index] = w0 * np.sqrt(d0_vals[src_index] / d0_vals[:src_index]) * (np.cos(w_phase) - 1.0j * np.sin(w_phase))
            u_spec[:src_index] = - w_spec[:src_index] * (k / kh**2) * m_src
            v_spec[:src_index] = - w_spec[:src_index] * (l / kh**2) * m_src
        else:
            # Airy function form in the case of trapped waves
            for zj in range(src_index):
                m_above = np.sqrt(abs(m_sqr(k, l, om_intr_vals[src_index:turn_ht_index], H_vals[src_index:turn_ht_index])))

                airy_arg = -((3.0 / 2.0) * simps(m_above, z_vals[src_index:turn_ht_index]) + m_src * abs(z_vals[zj] - z_vals[src_index]))**(2.0 / 3.0)

                refl_cnt = int(np.floor(t0 / (2.0 * refl_time)))
                refl_ph_shift = np.exp(-refl_cnt * (2.0 * refl_loss))
                refl_ph_shift = refl_ph_shift * (1.0 + np.sum(np.array([np.exp(-1.0j * (N - 1) * (2.0 * refl_phase - np.pi / 2.0)) for N in range(2, refl_cnt)])))

                airy_scaling = np.array([(-airy_arg)**0.25 * airy(airy_arg)[n] * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift for n in range(2)])

                u_spec[zj] = -2.0j * np.sqrt(np.pi) * w0 * (k / kh**2) * (np.sqrt(d0_vals[src_index] / d0_vals[zj])) * airy_scaling[1]
                v_spec[zj] = -2.0j * np.sqrt(np.pi) * w0 * (l / kh**2) * (np.sqrt(d0_vals[src_index] / d0_vals[zj])) * airy_scaling[1]
                w_spec[zj] = -2.0j * np.sqrt(np.pi) * w0 * (np.sqrt(d0_vals[src_index] / d0_vals[zj])) * airy_scaling[0]

        # Integrate solution upward
        for zj in range(src_index + 1, len(z_vals)):
            # compute propagation time to altitude
            if zj <= crit_lyr_index:
                prop_time = simps(1.0 / abs(cg(k, l, om_intr_vals[src_index:zj], H_vals[src_index:zj])), z_vals[src_index:zj])
            else:
                prop_time = np.inf 

            # check reference time vs. propagation time and that altitude isn't too far above turning height
            if prop_time < t0 and zj < min(len(z_vals) + 1, turn_ht_index + int(2.0 / dz)):
                # define attenuation losses
                if z_vals[zj] > 90.0:
                    m_imag_vals = m_imag(k, l, om_intr_vals[src_index:zj], z_vals[src_index:zj], H_vals[src_index:zj], T0_vals[src_index:zj], d0_vals[src_index:zj])
                    w_loss = simps(m_imag_vals, z_vals[src_index:zj])
                else:
                    w_loss = 0.0

                # define saturation spectra
                w_sat = np.sqrt(2.7e-2 * om_intr_vals[zj]**(1.0 / 3.0) / (m_sqr_vals[zj] * kh**2))

                # integrate using free solution if there is no turning point or if the 
                # evaluation time is less than the trapped propagation time
                m_curr = np.sqrt(m_sqr_vals[zj])
                d0_m_ratio = np.sqrt((d0_vals[src_index] / d0_vals[zj]) * (m_src / m_curr))

                if turn_ht_index > len(z_vals) or t0 < refl_time:
                    w_phase = simps(np.sqrt(m_sqr_vals[src_index:zj]), z_vals[src_index:zj])
                    w_spec[zj] = w0 * d0_m_ratio * (np.cos(w_phase) - 1.0j * np.sin(w_phase))

                    if abs(w_spec[zj]) > w_sat:
                        w_spec[zj] = w_spec[zj] * (w_sat / abs(w_spec[zj]))

                    u_spec[zj] = - w_spec[zj] * (k / kh**2) * m_curr * np.exp(-w_loss)
                    v_spec[zj] = - w_spec[zj] * (l / kh**2) * m_curr * np.exp(-w_loss)
                    w_spec[zj] = w_spec[zj] * np.exp(-w_loss)
                else:
                    if zj < turn_ht_index:
                        airy_arg = -((3.0 / 2.0) * simps(np.sqrt(m_sqr_vals[zj:turn_ht_index]), z_vals[zj:turn_ht_index]))**(2.0 / 3.0)
                    elif zj > turn_ht_index:
                        if np.all(m_sqr_vals[turn_ht_index:zj] < 0.0):
                            j_temp = zj
                        else:
                            j_temp = turn_ht_index + (np.where(m_sqr_vals > 0.0)[0][0] - 1)
                        airy_arg = ((3.0 / 2.0) * simps(np.sqrt(abs(m_sqr_vals[turn_ht_index:j_temp])), z_vals[turn_ht_index:j_temp]))**(2.0 / 3.0)
                    else:
                        airy_arg = 0.0

                    refl_cnt = int(np.floor(t0 / (2.0 * refl_time)))
                    refl_ph_shift = np.exp(-refl_cnt * (2.0 * refl_loss) - w_loss)
                    refl_ph_shift = refl_ph_shift * (1.0 + np.sum(np.array([np.exp(-1.0j * (N - 1) * (2.0 * refl_phase - np.pi / 2.0)) for N in range(2, refl_cnt)])))

                    airy_scaling = np.array([(-airy_arg)**0.25 * airy(airy_arg)[n] * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift for n in range(2)])

                    u_spec[zj] = -2.0j * np.sqrt(np.pi) * w0 * (k / kh**2) * d0_m_ratio * airy_scaling[1]
                    v_spec[zj] = -2.0j * np.sqrt(np.pi) * w0 * (l / kh**2) * d0_m_ratio * airy_scaling[1]
                    w_spec[zj] = -2.0j * np.sqrt(np.pi) * w0 * d0_m_ratio * airy_scaling[0]

                    if abs(w_spec[zj]) > w_sat:
                        u_spec[zj] = u_spec[zj] * (w_sat / abs(w_spec[zj]))
                        v_spec[zj] = v_spec[zj] * (w_sat / abs(w_spec[zj]))
                        w_spec[zj] = w_spec[zj] * (w_sat / abs(w_spec[zj]))
            else:
                u_spec[zj] = 0.0 + 0.0j
                v_spec[zj] = 0.0 + 0.0j
                w_spec[zj] = 0.0 + 0.0j

        f, a = plt.subplots(1, 3, sharey=True)
        a[0].semilogx(m_sqr_vals, z_vals, '--k')
        a[1].plot(om_intr_vals, z_vals, '-.k')
        a[2].plot(np.real(w_spec), z_vals, '-b')
        a[2].plot(np.imag(w_spec), z_vals, '-r')
        for n in range(3):
            a[n].axhline(z_vals[min(turn_ht_index, len(z_vals) - 1)], color="g")
            a[n].axhline(z_vals[min(crit_lyr_index, len(z_vals) - 1)], color="m")
        a[1].set_title("(k, l, om) = (" + str(k) + ", " + str(l) + ", " + str(om) + ")")
        plt.savefig("results-" + str(k) + "-" + str(l) + "-" + str(om) + ".png")

        return [u_spec, v_spec, w_spec]


def build_spec_v2(atmo_specification, t0=4.0 * 3600.0, dx=4.0, dz=0.2, Nk=128, N_om=5, ref_lat=30.0, random_phase=False, z_src=20.0, m_star=2.0*np.pi/2.5):
    """
        Version that integrates each Fourier component independently...

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

    # Set up atmosphere info
    z, Temp0, u0, v0, d0, p0 = np.loadtxt(atmo_specification, unpack=True)

    T0_interp = interp1d(z, Temp0)
    u0_interp = interp1d(z, u0 * 1.0e-3)
    v0_interp = interp1d(z, v0 * 1.0e-3)
    d0_interp = interp1d(z, d0 * 1.0e3)

    d0_finite_diffs = np.array([(d0[min(n + 1, len(z) - 1)] - d0[max(n - 1, 0)]) / (z[min(n + 1, len(z) - 1)] - z[max(n - 1, 0)]) for n in range(len(z))])
    sc_ht = interp1d(z, -d0 / d0_finite_diffs)

    z_vals = np.arange(0.0, z[-1], dz)
    atmo_info = [z_vals, sc_ht(z_vals), u0_interp(z_vals), v0_interp(z_vals), T0_interp(z_vals), d0_interp(z_vals)]
    src_index = np.argmin(abs(z_vals - z_src))

    # Define spectral component grid
    om_min = 2.0 * 7.292e-5 * np.sin(np.radians(ref_lat))
    om_max = np.max(BV_freq(sc_ht(z))) / np.sqrt(5)

    k_vals = np.fft.fftfreq(Nk, d=(2.0 * np.pi) / dx)
    l_vals = np.fft.fftfreq(Nk, d=(2.0 * np.pi) / dx)
    om_vals = np.linspace(om_min, om_max, N_om)

    # Define spectral information for each Fourier component
    u_spec = np.zeros((Nk, Nk, N_om, len(z_vals)), dtype=complex)
    v_spec = np.zeros_like(u_spec, dtype=complex)
    w_spec = np.zeros_like(u_spec, dtype=complex)
               
    # Compute each Fourier component
    for nk, k in enumerate(k_vals):
        for nl, l in enumerate(l_vals):
            for n_om, om in enumerate(om_vals):
                results = single_fourier_component(k, l, om, atmo_info, t0, src_index, m_star, om_vals[0], max(abs(k_vals)))

                u_spec[nk][nl][n_om] = results[0]
                v_spec[nk][nl][n_om] = results[1]
                w_spec[nk][nl][n_om] = results[2]

    du_vals = np.fft.ifft(u_spec, axis=0) * Nk / dx
    du_vals = np.fft.ifft(du_vals, axis=1) * Nk / dx
    du_vals = simps(du_vals, om_vals, axis=2)

    dv_vals = np.fft.ifft(v_spec, axis=0) * Nk / dx
    dv_vals = np.fft.ifft(dv_vals, axis=1) * Nk / dx
    dv_vals = simps(dv_vals, om_vals, axis=2)

    dw_vals = np.fft.ifft(w_spec, axis=0) * Nk / dx
    dw_vals = np.fft.ifft(dw_vals, axis=1) * Nk / dx
    dw_vals = simps(dw_vals, om_vals, axis=2)

    return z_vals, np.real(du_vals), np.real(dv_vals), np.real(dw_vals)

