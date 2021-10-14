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

from scipy.integrate import trapz
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
    return np.sqrt(0.0098 / H)


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


def m_imag(k, l, om_intr, H, T0, d0):
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
    return 3.563e-7 * (T0**0.69 / d0) * (abs(m_sqr(k, l, om_intr, H))**(3.0 / 2.0) / om_intr) * 1e-4


def prep(k, l, om, atmo_info, om_min, k_max):
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
            
        T0: float
            Ambient temperature in the atmosphere
        d0: float
            Ambient density in the atmosphere
        
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

    # reference vertical wavenumber for source definition 
    m_star = 2.0 * np.pi / 2.5
    
    kh = np.sqrt(k**2 + l**2)
    if 0.0 < kh and kh <= k_max:
        # check for critical layer
        om_intr = om - k * u0_vals - l * v0_vals

        if om_intr[0] <= 0.0:
            crit_lyr_index = 0
        elif np.all(om_intr > 0.0):
            crit_lyr_index = len(z_vals) + 1
        else:
            jz = np.where(om_intr < 0.0)[0][0]
            if jz > 2:
                crit_lyr_index = jz - 1

        # check for turning points and compute trapped solution info if necessary
        m_sqr_src = m_sqr(k, l, om_intr[0], H_vals[0])
        if m_sqr_src > 0.0 and crit_lyr_index > 0:
            w0 = 2.7e-2 * (m_sqr_src / (m_star**4 + m_sqr_src**2)) * (1.0 / (om_intr[0] * kh**2))
            w0 = w0 * om_min**(2.0 / 3.0) / (1.0 - (om_min / BV_freq(H_vals[0]))**(2.0 / 3.0))
            w0 = np.sqrt(abs(w0))

            m_sqr_vals = m_sqr(k, l, om_intr, H_vals)
            if np.all(m_sqr_vals > 0.0):
                turn_ht_index = len(z_vals) + 1
            else:
                jz = np.where(m_sqr_vals < 0.0)[0][0]
                if jz > crit_lyr_index:
                    turn_ht_index = len(z_vals) + 1
                elif jz > 2:
                    turn_ht_index = jz - 1
                    refl_phase = trapz(np.sqrt(m_sqr_vals[:turn_ht_index]), z_vals[:turn_ht_index])
                    refl_time = trapz(1.0 / abs(cg(k, l, om_intr[:turn_ht_index], H_vals[:turn_ht_index])), z_vals[:turn_ht_index])
                    refl_loss = trapz(m_imag(k, l, om_intr[:turn_ht_index], H_vals[:turn_ht_index], T0_vals[:turn_ht_index], d0_vals[:turn_ht_index]), z_vals[:turn_ht_index])

        # mask out waves for which C = N/m > C_{max} (~90 m/s)
        max_index = min(crit_lyr_index, turn_ht_index)
        if max_index > 0:
            velocity_check = np.max(BV_freq(H_vals[:max_index]) / np.sqrt(m_sqr(k, l, om_intr[:max_index], H_vals[:max_index])))
            if abs(velocity_check) > 0.09:
                turn_ht_index = 0
                crit_lyr_index = 0
    else:
        crit_lyr_index = 0
        turn_ht_index = 0
        
    return [crit_lyr_index, turn_ht_index, w0, refl_phase, refl_time, refl_loss]


def prop_upward(k, l, om, zj, t0, crit_lyr_index, turning_ht_index, w0, w_phase_prev, w_loss_prev, atmo_info, trapped_info):
    """
    Integrate the gravity wave spectrum upward using the freely propagating solution for Fourier
    combinations without turning points and the Airy function form for trapped components.



    """

    # extract atmospheric information
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
    w_phase, w_loss = w_phase_prev, w_loss_prev

    # Check that m was valid at the source, altitude is below critical layer, and no more than 2.0 km above a turning height
    if turning_ht_index > 0 and zj < crit_lyr_index and zj < min(len(z_vals) + 1, turning_ht_index + int(2.0 / dz)):
        kh = np.sqrt(k**2 + l**2)
        om_intr = om - k * u0_vals - l * v0_vals

        # Check reference time vs. propagation time to this altitude
        if zj <= crit_lyr_index:
            prop_time = trapz(1.0 / abs(cg(k, l, om_intr[:zj + 1], H_vals[:zj + 1])), z_vals[:zj + 1])
        else:
            prop_time = np.inf 
                            
        if prop_time < t0:
            if z_vals[zj] > 100.0:
                m_imag_prev = m_imag(k, l, om_intr[zj - 1], H_vals[zj - 1], T0_vals[zj - 1], d0_vals[zj - 1])
                m_imag_curr = m_imag(k, l, om_intr[zj], H_vals[zj], T0_vals[zj], d0_vals[zj])
                w_loss = w_loss_prev + dz * (m_imag_prev + m_imag_curr) / 2.0
            else:
                w_loss = 0.0

            # Define saturation spectra for this Fourier component
            w_sat_sqr = 2.7e-2 / (om_intr[zj] * m_sqr(k, l, om_intr[zj], H_vals[zj]) * kh**2)

            # integrate using free solution if there is no turning point or if the 
            # evaluation time is less than the trapped propagation time
            m_src = np.sqrt(abs(m_sqr(k, l, om_intr[0], H_vals[0])))
            m_curr = np.sqrt(abs(m_sqr(k, l, om_intr[zj], H_vals[zj])))
            d0_m_ratio = np.sqrt((d0_vals[0] / d0_vals[zj]) * (m_src / m_curr))

            if turning_ht_index > len(z_vals) or t0 < refl_time:
                m_prev = np.sqrt(abs(m_sqr(k, l, om_intr[zj - 1], H_vals[zj - 1])))
                w_phase = w_phase_prev + dz * (m_prev + m_curr) / 2.0
            
                w_spec = w0 * d0_m_ratio * (np.cos(w_phase) - 1.0j * np.sin(w_phase)) * np.exp(-w_loss)

                if abs(w_spec) > np.sqrt(abs(w_sat_sqr)):
                    w_spec = w_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))

                u_spec = - w_spec * (k / kh**2) * m_curr
                v_spec = - w_spec * (l / kh**2) * m_curr

            else:
                # Alternately, use the Airy function form in the case of trapped waves
                if zj < turning_ht_index:
                    m_vals = np.sqrt(abs(m_sqr(k, l, om_intr[zj:turning_ht_index], H_vals[zj:turning_ht_index])))
                    airy_arg = -((3.0 / 2.0) * trapz(m_vals, z_vals[zj:turning_ht_index]))**(2.0 / 3.0)
                elif zj > turning_ht_index:
                    # check that m(z) doesn't become real again at some altitude 
                    # above the first reflection point and stop integration if it does
                    m_sqr_vals = m_sqr(k, l, om_intr[turning_ht_index:zj], H_vals[turning_ht_index:zj])
                    if np.all(m_sqr_vals < 0.0):
                        j_temp = zj
                    else:
                        j_temp = turning_ht_index + (np.where(m_sqr_vals > 0.0)[0][0] - 1)

                    m_vals = np.sqrt(-m_sqr(k, l, om_intr[turning_ht_index:j_temp], H_vals[turning_ht_index:j_temp]))
                    airy_arg = ((3.0 / 2.0) * trapz(m_vals, z_vals[turning_ht_index:j_temp]))**(2.0 / 3.0)
                else:
                    airy_arg = 0.0

                refl_cnt = int(np.floor(t0 / (2.0 * refl_time)))
                refl_ph_shift = np.exp(-refl_cnt * (2.0 * refl_loss) - w_loss)
                refl_ph_shift = refl_ph_shift * (1.0 + np.sum(np.array([np.exp(-1.0j * (N - 1) * (2.0 * refl_phase - np.pi / 2.0)) for N in range(2, refl_cnt)])))

                airy_scaling0 = (-airy_arg)**0.25 * airy(airy_arg)[0] * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift
                airy_scaling1 = (-airy_arg)**0.25 * airy(airy_arg)[1] * np.exp(1.0j * np.pi / 4.0) * refl_ph_shift

                u_spec = -2.0j * np.sqrt(np.pi) * w0 * (k / kh**2) * d0_m_ratio * airy_scaling1
                v_spec = -2.0j * np.sqrt(np.pi) * w0 * (l / kh**2) * d0_m_ratio * airy_scaling1
                w_spec = -2.0j * np.sqrt(np.pi) * w0 * d0_m_ratio * airy_scaling0

                if abs(w_spec) > np.sqrt(abs(w_sat_sqr)):
                    u_spec = u_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))
                    v_spec = v_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))
                    w_spec = w_spec * (np.sqrt(w_sat_sqr) / abs(w_spec))

    return [u_spec, v_spec, w_spec, w_phase, w_loss]


# wrappers to use multiprocessing
def prep_wrapper(args):
    return prep(*args)

def prop_upward_wrapper(args):
    return prop_upward(*args)


def build_spec(atmo_specification, t0=4.0 * 3600.0, dx=4.0, dz=0.2, Nk=128, N_om=5, ref_lat=30.0, random_phase=False):
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
            Zonal (E/W) wind perturbations, du(x, y, z)
        dv_vals: 3darray
            Meridional (N/S) wind perturbations, dv(x, y, z)
        dw_vals: 3darray
            Vertical wind perturbations, dw(x, y, z)
            
    """

    # source height and vertical integration resolution
    z_src = 20.0

    # Set up atmosphere info
    z, Temp0, u0, v0, d0, _ = np.loadtxt(atmo_specification, unpack=True)

    T0_interp = interp1d(z, Temp0)
    u0_interp = interp1d(z, u0 * 1.0e-3)
    v0_interp = interp1d(z, v0 * 1.0e-3)
    d0_interp = interp1d(z, d0 * 1.0e3)

    d0_finite_diffs = np.array([(d0[min(n + 1, len(z) - 1)] - d0[max(n - 1, 0)]) / (z[min(n + 1, len(z) - 1)] - z[max(n - 1, 0)]) for n in range(len(z))])
    sc_ht = interp1d(z, -d0 / d0_finite_diffs)

    z_vals = np.arange(z_src, z[-1], dz)
    atmo_info = [z_vals, sc_ht(z_vals), u0_interp(z_vals), v0_interp(z_vals), T0_interp(z_vals), d0_interp(z_vals)]

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
    w_losses = np.zeros_like(w_spec, dtype=float)

    print("Running pre-analysis to identify turning heights and critical layers...")
    results = np.array([[[prep(k, l, om, atmo_info, om_vals[0], max(abs(k_vals))) for om in om_vals] for l in l_vals] for k in k_vals])
    for kn, ln in itertools.product(range(Nk), repeat=2):
        for om_n in range(N_om):
            crit_lyr_indices[kn][ln][om_n] = results[kn][ln][om_n][0]
            turn_ht_indices[kn][ln][om_n] = results[kn][ln][om_n][1]
            w0_vals[kn][ln][om_n] = results[kn][ln][om_n][2]
            refl_info[kn][ln][om_n] = results[kn][ln][om_n][3:]

    # Randomize initial phase and set initial solution
    if random_phase:
        w_phase[:,:,:,0] = np.random.rand((Nk, Nk, N_om)) * (2.0 * np.pi)
    w_spec[:,:,:,0] = w0_vals * (np.cos(w_phase[:,:,:,0]) - 1.0j * np.sin(w_phase[:,:,:,0]))   
    
    # integrate upward
    for j in range(1, len(z_vals)):
        if j % 5 == 0:
            print("Integrating upward.  Currently at z = " + str(np.round(z_vals[j], 2)) + " km...")
        results = np.array([[[prop_upward(k_vals[kn], l_vals[ln], om_vals[om_n], j, t0, crit_lyr_indices[kn][ln][om_n], turn_ht_indices[kn][ln][om_n], w0_vals[kn][ln][om_n], 
        w_phase[kn][ln][om_n][j - 1], w_losses[kn][ln][om_n][j - 1], atmo_info, refl_info[kn][ln][om_n]) for om_n in range(N_om)] for ln in range(Nk)] for kn in range(Nk)])
        
        for kn, ln in itertools.product(range(Nk), repeat=2):
            for om_n in range(N_om):
                u_spec[kn][ln][om_n][j] = results[kn][ln][om_n][0]
                v_spec[kn][ln][om_n][j] = results[kn][ln][om_n][1]
                w_spec[kn][ln][om_n][j] = results[kn][ln][om_n][2]

                w_phase[kn][ln][om_n][j] = np.real(results[kn][ln][om_n][3])
                w_losses[kn][ln][om_n][j] = np.real(results[kn][ln][om_n][4])

    du_vals = np.trapz(u_spec, om_vals, axis=2)
    du_vals = np.fft.ifft(du_vals, axis=0) * (Nk / dx)
    du_vals = np.fft.ifft(du_vals, axis=1) * (Nk / dx)

    dv_vals = np.trapz(v_spec, om_vals, axis=2)
    dv_vals = np.fft.ifft(dv_vals, axis=0) * (Nk / dx)
    dv_vals = np.fft.ifft(dv_vals, axis=1) * (Nk / dx)

    dw_vals = np.trapz(w_spec, om_vals, axis=2)
    dw_vals = np.fft.ifft(dw_vals, axis=0) * (Nk / dx)
    dw_vals = np.fft.ifft(dw_vals, axis=1) * (Nk / dx)

    return z_vals, np.real(du_vals), np.real(dv_vals), np.real(dw_vals)


