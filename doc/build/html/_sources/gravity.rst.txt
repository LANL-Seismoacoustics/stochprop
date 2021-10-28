.. _gravity:

==========================
Gravity Wave Perturbations
==========================
* Atmospheric specifications (e.g., G2S) available for a given location and time are averaged over some spatial and temporal scale so that sub-grid scale fluctuations can be estimated stochastically and applied in order to construct a suite of possible atmospheric states and the primary source of sub-grid fluctuations in the atmosphere is that of bouyancy or gravity waves.

* The methods available in :code:`stochprop` are based on the vertical ray tracing approach detailed in Drob et al. (2013) and are summarized below for reference.

++++++++++++++++++++++++++++++++++++++++++++
Freely Propagation and Trapped Gravity Waves
++++++++++++++++++++++++++++++++++++++++++++
* The gravity wave dynamics are governed by a pair relations describing the disperion and wave action conservation.  The dispersion relation describing the vertical wavenumber can be expressed as,

.. math::
	m^2 \left( k, l, \omega, z \right) = \frac{k_h^2}{\hat{\omega}^2} \left( N^2 - \hat{\omega}^2 \right) + \frac{1}{4H^2}
 
* In this relation :math:`k` and :math:`l` are the zonal and meridional wave numbers, :math:`k_h^2 = \sqrt{k^2 + l^2}` is the combined horizontal wavenumber, :math:`H = - \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}` is the density scale height, :math:`N^2 = \sqrt{-\frac{g}{\rho_0} \frac{\partial \rho_0}{\partial z}} = \sqrt{\frac{g}{H}}` is the atmospheric bouyancy frequency, and :math:`\hat{\omega}` is the intrinsic angular frequency (relative to the moving air) that is defined relative to the absolute angular frequency (relative to the ground), :math:`\omega`,

.. math::
	\hat{\omega} = \omega - k u_0 \left( z \right) - l v_0 \left( z \right)

* This dispersion relation can alternately be solved for :math:`\hat{\omega}` and used to define the vertical group velocity,

.. math::
	\hat{\omega} = \frac{k_h N \left( z \right)}{\sqrt{ k_h^2 + m^2 \left( z \right) + \frac{1}{4 H^2 \left( z \right)}}} \quad \rightarrow \quad 
	c_g \left(k, l, \omega, z \right) = \frac{\partial \hat{\omega}}{\partial m} = -\frac{m k_h N}{\left( k_h^2 + m^2 + \frac{1}{4 H^2} \right)^\frac{3}{2}} 

* The conservation of wave action leads to a condition on the vertical velocity perturbation spectrum that can be used to define a freely propagating solution,

.. math::
	\rho_0 m \left| \hat{w} \right|^2 = \text{constant} \; \rightarrow \;
	\hat{w} \left( k, l, \omega, z \right) = \hat{w}_0 e^{i \varphi_0} \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} e^{i \int_{z_0}^z{m \left( z^\prime \right) dz^\prime}}

* The above relation is valid in the case that :math:`m \left( k, l, \omega, z \right)` remains real through the integration upward in the exponential.  In the case that an altitude exists for which the vertical wavenumber becomes imaginary, the gravity wave energy reflects from this turning height and the above relation is not valid.  Instead, the solution is expressed in the form,

	.. math::
 		\hat{w} \left( k, l, \omega, z \right) = 2 i \sqrt{\pi} \hat{w}_0 \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} \times \left( - r \right)^\frac{1}{4} \text{Ai} \left( r \right) e^{-i \frac{\pi}{4}} S_n

	* The Airy function argument in the above is defined uniquely above and below the turning height :math:`z_t`,

	.. math::
		r = \left\{ \begin{matrix} - \left( \frac{3}{2} \int_z^{z_t} \left| m \left( z^\prime \right) \right| dz^\prime \right)^\frac{2}{3} & z < z_t \\ \left( \frac{3}{2} \int_{z_t}^z \left| m \left( z^\prime \right) \right| dz^\prime \right)^\frac{2}{3} & z > z_t \end{matrix} \right.

	* The reflection phase factor, :math:`S_n`, accounts for the caustic phase shifts from the :math:`n` reflections from the turning height,

	.. math::
   		S_n = \sum_{j = 1}^n{e^{i \left( j -1 \right) \left(2 \Phi - \frac{\pi}{2} \right)}}, \quad \Phi = \int_0^{z_t} m \left( z^\prime \right) d z^\prime

* The vertical velocity spectra defined here can be related to the horizontal velocity for the freely propagating and trapped scenarios through derivatives of the vertical velocity spectrum,

.. math::
	\hat{u}^\text{(free)} = - \frac{k m}{k_h^2} \hat{w}, \quad
   	\hat{u}^\text{(trapped)} = \frac{2 i \hat{w}_0 }{\sqrt{\pi}}\frac{k}{k_h^2} \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} \times \left( - r \right)^\frac{1}{4} \text{Ai}^\prime \left( r \right) e^{-i \frac{\pi}{4}} S_n

.. math::
	\hat{v}^\text{(free)} = - \frac{l m}{k_h^2} \hat{w}, \quad
   	\hat{v}^\text{(trapped)} = \frac{2 i \hat{w}_0 }{\sqrt{\pi}}\frac{l}{k_h^2} \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} \times \left( - r \right)^\frac{1}{4} \text{Ai}^\prime \left( r \right) e^{-i \frac{\pi}{4}} S_n

* Finally, once computed for the entire atmosphere, the spatial and temporal domain forms can be computed by an inverse Fourier transform,

.. math::
	w \left( x, y, z, t \right) = \int{e^{-i \omega t} \left( \iint{ \hat{w} \left( k, l, \omega, z \right) e^{i \left( kx + ly \right)} dk \, dl} \right) d \omega}



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Damping, Source and Saturation Spectra, and Critical Layers
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* At altitudes above about 100 km, gravity wave damping by molecular viscosity and thermal diffusion becomes increasingly important.  Following the methods developed by Drob et al. (2013), for altitudes above 100 km, an imaginary vertical wave number term can be defined, :math:`m \rightarrow m + m_i`, where,

	.. math::
		m_i \left(k, l, \omega, z \right) = -\nu \frac{m^3}{\hat{\omega}}, \quad \nu = 3.563 \times 10^{-7} \frac{T_0^{\, 0.69}}{\rho_0}

	* This produces a damping factor for the freely propagating solution that is integrated upward along with the phase,

	.. math::
		\hat{w} \left( k, l, \omega, z \right) = \hat{w}_0 e^{i \varphi_0} \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} e^{i \int_{z_0}^z{m \left( z^\prime \right) dz^\prime}} e^{-\int_{z_0}^{z}{m_i \left( z^\prime \right) dz^\prime}}

	* In the trapped solution, the reflection phase shift includes losses for each pass up to the turning height and back,

	.. math::
   		S_n = e^{-2 n \Psi} \sum_{j = 1}^n{e^{i \left( j -1 \right) \left(2 \Phi - \frac{\pi}{2} \right)}}, \quad \Phi = \int_0^{z_t} m \left( z^\prime \right) d z^\prime, \quad \Psi = \int_0^{z_t} m_i \left( z^\prime \right) d z^\prime,

	* Note that in both of these forms if :math:`z_t` is below 100 km there is no loss calculated and when it is above this altitude the losses are only computed from 100 km up to the turning height.

* The source spectra defined by Warner & McIntyre (1996) specifies the wave energy density for a source at 20 km altitude (note: :math:`\hat{\omega}` exponential corrected in publication errata),

	.. math::
		\mathcal{E}_\text{src} \left(m, \hat{\omega} \right) = 1.35 \times 10^{-2} \frac{m}{m_*^4 + m^4} \frac{N^2}{\hat{\omega}^\frac{5}{3}} \Omega, \quad \Omega = \frac{\hat{\omega}_\text{min}^\frac{2}{3}}{1 - \left( \frac{\hat{\omega}_\text{min}}{N} \right)^\frac{2}{3}}, \quad m_* = \frac{2 \pi}{2.5 \text{km}}
	
	* The wave energy density can be expressed in terms of spectral coordiantes using :math:`\mathcal{E} \left( k, l, \omega \right) = \mathcal{E} \left( m, \hat{\omega} \right) \frac{m}{k_h^2}` which can then be related to the vertical velocity spectrum producing the initial condition for starting the calculation, 

	.. math::
		\mathcal{E} \left(k, l, \omega \right) = \frac{1}{2} \frac{N^2}{\hat{\omega}^2} \left| \hat{w}_0 \right|^2 \quad \rightarrow \quad \left| \hat{w}_0 \right|^2 = 2.7 \times 10^{-2} \frac{m^2}{m^4_* + m^4}  \frac{\hat{\omega}^\frac{1}{3}}{k_h^2} \Omega.

* Gravity wave breaking in the atmosphere is included in analysis via a saturation limit following work by Warner & McIntyre (1996) where the spectral coordinate saturation spectrum is (note: the exponential for :math:`\hat{\omega}` is again corrected in publication errata),

	.. math::
		\mathcal{E}_\text{sat} \left(k, l, \omega \right) = 1.35 \times 10^{-2} \frac{N^2}{\hat{\omega}^\frac{5}{3} m^3}

	* Again using the relation between wave energy density and vertical velocity spectrum, this produces,

	.. math::
		\left| \hat{w}_\text{sat} \right|^2 = 2.7 \times 10^{-2} \frac{\hat{\omega}^\frac{1}{3}}{m^2 k_h^2}.
		
* Lastly, from the above definition for the vertical group velocity, :math:`c_g`, it is possible to have altitudes for which :math:`\hat{\omega} \rightarrow 0` and :math:`c_g` similarly goes to zero.  In such a location the wave eenrgy density becomes infinite; however, the propagation time to such an altitude is infinite and it is therefore considered a "critical layer" because the ray path will never reach the layer.  In computing gravity wave spectra using the methods here, a finite propagation time of several hours is defined and used to prevent inclusion of the critical layer effects and also quantify the number of reflections for trapped components.  Drob et al. included a damping factor for altitudes with propagation times more than 3 hours and that attenuation is included here as well.

++++++++++++++++++++++++++++++++++++++++
Gravity Wave implementation in stochprop
++++++++++++++++++++++++++++++++++++++++

* The implementation of the gravity wave analysis partially follows that summarized by Drob et al. (2013) and is sumamrized here

  * Atmospheric information is constructed from a provided atmospheric specification:

    * Interpolations of the ambient horizontal winds, :math:`u_0 \left( z \right)` and :math:`v_0 \left( z \right)`, density, :math:`\rho_0 \left( Z \right)`, and temperature, :math:`T_0 \left( z \right)` are defined.  

    * The density scale height, :math:`H \left( z \right) = - \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}`, is computed using finite differences of the ambient density.  
  
    * Atmospheric fields are then re-sampled on a set of altitudes with :math:`dz = 200` meters.
  
  * A grid of :math:`k`, :math:`l`, and :math:`\omega` values are defined:

    * The horizontal resolution, :math:`dx`, is set to 4 meters following Drob et al. (2013) with :math:`N_k = 128` (both of these quantities can be modified by the user, but default to the values from Drob et al.)

    * Five frequency values are defined for analysis covering a frequency band from :math:`\omega_\text{min} = 2 f_\text{Cor}` to :math:`\omega_\text{max} = \frac{N}{\sqrt{5}}` where :math:`f_\text{Cor}` is the Coriolis frequency,

     .. math::
	     f_\text{Cor} = 7.292 \times 10^{-5} \frac{\text{rad}}{\text{s}} \times \sin \left( \text{latitude} \right)

  * Unlike the implementation by Drob et al. (2013), the implementation in :code:`stochprop.gravity_waves` integrates each Fourier component individually and distributes calculations via :code:`multiprocessing`.  In preliminary evaluations, the implementation has able to compute the full gravity wave field in less than 2 minutes when using 10 CPUs.

  * For each Fourier component combination, :math:`k, l, \omega`, several checks are made and pre-analysis completed:

    * Those Fourier components for which :math:`k_h > k_\text{max}` are masked out of the calculation as well as those for which :math:`C = \frac{N}{m} > 90 \frac{\text{m}}{\text{s}}` 

    * Critical layers at which :math:`\hat{\omega} \rightarrow 0` are identified

    * Turning heights at which :math:`m^2 \left( z_t \right) \rightarrow 0` are identified and for each such Fourier combination the propagation time, phase shift, and attenuation factors are computed.

  * The relations above for :math:`\hat{w} \left( k, l, \omega, z \right)` are used to define the solution below the source height and to integrate the solution from the source height to the upper limit of the atmosphere using either the free or trapped form depending on whether a turning point exists

    * At each altitude, the propagation time to that point is computed and compared with a user specified propagation time that defaults to 4 hours to determine whether energy has reached that altitude.  
  
    * Similary, the number of reflections used in computing the trapped solution phase shift if determined by the ratio of the propagation time of the trapped solution with the specified time.

  * The gravity wave field in the spatial and time domain are obtained by inverting the spatial components using :code:`numpy.fft.ifft` on the appropriate axes and the :math:`\omega` integration is simplified by taking \(t=0\) in the solution which convert this Fourier inversion to a simple integration.

.. math::
	w \left( x, y, z, 0 \right) = \int{\left( \iint{ \hat{w} \left( k, l, \omega, z \right) e^{i \left( kx + ly \right)} dk \, dl} \right) d \omega}


* Use of the methods is summarized in the below example:

.. code-block:: python

	from stochprop import gravity_waves

	if __name__ == '__main__':
		atmo_spec = "profs/01/g2stxt_2010010100_39.7393_-104.9900.dat"
		output_path = "gw_perturb"

		t0 = 6.0 * 3600.0
		dx, Nk = 2.0, 128

		# Run gravity wave calculation
		gravity_waves.perturb_atmo(atmo_spec, "gw_perturb", t0=t0, dx=dx, Nk=Nk)

* Command line interface (CLI) options are also included for generating perturbed atmospheric specifications,

.. code-block:: bash

	stochprop gravity-waves --atmo-file profs/01/g2stxt_2010010100_39.7393_-104.9900.dat --out gw_perturb --dx 2.0 --nk 128 --cpu-cnt 8



