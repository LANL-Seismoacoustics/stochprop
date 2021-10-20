.. _gravity:

==========================
Gravity Wave Perturbations
==========================
* Atmospheric specifications (e.g., G2S) available for a given location and time are averaged over some spatial and temporal scale so that sub-grid scale fluctuations can be estimated stochastically and applied in order to construct a suite of possible atmospheric states and the primary source of sub-grid fluctuations in the atmosphere is that of bouyancy or gravity waves.

* The methods available in :code:`stochprop` are based on the vertical ray tracing approach detailed in Drob et al. (2013) and are summarized below for reference.

++++++++++++++++++++++++++++++++++++++++++++
Freely Propagation and Trapped Gravity Waves
++++++++++++++++++++++++++++++++++++++++++++
* The gravity wave dynamics are governed by a pair relations describing the disperion and wave action conservation.  The dispersion relation can be expressed as,

.. math::
	m^2 \left( k, l, \omega, z \right) = \frac{k_h^2}{\hat{\omega}^2} \left( N^2 - \hat{\omega}^2 \right) + \frac{1}{4H^2}
 
* In this relation :math:`k` and :math:`l` are the zonal and meridional wave numbers, :math:`k_h^2 = \sqrt{k^2 + l^2}` is the combined horizontal wavenumber, :math:`H = - \rho_0 \times \left( \frac{\partial \rho_0}{\partial z} \right)^{-1}` is the scale height, :math:`N^2 = \sqrt{-\frac{g}{\rho_0} \frac{\partial \rho_0}{\partial z}} = \sqrt{\frac{g}{H}}` is the atmospheric bouyancy frequency, and :math:`\hat{\omega}` is the intrinsic angular frequency (relative to the moving air) that is defined relative to the absolute angular frequency (relative to the ground), :math:`\omega`,

.. math::
	\hat{\omega} = \omega - k u_0 \left( z \right) - l v_0 \left( z \right)

* This dispersion relation can alternately be solved for :math:`\hat{\omega}` and used to define the vertical group velocity,

.. math::
	\hat{\omega} = \frac{k_h N \left( z \right)}{\sqrt{ k_h^2 + m^2 \left( z \right) + \frac{1}{4 H^2 \left( z \right)}}} \quad \rightarrow \quad 
	c_g = \frac{\partial \hat{\omega}}{\partial m} = -\frac{m k_h N}{\left( k_h^2 + m^2 + \frac{1}{4 H^2} \right)^\frac{3}{2}} 

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

* Lastly, the vertical velocity spectra defined here can be related to the horizontal velocity for the freely propagating and trapped scenarios,

.. math::
	\hat{u}^\text{(free)} = - \frac{k m}{k_h^2} \hat{w}, \quad
   	\hat{u}^\text{(trapped)} = \frac{2 i \hat{w}_0 }{\sqrt{\pi}}\frac{k}{k_h^2} \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} \times \left( - r \right)^\frac{1}{4} \text{Ai}^\prime \left( r \right) e^{-i \frac{\pi}{4}} S_n

.. math::
	\hat{v}^\text{(free)} = - \frac{l m}{k_h^2} \hat{w}, \quad
   	\hat{v}^\text{(trapped)} = \frac{2 i \hat{w}_0 }{\sqrt{\pi}}\frac{l}{k_h^2} \sqrt{ \frac{\rho_0 \left( z_0 \right)}{\rho_0 \left( z \right)} \frac{m \left( z_0 \right)}{m \left( z \right)}} \times \left( - r \right)^\frac{1}{4} \text{Ai}^\prime \left( r \right) e^{-i \frac{\pi}{4}} S_n


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

* The source spectra defined by Warner & McIntyre (1996) specifies the wave energy density as,

	.. math::
		\mathcal{E}_\text{src} \left(k, l, \omega \right) = 1.35 \times 10^{-2} \frac{m}{m_*^4 + m^4} \frac{N^2}{\hat{\omega}^3} \Omega, \quad \Omega = \frac{\hat{\omega}_\text{min}^\frac{2}{3}}{1 - \left( \frac{\hat{\omega}_\text{min}}{N} \right)^\frac{2}{3}}, \quad m_* = \frac{2 \pi}{2.5 \text{km}}
	
	* The wave energy density can be related to the vertical velocity spectrum producing the initial condition for starting the calculation, 

	.. math::
		\mathcal{E} \left(k, l, \omega \right) = \frac{1}{2} \frac{N^2}{\hat{\omega}^2} \left| \hat{w}_0 \right|^2 \quad \rightarrow \quad \left| \hat{w}_0 \right|^2 = 2.7 \times 10^{-2} \frac{m^2}{m^4_* + m^4}  \frac{\Omega}{\hat{\omega} k_h^2}.

* Gravity wave breaking in the atmosphere is included in analysis via a saturation limit following work by Warner & McIntyre (1996) where the spectral coordinate saturation spectrum is,

	.. math::
		\mathcal{E}_\text{sat} \left(k, l, \omega \right) = 1.35 \times 10^{-2} \frac{N^2}{\hat{\omega}^3 m^3}

	* Again using the relation between wave energy density and vertical velocity spectrum, this produces,

	.. math::
		\left| \hat{w}_\text{sat} \right|^2 = 2.7 \times 10^{-2} \frac{1}{\hat{\omega} m^2 k_h^2}.
		
* Lastly, from the above definition for the vertical group velocity, :math:`c_g`, it is possible to have altitudes for which :math:`\hat{\omega} \rightarrow 0` and :math:`c_g` similarly goes to zero.  In such a location the wave eenrgy density becomes infinite; however, the propagation time to such an altitude is infinite and it is therefore considered a "critical layer" because the ray path will never reach the layer.  In computing gravity wave spectra using the methods here, a propagation time of several hours is defined and used to prevent inclusion of the critical layer effects and also quantify the number of reflections for trapped components.

++++++++++++++++++++++++++++++++++++++++
Gravity Wave implementation in stochprop
++++++++++++++++++++++++++++++++++++++++

* The implementation here...
