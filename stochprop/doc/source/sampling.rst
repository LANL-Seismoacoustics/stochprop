.. _sampling:

===============================================
Atmospheric Fitting, Sampling, and Perturbation
===============================================

The Empirical Orthogonal Functions (EOFs) constructed using a suite of atmospheric specifications can be utilized in a number of different analyses of the atmospheric state.  In general, an atmospheric state can be constructed by defining a reference atmosphere, :math:`\vec{B}_0`, and a set of coefficients, :math:`\mathcal{C}_n`,

.. math::
	\hat{\vec{B}} = \vec{B}_0 + \sum_n{ \mathcal{C}_n \vec{\varepsilon}_n}
		

***********************************************
Fitting an Atmospheric Specification using EOFs
***********************************************

In the case that a specific state, :math:`\vec{B}`, is known, it can be approximated using the EOF basis functions by using the mean state pulled from the original SVD analysis, :math:`\bar{\vec{A}}`, and coefficients defined by projecting the atmospheric state difference from this mean onto each EOF,

.. math::
	\hat{\vec{B}} = \bar{\vec{A}} + \sum_n{ \mathcal{C}_n^{(\vec{B})} \vec{\varepsilon}_n}, \quad \quad \mathcal{C}_n^{(\vec{B})} = \vec{\varepsilon}_n \cdot \left( \vec{B} - \bar{\vec{A}} \right)

These coefficient calculations and construction of a new atmospheric specification can be completed using :code:`stochprop.eofs.fit_atmo` with the path to specific atmospheric state, a set of EOFs, and a specified number of coefficients to compute,


.. code:: Python

    prof_path = "profs/g2stxt_2011010118_39.1026_-84.5123.dat"
    eofs_path = "eofs/example"
	
    eofs.fit_atmo(prof_path, eofs_path, "eof_fit-N=30.met", eof_cnt=30)


This analysis is useful to determine how many coefficients are needed to accurately reproduce an atmospheric state from a set of EOFs.  Such an analysis is shown below for varying number of coefficients and convergence is found at 50 - 60 terms.

.. figure:: _static/_images/US_NE-fits.png
    :width: 700px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    Accuracy of fitting a specific atmospheric state (black) using varying numbers of EOF coefficients (red) shows convergence for approximately 50 - 60 terms in the summation

************************************************************
Sampling Specifications using EOF Coefficient Distributions
************************************************************

Samples can be generated that are representative of a given coefficient distributions as discussed in Blom et al. (2023) to re-sample a given vector space and provide some data reduction needed for propagation simulation campaigns.  Once an EOF basis function set has been defined, the coefficients defining the projection of the EOFs onto the original set is computed.  Given the set of coefficients for the set, a kernel density estimate (KDE) is applied to sample the coefficient values and produce representative atmospheric state vectors.  From the EOF result, :code:`stochprop.eofs.compute_coeffs` can be used to build coefficient values and :code:`stochprop.eofs.sample_atmo` used to generate samples.  

The mathematics involved are detailed in Blom et al. (2023),

.. math::
   	\hat{\vec{A}}_m = \bar{\vec{A}} + \sum_n{ \hat{\mathcal{C}}_n^{(m)} \vec{\varepsilon}_n}, \quad \quad \hat{\mathcal{C}}_n^{(m)}  \text{ from } \hat{\mathcal{P}}_n^{(A)} \left( \mathcal{C} \right) = \text{KDE} \left[ \mathcal{C}_n^{(\vec{A}_1)}, \mathcal{C}_n^{(\vec{A}_2)}, \ldots,  \mathcal{C}_n^{(\vec{A}_N)}  \right]
   	
Numerically, this can be accomplished via the :code:`stochprop.eofs.sample_atmo` function,

.. code:: Python

    coeffs = np.load("coeffs/example_winter-coeffs.npy")
    eofs.sample_atmo(coeffs, eofs_path, "samples/winter/example-winter", prof_cnt=25)

This analysis can be completed for each identified season to generate a suite of atmospheric specifications representative of the season as shown in the figure below.  This can often provide a significant amount of data reduction for propagation studies as multiple years of specifications (numbering in the 100's or 1,000's) can be used to construct a representative set of 10's of atmospheres that characterize the time period of interest as in the figure below.

.. figure:: _static/_images/US_RM-samples.png
    :width: 500px
    :align: center
    :alt: alternate text
    :figclass: align-center

    Samples for seasonal trends in the western US show the change in directionality of the stratospheric waveguide in summer and winter.


****************************************************
Perturbing Specifications to Account for Uncertainty
****************************************************

In most infrasonic analysis, propagation analysis through a specification for the approximate time and location of an event doesn't produce the exact arrivals observed due to the dynamic and sparsely sampled nature of the atmosphere. Because of this, it is useful to apply random perturbations to the estimated atmospheric state covering some confidence level and consider propagation through the entire suite of "possible" states.  In such a case, the reference atmosphere, :math:`\vec{A}_0` defines the initial states and coefficients are randomly generated from a normal distribution,

.. math::

    \tilde{\vec{A}}_m = \vec{A}_0 + \mathcal{W} \left( \varsigma \right) \sum_n{ w_n \mathcal{C}_n^{(m)} \vec{\varepsilon}_n}, \quad \mathcal{C}_n^{(m)} \text{ from } \mathcal{N} \left(0, 1 \right), 

the inter-EOF weighting, :math:`w_n`, and the overall perturbation scaling, :math:`\mathcal{W}`, along with application to localization and characterization analyses are ongoing areas of R&D...
    

.. 
    COMMENTED OUT SECTION

    where the weighting of each EOF is 

    .. math::

        w_n \left( \gamma, \eta \right) = \sigma_n^{\gamma} \; \bar{z}_n^{\eta}


     and this process is likely to evolve as methods are evaluated.  The current implementation uses the singular values, :math:`\sigma_n` and the mean altitude of the perturbation, :math:`\bar{z}_n = \sum_j{z_j \varepsilon_n \left( z_j \right)` in order to avoid rapidly oscillating EOFs from contributing too much noise and to focus perturbations at higher altitudes where uncertainties are larger, respectively.  The exponential coefficients have default values of :math:`\gamma = 0.25` and :math:`\eta=2`, but can be modified in the function call (again, these default values are an area of ongoing R&D).  Once computed, the set of perturbations is scaled to match the specified standard deviation, :math:`\mathcal{W} \left( \varsigma \right)` (averaged over the entire set of altitudes).  This perturbation analysis can be completed using :code:`stochprop.eofs.perturb_atmo` with a specified starting atmosphere, set of EOFs, output path, uncertainty measure in meters-per-second, and number of samples needed,

    .. code:: Python

        eofs.perturb_atmo(prof_path, eofs_path, "eof_perturb", uncertainty=5.0, sample_cnt=10)

    * The below figure shows a sampling of results using uncertainties of 5.0, 10.0, and 15.0 meters-per-second.  The black curve is input as the estimated atmospheric state and the red curves are generated by the perturbations.

    .. figure:: _static/_images/atmo_perturb.png
        :width: 500px
        :align: center
        :alt: alternate text
        :figclass: align-center

        Perturbations to a reference atmospheric state can be computed using randomly generated coefficients for a suite of EOFs with specified standard deviation


**********************
Command Line interface
**********************

The command line interface (CLI) for the fitting and sampling via EOFs can be utilized more easily as summarized in the :ref:`quickstart`.  Fitting of an atmospheric specification using EOFs can be done through a visualization method and the sampling methods are included in the statistics methods, :code:`stochprop stats`. Usage info can be summarized using the :code:`--help` (or :code:`-h`) option:

    .. code-block:: none 

        Usage: stochprop plot eof-fit [OPTIONS]

          stochprop plot eof-fit
          -----------------------
  
          Example Usage:
               stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 25

        Options:
          --atmo-file TEXT    Reference atmospheric specification (required)
          --eofs-path TEXT    EOF output path and prefix (required)
          --eof-cnt INTEGER   Number of EOFs to visualize (default: 5)
          --output-file TEXT  Output file to save fit (optional)
          -h, --help          Show this message and exit.

Fitting an atmospheric specification can be done via the visualization methods,

    .. code-block:: none 

        stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 25

The coefficient calculations and sampling can be accomplished using the :code:`stochprop stats` methods as summarized in the :ref:`quickstart`,

    .. code-block:: none 

        stochprop stats eof-coeffs --atmo-dir profs/ --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --week-selection '38:52,1:15'

        stochprop stats sample-eofs --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --sample-path samples/winter/example_winter --sample-cnt 50


