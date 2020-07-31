.. _analysis:

=====================================
Stochastic Propagation Analysis
=====================================

* The atmospheric state at a given time and location is uncertain due to its dynamic and sparsely sampled nature
* Propagation effects for infrasonic signals must account for this uncertainty in order to properly quantify uncertainty in analysis results
* A methodology of constructing propagation statistics has been developed that identifies a suite of atmospheric states that characterize the possible space of scenarios, runs propagation simulations through each possible state, and builds statistical distributions for propagation effects

.. figure:: _static/_images/stochprop_fig1.jpg
    :width: 500px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    Stochastic propagation models are constructing using a suite of possible atmospheric states, propagation modeling applied to each, and a statistical model describing the variability in the resulting set of predicted effects
    
* The tools included here provide a framework for constructing such models as well as perform a number of other analyses related to atmospheric variability and uncertainty
    
_______________________________________
:ref:`eofs`
_______________________________________
* Empirical Orthogonal Functions (EOFs) provide a mathematical means of measuring variations in the atmospheric state
* Methods measure EOF statistics to reduce the number of atmospheric samples necessary to characterize the atmosphere at a given location during a specified time period


_________________________________
:ref:`sampling`
_________________________________
* EOFs can be used to fit a specified atmosphere by computing coefficients for each EOF
* Statistics of the coefficients for a suite of atmospheric states can be used to generate a set of characteristics samples
* Randomly generated EOF coefficients can be used to generate perturbations to an initial atmospheric specification and construct a suite of atmospheric states that fall within expected uncertainty


_________________________________
:ref:`propagation`
_________________________________
* InfraGA/GeoAc ray tracing analysis can be applied to a suite of atmospheric states to predict geometric propagation characteristics such as arrival location, travel time, and direction of arrival needed to estimate the source location
* NCPAprop modal simulations can be applied to a suite of atmospheric states to predict finite frequency transmission loss needed to characterize the infrasonic source


    .. toctree::
        :hidden:

        eofs
        sampling
        propagation