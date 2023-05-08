.. stochprop documentation master file, created by
   sphinx-quickstart on Wed July 22 08:00:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index:

=====================================
stochprop
=====================================

Simulations of infrasonic propagation in the atmosphere typically utilize a single atmospheric specification describing the acoustic sound speed, ambient winds, and density as a function of altitude.  Due to the dynamic and sparsely sampled nature of the atmosphere, there is a notable amount of uncertainty in the atmospheric state at a given location and time so that a more robust analysis of infrasonic propagation requires inclusion of this uncertainty.  This Python library, stochprop, has been implemented using methods developed jointly by infrasound scientists at Los Alamos National Laboratory (LANL) and the University of Mississippi's National Center for Physical Acoustics (NCPA).  This software library includes methods to quantify variability in the atmospheric state, identify typical seasonal variability in the atmospheric state and generate suites of representative atmospheric states during a given season, as well as perform uncertainty analysis on a specified atmospheric state given some level of uncertainty.  These methods have been designed to interface between propagation modeling capabilities such as InfraGA/GeoAc and NCPAprop and signal analysis methods in the LANL InfraPy tool.  

_____________________________________
Contents
_____________________________________
.. automodule:: stochprop
    :members:

.. toctree::
    :maxdepth: 3
    :titlesonly:

    authorship
    installation
    quickstart
    analysis
    stochprop
    references
