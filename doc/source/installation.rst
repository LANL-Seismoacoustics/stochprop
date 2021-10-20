.. _installation:

=====================================
Installation
=====================================

-------------------------------------
Anaconda
-------------------------------------

The installation of stochprop is ideally completed using pip through Anaconda to resolve and download the correct python libraries. If you don't currently have anaconda installed
on your system, please do that first.  Anaconda can be downloaded from https://www.anaconda.com/distribution/.


----------------------------------------
Installing Dependencies
----------------------------------------

****************************************
Propagation Modeling Methods
****************************************

A subset of the stochprop methods require access to the  LANL InfraGA/GeoAc ray tracing methods as well as the NCPAprop normal mode methods.  Many of the 
empirical orthogonal function (EOF) based atmospheric statistics and gravity wave pertorbation methods can be used without these propagation tools, but full usage of stochprop requires them.

* InfraGA/GeoAc: https://github.com/LANL-Seismoacoustics/infraGA
* NCPAprop: https://github.com/chetzer-ncpa/ncpaprop

****************************************
InfraPy Signal Analysis Methods
****************************************

The propagation models constructed in stochprop are intended for use in the Bayesian Infrasonic Source Localization (BISL) and Spectral Yield Estimation (SpYE)
methods in the LANL InfraPy signal analysis software suite.  As with the InfraGA/GeoAc and NCPAprop linkages, many of the EOF-based atmospheric statistics methods
can be utilized without InfraPy, but full usage will require installation of InfraPy (https://github.com/LANL-Seismoacoustics/infrapy).

-------------------------------------
Installing stochprop
-------------------------------------

Once Anaconda is installed, you can install stochprop using pip by navigating to the base directory of the package (there will be a file there
named setup.py).  Assuming InfraPy has been installed within a conda environment called infrapy_env, it is recommended to install stochprop in the same environment using:

.. code-block:: none

    >> conda activate infrapy_env
    >> pip install -e .

Otherwise, a new conda environment should be created with the underlying dependencies and pip should be used to install there (work on this later):

.. code-block:: none

    >> conda env create -f stochprop_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on how to activate and deactivate the new environment:

To activate the environment, use:

.. code-block:: none

    >> conda activate stochprop_env

To deactivate an active environment, use

.. code-block:: none

    >> conda deactivate


-------------------------------------
Testing stochprop
-------------------------------------

Once the installation is complete, you can test the methods by navigating to the /examples directory located in the base directory, and running:

.. code-block:: none

    >> python eof_analysis.py
    >> python atmo_analysis.py

A set of propagation analyses are included, but require installation of infraGA/GeoAc and NCPAprop.  These analysis can be run to ensure linkages are
working between stochprop and the propagation libraries, but note that the simulation of propagation through even the example suite of atmosphere
takes a significant amount of time.

