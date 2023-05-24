.. _installation:

=====================================
Installation
=====================================

-------------------------------------
Anaconda
-------------------------------------

The installation of stochprop is ideally completed using pip through `Anaconda <https://www.anaconda.com/distribution/>`_ to resolve and download the correct python libraries. If you don't currently have anaconda installed on your system, please do that first.


----------------------------------------
Installing Dependencies
----------------------------------------

****************************************
Propagation Modeling Methods
****************************************

A subset of the stochprop methods require access to the  LANL `InfraGA/GeoAc <https://github.com/LANL-Seismoacoustics/infraGA>`_ ray tracing methods as well as the `NCPAprop <https://github.com/chetzer-ncpa/ncpaprop>`_ normal mode methods.  Many of the empirical orthogonal function (EOF) based atmospheric statistics and gravity wave pertorbation methods can be used without these propagation tools, but full usage of stochprop requires them.


****************************************
InfraPy Signal Analysis Methods
****************************************

The propagation models constructed in stochprop are intended for use in the Bayesian Infrasonic Source Localization (BISL) and Spectral Yield Estimation (SpYE)
methods in the LANL `InfraPy <https://github.com/LANL-Seismoacoustics/infrapy>`_ signal analysis software suite.  As with the InfraGA/GeoAc and NCPAprop linkages, many of the EOF-based atmospheric statistics methods
can be utilized without InfraPy, but full usage will require installation of InfraPy.

-------------------------------------
Installing stochprop (stand alone)
-------------------------------------

A stand alone conda environment can be created with only the stochprop using pip:

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
Installing stochprop (PyGS)
-------------------------------------

*Add some explanation of the Python Geophysics Suite (PyGS)*

Once Anaconda is installed, you can install stochprop using pip by navigating to the base directory of the package (there will be a file there
named setup.py).  Assuming InfraPy has been installed within a conda environment called infrapy_env, it is possible to install stochprop in the same environment buy first creating a clone of the infrapy environment,


.. code-block:: none

    conda create --name pygs_env --clone infrapy_env

then install the stochprop methods into the new Python Geophysics Suite environment (pygs_env),

.. code-block:: none

    >> conda activate pygs_env
    >> pip install -e .

A similar installation can be done for the `InfraGA/GeoAc <https://github.com/LANL-Seismoacoustics/infraGA>`_ anaconda environment so that all methods are available in a single work space.  Note: this combined installation is still in development and there are sone package conflicts that may break as the various libraries advance.

-------------------------------------
Testing stochprop
-------------------------------------

Once the installation is complete, you can test the methods by running the command line interface help.  Firstly, activate either the :code:`stochprop_env` or :code:`pygs_env`, then run the :code:`--help` option for stochprop.

.. code-block:: none

    stochprop --help

This command will show the general usage of the stochprop package:

.. code-block:: none

    Usage: stochprop [OPTIONS] COMMAND [ARGS]...

      stochprop
      ---------

      Python-based tools for quantifying infrasonic propagation uncertainty via
      stochastic analyses

    Options:
      -h, --help  Show this message and exit.

    Commands:
      eof      Empirical Orthogonal Function (EOF) methods
      perturb  Atmospheric specification perturbing methods
      prop     Propagation model construction methods

Usage of the indivitual packages and sub-commands can be similarly displayed with the :code:`--help` flag (e.g., :code:`stochprop eof build --help`).
