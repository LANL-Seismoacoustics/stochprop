.. _installation:

=====================================
Installation
=====================================

-----------------
Operating Systems
-----------------

stochprop has been installed on machines running newer versions of Linux or Apple OS X.  Installation on a Windows system has not been tested, but requires an Anaconda Python installation, so it should be reasonably straightforward.  Installation of propagation and signal analysis tools on such a system might be more challenging however. 

----------------------------------------
Dependencies
----------------------------------------


**Anaconda**

The installation of stochprop is ideally completed using pip through `Anaconda <https://docs.anaconda.com/free/anaconda/install/index.html>`_ to resolve and download the correct python libraries. If you don't currently have Anaconda installed on your system, you will need to do that first.


**Propagation Modeling Methods**

A subset of the stochprop methods require access to the  LANL `InfraGA/GeoAc <https://github.com/LANL-Seismoacoustics/infraGA>`_ ray tracing methods as well as the `NCPAprop <https://github.com/chetzer-ncpa/ncpaprop-release>`_ normal mode methods.  Many of the empirical orthogonal function (EOF) based atmospheric statistics and gravity wave pertorbation methods can be used without these propagation tools, but full usage of stochprop requires them.


**InfraPy Signal Analysis Methods**

The propagation models constructed in stochprop are intended for use in the Bayesian Infrasonic Source Localization (BISL) and Spectral Yield Estimation (SpYE)
methods in the LANL `InfraPy <https://github.com/LANL-Seismoacoustics/infrapy>`_ signal analysis software suite.  As with the InfraGA/GeoAc and NCPAprop linkages, many of the EOF-based atmospheric statistics methods
can be utilized without InfraPy, but full usage will require installation of InfraPy.

-----------------------------
stochprop Installation
-----------------------------

**Stand Alone Install**

A stand alone Anaconda environment can be created with the stochprop YML file,

.. code-block:: none

    conda env create -f stochprop_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on how to activate and deactivate the new environment:

To activate the environment, use:

.. code-block:: none

    conda activate stochprop_env

To deactivate an active environment, use

.. code-block:: none

    conda deactivate

**Installing via PyGS**

A number of LANL infrasound software tools have been developed and made available to the community through the `LANL Seismoacoustics github page <https://github.com/LANL-Seismoacoustics/infrapy>`_.  Collectively, these tools are referred to as the Python Geophysics Suite (PyGS).  As noted above, the methods in *stochprop* function as a linkage between propagation modeling methods in *InfraGA/GeoAc* and *NCPAprop* to localization and yield estimation methods in *InfraPy*.  Because of this, it's useful to install these various packages into a single Anaconda environment; however, package conflicts exist that make a full installation of these tools difficult.  Resolving these conflicts is an ongoing effort.

For now, there are minimal conflicts between *stochprop* and the ray tracing methods in *infraGA/GeoAc*.  One can install *stochprop* directly using pip by navigating to the base directory of the package (there will be a file there named 'setup.py').  Assuming *infraga* has been installed within a conda environment called infraga_env, it is possible to install *stochprop* in the same environment buy first creating a clone of the *infraga* environment,

.. code-block:: none

    conda create --name pygs_env --clone infraga_env

then install the stochprop methods into the new Python Geophysics Suite environment (pygs_env),

.. code-block:: none

    conda activate pygs_env
    pip install -e .

Once conflicts have been resolved, a single 'pygs.yml' file will be made available that installs the various LANL geophysics python tools into a single Anaconda environment.

**Testing stochprop**

Once the installation is complete, you can test the methods by running the command line interface help.  Firstly, activate either the :code:`stochprop_env` or :code:`pygs_env`, then run the :code:`--help` or :code:`-h` option for stochprop.

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
      plot   Visualization methods
      prop   Propagation model construction methods
      stats  Atmosphere statistics methods

Usage of the individual packages and sub-commands can be similarly displayed with the help flag (e.g., :code:`stochprop stats build-eofs -h`).
