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


-------------------------------------------
Python Geophysics Suite (PyGS) Installation
-------------------------------------------

Infrasound software tools developed by LANL SMEs have become increasing coupled in usage so that having them in a common Python environment is useful.  An in-development Python Geophysics Suite (PyGS) YML file is included in the InfraPy repository that will build an environment and install InfraPy, infraGA/GeoAc, and stochprop from GitHub.  It can be run using the same syntax as above,

.. code-block:: bash

    >> conda env create -f pygs_env.yml

All dependencies will be installed and the LANL Python libraries pulled from GitHub to complete the environment.  To finish setting up, activate the environment and compile the infraGA/GeoAc software,

.. code-block:: bash

    >> conda activate pygs
    >> infraga compile 

---------------------------------------------------------
Python Geophysics Suite (PyGS) Installation - Dev Version
---------------------------------------------------------

Because the PyGS YML file installs via GitHub cloning, it doesn't copy the examples/ directories from the various libraries for demonstration and also doesn't leave the source code easily accessible for any de-bugging or customization.  A separate developer version is also included that requires a few more steps.  Build an instance of the environment with just stochprop included using the included YML file,

.. code-block:: bash

    >> conda env create -f pygs-dev_env.yml

Next, clone the other repositories if you don't have them,

.. code-block:: bash

    >> git clone https://github.com/LANL-Seismoacoustics/infraga.git
    >> git clone https://github.com/LANL-Seismoacoustics/infrapy.git

If you have SSH keys set up for GitHub, you can alternately clone as,

.. code-block:: bash
	
    >> git clone git@github.com:LANL-Seismoacoustics/infraga.git
    >> git clone git@github.com:LANL-Seismoacoustics/infrapy.git

Once the PyGS development environment is built, activate it using :code:`conda activate pygs_dev` and then use pip with the :code:`-e` flag to install infraGA/GeoAc and InfraPy and compile the infraGA/GeoAc ray tracing methods,

.. code-block:: bash

    >> cd /path/to/infraga
    >> pip install -e .

    >> infraga compile 

    >> cd /path/to/infrapy
    >> pip install -e .

This installation will clone the example directories with all relevant data and also allow you to interact with other :code:`git` branches for customization.

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
