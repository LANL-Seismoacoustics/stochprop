# -*- coding: utf-8 -*-
try:
    import setuptools
except ImportError:
    pass

import os
import glob
from distutils.core import setup

setup(name = "stochprop",
    version = '0.1.3',
    description = "A set of tools for analyzing atmospheric variability and uncertainty in infrasound propagation studies",
    author_email = 'pblom at lanl dot gov',
    packages = ['stochprop',
                'stochprop.cli'],
    scripts=[],

    entry_points = {'console_scripts':['stochprop=stochprop.cli.__main__:main']},

    install_requires = ['click',
                        'matplotlib',
                        'numpy',
                        'numpydoc',
                        'scipy',
                        'netcdf4',
                        'ipython',
                        'cartopy',
                        'pyproj',
                        'pathos',
                        'sphinx',
                        'sphinx_rtd_theme']
     )
