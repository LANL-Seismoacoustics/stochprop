# -*- coding: utf-8 -*-
try:
    import setuptools
except ImportError:
    pass

import os
import glob
from distutils.core import setup

setup(name = "stochprop",
      version = '0.1',
      description = "A set of tools for analyzing atmospheric variability and uncertainty for infrasound propagation studies",
      author_email = 'pblom at lanl dot gov',
      packages = ['stochprop'],
      scripts=[],
      install_requires = ['numpy',
                          'scipy',
                          'matplotlib',
                          'pyproj',
                          'pathos',
                          'sphinx']
     )
