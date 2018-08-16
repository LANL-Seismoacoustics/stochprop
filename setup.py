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
      description = "A set of tools for computing representative atmospheres from statistical analysis of a suite of G2S specifications and then run propagation simulations through the resulting set of atmospheres to build stochastic propagation models.",
      author_email = 'pblom at lanl dot gov',
      packages = ['stochprop'],
      scripts=[],
      install_requires = ['numpy', 'scipy', 'obspy', 'IPython', 'infrapy']
     )
