# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Concentrate on deep neutral lines in 30 Dor
#
# These lines need a more careful extraction for the following reasons:
#
# 1. They are very weak
# 2. They are affected by inaccuracies in the global continuum subtraction, especially near the Paschen jump, so need some hyperlocal continuum subtraction instead
# 3. Some of them are blended with other lines, which we could maybe subtract using a Gaussian fit to the integrated profile.
#

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import cmasher as cm
from mpdaf.obj import Cube
import regions
import sys
import warnings
from whispy.linetools import EmissionLine, SpectralRange, VelocityScale
from whispy import moments
from whispy import extract
from whispy import sky

import astropy.units as u
import astropy.constants as const

sns.set_context("talk")
sns.set_color_codes()
# -


