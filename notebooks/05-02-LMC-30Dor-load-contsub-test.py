# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,md,py:light
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

# Test of the 30 Dor continuum-subtracted cube containing Hα.  Mabel was having trouble with it.

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube

cubeA = Cube("../big-data/lmc-30dor-A-subcube-62-71-contsub.fits")

# Look at average spectrum of cube. 

cubeA.mean(axis=(1, 2)).plot()

# Looks bad.  Try setting the y limits to something more sensible:

fig, ax = plt.subplots()
cubeA.mean(axis=(1, 2)).plot()
ax.set(
    ylim=[-1e5, 1e5],
);

# So that still looks bad - must be some pixels with spurious data.
#
# Look at the summed image:

cubeA.sum(axis=0).plot(vmin=-10, vmax=1e6);

# So this shows noise in the middle parts that are due to the saturation of Hα, but that is not the problem.  Also has suspicious-looking pixels around the edges. So we will trim these off by expanding the mask by 1 pixel:

cubeA.spatial_erosion(1).mean(axis=(1, 2)).plot()

# That looks a lot better, so we will make the change to the cube

cubeA.spatial_erosion(1, inplace=True)

# Check the image again

cubeA.sum(axis=0).plot(vmin=-10, vmax=1e6);

# Replot the average spectrum, but zoom the y acis so we can see the faint parts:

fig, ax = plt.subplots(figsize=(12, 6))
cubeA.mean(axis=(1, 2)).plot()
ax.axhline(0.0, linestyle="dashed", color="k")
ax.set(
    ylim=[-30, 200],
);


