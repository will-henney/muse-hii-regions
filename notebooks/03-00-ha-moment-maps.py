# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
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

# # Estimate Hα velocity moment maps for NGC 346

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import astropy.units as u
sns.set_context("talk")

# Load the data cube from the FITS file. 

datapath = Path("/Users/will/Work/Muse-Hii-Data/SMC-NGC-346/")
fitsfilepath = datapath / "ADP.2017-10-16T11_04_19.247.fits"
cube = Cube(str(fitsfilepath))

# ## Select the wavelength range
#
# Then take a narrow band about the H alpha line that avoids the flanking [N II] lines:

hacube = cube.select_lambda(6558.0, 6576.0)

# Plot various aggregate measures over all the spatial pixels (mean, median, max, min):

fig, ax = plt.subplots(figsize=(10, 5))
hacube.mean(axis=(1, 2)).plot(label="mean")
(0.01*hacube.max(axis=(1, 2))).plot(label="max x 0.01")
hacube.median(axis=(1, 2)).plot(label="median")
(0.1*hacube.min(axis=(1, 2))).plot(label="min x 0.1")
ax.legend()

# ## Subtract the continuum
#
# On the plus side, it looks like there is never any significant continuum slope, so we can just calculate an average continuum level to subtract:

left_continuum = cube.select_lambda(6558.0, 6562.0).mean(axis=0)
right_continuum = cube.select_lambda(6572.0, 6576.0).mean(axis=0)
ha_cont = 0.5 * (left_continuum + right_continuum)

ha_cont.plot(colorbar="v")

# Promote the continuum image to a cube by taking a copy of the the line cube and pasting in new data and variance: 

contcube = hacube.copy()
contcube.data = np.ones_like(contcube.data) * ha_cont.data[None, :, :]
contcube.var = np.ones_like(contcube.var) * ha_cont.var[None, :, :]

# Subtract the continuum.

hacube_contsub = hacube - contcube 

# Another way of doing it would have been to just subtract `ha_cont.data[None, :, :]` from the `.data` array of `hacube`.  But in that case, the `.var` would not have been propagated automatically and we would have to calculate it by hand. 
#
# Now we do the same plot:

fig, ax = plt.subplots(figsize=(10, 5))
hacube_contsub.mean(axis=(1, 2)).plot(label="mean")
(0.03*hacube_contsub.max(axis=(1, 2))).plot(label="max x 0.03")
hacube_contsub.median(axis=(1, 2)).plot(label="median")
(0.1*hacube_contsub.min(axis=(1, 2))).plot(label="min x 0.1")
ax.legend()

# This looks good – the continuum has been eliminated in all but the min aggregation, which is presumably dominated by a small number of pixels that have strange profiles.

# ## Make a cube of wavelengths

# It is is easy to get a 1D array of the wavelengths:

waves = hacube_contsub.wave.coord()
waves

# But there seems no built-in way to get a cube of them, so we will construct it ourselves:

wavcube = hacube_contsub.clone(np.ones, np.zeros)
wavcube.data *= waves[:, None, None]
wavcube.unit = u.angstrom
wavcube.info()

# Note that we used the `.clone` method to get a new empty cube with the same coordinates, which we initialized with all ones as data and all zeros as variances (because there are no uncertainties in the wavelengths).  Then we used broadcasting to multiply the data by the 1D wavelength array and set the unit correctly.

# ## Calculate the moments
#
# To reduce effects noise, especially in the higher moments, we use the core window of 6562 to 6572 in between the continuum-fitting ranges

hacore = hacube_contsub.select_lambda(6562.0, 6572.0)
wavcore = wavcube.select_lambda(6562.0, 6572.0)

# The velocity moments are now trivial sums over this window.  

mom0 = hacore.sum(axis=0)
wav0 = 6566.5
mom1 = mom0.copy()
mom1.data = np.sum(hacore.data * (wavcore.data - wav0), axis=0) / mom0.data
# mom1.mask = mom1.mask | (mom0.data < 0.0)

fig, ax = plt.subplots(figsize=(8, 8))
mom1.plot(
    cmap="seismic",
    vmin=-1.0, 
    vmax=1.0,
    colorbar="v",
)

# ## Try and deal with sky over-subtraction

# Look at where Ha brightness is negative, and overlay with positions of stars (continuum is high)

fig, ax = plt.subplots(figsize=(8, 8))
mom0.plot(
#    vmin=0.5*mom0.data.min(),
#    vmax=0.0,
    vmin=-10000,
    vmax=10000,
    cmap="viridis",
    colorbar="v",
)
ax.contour(ha_cont.data, levels=[5000], colors="r")

# Look at some of the line profiles in the negative zones:

hacube[:, 300, 60].plot()
hacube[:, 10, 300].plot()
hacube[:, 150, 150].plot()
hacube[:, 260, 160].plot()

# So what happens if we take the deepest one of those (in the lower right corner) and subtract it from each of the others:

# skyspec = hacube[:, 10, 300]
starmask2d = ha_cont.data > 1000.0
bmask2d = hacore.sum(axis=0).data > -20000
skycube = hacube.copy()
skycube.mask = skycube.mask | starmask2d[None, :, :] | bmask2d[None, :, :]
skyspec = skycube.mean(axis=(1, 2))
(hacube[:, 300, 60] - skyspec).plot()
(hacube[:, 150, 150] - skyspec).plot()
(hacube[:, 260, 160] - skyspec).plot()

# That  doesn't quite work because it is from a region that is particularly blue-shifted.  But we will try it anyhow.

(hacube - skyspec).sum(axis=0).plot(vmin=0.0, vmax=1e5, colorbar="v")

# Redo the moments, but subtracting skyspec first:

_hacube_contsub = hacube - contcube - skyspec
_hacore = _hacube_contsub.select_lambda(6562.0, 6572.0)
_mom0 = _hacore.sum(axis=0)
_mom1 = _mom0.copy()
_mom1.data = np.sum(_hacore.data * (wavcore.data - wav0), axis=0) / _mom0.data

_mom0.plot(vmin=0.0, vmax=1e5, colorbar="v")

fig, ax = plt.subplots(figsize=(8, 8))
_mom1.plot(
    cmap="seismic",
    vmin=-0.75, 
    vmax=0.25,
    colorbar="v",
)

hacore[:, 70, 250].plot()
_hacore[:, 70, 250].plot()

hacore[:, 75, 200].plot()
_hacore[:, 75, 200].plot()

hacore[:, 250, 160].plot()
_hacore[:, 250, 160].plot()

hacore[:, 100, 30].plot()
_hacore[:, 100, 30].plot()

hacore[:, 150, 150].plot()
_hacore[:, 150, 150].plot()

# ### Try looking at joint distribution of unnormalized moments
#
# The unnormalized moments should behave better. 

mom0 = hacore.sum(axis=0)
mom1 = mom0.copy()
#wav0 = 6566.5
wav0 = 6566.6
mom1.data = np.sum(hacore.data * (wavcore.data - wav0), axis=0)
mom2 = mom0.copy()
mom2.data = np.sum(hacore.data * (wavcore.data - wav0)**2, axis=0)

# Note that we do not divide `mom1` by `mom0`

import pandas as pd
sns.set_color_codes()

starmask = ha_cont.data > 1e4
m = starmask | mom0.mask | (mom0.data > 3e4) | (np.abs(mom1.data) > 1e4)
df = pd.DataFrame({
    "mom0": mom0.data[~m],
    "mom1": mom1.data[~m],
    "mom2": mom2.data[~m],
})
df.describe()

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.axes[1, 0].axvline(0.0, color="r", linestyle="dashed")
g.axes[2, 0].axvline(0.0, color="r", linestyle="dashed")

_mom0 = _hacore.sum(axis=0)
_mom1 = _mom0.copy()
_mom1.data = np.sum(_hacore.data * (wavcore.data - wav0), axis=0)
_mom2 = _mom0.copy()
_mom2.data = np.sum(_hacore.data * (wavcore.data - wav0)**2, axis=0)

m = (
    starmask 
    | _mom0.mask 
    | (_mom0.data > 10e4) 
    | (_mom1.data < -2.5e4) | (_mom1.data > 0.5e4)
    | (_mom2.data < -1e4)  | (_mom2.data > 1.5e5) 
)
df = pd.DataFrame({
    "mom0": _mom0.data[~m],
    "mom1": _mom1.data[~m],
    "mom2": _mom2.data[~m],
})
df.describe()

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.axes[1, 0].axvline(0.0, color="r", linestyle="dashed")
g.axes[2, 0].axvline(0.0, color="r", linestyle="dashed")

# Back to normalized moments

_mom0 = _hacore.sum(axis=0)
_mom1 = _mom0.copy()
_mom1.data = np.sum(_hacore.data * (wavcore.data - wav0), axis=0) / _mom0.data
_mom2 = _mom0.copy()
_mom2.data = np.sum(_hacore.data * (wavcore.data - wav0 - _mom1.data)**2, axis=0) / _mom0.data

m = (
    starmask 
    | _mom0.mask 
    | (_mom0.data > 10e4) 
    | (_mom1.data < -0.5) | (_mom1.data > 0.25)
    | (_mom2.data > 1.75) | (_mom2.data < 0.5)
)
df2 = pd.DataFrame({
    "log10 mom0": np.log10(_mom0.data[~m]),
    "mom1": _mom1.data[~m],
    "mom2": _mom2.data[~m],
})
df2.describe()

g = sns.pairplot(
    df2,
    kind="hist",
    height=4,
    corner=True,
)

fig, ax = plt.subplots(figsize=(8, 8))
_mom2.mask = _mom2.mask | m
_mom2.plot(
    cmap="gray",
    vmin=0.5, 
    vmax=1.75,
    colorbar="v",
)
ax.contour(_mom1.data, levels=[-0.3, -0.2, -0.1, 0.0, 0.1], cmap="bwr")


