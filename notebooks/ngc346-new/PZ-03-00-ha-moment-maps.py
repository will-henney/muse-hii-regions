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
fitsfilepath = datapath / "PeterZeidler" / "DATACUBE_FINAL_fwhm_cor.fits"
cube = Cube(str(fitsfilepath))

import sys

from whispy import moments

moments.FIGPATH = Path("../../figs")
moments.SAVEPATH = Path("../../data")

# ## Select the wavelength range
#
# Then take a narrow band about the H alpha line that avoids the flanking [N II] lines:

hacube = cube.select_lambda(6558.0, 6576.0)

# Plot various aggregate measures over all the spatial pixels (mean, median, max, min):

fig, ax = plt.subplots(figsize=(10, 5))
hacube.mean(axis=(1, 2)).plot(label="mean")
(0.01 * hacube.max(axis=(1, 2))).plot(label="max x 0.01")
hacube.median(axis=(1, 2)).plot(label="median")
(hacube.min(axis=(1, 2))).plot(label="min")
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
(0.03 * hacube_contsub.max(axis=(1, 2))).plot(label="max x 0.03")
hacube_contsub.median(axis=(1, 2)).plot(label="median")
(0.1 * hacube_contsub.min(axis=(1, 2))).plot(label="min x 0.1")
ax.legend()

# This looks good – the continuum has been eliminated in all but the min aggregation, where it goes slightly negative on the blue side, but this is probably real evidence of the underlying stellar absorption

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
    vmin=0,
    vmax=200000,
    cmap="viridis",
    colorbar="v",
)
ax.contour(ha_cont.data, levels=[5000], colors="r")

# We have no negative zones, so we do not need to do what we did previously. So we can ditch all this.

# ### Try looking at joint distribution of unnormalized moments
#
# The unnormalized moments should behave better since they can pass through zero without having the 0/0 problem. .

# #### Before correcting the sky

mom0 = hacore.sum(axis=0)
mom1 = mom0.copy()
# wav0 = 6566.5
wav0 = 6566.4
mom1.data = np.sum(hacore.data * (wavcore.data - wav0), axis=0)
mom2 = mom0.copy()
mom2.data = np.sum(hacore.data * (wavcore.data - wav0) ** 2, axis=0)

# Note that we do not divide `mom1` by `mom0`

import pandas as pd

sns.set_color_codes()

# We mask out the brightest parts and also the outliers in velocity so that we can see the part that passes through zero better.

starmask = ha_cont.data > 1e4
m = starmask | mom0.mask | (mom0.data > 2e5) | (np.abs(mom1.data) > 3e4)
df = pd.DataFrame(
    {
        "mom0": mom0.data[~m],
        "mom1": mom1.data[~m],
        "mom2": mom2.data[~m],
    }
)
df.describe()

# Make a corner plot of the moment distributions:

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.axes[1, 0].axvline(0.0, color="r", linestyle="dashed")
g.axes[2, 0].axvline(0.0, color="r", linestyle="dashed")
g.axes[2, 0].axhline(0.0, color="r", linestyle="dashed")
g.axes[2, 1].axhline(0.0, color="r", linestyle="dashed")
g.fig.suptitle("Uncorrected, unnormalized moments")
g.tight_layout()

# Note that for an emission line, mom0 and mom2 should be positive definite.  They are now, so that is good.
#

# #### (What used to be) After correcting the sky
#
# **We do not need to do this any more really. **
#
# I have removed everything that creates the underscore-prefixed variables, and then swap back the original variables in the plots

# We take a less restrictive mask, only masking out the outliers in any of the moments:

m = (
    starmask
    | mom0.mask
    | (mom0.data > 45e4)
    | (mom1.data < -3e4)
    | (mom1.data > 3e4)
    | (mom2.data < 0)
    | (mom2.data > 4.5e5)
)
df = pd.DataFrame(
    {
        "mom0": mom0.data[~m],
        "mom1": mom1.data[~m],
        "mom2": mom2.data[~m],
    }
)
df.describe()

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(color="g"),
    diag_kws=dict(color="g"),
)
g.fig.suptitle("Corrected, unnormalized moments")
g.tight_layout(pad=0)

# So that looks like it has worked. Note that mom1 fans out as a triangle as mom0 increases, but it stays roughly centered on zero. This is because of the judicious adjustment to wav0 above.

# #### Back to normalized moments again
#
# If that is sorted, then we should be able to look at the normalized moments:

mom0 = hacore.sum(axis=0)
mom1 = mom0.copy()
mom1.data = np.sum(hacore.data * (wavcore.data - wav0), axis=0) / mom0.data
mom2 = mom0.copy()
mom2.data = (
    np.sum(hacore.data * (wavcore.data - wav0 - mom1.data) ** 2, axis=0) / mom0.data
)

m = (
    starmask
    | mom0.mask
    | (mom0.data < 1.6e4)
    | (mom0.data > 45e4)
    | (mom1.data < -0.4)
    | (mom1.data > 0.4)
    | (mom2.data > 1.75)
    | (mom2.data < 0.5)
)
KMS_PER_ANGSTROM = 3e5 / 6563.0
df2 = pd.DataFrame(
    {
        "log10 I(Ha)": np.log10(mom0.data[~m]),
        "V(Ha)": mom1.data[~m] * KMS_PER_ANGSTROM,
        "sig(Ha)": np.sqrt(mom2.data[~m]) * KMS_PER_ANGSTROM,
    }
)
df2.describe()


fig, ax = plt.subplots(figsize=(8, 8))
mom2.mask = mom2.mask | m
mom2.plot(
    cmap="gray",
    vmin=0.5,
    vmax=1.75,
    colorbar="v",
)
ax.contour(KMS_PER_ANGSTROM * mom1.data, levels=[-10, -5, 0, 5, 10], cmap="bwr")

# ## Redo it all using my moments library
#
# Lots of the above can be done more easily now:

mom6563 = moments.find_moments(hacore)

mom6563[0].info()

mom_pars_6563 = dict(
    restwav=6562.79,
    irange=[3.0e4, 1.0e6],
    vrange=[150, 180],
    srange=[40, 60],
)

moments.save_moments_to_fits(
    mom6563,
    label="6563",
    flabel="ngc346-PZ-hi",
    **mom_pars_6563,
)

plot_pars_6563 = dict(
    ilabel="H I",
    label="6563",
    flabel="ngc346-PZ-hi",
    **mom_pars_6563,
)
g = moments.moments_corner_plot(mom6563, rebin=1, **plot_pars_6563)

g = moments.moments_corner_plot(mom6563, rebin=8, **plot_pars_6563)

# So, these results look a bit more "boring" than the ones from the ESO cube, but they are probably more reliable.
#
# The range of velocities is reduced and the width does not show any significant change with velocity.
