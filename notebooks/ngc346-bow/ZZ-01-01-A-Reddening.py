# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Notebook A of line ratios from saved maps from the PZ cube
#
# This is from the splitting up of ZZ-01-01-more-line-ratios.ipynb.  It includes calculation of the reddening map

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
import regions
import sys
import pandas as pd
import cmasher as cmr
import pyneb as pn

sns.set_context("talk")
sns.set_color_codes()
# -
# ## Path to the root of this repo

ROOT = Path.cwd().parent.parent

mapspath = ROOT / "data/n346-bow-lines/refined"


# +
def get_image(lineid, combo="P-007", ext="SUM"):
    im = Image(str(mapspath / f"map-{lineid}-{combo}.fits"), ext=ext)
    im.data = im.data.astype(float)
    return im


def get_image_raw(lineid, variant="csub", channel="ABC"):
    p = ROOT / "data" / "n346-bow-lines" / f"maps-n346-PZ-2pass-{variant}-007"
    candidates = list(p.glob(f"**/*-{lineid}-{channel}.fits"))
    assert candidates
    imfile = candidates[0]
    im = Image(str(imfile))
    im.data = im.data.astype(float)
    return im


# -

# ## Calculate reddening from Balmer decrement

# Load the Hα and Hβ maps
imha = get_image_raw("h-i-6562-79")
imhb = get_image_raw("h-i-4861-32")

# ### Look at the raw Hα/Hβ ratio:

fig, ax = plt.subplots(figsize=(12, 12))
(imha / imhb).plot(vmin=2.7, vmax=3.7, cmap="gray", colorbar="v")

# So we see very little structure there, compared with in the ESO pipeline cube case.  In priniciple, lighter means more extinction.  There might be a hint of this at the bottom of the image, where we see possible signs of the foreground filament.
#
# But in other parts, we just see the stars (which have a different ratio because of photospheric absorption)

# ### PyNeb calculation of intrinsic Balmer decrement

hi = pn.RecAtom("H", 1)

# Calculate the theoretical Balmer decrement from PyNeb. Density and temperature from Valerdi:2019a

tem, den = 12500, 100
R0 = hi.getEmissivity(tem, den, wave=6563) / hi.getEmissivity(tem, den, wave=4861)
R0

# ### Look at correlation between Hα and Hβ in the faint limit
#
# To make thinks easier, I multiply the Hb values by R0 so we have a square plot.  I zoom in on the faint parts:

imax = 100000
m = imha.data < imax
m = m & (R0 * imhb.data < imax)
m = m & ~imha.mask & ~imhb.mask
df = pd.DataFrame(
    {
        "ha": imha.data[m],
        "hb": R0 * imhb.data[m],
    }
)

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, 0.9 * imax], "--", color="r", linewidth=1.4)
g.axes[1, 0].plot([0, imax], [0, 0.8 * imax], "--", color="r", linewidth=0.7)
g.fig.suptitle("Correlation between Ha and Hb brightness")

# So, the slope is not unity, meaning the extinction is not zero.  But the intercept is zero, which is great. So nothing more to do,

# Now define some regions to take averages

boxes = {
    "sw filament": regions.RegionBoundingBox(
        iymin=20,
        iymax=40,
        ixmin=200,
        ixmax=310,
    ),
    "bow shock": regions.RegionBoundingBox(
        iymin=165,
        iymax=205,
        ixmin=240,
        ixmax=290,
    ),
    "w filament": regions.RegionBoundingBox(
        iymin=100,
        iymax=130,
        ixmin=25,
        ixmax=55,
    ),
    "c filament": regions.RegionBoundingBox(
        iymin=195,
        iymax=210,
        ixmin=155,
        ixmax=195,
    ),
}

# Plot on a better scale and show the regions:

# +
fig, ax = plt.subplots(figsize=(12, 12))
(imha / (imhb)).plot(
    vmin=R0,
    vmax=1.5 * R0,
    scale="linear",
    cmap="gray_r",
    colorbar="v",
)

for box, c in zip(boxes.values(), "yrmgc"):
    box.plot(
        ax=ax,
        lw=3,
        edgecolor=c,
        linestyle="dashed",
        #        facecolor=(1.0, 1.0, 1.0, 0.4),
        fill=False,
    )
# -

# We can see some very high extinction at in the S filaments.  And some small increase in extinction in the main diagonal filament.  This is probably limited having foreground emission to some extent.
#
# Look at average values in the sample boxes

for label, box in boxes.items():
    (yslice, xslice), _ = box.get_overlap_slices(imha.shape)
    ha = np.median(imha[yslice, xslice].data.data)
    hb = np.median(imhb[yslice, xslice].data.data)
    print(f"{label}: {ha / hb:.3f}")

# For reference, the values from the PZ notebook
# ```
# sw filament: 3.318
# bow shock: 3.156
# w filament: 3.226
# c filament: 3.224
# ```

# I tried mean and median, and it made very little difference.  Lowest in the bow shock region; slightly higher in the west and central filaments.  Even higher in the southwest filament.

# ### Equivalent widths
# *New 2024-04-23*

# imnii = Image(str(ROOT / "data/ngc346-PZ-nii-6583-bin01-sum.fits"))
imcont_ha = get_image_raw("h-i-6562-79", variant="cont")
imcont_hb = get_image_raw("h-i-4861-32", variant="cont")

imcont_ha.wcs

imha.wcs

# We now have a local continuum for all the lines

ewha = 3 * 1.25 * imha / imcont_ha

# +
fig, ax = plt.subplots(figsize=(12, 12))
(ewha).plot(
    vmin=0,
    vmax=1000,
    scale="linear",
    cmap="gray_r",
    colorbar="v",
)

for box, c in zip(boxes.values(), "yrmgc"):
    box.plot(
        ax=ax,
        lw=3,
        edgecolor=c,
        linestyle="dashed",
        #        facecolor=(1.0, 1.0, 1.0, 0.4),
        fill=False,
    )
# -

for label, box in boxes.items():
    (yslice, xslice), _ = box.get_overlap_slices(imha.shape)
    ew = np.mean(ewha[yslice, xslice].data.data)
    dew = np.std(ewha[yslice, xslice].data.data)
    print(f"{label}: {ew:.3f} +/- {dew:.3f}")


# ### The reddening law

pn.RedCorr().getLaws()


# PyNeb does not seem to have anything specifically tailored to the SMC.  The average SMC extinction law is supposedly simply $1/\lambda$.
#
# But, it is possible to get a SMC curve by using the "F99-like" option, which uses the curve of Fitzpatrick & Massa 1990, ApJS, 72, 163. This depends on $R_V$ and 6 other parameters (!!!).  Most of the parameters only affect the UV part of the curve, which does not concern us.
#
# Then, we can use the average values of $R_V$ and the other parameters, which were fit by Gordon:2003l to SMC stars. This is $R_V = 2.74 \pm 0.13$.
#
# So here I compare that SMC curve with $1/\lambda$ and with the Clayton curve for Milky Way (but also adjusted to $R_V = 2.74$):


# +
def A_lam(wave):
    return 4861.32 / wave


def my_X(wave, params=[]):
    """A_lam / E(B - V) ~ lam^-1"""
    return A_lam(wave) / (A_lam(4400) - A_lam(5500))


rc = pn.RedCorr()
rc.UserFunction = my_X
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
f, ax = plt.subplots(figsize=(10, 10))
rc.plot(laws=["user", "F99", "CCM89"], ax=ax)

ax.set(
    xlim=[4000, 9300],
    ylim=[-2.5, 1.0],
    #    xlim=[4000, 7000],
    #    ylim=[-1, 1],
)
# -

# So the Gordon curve is flatter in the blue, steeper in green, and flatter in red, as compared to $1/\lambda$.

rc = pn.RedCorr()
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
rc.law = "F99"

# Test it out for the bow shock region:

rc.setCorr(obs_over_theo=3.123 / R0, wave1=6563.0, wave2=4861.0)
rc.E_BV, rc.cHbeta

# And for the highest extinction region

rc.setCorr(obs_over_theo=4.165 / R0, wave1=6563.0, wave2=4861.0)
rc.E_BV, rc.cHbeta

# So $E(B - V)$ varies from about 0.1 to about 0.35. This is similar to what is found for the stars.

# ### The reddening map
#
# We can now make a map of $E(B - V)$

R = imha / (imhb)
rc.setCorr(obs_over_theo=R.data / R0, wave1=6563.0, wave2=4861.0)
imEBV = R.copy()
imEBV.data = rc.E_BV
imEBV.mask = imha.mask | imhb.mask

fig, ax = plt.subplots(figsize=(12, 12))
imEBV.rebin(4).plot(
    vmin=0.0,
    vmax=1.0,
    scale="linear",
    cmap="magma_r",
    colorbar="v",
)

# Looks like I would expect. Check values in the boxes:

for label, box in boxes.items():
    (yslice, xslice), _ = box.get_overlap_slices(imha.shape)
    ebv = np.median(imEBV[yslice, xslice].data.data)
    print(f"{label}: {ebv:.3f}")

# These seem the same as before.  But we want to eliminate extreme values.

badpix = (imEBV.data > 1.0) | (imEBV.data < 0.0)
imEBV.mask = imEBV.mask | badpix

# Save it to a file:

imEBV.write(str(ROOT / "data/ngc346-ZZ-reddening-E_BV.fits"), savemask="nan")

# Lots of regions are affected by the stellar absorption.  There are apparent increases in reddening at the position of each star.  This is not real, but is due to the photospheric absorption having more of an effect on Hb (mainly because the emission line is weaker).
#
# At some point, I am going to have to deal with that. But it is not an issue for the bow shock emission, since this is in an area free of stars.  We should just use the median bow shock reddening of $E(B-V) = 0.139$ so that we don't introduce any extra noise.

# +
rc.E_BV = 0.139
wavs = np.arange(4600, 9300)
Alam = rc.E_BV * rc.X(wavs)

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(wavs, Alam)
ax.set(
    xlabel="Wavelength, $\lambda$, Å",
    ylabel="Extinction, $A_\lambda$, mag",
    title=f"Extinction curve for $E(B - V) = {rc.E_BV:.3f}$",
)
sns.despine()
