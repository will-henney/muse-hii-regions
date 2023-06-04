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

# # Kinematics of [O III] with PZ cube
#
# We will calculate the velocity moments for [O III] from the Peter Zeidler cube to compare with the ones we already calculated for H alpha, and also with the [O III] measurements from the pipeline ESO cube.  There are some advantages in using [O III] over H alpha:
#
# * There are two lines with a known ratio, which allows checking the intensity zero-point easily
# * Thermal broadening is less
# * Photospheric absorption is less
#
# On the other hand, there is at least one disadvantage:
#
# * The spectral resolution is lower (in velocity units)
#

# Load libraries and data cube (identical to earlier notebooks):

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import astropy.units as u

sns.set_context("talk")
sns.set_color_codes()
# -

datapath = Path("/Users/will/Work/Muse-Hii-Data/SMC-NGC-346/")
fitsfilepath = datapath / "PeterZeidler" / "DATACUBE_FINAL_fwhm_cor.fits"
cube = Cube(str(fitsfilepath))

# I have now moved functions for dealing with the velocity moments to a separate library:

import sys
from whispy import moments

# Where to save figures and FITS images:

moments.FIGPATH = Path("../../figs")
moments.SAVEPATH = Path("../../data")

# ## Choose the wavelength range
#
# ### Broad overview
#
# First we inspect the spectrum in 6 broad horizontal strips across the image, from south to north.

jstrips = [
    [0, 50],
    [50, 100],
    [100, 150],
    [150, 200],
    [200, 250],
    [250, -1],
]

wide_band = cube.select_lambda(4850, 5050)
fig, ax = plt.subplots(figsize=(10, 8))
for j1, j2 in jstrips:
    (wide_band[:, j1:j2, :].mean(axis=(1, 2)).plot(label=f"j = {j1}:{j2}"))
ax.legend(ncol=2)
ax.set(yscale="log")
sns.despine()

# So, this makes it look like all the lines are positive (and hopefully for the PZ cube this is the case).  
# * The other strong line we see is H beta 4861. 
# * The other weak lines that we see are He I 4922 and 5016 Å (and marginally [Fe III] 4987)
# * If we compare with the original version, the continuum is not so flat, but instead shows lots of undulations

# ### Narrow in on the [O III] doublet
#
# Set some limits for the continuum and line extraction:

wlim = {
    "4959": {
        "core": [4953.0, 4969.0],
        "blue": [4947.0, 4953.0],
        "red": [4969.0, 4975.0],
    },
    "5007": {
        "core": [5002.0, 5016.0],
        "blue": [4995.0, 5001.0],
        "red": [5022.0, 5028.0],
    },
}
rangecolors = {"core": "g", "blue": "b", "red": "r"}

# Plot these limits zoomed in on the two [O III] lines:

# +
medium_band = cube.select_lambda(4940, 5040)
fig, ax = plt.subplots(figsize=(10, 8))
for j1, j2 in jstrips:
    (
        medium_band[:, j1:j2, :]
        .mean(axis=(1, 2))
        .plot(label=f"j = {j1}:{j2}", linewidth=2)
    )

for line, linedata in wlim.items():
    for span, spandata in linedata.items():
        ax.axvspan(
            spandata[0], spandata[1], alpha=0.2, color=rangecolors[span], zorder=-100
        )
ax.legend(ncol=2)
ax.set(yscale="linear", ylim=[0.0, 1800])
sns.despine()


# -

# That seems to look fine.  Note that the red continuum band for 5007 is separated a bit from the core to allow space for the He I line (although it is so weak that this probably doesn't matter).
#
# ### Do continuum subtraction
#
# Now, we can use the same `wlim` data to extract the line and continuum:


def extract_core_and_cont(cube, spandata):
    """Return continuum-subtracted line core and continuum map

    The line core is a 3D cube over the narrow core wavelengths
    """
    cblue = cube.select_lambda(*spandata["blue"]).mean(axis=0)
    cred = cube.select_lambda(*spandata["red"]).mean(axis=0)
    cont = 0.5 * (cblue + cred)
    core = cube.select_lambda(*spandata["core"]) - cont
    return core, cont


core5007, cont5007 = extract_core_and_cont(medium_band, wlim["5007"])
core4959, cont4959 = extract_core_and_cont(medium_band, wlim["4959"])

# Plot the line core and the continuum for 5007 and 4959

fig, axes = plt.subplots(
    2,
    2,
    figsize=(10, 10),
    sharex=True,
    sharey=True,
)
core5007.sum(axis=0).plot(ax=axes[0, 0])
core4959.sum(axis=0).plot(ax=axes[0, 1])
cont5007.plot(ax=axes[1, 0], scale="sqrt")
cont4959.plot(ax=axes[1, 1], scale="sqrt")
axes[0, 0].contour(core5007.sum(axis=0).data, levels=[0.0], colors="r")
axes[0, 1].contour(core4959.sum(axis=0).data, levels=[0.0], colors="r")
axes[1, 0].contour(cont5007.data, levels=[0.0], colors="r")
axes[1, 1].contour(cont4959.data, levels=[0.0], colors="r")
fig.tight_layout(pad=0)

# The red contours show the zero level. Unlike with the ESO cube, there are no negative regions in either line or continuum. This is good.

# ## Wavelength moments
#
# We now use the functions defined in my `moments` module (part of whispy package).

# ### First look without having fixed the sky
#
# This comes out fine

mom5007 = moments.find_moments(core5007)
fig, ax = plt.subplots(figsize=(8, 8))
mom5007[1].plot(
    cmap="seismic",
    vmin=5009.2,
    vmax=5009.85,
    colorbar="v",
)

# This looks similar to the sky-corrected version that I obtained from the ESO cube, although the amplitude of artifacts (cross-hatched patterns) seems larger

# So we skip entirely the correction steps we did lat time.
#

# But we do do the side-by-side comparison of the first moment for 5007 and 4959:

mom5007 = moments.find_moments(core5007)
mom4959 = moments.find_moments(core4959)
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
mom5007[1].plot(
    cmap="seismic",
    vmin=5009.2,
    vmax=5009.85,
    colorbar="h",
    ax=axes[0],
)
mom4959[1].plot(
    cmap="seismic",
    vmin=4961.25,
    vmax=4961.9,
    colorbar="h",
    ax=axes[1],
)
fig.tight_layout(h_pad=10)

# These look very similar to the H alpha map in general.  They are also extremely similar to one another.  Both of these facts are reassuring.
#
# On the negative side, as mentioned above, we see a lot of small scale patterming, which seem to be due to rebinning onto the uniform wavelength grid of the cube

fig, axes = plt.subplots(
    1,
    2,
    figsize=(10, 5),
    sharey=True,
)
mom5007[2].plot(
    cmap="gray",
    vmin=0.5,
    vmax=1.5,
    colorbar="v",
    ax=axes[0],
)
mom4959[2].plot(
    cmap="gray",
    vmin=0.5,
    vmax=1.5,
    colorbar="v",
    ax=axes[1],
)
fig.tight_layout()

# These maps of sigma seem to show less artifacts than the equivalent ones from the ESO cube. 

fig, ax = plt.subplots()
wav5007 = np.median(mom5007[1].data.data)
wav4959 = np.median(mom4959[1].data.data)
# ax.plot(mom5007[1].data.mean(axis=1) - wav5007)
# ax.plot(mom4959[1].data.mean(axis=1) - wav4959)
ax.plot(np.median(mom5007[1].data.data, axis=1) - wav5007)
ax.plot(np.median(mom4959[1].data.data, axis=1) - wav4959)
ax.set(ylim=[-0.12, 0.12])
for i0 in [80, 160, 240]:
    ax.axvline(i0, lw=0.5, color="k", alpha=0.2)
wav5007, wav4959

# This is the median N-S (vertical) profile of the map, which nicely shows the blue-red-blue pattern as one goes from south to north. The agreement between the two doublet components is excellent. The pattern is similar to what I got with the ESO cube, except the total variation is a bit smaller. 

fig, ax = plt.subplots()
ax.plot(np.median(mom5007[1].data.data, axis=0) - wav5007)
ax.plot(np.median(mom4959[1].data.data, axis=0) - wav4959)
# ax.plot(mom5007[1].data.mean(axis=0) - wav5007)
# ax.plot(mom4959[1].data.mean(axis=0) - wav4959)
ax.set(ylim=[-0.12, 0.12])
# for i0 in [80, 160, 240]:
#    ax.axvline(i0, lw=0.5, color="k", alpha=0.2)

# And this is the same in the E-W (horizontal) direction. There is less horizontal variation than vertical variation, as is obvious from just looking at the maps. 

import pandas as pd

save_pars_5007 = dict(
    label="5007",
    flabel="ngc346-oiii",
    restwav=5006.84,
    irange=[5.0e4, 3.0e5],
    vrange=[150, 175],
    srange=[60, 90],
)
plot_pars_5007 = dict(
    ilabel="[O III]",
    **save_pars_5007,
)

moments.save_moments_to_fits(
    mom5007,
    **save_pars_5007,
)

g = moments.moments_corner_plot(mom5007, rebin=1, **plot_pars_5007)

g = moments.moments_corner_plot(mom5007, rebin=4, **plot_pars_5007)

g = moments.moments_corner_plot(
    mom5007,
    rebin=16,
    **plot_pars_5007,
    hist_bins=40,
    image_bins=20,
)

m = (
    mom5007[0].mask
    | (mom5007[0].data < 2e4)
    | (mom5007[1].data < 5009.2)
    | (mom5007[1].data > 5009.8)
    | (mom5007[2].data < 0.8)
    | (mom5007[2].data > 1.5)
)
df = pd.DataFrame(
    {
        "log10 I(5007)": np.log10(mom5007[0].data[~m]),
        "V(5007)": 3e5 * (mom5007[1].data[~m] - 5006.84) / 5006.84,
        "sig(5007)": 3e5 * mom5007[2].data[~m] / 5006.84,
    }
)
df.describe()

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(color="r"),
    diag_kws=dict(color="r"),
)
g.fig.suptitle("[O III] 5007 corrected, normalized moments")
g.tight_layout(pad=0)

m = (
    mom4959[0].mask
    | (mom4959[0].data < 0.67e4)
    | (mom4959[1].data < 4961.2)
    | (mom4959[1].data > 4961.8)
    | (mom4959[2].data < 0.8)
    | (mom4959[2].data > 1.5)
)
df4959 = pd.DataFrame(
    {
        "log10 I(4959)": np.log10(mom4959[0].data[~m]),
        "V(4959)": 3e5 * (mom4959[1].data[~m] - 4958.91) / 4958.91,
        "sig(4959)": 3e5 * mom4959[2].data[~m] / 4958.91,
    }
)
df4959.describe()

g = sns.pairplot(
    df4959,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(color="m"),
    diag_kws=dict(color="m"),
)
g.fig.suptitle("[O III] 4959 corrected, normalized moments")
g.tight_layout(pad=0)

# The distributions are very consistent between the two [O III] lines.  They are also quite similar to H alpha, especially in the mid-range of intensity.
#
# Unlike in the ESO cube, there is not really a bimodal distribution in the velocities of either line. The marginal histogram is single-peaked, with a maximum at about 162 km/s. 
#
# However, on the joint I-v distribution, we can see evidence of a multi-peaked structure. The 162 km/s peak corresponds to the intermediate brightness values. Lower brightness values have a distinct peak at 159 km/s. And then the brightest pixels also seem to move twoards slightly bluer velocities (160 km/s). Finally, the intermediate brightnesses show a "spike" towards higher velocities: 165 km/s.  
#
# The line width does not vary with velocity, but there is some evidence for structure in the I-sigma plot.

# Now look at the cross-correlations between the two lines:

m = (
    mom5007[0].mask
    | (mom5007[0].data < 3e4)
    | (mom5007[1].data < 5009.3)
    | (mom5007[1].data > 5009.8)
    | (mom5007[2].data < 1.0)
    | (mom5007[2].data > 1.5)
    | mom4959[0].mask
    | (mom4959[0].data < 1e4)
    | (mom4959[1].data < 4961.3)
    | (mom4959[1].data > 4961.8)
    | (mom4959[2].data < 1.0)
    | (mom4959[2].data > 1.5)
)
df2 = pd.DataFrame(
    {
        "log10 I(5007)": np.log10(mom5007[0].data[~m]),
        "V(5007)": 3e5 * (mom5007[1].data[~m] - 5006.84) / 5006.84,
        "sig(5007)": 3e5 * mom5007[2].data[~m] / 5006.84,
        "log10 I(4959)": np.log10(mom4959[0].data[~m]),
        "V(4959)": 3e5 * (mom4959[1].data[~m] - 4958.91) / 4958.91,
        "sig(4959)": 3e5 * mom4959[2].data[~m] / 4958.91,
    }
)
df2.corr()

xvars = [_ for _ in df2.columns if "5007" in _]
yvars = [_ for _ in df2.columns if "4959" in _]
xvars, yvars

g = sns.pairplot(
    df2,
    kind="hist",
    height=4,
    x_vars=xvars,
    y_vars=yvars,
    plot_kws=dict(color="b"),
)
g.axes[1, 1].plot([150, 175], [150, 175], c="r")
g.axes[2, 2].plot([60, 90], [60, 90], c="r")
g.fig.suptitle("Correlations between 5007 and 4959")
g.tight_layout(pad=0)

# That all looks good, although the correlation coefficients for V-V and sig-sig are a bit lower than for the ESO cube, presumably due to the pattern noise. 

df3 = df2[["log10 I(5007)"]].copy()
df3["5007 / 4959"] = 10 ** (df2["log10 I(5007)"] - df2["log10 I(4959)"])
df3["dV"] = df2["V(5007)"] - df2["V(4959)"]
df3["sig ratio"] = df2["sig(5007)"] / df2["sig(4959)"]
df3.describe()

# +
m = (
    (df3["5007 / 4959"] < 2.9)
    | (df3["5007 / 4959"] > 3.1)
    | (np.abs(df3["dV"]) > 5.0)
    | (df3["sig ratio"] < 0.8)
    | (df3["sig ratio"] > 1.1)
)

df3 = df3[~m]
df3.corr()
# -

g = sns.pairplot(
    df3,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(color="r"),
    diag_kws=dict(color="r"),
)
g.fig.suptitle("[O III] 5007 vs 4959 ratios and differences")
g.tight_layout(pad=0)

# No meaningful correlations or trends here, which is good. And better than the ESO cube, where there were several. There is a slight tendency for sigma ratio to increase with intensity ratio, which could be evidence that there is an additional line that blends with 5007. 


