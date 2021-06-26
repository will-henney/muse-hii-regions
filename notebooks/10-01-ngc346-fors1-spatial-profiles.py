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

# # Spatial profiles along slit from FORS1 spectra
#
# I would like to achieve the following:
#
# 1. [X] Identify W 3 and the bow shock in the spectra
# 2. [ ] Trace the full western extent of the bow shock in [Ar IV]
#     - the MUSE field extends only 22 arcsec from the star
#     - the FORS1 spectrum shows a roughly linear ramp that extends about 34 arcsec (see image below)
# 3. [ ] Calculate the [O III] 4363/5007 temperature profile of the bow shock
#     - See if there is any evidence for temperatures as high as the [Ar IV] temperature
# 4. [ ] Measure the [Ne III] 3869 profile.  Ne$^+$ has the same ionization potential as Ar$^{+2}$ (40 eV), but [Ne III] seems to be much more broadly distributed than [Ar IV]
#
#
#
# ![FORS1 profile of Ar IV](assets/ngc346-screenshot-fors1-ariv-profile.png)
#
#

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import regions
from astropy.io import fits
from astropy.wcs import WCS
import sys

sys.path.append("../lib")
import extract

sns.set_context("talk")
sns.set_color_codes()
# -

# Load line and continuum images for the two blue spectral ranges.  We may use the continuum for masking out regions that have excessive stellar contamination.

hdu1 = fits.open("../data/ngc346-fors1-A-3500-4600-contsub.fits")[0]
hdu2 = fits.open("../data/ngc346-fors1-A-4400-6000-contsub.fits")[0]
hdu1c = fits.open("../data/ngc346-fors1-A-3500-4600-cont.fits")[0]
hdu2c = fits.open("../data/ngc346-fors1-A-4400-6000-cont.fits")[0]

hdu1B = fits.open("../data/ngc346-fors1-B-3562-4600-contsub.fits")[0]
hdu2B = fits.open("../data/ngc346-fors1-B-4400-6000-contsub.fits")[0]
hdu1cB = fits.open("../data/ngc346-fors1-B-3562-4600-cont.fits")[0]
hdu2cB = fits.open("../data/ngc346-fors1-B-4400-6000-cont.fits")[0]

# Set up the world coords for the two ranges

wcs1 = WCS(hdu1)
wcs2 = WCS(hdu2)
wcs1B = WCS(hdu1B)
wcs2B = WCS(hdu2B)

# Set up the lines that we want with wavelength limits:

restwav = {
    "He II 4686": 4685.68,
    "[Ar IV] 4740": 4740.17,
    "[O III] 5007": 5006.84,
    "[O III] 4363": 4363.209,
    "[Ne III] 3869": 3869.07,
    "H I 4861": 4861.32,
    "H I 4340": 4340.463,
    "H I 4102": 4101.735,
    "He I 5876": 5875.62,
    "[Cl III] 5518": 5517.71,
    "[Cl III] 5538": 5537.88,
    "[N II] 5755": 5755.08,
    "[Fe III] 4658": 4658.10,
    "[Fe III] 4987": 4987.20,
    "[Fe III] 5270": 5270.4,
    "O II 4650": 4650.00,
}

vsys = 160.0
vlim = np.array([50.0, 400.0])
restwav["He II 4686"] * (1.0 + vlim / 3e5)
restwav["[O III] 5007"] * (1.0 + vlim / 3e5)

dwav = wcs2.wcs.cdelt[0]


# +
def extract_line(wav0, pvim, wcs, vlim=[-50, 500]):
    """
    Integrated line flux along the slit
    """
    wavrange = wav0 * (1.0 + np.array(vlim) / 3e5)
    imwin, ww = extract.pvslice(pvim, wcs, wavrange, None)
    return imwin.sum(axis=1)

def extract_ew(wav0, pvim, pvcont, wcs, vlim=[-50, 500]):
    """
    Line equivalent width (angstrom units) along the slit
    """
    wavrange = wav0 * (1.0 + np.array(vlim) / 3e5)
    imwin, ww = extract.pvslice(pvim, wcs, wavrange, None)
    imcont, _ = extract.pvslice(pvcont, wcs, wavrange, None)
    dwav = wcs.wcs.cdelt[0]
    return dwav * (imwin / imcont).sum(axis=1)


# -

sA = {}
ewA = {}
for label, wav0 in restwav.items():
    if wav0 > 4400.0:
        pvim, pvcont, wcs = hdu2.data, hdu2c.data, wcs2
    else:
        pvim, pvcont, wcs = hdu1.data, hdu1c.data, wcs1
    sA[label] = extract_line(wav0, pvim, wcs)
    ewA[label] = extract_ew(wav0, pvim, pvcont, wcs)

ewA

sB = {}
ewB = {}
for label, wav0 in restwav.items():
    if wav0 > 4400.0:
        pvim, pvcont, wcs = hdu2B.data, hdu2cB.data, wcs2B
    else:
        pvim, pvcont, wcs = hdu1B.data, hdu1cB.data, wcs1B
    sB[label] = extract_line(wav0, pvim, wcs)
    ewB[label] = extract_ew(wav0, pvim, pvcont, wcs)    

ny = len(sA["He II 4686"])
_, posA = wcs1.pixel_to_world_values(
    [0] * ny,
    np.arange(ny),
)
ny = len(sB["He II 4686"])
_, posB = wcs1.pixel_to_world_values(
    [0] * ny,
    np.arange(ny),
)
posA *= -1.0
posB *= -1.0

# +
fig, axes = plt.subplots(
    2,
    1,
    figsize=(12, 12),
    sharex=True,
)
axes[0].plot(posA, sA["[Ar IV] 4740"], linewidth=1.0, alpha=1.0)
axes[0].plot(posA, sA["He II 4686"], linewidth=2.0, alpha=1.0)
axes[0].axhline(0.0, linestyle="dashed", color="k")
axes[0].axvline(0.0, linestyle="dotted", color="k")

axes[1].plot(posB, sB["[Ar IV] 4740"], linewidth=1.0, alpha=1.0)
axes[1].plot(posB, sB["He II 4686"], linewidth=2.0, alpha=1.0)
axes[1].axhline(0.0, linestyle="dashed", color="k")
axes[1].axvline(0.0, linestyle="dotted", color="k")

for ax in axes:
    ax.set(
        xlim=[-25, 75],
        ylim=[-30, 60],
    )

# +
fig, axes = plt.subplots(2, 1, figsize=(12, 12))
axes[0].plot(posA, 1e-3 * sA["[O III] 5007"], linewidth=1.0, alpha=1.0)
axes[0].plot(posA, 0.075 * sA["[O III] 4363"], linewidth=1.0, alpha=1.0)
axes[0].plot(posA, sA["He II 4686"], linewidth=2.0, alpha=1.0)
axes[1].plot(posB, 1e-3 * sB["[O III] 5007"], linewidth=1.0, alpha=1.0)
axes[1].plot(posB, 0.075 * sB["[O III] 4363"], linewidth=1.0, alpha=1.0)
axes[1].plot(posB, sB["He II 4686"], linewidth=2.0, alpha=1.0)

for ax in axes:
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.set(
        xlim=[-25, 75],
        ylim=[-30, 60],
    )
# -

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, s in zip(axes, [posA, posB], [sA, sB]):
    ax.plot(pos, s["[O III] 4363"] / s["[O III] 5007"], linewidth=3.0, alpha=1.0)
    ax.plot(pos, 0.0003 * s["[O III] 5007"], linewidth=2.0, alpha=1.0)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.set(
        xlim=[-10, 50],
        ylim=[-0.003, 0.02],
    )

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, s in zip(axes, [posA, posB], [sA, sB]):
    ax.plot(
        pos,
        s["[O III] 4363"]
        / s["[O III] 5007"]
        / np.median(s["[O III] 4363"] / s["[O III] 5007"]),
        linewidth=0.3,
        alpha=1.0,
        label="[O III] 4363 / 5007",
    )
    ax.plot(
        pos,
        s["[Ne III] 3869"]
        / s["[O III] 5007"]
        / np.median(s["[Ne III] 3869"] / s["[O III] 5007"]),
        linewidth=0.3,
        alpha=1.0,
        label="[Ne III] 3869 / [O III] 5007",
    )
    ax.plot(
        pos,
        0.1 * s["[O III] 5007"] / np.median(s["[O III] 5007"]),
        #    0.1* s4341 / np.median(s4341),
        linewidth=1.0,
        alpha=1.0,
        label="[O III] 5007",
    )
    ax.plot(
        pos,
        s["H I 4340"] / s["H I 4861"] / np.median(s["H I 4340"] / s["H I 4861"]),
        linewidth=0.5,
        alpha=1.0,
        label="Hγ / Hβ",
    )
    # ax.plot(positions, 0.0003 * s4686, linewidth=2.0, alpha=1.0)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")

    ax.legend(fontsize="x-small")
    ax.set(
        # xlim=[-10, 50],
        xlim=[-230, 200],
        ylim=[-0.01, 1.5],
    )

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line in "[O III] 5007", "H I 4861":
        ax.plot(pos, ew[line], linewidth=3.0, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        #xlim=[-10, 50],
        #ylim=[-0.003, 0.02],
    );

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line in "He I 5876", "[Ne III] 3869":
        ax.plot(pos, ew[line], linewidth=3.0, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        #xlim=[-10, 50],
        #ylim=[-0.003, 0.02],
    );

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line in "[O III] 4363", "[N II] 5755":
        ax.plot(pos, ew[line], linewidth=3.0, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        #xlim=[-10, 50],
        #ylim=[-0.003, 0.02],
    );

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    m = ew["[O III] 5007"] > 250.0
    for line in "[Ar IV] 4740", "He II 4686":
        y = ew[line][:]
        y[~m] = np.nan
        ax.plot(pos, y, linewidth=1.5, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        xlim=[-10, 50],
        ylim=[-2, 3],
    );

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, s, ew in zip(axes, [posA, posB], [sA, sB], [ewA, ewB]):
    m = ew["[O III] 5007"] > 250.0
    for line in "[Fe III] 4658", "O II 4650":
        y = s[line][:]
        y[~m] = np.nan
        ax.plot(pos, y, linewidth=0.5, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        xlim=[-100, 50],
        ylim=[-20, 40],
    );

a = np.arange(4)
np.kron(a, np.ones((4,)))

sys.path.append("../../multibin-maps")
import rebin_utils

fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line, color in zip(["[Ar IV] 4740", "He II 4686"], "cr"):
        m = ew["[O III] 5007"] > 250.0
        x, y = pos[:], ew[line][:]
        w = np.ones_like(y)
        lw = 0.5
        for n in [1, 2, 4, 8]:
            ax.plot(
                x, y, 
                linewidth=lw, alpha=1.0, label=line, color=color, drawstyle="steps-mid",
            )
            [x, y], m, w = rebin_utils.downsample1d([x, y], m, weights=w)
            y[~m] = np.nan
            lw += 0.5
            line="_nolabel"
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend(fontsize="x-small")
    ax.set(
        xlim=[-20, 50],
        #xlim=[-100, 50],
        ylim=[-2, 3],
    );









# +
fig, ax = plt.subplots(figsize=(12, 12))
hb_hg = s4861 / s4341
hg_hd = s4341 / s4102

ax.scatter(hb_hg, hg_hd, s=300, linewidth=0, alpha=0.1, marker=".")
ax.set(
    xlim=[1.5, 2.5],
    ylim=[1.0, 2.0],
)
