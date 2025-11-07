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

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# # Notebook F of more line ratios from saved maps from the PZ cubes
#
# This is adapted from 01-02-yet-more-line-ratios. It includes the He II analysis
#


# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
from zz_utils import get_image_raw, ROOT
import regions
import sys
import pandas as pd
import cmasher as cmr
import pyneb as pn

sns.set_context("talk")
sns.set_color_codes()


# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# ## He II emission measure

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
he2 = pn.RecAtom("He", 2)
he1 = pn.RecAtom("He", 1)
h1 = pn.RecAtom("H", 1)


# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
temperatures = [12500, 13800, 15500]
e4686 = he2.getEmissivity(tem=temperatures, den=1, wave=4686)
e4861 = h1.getEmissivity(tem=temperatures, den=1, wave=4861)
e5875 = he1.getEmissivity(tem=temperatures, den=1, wave=5876)
e4686, e4861, e5875
# -

np.mean(e4861[1])

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# From the other notebook, we measure 5875 / 4861 = 0.108 +/- 0.001
#
# Whereas the ratio of emissivities is

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e5875 / e4861

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# This implies a He abundance of

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
np.array([0.107, 0.108, 0.109])[None, :] / (e5875 / e4861)[:, None]

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# From Table 5 of Valerdi+ (2019), they have a He+ abundance of 10.915 or 10.917, with error of +/- 0.004 , so:

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
10 ** (np.array([10.915, 10.917]) - 12.0)

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# They find a small He++ abundance of

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
f"He++/H+ = {10 ** (8.30 - 12.0):.3e}; He++/He+ = {10 ** (8.30 - 10.915):.3e}"

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# In the bow shock I measure 4686 / 4861 = 0.015 or so, which would be 1.5 on a scale of Hβ = 100.  Mabel's Table 2 gives 0.24, which is 6 times less. This is consistent with the general value we find away from the bow's inner edge.
#
# So, I can work out my own He++/H+ abundance:

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
#

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
y_heiii_hii = 0.015 / (e4686 / e4861)
y_heiii_heii = y_heiii_hii / 0.08167471
f"Bow shock He++/H+ = {np.round(y_heiii_hii, 4)}; He++/He+ = {np.round(y_heiii_heii, 4)}"

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# In other words, the He++ emission measure, $\int n(\mathrm{He^{++}})\, n_\mathrm{e}\, dz$, is ony 2.5% of the total $\int n(\mathrm{He^{+}})\, n_\mathrm{e}\, dz$.
#
# *This is totally consistent, with the observed small change in 5875/4861, which is also about 2%*
#
# **TODO**
#
# - [ ] By assuming a path length, we can derive a density
# - [ ] We can do the same with [Ar IV] and [Ar III]
# - [ ] Ask Mabel to look at her slit A to see if she sees the He II and [Ar IV] signatures of the bow shock.  Also, maybe measure [O III] temperature to see if it goes up

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# ### Conversion to a density
#
# We want the surface brightness of He++ in physical units.
#
# MUSE flux units are $10^{-20}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ Å^{-1}\ pix^{-1}}$ in the cube, but we have summed over wavelength pixels, which are 1.4 Å.  The spatial pixels are 0.2 arcsec.

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import astropy.units as u

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
muse_bright_unit = 1e-20 * 1.4 * u.erg / u.s / u.cm**2 / (0.2 * u.arcsec) ** 2
muse_bright_unit.to(u.erg / u.s / u.cm**2 / u.sr)

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
datadir = ROOT / "data"
figdir = ROOT / "figs"

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
im4686 = Image(str(datadir / "ngc346-ZZ-heii-4686-correct.fits"))
im5875 = Image(str(datadir / "ngc346-ZZ-hei-5875-correct.fits"))
imhb = Image(str(datadir / "ngc346-ZZ-hi-4861-correct.fits"))
imariv = Image(str(datadir / "ngc346-ZZ-ariv-4711-plus-4740-correct.fits"))
imariii = Image(str(datadir / "ngc346-ZZ-ariii-7136-correct.fits"))
imoiii = Image(str(datadir / "ngc346-ZZ-oiii-5007-correct.fits"))
# -

# Use the same bow shock window as in the E notebook

xslice, yslice = slice(230, 300), slice(144, 245)

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 10))
(im4686 / imhb)[yslice, xslice].plot(vmin=0.0, vmax=0.01, colorbar="v")
ax_x = ax.secondary_xaxis(-0.06, functions=(lambda x: x + xslice.start, lambda x: x - xslice.start), color="r")
ax_y = ax.secondary_yaxis(-0.12, functions=(lambda y: y + yslice.start, lambda y: y - yslice.start), color="r")

                          

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"} active=""
# yslice

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
xxslice = slice(None, None)
# yyslice = slice(164, 204) # original
# yyslice = slice(160, 210) # broader
# yyslice = slice(170, 200) # narrower
yyslice = slice(180, 200)  # top half ultra narrow


def make_profile(im):
    # return np.make(im[yyslice, xxslice].data, axis=0)
    return np.mean(im[yyslice, xxslice].data, axis=0)


heii_profile = make_profile(im4686)
hei_profile = make_profile(im5875)
hb_profile = make_profile(imhb)
ariv_profile = make_profile(imariv)
ariii_profile = make_profile(imariii)
oiii_profile = make_profile(imoiii)

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 4))
ix0 = 227.5
nx = len(heii_profile)
pos = (np.arange(nx) - ix0) * 0.2

def _norm(y, x, x1=-12, x2=-7):
    m = (x >= x1) & (x <= x2)
    return y / np.mean(y[m])


ax.plot(pos, _norm(heii_profile, pos, x1=2, x2=5) / 2, label="He II", lw=4)
ax.plot(pos, _norm(ariv_profile, pos) / 4, label="[Ar IV]", lw=3)
fac = 3 * 0.0014
ax.plot(pos, _norm(oiii_profile, pos), label="[O III]", lw=2.0)
fac = 3 * 0.145
ax.plot(pos, _norm(ariii_profile, pos), label="[Ar III]", lw=1.5)
fac = 3 * 0.100
ax.plot(pos, _norm(hei_profile, pos), label="He I", lw=1.0)
fac = 3 * 0.0105
ax.plot(pos, _norm(hb_profile, pos), label="Hβ", lw=0.5)


ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=3, loc="upper left", fontsize="x-small")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Surface brightness",
    xlim=[-12, 22],
    ylim=[-0.3, 2.5],
)
sns.despine()
fig.savefig(figdir / "ngc346-ZZ-bow-shock-brightness-cuts.pdf")

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# From the profile graph above, the peak He II brightness is about 400 MUSE brightness units
# -

peak_heii_muse = np.max(heii_profile[(pos > 0) & (pos < 10)])
peak_heii_muse

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
peak_heii = muse_bright_unit * peak_heii_muse
peak_heii
# -

mshell = (abs(pos) < 10) & (heii_profile > 100)
heii_profile[mshell]

pos[mshell]

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# From the image, the chord length through the bow is about 60 pixels.  We can assume that the line-of-sight depth is similar:

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
depth_heii = 60 * 0.2 * 61700 * u.au
depth_heii.to(u.pc)

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# Surface brightness assuming optically thin emission with isotropic line emissivity, $e(\lambda)$, is given by
# $$
# I(\lambda) = \int \frac{e(\lambda)\, n_e\, n_i}{4 \pi} \, dz
# $$
# where $e(\lambda)$ is in the units given by pyneb: `u.erg * u.cm**3 / u.s` and $n_e$, $n_i$ are the electron and ion densities.
#
# Assume neutral fractions of He and H are negligible and hydrogen number density is $n$.  If the He abundance is $y = n(\mathrm{He}) / n(\mathrm{H})$ and the He++ ion fraction is $x_{++}$, then we have:
# $$
# n(\mathrm{He^{++}}) = y\, x_{++}\, n \quad \text{and} \quad n_e = [1 + y\, (1 + x_{++})]\, n
# $$
# implying that
# $$
# n(\mathrm{He^{++}}) \, n_e = n_e^2 \, \frac{y\, x_{++}}{1 + y\, (1 + x_{++})}
# $$
#
# So, with homogeneous conditions, we have
# $$
# I(4686) =
# \frac{e(4686)}{4\pi}\,
# \frac{y\, x_{++}}{1 + (1 + x_{++}) y}
# \, n_e^2 \, \delta z
# $$
# which can be solved for density to yield
# $$
# n_e = \left[
# \frac{4\pi\, I(4686)}{\delta z\, e(4686)} \, \frac{1 + (1 + x_{++}) y}{x_{++}\,y}
# \right]^{1/2}
# $$

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# We can take the helium abundance from above and get ...

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
pn_e_units = u.erg * u.cm**3 / u.s
yHe = 0.0824
ne = np.sqrt(
    4
    * np.pi
    * u.sr
    * peak_heii
    / (depth_heii.cgs * e4686 * pn_e_units)
    * (1 + 2 * yHe)
    / yHe
).to(u.cm**-3)
ne

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# **So electron density of 11 pcc!**
#
# Note, however that this assumes that the helium is 100% doubly ionized in the 4686 emitting region. If it is only partially ionized, then this is a lower limit (density would scale approximately as $x_{++}^{-1/2}$).
#
# The He II / Hb ratio is so low that the the ion fraction must be small. 

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# #### Find the He II flux and the He++ ionizing luminosity
#
#

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
im4686[160:210, 235:255].plot(vmin=0, vmax=400)

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
muse_flux_unit = 1e-20 * 1.4 * u.erg / u.s / u.cm**2
# -

# Calculate the intrinsic flux of the He II line, taking into account the foreground extinction (calculated from the Balmer decrement in the previous notebook)

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
cutout = im4686[160:210, 235:255]
m = (cutout.data > 0.0) & ~cutout.mask
A_heii = 0.34
F_heii = muse_flux_unit * np.sum(cutout.data[m]) * 10 ** (0.4 * A_heii)
F_heii

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
D_lmc = 61700 * u.pc

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
L_heii = 4 * np.pi * D_lmc.cgs**2 * F_heii
L_heii

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
L_heii.to(u.solLum)

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# Effective recomb rate:

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import astropy.constants as constants

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
hnu4686 = (constants.h * constants.c / (4686 * u.Angstrom)).cgs
hnu4686

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
alpha_eff_4686 = e4686 * pn_e_units / hnu4686
alpha_eff_4686

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
pn.atomicData.getAllAvailableFiles("He2")

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# It only works as follows:

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
pn.atomicData.setDataFile("he_ii_trc_SH95-caseB.dat")
alphaB_He_plus = pn.RecAtom("He", 2).getTotRecombination(tem=temperatures, den=100)
alphaB_He_plus *= u.cm**3 / u.s
alphaB_He_plus

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# Solid angle:
#
# Originally, I had assumed a +
# /- 75 degree end cap. But after doing some estimates in the org file, I find that 54 +/- 16 deg is a better estimate

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
Omega_over_4pi = (1 - np.cos(54 * u.deg)) / 2
Omega_over_4pi

# + [markdown] jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"}
# $$
# Q \frac{\Omega}{4\pi} = \int_{\mathcal{V}} n_e \, n_i \, \alpha_B\, d\mathcal{V}
# $$
# and
# $$
# L(4686) = \int_{\mathcal{V}} n_e \, n_i \, e(4686)\, d\mathcal{V}
# $$
# so that
# $$
# Q = \frac{\alpha_B \, L(4686)} {e(4686)\, (\Omega/4\pi)}
# $$

# + jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
Q2 = alphaB_He_plus * L_heii / (e4686 * pn_e_units) / Omega_over_4pi
Q2
# -

# This is way higher than what we got last time, which is due to including the effects of foreground extinction plus a smaller esitimate of the covering fraction. The uncertainty is +/- 50%, which is dominated by the uncertainty in the covering fraction.
