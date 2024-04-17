# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# # Further line analysis from saved moment images
#
# The previous notebook has grown too long, so I am starting a new one.  What I plan to do here is:
#
# 1. [X] [Cl III] density ratio. *Conclusion:* [Cl III] is useless
# 2. [X] [S II] density - this is the best
# 3. [X] He II analysis – ionizing luminosity and density
#     *This has a mistake at the moment, which I need to track down* Now fixed in Org file
# 4. [X] [Ar IV] and [Ar III] ionization balance
# 2. [ ] Maybe tetrablok binning and make figure
# 3. [X] Profile cuts across the bow shock
# 4. [ ] Analysis of average line velocities
#

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
import regions
import pandas as pd
import cmasher as cmr
import pyneb as pn

sns.set_context("talk")
sns.set_color_codes()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
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


# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
def trim_edges(im, m):
    """Trim m pixels of each edge of image in place by setting mask"""
    im.mask[:m, :] = True
    im.mask[-m:, :] = True
    im.mask[:, :m] = True
    im.mask[:, -m:] = True
    return None


# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
datadir = Path.cwd().parent.parent / "data"
figdir = Path.cwd().parent.parent / "figs"

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im6716 = Image(str(datadir / "ngc346-sii-6716-bin01-sum.fits"))
im6731 = Image(str(datadir / "ngc346-sii-6731-bin01-sum.fits"))
imha = Image(str(datadir / "ngc346-hi-6563-bin01-sum.fits"))
imcont = Image(str(datadir / "ngc346-cont-6312-mean.fits"))
imha

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
s2 = pn.Atom("S", 2)
e6716 = s2.getEmissivity(tem=13000, den=[3, 10, 30, 100, 300, 1000], wave=6716)
e6731 = s2.getEmissivity(tem=13000, den=[3, 10, 30, 100, 300, 1000], wave=6731)
Rgrid = e6716 / e6731
Rgrid

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
m = imcont.data > 5e3
im6716.mask = im6716.mask | m
im6731.mask = im6731.mask | m
trim_edges(im6716, 10)
trim_edges(im6731, 10)
trim_edges(imcont, 10)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 1
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
im6716.rebin(n).plot(vmin=-10, vmax=12000, ax=axes[0, 0], colorbar="v")
im6731.rebin(n).plot(vmin=-10, vmax=12000, ax=axes[0, 1], colorbar="v")
imcont.rebin(n).plot(vmin=0, vmax=1e4, ax=axes[1, 0], colorbar="v")
(im6716.rebin(n) / im6731.rebin(n)).plot(
    ax=axes[1, 1],
    vmin=0.0,
    vmax=2.0,
    cmap="gray",
    colorbar="v",
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 1

imx = im6731.rebin(n)
imy = (im6716 / 1.44).rebin(n)

imin, imax = -10, 30000
m = (imx.data < imax) & (imy.data < imax)
m = m & (imx.data > imin) & (imy.data > imin)
m = m & (imy.data / imx.data > 0.5) & (imy.data / imx.data < 1.5)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6731": imx.data[m],
        "6716": imy.data[m],
        "6716 / 6731": imy.data[m] / imx.data[m],
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
g.axes[1, 0].plot([imin, imax], [imin, imax], "--", color="r")
g.axes[1, 0].plot([imin, imax], [imin, imax], "--", color="r")
g.fig.suptitle("Correlation between [S II] 6731 and 6716 brightness")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 1
imx = imha.rebin(n)
imy = im6716.rebin(n) / im6731.rebin(n)

m = (imy.data > 0.8) & (imy.data < 1.6)
m = m & (imx.data < 1.5e6)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": np.log10(imx.data[m]),
        "6716 / 6731": np.log10(imy.data[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=150,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=150,
    ),
)
for R in Rgrid:
    g.axes[1, 0].axhline(np.log10(R), color="r")
g.fig.suptitle("Correlation between [S II] 6716 / 6731 sum and Ha brightness")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 1
imx = imha.rebin(n)
imy = im6716.rebin(n) / im6731.rebin(n)

m = (imy.data > 0.8) & (imy.data < 1.6)
m = m & (imx.data < 1.5e6)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": np.log10(imx.data[m]),
        "6716 / 6731": np.log10(imy.data[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=150,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=150,
    ),
)
for R in Rgrid:
    g.axes[1, 0].axhline(np.log10(R), color="r")
g.fig.suptitle("Correlation between [S II] 6716 / 6731 sum and Ha brightness")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
xslice, yslice = slice(230, 300), slice(144, 245)

n = 1

imx = imha[yslice, xslice].rebin(n)
imy = im6716[yslice, xslice].rebin(n) / im6731[yslice, xslice].rebin(n)

m = (imy.data > 1.1) & (imy.data < 2)
m = m & (imx.data < 1e6)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": np.log10(imx.data[m]),
        "6716 / 6731": np.log10(imy.data[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=100,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=100,
    ),
)
g.axes[1, 1].axvline(np.log10(Rgrid[0]), color="r")
for R in Rgrid[:-1]:
    g.axes[1, 0].axhline(np.log10(R), color="r")

g.fig.suptitle("Correlation between [S II] 6716 / 6731 ratio and Ha brightness")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
xslice, yslice = slice(230, 300), slice(144, 245)

n = 1

imx = imha[yslice, xslice].rebin(n)
imy = im6716[yslice, xslice].rebin(n) / im6731[yslice, xslice].rebin(n)

m = (imy.data > 1.1) & (imy.data < 2)
m = m & (imx.data < 80000)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": imx.data[m],
        "6716 / 6731": imy.data[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=100,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=100,
    ),
)
g.axes[1, 1].axvline(Rgrid[0], color="r")
for R in Rgrid[:-1]:
    g.axes[1, 0].axhline(R, color="r")

g.fig.suptitle("Correlation between [S II] 6716 / 6731 ratio and Ha brightness")

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# So, in the bow shock region, we see ratios as low as 1.3 in the brightest parts, but these are globule surfaces.  The bulk of the emission has around 1.4

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
s2.getTemDen([1.4, 1.3, 0.8], tem=12000, wave1=6716, wave2=6731)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# So the density is about 50 +/- 30 pcc in the diffuse gas.  We get ten times higher density in the case of Source E, which has a ratio as low as 0.8

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# ### Make a map of [S II] density

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
r_s2_grid = np.linspace(0.5, 1.44, 1001)
n_s2_grid = s2.getTemDen(r_s2_grid, tem=12000.0, wave1=6716, wave2=6731)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n_s2_grid

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
iew6716 = imcont / (1.25 * im6716)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(10, 10))
iew6716.plot(colorbar="v", cmap="gray_r", scale="linear", vmin=-1.0, vmax=10.0)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fixmask = im6716.mask | (iew6716.data > 10.0) | (iew6716.data < -0.2)
fixmask[90:97, 147:152] = True
fixmask[79:86, 191:197] = True


im_n_sii = im6716.clone(data_init=np.empty)
im_n_sii.mask = fixmask
trim_edges(im_n_sii, 10)
im_n_sii.data[~fixmask] = np.interp(
    im6716.data[~fixmask] / im6731.data[~fixmask],
    r_s2_grid,
    n_s2_grid,
    left=np.nan,
    right=np.nan,
)
im_n_sii.mask = im_n_sii.mask | ~np.isfinite(im_n_sii.data)

# + jupyter={"source_hidden": true, "outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 12))
im_n_sii.rebin(2).plot(colorbar="v", cmap="gray_r", scale="sqrt", vmin=0.0, vmax=3000.0)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# This seems to be good enough in some of the diffuse regions. although it is way to noisy in the faint parts.

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im_n_sii.write(str(datadir) + "/ngc346-N-sii.fits", savemask="nan")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im_T_siii = Image(str(datadir / "ngc346-T-siii.fits"))
imhb = Image(str(datadir / "ngc346-hi-4861-correct.fits"))

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 16
imx = im_T_siii.rebin(n)
imy = im_n_sii.rebin(n)
imz = imhb.rebin(n)

m = ~imx.mask & ~imy.mask
m = m & (imx.data > 5000) & (imx.data < 20000)
m = m & (imy.data > 0) & (imy.data < 300)
df = pd.DataFrame(
    {
        "T([S III])": imx.data[m],
        "n([S II])": imy.data[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imz.data[m],
        bins=50,
    ),
    diag_kws=dict(
        weights=imz.data[m],
        bins=50,
    ),
)
g.fig.suptitle("Temperature vs Density")
# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
#     ## He II emission measure

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
he2 = pn.RecAtom("He", 2)
he1 = pn.RecAtom("He", 1)
h1 = pn.RecAtom("H", 1)


# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
temperatures = [12500, 13800, 15500]
e4686 = he2.getEmissivity(tem=temperatures, den=1, wave=4686)
e4861 = h1.getEmissivity(tem=temperatures, den=1, wave=4861)
e5875 = he1.getEmissivity(tem=temperatures, den=1, wave=5876)
e4686, e4861, e5875
# -

np.mean(e4861[1])

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# From the other notebook, we measure 5875 / 4861 = 0.108 +/- 0.001
#
# Whereas the ratio of emissivities is

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e5875 / e4861

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# This implies a He abundance of

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
np.array([0.107, 0.108, 0.109])[None, :] / (e5875 / e4861)[:, None]

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# From Table 5 of Valerdi+ (2019), they have a He+ abundance of 10.915 or 10.917, with error of +/- 0.004 , so:

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
10 ** (np.array([10.915, 10.917]) - 12.0)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# They find a small He++ abundance of

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
f"He++/H+ = {10**(8.30 - 12.0):.3e}; He++/He+ = {10**(8.30 - 10.915):.3e}"

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# In the bow shock I measure 4686 / 4861 = 0.015 or so, which would be 1.5 on a scale of Hβ = 100.  Mabel's Table 2 gives 0.24, which is 6 times less. This is consistent with the general value we find away from the bow's inner edge.
#
# So, I can work out my own He++/H+ abundance:

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
#

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
y_heiii_hii = 0.015 / (e4686 / e4861)
y_heiii_heii = y_heiii_hii / 0.08167471
f"Bow shock He++/H+ = {np.round(y_heiii_hii, 4)}; He++/He+ = {np.round(y_heiii_heii, 4)}"

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# In other words, the He++ emission measure, $\int n(\mathrm{He^{++}})\, n_\mathrm{e}\, dz$, is ony 2.5% of the total $\int n(\mathrm{He^{+}})\, n_\mathrm{e}\, dz$.
#
# *This is totally consistent, with the observed small change in 5875/4861, which is also about 2%*
#
# **TODO**
#
# - [ ] By assuming a path length, we can derive a density
# - [ ] We can do the same with [Ar IV] and [Ar III]
# - [ ] Ask Mabel to look at her slit A to see if she sees the He II and [Ar IV] signatures of the bow shock.  Also, maybe measure [O III] temperature to see if it goes up

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# ### Conversion to a density
#
# We want the surface brightness of He++ in physical units.
#
# MUSE flux units are $10^{-20}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ Å^{-1}\ pix^{-1}}$ in the cube, but we have summed over wavelength pixels, which are 1.4 Å.  The spatial pixels are 0.2 arcsec.

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
import astropy.units as u

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
muse_bright_unit = 1e-20 * 1.4 * u.erg / u.s / u.cm ** 2 / (0.2 * u.arcsec) ** 2
muse_bright_unit.to(u.erg / u.s / u.cm ** 2 / u.sr)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im4686 = Image(str(datadir / "ngc346-heii-4686-correct.fits"))
im5875 = Image(str(datadir / "ngc346-hei-5875-correct.fits"))
imhb = Image(str(datadir / "ngc346-hi-4861-correct.fits"))
imariv = Image(str(datadir / "ngc346-ariv-4711-plus-4740-correct.fits"))
imariii = Image(str(datadir / "ngc346-ariii-7136-correct.fits"))
imoiii = Image(str(datadir / "ngc346-oiii-5007-bin01-sum.fits"))

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(10, 10))
(im4686 / imhb)[yslice, xslice].plot(vmin=0.0, vmax=0.016, colorbar="v")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
yslice

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
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

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(heii_profile)
pos = (np.arange(nx) - ix0) * 0.2
pos2 = (np.arange(len(oiii_profile)) - ix0) * 0.2

ax.plot(pos, heii_profile, label="He II", lw=4)
ax.plot(pos, 1.00 * ariv_profile, label="[Ar IV]", lw=3)
fac = 3 * 0.0014
ax.plot(pos2, fac * oiii_profile, label=f"[O III] / {1/fac:.1f}", lw=2.0)
fac = 3 * 0.145
ax.plot(pos, fac * ariii_profile, label=f"[Ar III] / {1/fac:.1f}", lw=1.5)
fac = 3 * 0.100
ax.plot(pos, fac * hei_profile, label=f"He I / {1/fac:.1f}", lw=1.0)
fac = 3 * 0.0105
ax.plot(pos, fac * hb_profile, label=f"Hβ / {1/fac:.1f}", lw=0.5)


ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=3, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Surface brightness",
    xlim=[-12, 22],
    ylim=[-199, 1400],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-brightness-cuts.pdf")

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# From the profile graph above, the peak He II brightness is about 400 MUSE brightness units

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
peak_heii = muse_bright_unit * 400
peak_heii
# -

mshell = (abs(pos) < 10) & (heii_profile > 200)
heii_profile[mshell]

pos[mshell]

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# From the image, the chord length through the bow is about 60 pixels.  We can assume that the line-of-sight depth is similar:

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
depth_heii = 60 * 0.2 * 61700 * u.au
depth_heii.to(u.pc)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
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
# \approx n_e^2 \, \frac{y}{1 + 2 y}
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

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# We can take the helium abundance from above and get ...

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
pn_e_units = u.erg * u.cm ** 3 / u.s
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

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# **So electron density of 11 pcc!**
#
# Note, however that this assumes that the helium is 100% doubly ionized in the 4686 emitting region. If it is only partially ionized, then this is a lower limit (density would scale approximately as $x_{++}^{-1/2}$).

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# #### Find the He II flux and the He++ ionizing luminosity
#
#

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im4686[160:210, 235:255].plot(vmin=0, vmax=400)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
muse_flux_unit = 1e-20 * 1.4 * u.erg / u.s / u.cm ** 2
# -

# Calculate the intrinsic flux of the He II line, taking into account the foreground extinction (calculated from the Balmer decrement in the previous notebook)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
cutout = im4686[160:210, 235:255]
m = (cutout.data > 0.0) & ~cutout.mask
A_heii = 0.34
F_heii = muse_flux_unit * np.sum(cutout.data[m]) * 10**(0.4*A_heii)
F_heii

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
D_lmc = 61700 * u.pc

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
L_heii = 4 * np.pi * D_lmc.cgs ** 2 * F_heii
L_heii

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
L_heii.to(u.solLum)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# Effective recomb rate:

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
import astropy.constants as constants

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
hnu4686 = (constants.h * constants.c / (4686 * u.Angstrom)).cgs
hnu4686

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
alpha_eff_4686 = e4686 * pn_e_units / hnu4686
alpha_eff_4686

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
pn.atomicData.getAllAvailableFiles("He2")

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# It only works as follows:

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
pn.atomicData.setDataFile("he_ii_trc_SH95-caseB.dat")
alphaB_He_plus = pn.RecAtom("He", 2).getTotRecombination(tem=temperatures, den=100)
alphaB_He_plus *= u.cm ** 3 / u.s
alphaB_He_plus

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# Solid angle: 
#
# Originally, I had assumed a +
# /- 75 degree end cap. But after doing some estimates in the org file, I find that 54 +/- 16 deg is a better estimate

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
Omega_over_4pi = (1 - np.cos(54 * u.deg)) / 2
Omega_over_4pi

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
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

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
Q2 = alphaB_He_plus * L_heii / (e4686 * pn_e_units) / Omega_over_4pi
Q2
# -

# This is way higher than what we got last time, which is due to including the effects of foreground extinction plus a smaller esitimate of the covering fraction. The uncertainty is +/- 50%, which is dominated by the uncertainty in the covering fraction. 

# ## Alternative route to density from Balmer lines
#
# We can estimate the bow shock shell emission measure from the jump in the surface brightness of H beta, say. This has the advantage that we know that the hydrogen is fully ionized.

A_hb = 0.33
jump_hb = 31.7 * 300 * muse_bright_unit * 10**(0.4*A_hb)
f"{jump_hb=}"

# The 300 is my estimate from the above figure of the Hb jump at the bow shock inner edge, while the 31.7 is what we divided Hb by before plotting it. 
#
# Thhen we also correct for the dust extinction `A_hb`. 
#
# So this should be the intrinsic surface brightness of Hb, which should be EM times emission coefficient / 4 pi
#

EM = jump_hb * (4 * np.pi * u.sr).to(u.arcsec**2) / (e4861 * pn_e_units)
EM

# There are three answers, corrending to assumed temperature of 12,500, 13,800, or 15,500 K

np.sqrt(EM / depth_heii).to(u.cm**-3)

# So that gives a bigger density than I was getting before. On the other hand, it might be a slight over-estimate since the H+ path length may be largerby a factor of 2 (see next section). There is also a small correction for the ionized helium contribution to the electron density, but that is much smaller than the uncertainty in the path length. 
#
# So, we would finally get
#

np.sqrt(EM / (2 * depth_heii)).to(u.cm**-3)

# So, taking into account the uncertainty in the path length would give us a hydrogen density of 35 +/- 10 pcc

# ### Sanity check from He II / H I
#
# In principle, we can come to the same conclusion starting from the 4686/4861 ratio, so long as we account for the fraction of the H beta emission that comes from the bow shock.

# +
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(heii_profile)
pos = (np.arange(nx) - ix0) * 0.2
imargin = 10
ratio = heii_profile / hb_profile
ax.plot(pos[imargin:], ratio[imargin:], label="He II / H beta", lw=4)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=3, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Line ratio",
#    xlim=[-12, 22],
    ylim=[-0.005, 0.018],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-heii-hbeta-cuts.pdf")
# -

# So the peak ratio is about 0.015. Note that we already did this once above in the *He II emission measure* section. We got exactly the same value. 
#
# We now have to divide this by the fraction of H beta that comes frm the bow shock.

hb_bow_frac = 0.3
heii_hb_bow_max = np.max(ratio[imargin:]) / hb_bow_frac
heii_hb_bow_max

# Then we should have
# $$
# \texttt{heii\_hb\_bow\_max} = 
# \frac{I(4686)}{I(4681)} =
# \frac{
# \int e(4686)\, n(\mathrm{He^{++}})\, n_e\, dz
# }{
# \int e(4861)\, n(\mathrm{H^+})\, n_e\, dz
# }
# \approx
# \frac{e(4686)}{e(4861)}\,
# y(\mathrm{He})\,
# x(\mathrm{He}^{++})
# $$

e4686, e4861, e4686 / e4861

# So the He II emission coefficient is about 12 times larger than H beta, which almost exactly cancels out the He abundance factor. And the temperature dependence is very slight.

yhe = 0.08325295
xheiii = (heii_hb_bow_max * e4861) / (yhe * e4686)
xheiii

# Note that this assumes that the line-of-sight path length is the same for He++ and for H+. In reality it may be longer for H+ because the thickness of the emitting shell is greater. In the thin shell approximation, the ratio of path lengths is the sqrt of the ratio of thickness. 
#
# The ratio of thicknesses is between 2 and 5, depending on just how we measure it. So this would give a ratio of path lengths of about (2 +/- 1), which is also what you would have guessed by looking at the chord lengths on the plane of the sky (for [Ar IV] and He II). 
#
# This would give a final value of $x(\mathrm{He^{++}}) \approx 0.1$





# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# ## Also do profiles of low ionization lines

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im6300 = Image(str(datadir / "ngc346-oi-6300-bin01-sum.fits"))
im5518 = Image(str(datadir / "ngc346-cliii-5518-bin01-sum.fits"))
im5538 = Image(str(datadir / "ngc346-cliii-5538-bin01-sum.fits"))
im9069 = Image(str(datadir / "ngc346-siii-9069-bin01-sum.fits"))

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
oi_profile = make_profile(im6300)
cliiis_profile = make_profile(im5518)
cliiil_profile = make_profile(im5538)
siis_profile = make_profile(im6716)
siil_profile = make_profile(im6731)
siii_profile = make_profile(im9069)

sii_profile = 0.5 * (siis_profile + siil_profile)
cliii_profile = 0.5 * (cliiis_profile + cliiil_profile)


# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(oi_profile)
pos = (np.arange(nx) - ix0) * 0.2
pos2 = (np.arange(len(siil_profile)) - ix0) * 0.2

ax.plot(pos, 1.0 * oi_profile / np.median(oi_profile), label="[O I] brightness", lw=3)
ax.plot(
    pos2, 1.0 * sii_profile / np.median(sii_profile), label="[S II] brightness", lw=3
)
ax.plot(
    pos,
    1.7 * cliii_profile / np.median(cliii_profile),
    label="[Cl III] brightness",
    lw=3,
)
ax.plot(
    pos, 1.7 * siii_profile / np.median(siii_profile), label="[S III] brightness", lw=3
)


ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=2, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[0, 3.9],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-lowion-cuts.pdf")
# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}



# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# ## Ar ionization balance
#
# We will calculate the conditions at the very peak of the [Ar IV] emission.
#
# We need a very careful slection of the background (BG) and bow shock (BS) samples, since we want to make sure we are in the little triangle window where the intermediate ionization lines are not too contaminated by globule i-fronts and unrelated filaments.

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
i0, j0, w, h = 234, 200, 12, 8
bgbox = regions.RegionBoundingBox(
    iymin=j0 - h // 2,
    iymax=j0 + h // 2,
    ixmin=i0 - w // 2,
    ixmax=i0 + w // 2,
)
i0, j0, w, h = 250, 193, 8, 8
bsbox = regions.RegionBoundingBox(
    iymin=j0 - h // 2,
    iymax=j0 + h // 2,
    ixmin=i0 - w // 2,
    ixmax=i0 + w // 2,
)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
imariv.plot(ax=axes[0], vmin=0, vmax=700, colorbar="v")
imariii.plot(ax=axes[1], vmin=0, vmax=2500, colorbar="v")
(imariv / imariii).plot(ax=axes[2], vmin=0, vmax=0.45, cmap="magma_r", colorbar="v")
for ax in axes:
    bsbox.plot(ax=ax, color="w")
    bgbox.plot(ax=ax, color="w")
    ax.set(
        xlim=[200, 300],
        ylim=[100, 250],
    )
fig.tight_layout()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
bs_slices, _ = bsbox.get_overlap_slices(imariv.shape)
bg_slices, _ = bgbox.get_overlap_slices(imariv.shape)

bs_ariv = imariv[bs_slices].data.mean()
bg_ariv = imariv[bg_slices].data.mean()
sbs_ariv = imariv[bs_slices].data.std()
sbg_ariv = imariv[bg_slices].data.std()

bs_slices, _ = bsbox.get_overlap_slices(imariii.shape)
bg_slices, _ = bgbox.get_overlap_slices(imariii.shape)
bs_ariii = imariii[bs_slices].data.mean()
bg_ariii = imariii[bg_slices].data.mean()
sbs_ariii = imariii[bs_slices].data.std()
sbg_ariii = imariii[bg_slices].data.std()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
f"[Ar IV]: BS = {bs_ariv:.2f} +/- {sbs_ariv:.2f}, BG = {bg_ariv:.2f} +/- {sbg_ariv:.2f}"

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
f"[Ar III]: BS = {bs_ariii:.2f} +/- {sbs_ariii:.2f}, BG = {bg_ariii:.2f} +/- {sbg_ariii:.2f}"

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
BS_ariv = bs_ariv - bg_ariv
sBS_ariv = np.hypot(sbs_ariv, sbg_ariv)

BS_ariii = bs_ariii - bg_ariii
sBS_ariii = np.hypot(sbs_ariii, sbg_ariii)

BS_ar_iv_iii = BS_ariv / BS_ariii
sBS_ar_iv_iii = BS_ar_iv_iii * np.hypot(sBS_ariii / BS_ariii, sBS_ariv / BS_ariv)

f"BG-subtracted BS: [Ar IV] / [Ar III] = {BS_ar_iv_iii:.3f} +/- {sBS_ar_iv_iii:.3f}"

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
ar4 = pn.Atom("Ar", 4)
ar3 = pn.Atom("Ar", 3)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# First T is for BG nebula. Second and third are lower and upper limits for bow shock.  See analysis of [Ar IV] temperature in `10-01` notebook

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
Ts = [12500, 14000, 16000]
e4711 = ar4.getEmissivity(tem=Ts, den=10.0, wave=4711)
e4740 = ar4.getEmissivity(tem=Ts, den=10.0, wave=4740)
e7136 = ar3.getEmissivity(tem=Ts, den=10.0, wave=7136)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e4711 + e4740

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e7136

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# So the emissivity is very similar, strangely.
#
# Anyway, we should have:
# $$
# \frac {n(\mathrm{Ar^{3+}})} {n(\mathrm{Ar^{2+}})}
# =
# \frac{I(4711 + 4740)}{I(7136)}
# \, \frac{e(7136)}{e(4711 + 4740)}
# $$

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e7136 / (e4711 + e4740)

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# So the T uncertainty of +/- 1000 K would give +/- 10% uncertainty in the emissivity ratio. For the time being we take the middle value:

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
Ar3p_over_Ar2p = BS_ar_iv_iii * (e7136 / (e4711 + e4740))[1]
Ar3p_over_Ar2p

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# Or, close enough to unity.  So, at the inner edge of the bow shock we have 50% Ar++ and 50% Ar+++. We can use this to constrain the stellar spectrum if we run some Cloudy models
#
# Now do the same, but for the BG nebula

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
Ar3p_over_Ar2p_BG = (bg_ariv / bg_ariii) * (e7136 / (e4711 + e4740))[0]
Ar3p_over_Ar2p_BG

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# So this implies 20% Ar+++ in the BG region.
#
# However, these are both still integrals along the line of sight, so the actual variation in ionization might be larger.
#
# For instance, the Ar+++ fraction might be higher than 50% at the inner edge.  We could do a LOS integration on the Cloudy emissivities to investigate this.
#
# Likewise, the Ar+++ fraction in the BG might be lower than 0.2 since that could be due to contamination by the bow shock wing emission.

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# ## Can we get a He I density?

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
dgrid = [1.0, 10.0, 100.0, 1000.0]
T0 = [11000, 13000, 18000]
he1.getEmissivity(tem=T0, den=dgrid, wave=5876) / he1.getEmissivity(
    tem=T0, den=dgrid, wave=6678
)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
he1.getEmissivity(tem=T0, den=dgrid, wave=4922) / he1.getEmissivity(
    tem=T0, den=dgrid, wave=5876
)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
he1.getEmissivity(tem=T0, den=dgrid, wave=5048) / he1.getEmissivity(
    tem=T0, den=dgrid, wave=5876
)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im5048 = Image(str(datadir / "ngc346-hei-5048-bin01-sum.fits"))
im5876 = Image(str(datadir / "ngc346-hei-5875-bin01-sum.fits"))

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 1
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
im5048.rebin(n).plot(vmin=-10, vmax=60, ax=axes[0, 0], colorbar="v")
im5876.rebin(n).plot(vmin=-10, vmax=3000, ax=axes[0, 1], colorbar="v")
imcont.rebin(n).plot(vmin=0, vmax=1e4, ax=axes[1, 0], colorbar="v")
(im5048.rebin(n) / im5876.rebin(n)).plot(
    ax=axes[1, 1],
    vmin=0.0,
    vmax=0.03,
    cmap="magma",
    colorbar="v",
)
for ax in axes.flat:
    bsbox.plot(ax=ax, color="w")
    bgbox.plot(ax=ax, color="w")
    ax.set(
        xlim=[200, 300],
        ylim=[100, 250],
    )
fig.tight_layout()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
yyslice = slice(164, 204)  # original
# yyslice = slice(160, 210) # broader
# yyslice = slice(170, 200) # narrower
# yyslice = slice(180, 200) # top half ultra narrow
hei_5048_profile = make_profile(im5048)
hei_5876_profile = make_profile(im5876)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(hei_profile)
pos = (np.arange(nx) - ix0) * 0.2

ax.plot(pos, 0.01 * hei_5876_profile / np.median(hei_5876_profile), label="He I", lw=4)
ax.plot(pos, hei_5048_profile / hei_5876_profile, label="5048 / 5875", lw=3)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=3, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Surface brightness",
    xlim=[-12, 22],
    ylim=[0, 0.03],
)
sns.despine()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
bs_5048 = im5048[bs_slices].data.mean()
bs_5876 = im5876[bs_slices].data.mean()
bs_5048 / bs_5876

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
bg_5048 = im5048[bg_slices].data.mean()
bg_5876 = im5876[bg_slices].data.mean()
bg_5048 / bg_5876

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# These are all way lower than the theoretical values for reasonable temperatures.  Maybe the 5048 line is affected by underlying stellar absorption.

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
# yyslice = slice(164, 204) # original
# yyslice = slice(160, 210) # broader
# yyslice = slice(170, 200) # narrower
yyslice = slice(180, 200)  # top half ultra narrow
n_sii_profile = make_profile(im_n_sii)
T_siii_profile = make_profile(im_T_siii)
sii_profile = make_profile(im6731)
im9069 = Image(str(datadir / "ngc346-siii-9069-bin01-sum.fits"))
siii_profile = make_profile(im9069)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(sii_profile)
pos = (np.arange(nx) - ix0) * 0.2
pos2 = (np.arange(len(siii_profile)) - ix0) * 0.2

ax.plot(
    pos, 0.01 * n_sii_profile, ds="steps-mid", label="$n$([S II]) / 100 cm$^{-3}$", lw=2
)
ax.plot(
    pos2, 0.0001 * T_siii_profile, ds="steps-mid", label="$T$([S III]) / 10,000 K", lw=2
)
ax.plot(
    pos, 1.0 * sii_profile / np.median(sii_profile), label="[S II] brightness", lw=3
)
ax.plot(
    pos2, 1.8 * siii_profile / np.median(siii_profile), label="[S III] brightness", lw=3
)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=2, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[0, 3.45],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-sii-siii-ne-te.pdf")

# + [markdown] pycharm={"name": "#%% md\n"} tags=["temperature"] jupyter={"outputs_hidden": false}
# So the [S III] temperature is significantly larger than the [O III] temperature from the nebula. It is around 14000 K and has a drop towards W 3 coming from the east side, and then a step and a very constant region.  But the step occurs before the bow shock, so is probably unrelated.

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(sii_profile)
pos = (np.arange(nx) - ix0) * 0.2
pos2 = (np.arange(len(siii_profile)) - ix0) * 0.2

ax.plot(pos, 0.01 * n_sii_profile, label="$n$([S II]) / 100 cm$^{-3}$", lw=4)
ax.plot(pos2, 0.0001 * T_siii_profile, label="$T$([S III]) / 10,000 K", lw=3)
ax.plot(pos2, 1.0 * sii_profile[1:] / siii_profile, label="[S II] / [S III]", lw=3)
ax.plot(pos2, 4.0 * oi_profile / sii_profile[1:], label="[O I] / [S II]", lw=3)


ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=2, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[0, 2.1],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-sii-siii-ratio-ne-te.pdf")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
imhei_c = Image(str(datadir / "ngc346-hei-5875-correct.fits"))
imhi_c = Image(str(datadir / "ngc346-hi-4861-correct.fits"))
imheii_c = Image(str(datadir / "ngc346-heii-4686-correct.fits"))
imcont2 = Image(str(datadir / "ngc346-cont-4686-mean.fits"))


# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
m = imcont2.data > 500.0
for im in imhei_c, imheii_c, imhi_c:
    im.mask = im.mask | m

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
# yyslice = slice(164, 204) # original
yyslice = slice(160, 210)  # broader
# yyslice = slice(170, 200) # narrower
# yyslice = slice(180, 200) # top half ultra narrow
hei_c_profile = make_profile(imhei_c)
hi_c_profile = make_profile(imhi_c)
heii_c_profile = make_profile(imheii_c)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(hei_profile)
pos = (np.arange(nx) - ix0) * 0.2

ax.plot(
    pos,
    heii_c_profile / hi_c_profile,
    ds="steps-mid",
    label="He II λ4686 / H I λ4861",
    lw=2,
)
ax.plot(
    pos,
    hei_c_profile / hi_c_profile - 0.10,
    ds="steps-mid",
    label="(He I λ5875 / H I λ4861) – 0.1",
    lw=2,
)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=1, loc="upper right")
ax.set_yticks([0.0, 0.005, 0.010])
ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[-0.005, 0.015],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-he-ratios.pdf")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
im5518 = Image(str(datadir / "ngc346-cliii-5518-bin01-sum.fits"))
im5538 = Image(str(datadir / "ngc346-cliii-5538-bin01-sum.fits"))
imha = Image(str(datadir / "ngc346-hi-6563-bin01-sum.fits"))
imcont = Image(str(datadir / "ngc346-cont-4686-mean.fits"))

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
cl3 = pn.Atom("Cl", 3)
Rlo = cl3.getLowDensRatio(wave1=5518, wave2=5538)
Rhi = cl3.getHighDensRatio(wave1=5518, wave2=5538)
Rlo, Rhi

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
m = imcont.data > 3e2
im5518.mask = im5518.mask | m
im5538.mask = im5538.mask | m
trim_edges(im5518, 20)
trim_edges(im5538, 20)
trim_edges(imcont, 20)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
shift5538 = 15.0
shift5518 = 23.0
im5538.data += shift5538
im5518.data += shift5518

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 16
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
im5538.rebin(n).plot(vmin=-10, vmax=120, ax=axes[0, 0], colorbar="v")
im5518.rebin(n).plot(vmin=-10, vmax=120, ax=axes[0, 1], colorbar="v")
imcont.rebin(n).plot(vmin=0, vmax=1e4, ax=axes[1, 0], colorbar="v")
(im5518.rebin(n) / im5538.rebin(n)).plot(
    ax=axes[1, 1],
    vmin=0.0,
    vmax=2.0,
    cmap="gray",
    colorbar="v",
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 8

imx = im5538.rebin(n)
imy = (im5518 / Rlo).rebin(n)

imin, imax = -10, 100
m = (imx.data < imax) & (imy.data < imax)
m = m & (imx.data > imin) & (imy.data > imin)
m = m & (imy.data / imx.data > -2) & (imy.data / imx.data < 4.0)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "5538": imx.data[m],
        "5518": imy.data[m],
        "5518 / 5538": imy.data[m] / imx.data[m],
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
g.axes[1, 0].plot([imin, imax], [imin, imax], "--", color="r")
g.axes[1, 0].plot([imin, imax], [imin, imax], "--", color="r")
g.fig.suptitle("Correlation between [Cl III] 5538 and 5518 brightness")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
n = 32

imx = imha.rebin(n)
imy = im5518.rebin(n) / im5538.rebin(n)

imin, imax = -10, 100
m = (imy.data > 0.5) & (imy.data < 2.0)
m = m & (imx.data < 1.5e5)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": imx.data[m],
        "5518 / 5538": imy.data[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=50,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=50,
    ),
)
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(1.5, color="r")
g.fig.suptitle("Correlation between [Cl III] 5538, 5518 sum and ratio")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
xslice, yslice = slice(230, 300), slice(144, 245)

n = 2

imx = imha[yslice, xslice].rebin(n)
imy = im5518[yslice, xslice].rebin(n) / im5538[yslice, xslice].rebin(n)

imin, imax = -10, 100
m = (imy.data > 0.15) & (imy.data < 15.0)
m = m & (imx.data < 1.5e5)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": np.log10(imx.data[m]),
        "5518 / 5538": np.log10(imy.data[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=30,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=30,
    ),
)
g.axes[1, 0].axhline(np.log10(1.5), color="r")
g.axes[1, 1].axvline(np.log10(1.5), color="r")


g.fig.suptitle("Correlation between [Cl III] 5538 / 5518 ratio and Ha brightness")

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
m = ~im5518[yslice, xslice].mask
npix = m.sum()
y = im5518[yslice, xslice].data[m]
x = im5538[yslice, xslice].data[m]
w = imha[yslice, xslice].data[m]

# unweighted
ym = y.mean()
xm = x.mean()
R1 = ym / xm
dR1 = np.sqrt((y.var() / ym ** 2 + x.var() / xm ** 2) / (npix - 1))

# weighted
ymw = np.average(y, weights=w)
xmw = np.average(x, weights=w)
R2 = ymw / xmw
dR2 = np.sqrt(
    (
        np.average(((y - ymw) / ymw) ** 2, weights=w)
        + np.average(((x - ymw) / xmw) ** 2, weights=w)
    )
    / (npix - 1)
)
f"Unweighted R = {R1:.4f} +/- {dR1:.4f}; Weighted R = {R2:.4f} +/- {dR2:.4f}"

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
cl3.getTemDen(1.44, tem=12000, wave1=5518, wave2=5538)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e5518 = cl3.getEmissivity(tem=12000, den=[1, 10, 100, 1000], wave=5518)
e5538 = cl3.getEmissivity(tem=12000, den=[1, 10, 100, 1000], wave=5538)
e5518 / e5538

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e5518 = cl3.getEmissivity(tem=15000, den=[1, 10, 100, 1000], wave=5518)
e5538 = cl3.getEmissivity(tem=15000, den=[1, 10, 100, 1000], wave=5538)
e5518 / e5538

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
e5518 = cl3.getEmissivity(tem=8000, den=[1, 10, 100, 1000], wave=5518)
e5538 = cl3.getEmissivity(tem=8000, den=[1, 10, 100, 1000], wave=5538)
e5518 / e5538

# + [markdown] pycharm={"name": "#%% md\n"} jupyter={"outputs_hidden": false}
# Check the 2 sigma lower limit

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
cl3.getTemDen(R2 - 5 * dR2, tem=15000, wave1=5518, wave2=5538)

# + pycharm={"name": "#%%\n"} jupyter={"outputs_hidden": false}
cl3.getTemDen(R2 - 5 * dR2, tem=1000, wave1=5518, wave2=5538)
# -


