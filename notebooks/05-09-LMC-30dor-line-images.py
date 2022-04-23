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

# # Line images and ratios from 30 Dor MUSE data

# +
from pathlib import Path
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
import cmasher as cm
from mpdaf.obj import Image
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd

sns.set_context("talk")
sns.set_color_codes()
# -

DATADIR = Path("../data")
def p(s):
    """Return path to image file"""
    return str(DATADIR / f"lmc-30dor-ABCD-{s}-bin01-sum.fits")


imcont = Image("../data/lmc-30dor-ABCD-46-55-avcont.fits")

# ## The misterious molecular lines

# +
im9114 = Image(p("xxx-9114"))
im9112c = Image(p("CONT-9112"))

im9114.data -= im9112c.data
im9114.plot(vmin=-10, vmax=50, cmap="gray_r")

# +
im8152 = Image(p("xxx-8152"))

im8152.plot(vmin=-10, vmax=50, cmap="gray_r")
# -

# And this is the other unidentified line, which is a bit more diffuse. 

# +
im9145 = Image(p("xxx-9145"))
im9145c = Image(p("CONT-9145"))

im9145.data -= im9145c.data
im9145.plot(vmin=-30, vmax=50, cmap="gray_r")
# -

# Although now that I look at it, it really isn't *so* different from the other two. 

# The forbidden [C I] line.

im8727 = Image(p("ci-8727"))
im8727.plot(vmin=-10, vmax=50, cmap="gray_r")

# And improve the map by subtracting off the hyperlocal continuum. 

# +
im8730c = Image(p("CONT-8730"))
im8727.data -= im8730c.data

im8727.plot(vmin=-150, vmax=150, cmap="gray_r")
# -
# Look at the fluorescent O I emission for comparison.  Also plot the continuum contours on top. 

# + tags=[]
im8446 = Image(p("oi-8446"))
fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(im8446.data, origin="lower", vmin=0, vmax=0.3e4, cmap="inferno")
ax.contour(imcont.data, levels=[5e3, 1e4, 1e5], colors="c")
...;
# -

# Now plot a color image including the 3 lines: 
# * C I 8727
# * Si I 8152 - pretty certain of ID
# * C I 9114 - unsure of ID
#

fig, ax = plt.subplots(figsize=(10, 10))
rgb = np.stack(
    [
        im8727.data / 100,
        (im8152.data + 10) / 40,
        (im9114.data + 10) / 40,
    ],
    axis=-1,
)
ax.imshow(rgb, origin="lower")
ax.contour(
    im8446.data,
    levels=[0.5e3, 1e3, 2e3, 4e3, 6e3],
    linewidths=[0.4, 0.8, 1.2, 1.6, 2.0],
    colors="y",
)
...;

# There is an anticorrelation between mystery lines (green/blue) and both the O I (contours) and the [C I] (red). The only exception is "Cloud 2" (magenta colored in above image), which shows both 8727 and 9114, but not 8152!
#
# Also the [C I] 8727 looks like it only comes from clouds that have a reasonable amount of CO. 

# Now look at the Ca I 9095 line, which Jorge García Rojas had down as a C I, which cannot be true.

im9095 = Image(p("cai-9095"))
im9095.data -= im9112c.data
starmask = imcont.data > 5e3
im9095.data[starmask] = 0.0

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(im9095.data, origin="lower", vmin=-30, vmax=150, cmap="inferno")
ax.contour(
    im8446.data,
    levels=[0.5e3, 1e3, 2e3, 4e3, 6e3],
    linewidths=[0.4, 0.8, 1.2, 1.6, 2.0],
    colors="b",
)

# +
im9204 = Image(p("cai-9204"))
im9219 = Image(p("cai-9219"))
im9244 = Image(p("cai-9244"))

starmask = imcont.data > 5e3
for im in im9204, im9219, im9244:
    im.data[starmask] = 0.0
# -

fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)
ims = im9095, im9204, im9219, im9244
vmaxes = 150, 200, 400, 100
for im, ax, vmax in zip(ims, axes.flat, vmaxes):
    ax.imshow(im.data, origin="lower", vmin=-vmax/5, vmax=vmax, cmap="inferno")
    ax.set_title(im.filename.split("-")[4])
fig.suptitle("Supposed Ca I lines")
fig.tight_layout()
...;

# So, of these four, 9095 and 9244 look kind of similar, and more or less like I would expect from a Ca I line.  Whereas 9219 looks more like an ionized gas tracer and 9204 too but with some added knots.
#
# I had already noted in notebook 05-04 that 9219 is likely not Ca I.
#

# +
im6104 = Image(p("ca-i-6104"))
im7890 = Image(p("cai-7890"))
im8125 = Image(p("cai-8125"))

starmask = imcont.data > 5e3
for im in im6104, im7890, im8125:
    im.data[starmask] = 0.0
# -

fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)
ims = im6104, im7890, im8125, im9244
iis = 5, 4, 4, 4
vmaxes = 80, 50, 200, 100
for im, ax, vmax, ii in zip(ims, axes.flat, vmaxes, iis):
    ax.imshow(im.data, origin="lower", vmin=-vmax/2, vmax=vmax, cmap="magma")
    ax.set_title(im.filename.split("-")[ii])
fig.suptitle("More supposed Ca I lines")
fig.tight_layout()
...;

# So the 6104 line is clearly not Ca I since its distribution looks like a higly ionized line. It has highest value around the WR star in the north.
#
# The 7890 line might possibly be Ca I but its distribution is very different form the other lines. It looks more like Fe III maybe. *Now identified as [Ni III] nebular line* – see `50-01-astroquery-lines` notebook.
#
# The 8125 line does look exactly like 9244 and 9095, except that there seems to be a problem with the sky subtraction in panels ABC.







# ## Remap Brackett gamma IR line
#
# This comes from Sherry Weh


BIGDATAPATH = Path.cwd().parent / "big-data" / "30-Dor-Radio"

hdu_brg = fits.open(BIGDATAPATH / "30DOR_BRG_FL_SCCY.fits")[0]

# Load a random Paschen line to get the MUSE wcs

hdu_pa = fits.open(p("hi-9229"))[0]

from reproject import reproject_interp

# Reproject the near-IR map to the MUSE grid

brg_ABCD = reproject_interp(hdu_brg, hdu_pa.header, return_footprint=False)

# And save to a file

fits.PrimaryHDU(
    data=brg_ABCD, 
    header=hdu_pa.header,
).writeto(
    DATADIR / "lmc-30dor-ABCD-hi-21661-reproject.fits",
    overwrite=True,
)

hdu_brg.header["BUNIT"]

# Is that really the brightness units of the image?

# ## H I line ratios for reddening

hifiles = DATADIR.glob(f"lmc-30dor-ABCD-hi-*-bin01-sum.fits")
sorted(hifiles)

im4861 = Image(p("hi-4861"))
im8467 = Image(p("hi-8467"))
im8545 = Image(p("hi-8545"))
im8598 = Image(p("hi-8598"))
im8665 = Image(p("hi-8665"))
im8750 = Image(p("hi-8750"))
im9229 = Image(p("hi-9229"))

im21661 = 1e20 * Image(str(DATADIR / "lmc-30dor-ABCD-hi-21661-reproject.fits"))

im6563 = Image("../data/MUSE_R136toWill/GAUS_Ha6562.8_060_Will.fits")
saveunit = im6563.unit
sig6563 = Image("../data/MUSE_R136toWill/GAUS_Ha6562.8_060_Will.fits", ext=3)
im6563.data *= sig6563.data / 1.25
im6563.unit = saveunit

im6563.unit

im7751 = Image(p("ariii-7751"))
im7136 = Image(p("ariii-7136"))
(im7751 / im7136).data.mean()

np.nanmedian(im21661.data), np.nanmedian(im9229.data.data)

from astropy.convolution import convolve, Gaussian2DKernel
kernel = Gaussian2DKernel(x_stddev=3.0)

# + tags=[]
fig, axes = plt.subplots(3, 2, figsize=(12, 15), sharex=True, sharey=True)
(im9229 / im4861).plot(
    vmin=0.025, vmax=0.25, 
    cmap=cm.ghostlight_r,
    ax=axes[0, 0],
    colorbar="v",
)
axes[0, 0].set_title("9229 / 4861")
(im6563 / im4861).plot(
    vmin=1.2, vmax=6.0, 
    cmap=cm.ghostlight_r,
    ax=axes[0, 1],
    colorbar="v",
)
axes[0, 1].set_title("6563 / 4861")
(im9229 / im6563).plot(
    vmin=0.022, vmax=0.066, 
    cmap=cm.ghostlight_r,
    ax=axes[1, 0],
    colorbar="v",
)
axes[1, 0].set_title("9229 / 6563")
(im7751 / im7136).plot(
    vmin=0.25, vmax=0.4, 
    cmap=cm.amethyst_r,
    ax=axes[1, 1],
    colorbar="v",
)

rsmooth = convolve(im7751.data, kernel) / convolve(im7136.data, kernel)
axes[1, 1].contour(rsmooth, levels=[0.266], colors="r")
axes[1, 1].set_title("7751 / 7136")

(im21661 / im9229).plot(
    vmin=0.035, vmax=0.15, 
    cmap=cm.ghostlight_r,
    ax=axes[2, 0],
    colorbar="v",
)
axes[2, 0].set_title("21661 / 9229")

np.sqrt(im9229).plot(
    #vmin=2.5, vmax=4.0, 
    vmin=0.0,
    vmax=2e2,
    cmap=cm.cosmic,
    ax=axes[2, 1],
    colorbar="v",
)
axes[2, 1].set_title("I(9229)")

fig.tight_layout()
...;

# + [markdown] tags=[]
# In general, there is a good correlation between all the reddening indicators
#
#
# There are some cases, however, where the blue–red reddening and the red–infared are different.  For instance, the cloud at (x, y) = (400, 10) has a high 6563/4861 and high 7751/7136, but is not particularly visible in 9229 / 6563
#
# But note that this same cloud is prominent in Br-gamma / 9229, which is odd.
#
# Also, we need to improve the sky subtraction for the [Ar III] lines, since there is a clear variation between the 4 fields.
#
# And another issue is the zero point in the Br-gamma map.  Presumably, that is what gives the dark Brg/9229 regions in the low brightness region to the West (where all other ratios suggest low extinction).

# +
n = 4
r_hi_hb = im9229.rebin(n) / im4861.rebin(n)
r_ha_hb = im6563.rebin(n) / im4861.rebin(n)
r_hi_ha = im9229.rebin(n) / im6563.rebin(n)
r_ariii = im7751.rebin(n) / im7136.rebin(n)
r_bg_hi = im21661.rebin(n) / im9229.rebin(n)
i_ha = im6563.rebin(n)
m = ~r_hi_ha.data.mask & (imcont.rebin(n).data < 1e3)
m = m & (r_hi_ha.data > 0.02) & (r_hi_ha.data < 0.1)
m = m & (r_ha_hb.data > 1.3) & (r_ha_hb.data < 7.0)
m = m & (r_ariii.data > 0.24) & (r_ariii.data < 0.4)
m = m & (r_bg_hi.data > 0.025) & (r_bg_hi.data < 0.25)
m
df = pd.DataFrame(
    {
        "I(Ha)": np.log10(i_ha.data[m]),
        "Ha / Hb": np.log10(r_ha_hb.data[m]),
        "H9229 / Ha": np.log10(r_hi_ha.data[m]),
        "Br g / H9229": np.log10(r_bg_hi.data[m]),
        "7751 / 7136": np.log10(r_ariii.data[m]),
    }
)

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
#       weights=im21661.rebin(n).data[m],
        weights=i_ha.data[m],
#        weights=r_hi_hb.data[m] - 0.025,
        bins=100,
    ),
)
...;

# +
n = 8
r_hi_hb = im9229.rebin(n) / im4861.rebin(n)
r_hi_ha = im9229.rebin(n) / im6563.rebin(n)
r_ha_hb = im6563.rebin(n) / im4861.rebin(n)
r_ariii = im7751.rebin(n) / im7136.rebin(n)
r_bg_ha = im21661.rebin(n) / im6563.rebin(n)
r_bg_hi = im21661.rebin(n) / im9229.rebin(n)
i_ha = im6563.rebin(n)
i_bg = im21661.rebin(n)
m = ~r_hi_hb.data.mask & (imcont.rebin(n).data < 1e3)
m = m & (r_hi_hb.data > 0.026) & (r_hi_hb.data < 0.7)
m = m & (r_ha_hb.data > 1.3) & (r_ha_hb.data < 7.0)
m = m & (r_ariii.data > 0.24) & (r_ariii.data < 0.4)
m = m & (r_bg_ha.data > 0.0003) & (r_bg_ha.data < 0.01)
m
df = pd.DataFrame(
    {
        "I(Bg)": np.log10(i_bg.data[m]),
        "Ha / Hb": np.log10(r_ha_hb.data[m]),
        "H9229 / Ha": np.log10(r_hi_ha.data[m]),
        "Br g / H9229": np.log10(r_bg_hi.data[m]),
        "7751 / 7136": np.log10(r_ariii.data[m]),
    }
)

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
#       weights=im21661.rebin(n).data[m],
        weights=i_ha.data[m],
#        weights=r_hi_hb.data[m] - 0.025,
        bins=50,
    ),
)
...;
# -

# These show very high correlations between the reddening indicators.  There is a subtle curvature in some of them though.
#
# In particular, it looks like the infrared /optical ratios increase slower vs red/blue at higher extinctions.  *What could that mean?* *Does it make sense in terms of higher R_V?*
#
# Also, maybe the underlying photospheric lines might be affecting things.
#
# There seems to be a change in behavior around 0.3 in log10(Ha/Hb). Note that I am not convinced by the normalization since the Ha/Hb ratio seems too small in the low extinction parts.

# ## The O++ lines

im4959 = Image(p("oiii-4959"))
im4931 = Image(p("oiii-4931"))
im4642 = Image(p("oii-4642"))
im4650 = Image(p("oii-4650"))

r_orl_cel = ((im4650 + im4642)/ im4959)
r_orl_cel.mask = r_orl_cel.mask | (im4959.data < 3e4) | (imcont.data > 1e3)
fig, ax = plt.subplots(figsize=(10, 10))
r_orl_cel.plot(vmin=0, vmax=0.003, cmap=cm.arctic_r, colorbar="v")
ax.contour(im4959.data, levels=[1e5, 2e5, 4e5], linewidths=[1.0, 2.0, 3.0], colors="k")

# This shows that the collisional lines are relatively stronger where the intensity is highest.  And that teh permitted lines are relatively higher in the inner parts. 

n = 4
im1, im2 = im4642.copy(), im4650.copy()
im1.mask = im1.mask | (im4959.data < 3e4) | (imcont.data > 1e3)
im2.mask = im1.mask
r = (im1.rebin(n) / im2.rebin(n))
fig, ax = plt.subplots(figsize=(10, 10))
r.plot(vmin=0.6, vmax=1.4, cmap=cm.fusion, colorbar="v")

im4740 = Image(p("ariv-4740"))
im7330 = Image(p("oii-7330"))

r_ariv_oiii = (im4740 / im4959)
fig, ax = plt.subplots(figsize=(10, 10))
r_ariv_oiii.plot(vmin=0.0, vmax=0.004, cmap=cm.neutral_r, colorbar="v")
ax.contour(im4959.data, levels=[1e5, 2e5, 4e5], linewidths=[1.0, 2.0, 3.0], colors="r")

r_oiii_oii = (im4959 / im7330)
fig, ax = plt.subplots(figsize=(10, 10))
r_oiii_oii.plot(vmin=0.0, vmax=300.0, cmap=cm.neutral_r, colorbar="v")
ax.contour(im4959.data, levels=[1e5, 2e5, 4e5], linewidths=[1.0, 2.0, 3.0], colors="r")

# +
n = 4
m = ~r_orl_cel.rebin(n).data.mask
m = m & (r_oiii_oii.rebin(n).data > 3.0) & (r_oiii_oii.rebin(n).data < 600.0)
m = m & (r_ariv_oiii.rebin(n).data > 3e-4) & (r_ariv_oiii.rebin(n).data < 0.007)
m = m & (r_orl_cel.rebin(n).data > 6e-4) & (r_orl_cel.rebin(n).data < 0.007)
df = pd.DataFrame(
    {
        "R(V1/4959)": np.log10(r_orl_cel.rebin(n).data[m]),
        "R([O III] / [O II])": np.log10(r_oiii_oii.rebin(n).data[m]),
        "R([Ar IV] / [O III])": np.log10(r_ariv_oiii.rebin(n).data[m]),
    }
)

g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=im4740.rebin(n).data[m],
        bins=50,
    ),
);

# -

# ## The C++ lines
#
# This would be interesting to compare with my Orion project on the same

im6462 = Image(p("cii-6462"))
im7231 = Image(p("cii-7231"))
im7236 = Image(p("cii-7236"))

(im7236 + im7231 - 300).plot(vmin=0, vmax=1500)

im5820 = Image(p("ni-iv-5820"))

# + tags=[]
im5820.plot(vmin=0, vmax=300)
# -

im8046 = Image(p("cliv-8046"))

# + jupyter={"source_hidden": true} tags=[]
im8046.plot(vmin=0, vmax=1000)
# -

im4740.plot(vmin=0, vmax=1500)

r_cii_oiii = ((im7231 + im7236 - 300) / im6563)
fig, ax = plt.subplots(figsize=(10, 10))
r_cii_oiii.plot(vmin=0.0, vmax=0.004, cmap=cm.arctic_r, colorbar="v")
ax.contour(im4959.data, levels=[1e5, 2e5, 4e5], linewidths=[1.0, 2.0, 3.0], colors="r")


# So there is slight evidence for an excess of C II over Ha on the inside rim, but we would have to do a much better job of the sky subtraction to be sure. 

# ## Look at the Raman wings

# Where are the files of these?

def rb(s):
    """Return path to raman band file"""
    return str(DATADIR / f"30D_ABCD-raman-{s}.fits")


imR1 = Image(rb("R040"))
imR2 = Image(rb("R058"))
imR3 = Image(rb("R087"))
imR4 = Image(rb("R136"))

imB1 = Image(rb("B033"))
imB2 = Image(rb("B054"))
imB3 = Image(rb("B080"))
imB4 = Image(rb("B133"))

fig, [axR, axB] = plt.subplots(1, 2, figsize=(12, 8))
(imR1 + imR2 + imR3 + imR4).plot(ax=axR, vmin=0, vmax=2000, cmap=cm.neutral)
(imB1 + imB2 + imB3 + imB4).plot(ax=axB, vmin=0, vmax=2000, cmap=cm.neutral)
...;

# +
fig, [axR, axB] = plt.subplots(1, 2, figsize=(14, 8))
rgbR =np.stack(
    [
        imR1.data / 600,
        imR2.data / 600,
        (imR3).data / 200,
    ],
    axis=-1,
)
axR.imshow(rgbR[:, :, ::-1], origin="lower")
axR.set_title("Red wing")


rgbB =np.stack(
    [
        imB1.data / 1000,
        imB2.data / 500,
        (imB3).data / 200,
    ],
    axis=-1,
)
axB.imshow(rgbB[:, :, ::-1], origin="lower")
axB.set_title("Blue wing")
fig.suptitle("Individual Raman wings")
fig.tight_layout()

...;

# +
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
rgb =np.stack(
    [
        (imR1 + imB1).data / 1800,
        (imR2 + imB2).data / 1000,
        (imR3 + imB3).data / 500,
    ],
    axis=-1,
)
starmask = imcont.data > 2e3
rgb *= (~starmask).astype(float)[:, :, None]

ax.imshow(rgb[:, :, ::-1], origin="lower")
oi = im8446.data.copy()
oi[starmask] = np.nan

ax.contour(oi, levels=[1e3, 2e3, 4e3], linewidths=np.arange(1, 4), colors=["y", "g", "r"])
fig.suptitle("Combined Raman wings with O I contours and masked out stars")
fig.tight_layout()
...;

# -

# So in general the Raman and O I are very well correlated.  Better than for anyy other line I would imagine, although we need to check that. 
#
# The brightest Raman emission coincides with the brightest O I 8446 in MC 10 in quadrant A.
#
# However, in quadrant D the brightest O I is in a globule whereas the brightest Raman is at the Bar feature. 
#
# Maybe the globule has some internal excitation, since it does not have any CO emission associated with it. 


