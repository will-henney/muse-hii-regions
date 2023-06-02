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

# # Extract [O I] and O I lines

# ## Library imports and data loading
#
# Copied from the `02-01-Raman-Wings-IR-Source-C` notebook:

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import regions
import sys
sys.path.append("../lib")
import moments
import extract

sns.set_context("talk")
sns.set_color_codes()
# -

moments.FIGPATH = Path("../figs")
moments.SAVEPATH = Path("../data")

# *Update 2021-05-4:* Do the co-added cube as well

datapath = Path("/Users/will/Work/Muse-Hii-Data/SMC-NGC-346/")
file1 = "ADP.2017-10-16T11_04_19.247.fits"
file2 = "ADP.2016-10-12T05_43_23.882.fits"
cube1 = Cube(str(datapath / file1))
cube2 = Cube(str(datapath / file2))

# ## Use moments to extract the [O I] 6300 line

# +
b1, b2 = 6265, 6275
r1, r2 = 6332, 6345
w1, w2 = 6298, 6308

oi6300cube = cube2.select_lambda(b1-200, r1+200)
moi6300cube = cube1.select_lambda(b1-200, r1+200)

# +
# Position of mYSO
x, y = 116, 147
# position of B0 star SSN 43
xb, yb = 95, 163
# position of B0 star SSN 40
xbb, ybb = 121, 133
m = 3

boxes = {
    "SSN 152+168": regions.BoundingBox(
        iymin=y-m, iymax=y+m, ixmin=x-m, ixmax=x+m,
    ),
#    "BG 152+168": regions.BoundingBox(
#        iymin=y-2*m, iymax=y+2*m, ixmin=x-2*m, ixmax=x+2*m,
#    ),
    "SSN 43": regions.BoundingBox(
        iymin=yb-m, iymax=yb+m, ixmin=xb-m, ixmax=xb+m,
    ),
#    "BG 43": regions.BoundingBox(
#        iymin=yb-2*m, iymax=yb+2*m, ixmin=xb-2*m, ixmax=xb+2*m,
#    ),
    "SSN 40": regions.BoundingBox(
        iymin=ybb-m, iymax=ybb+m, ixmin=xbb-m, ixmax=xbb+m,
    ),
#    "BG 40": regions.BoundingBox(
#        iymin=ybb-2*m, iymax=ybb+2*m, ixmin=xbb-2*m, ixmax=xbb+2*m,
#    ),
    "blue": regions.BoundingBox(
        iymin=75, iymax=140, ixmin=15, ixmax=40
    ),
    "red": regions.BoundingBox(
        iymin=200, iymax=250, ixmin=210, ixmax=300
    ),
    "magenta": regions.BoundingBox(
        iymin=10, iymax=50, ixmin=100, ixmax=150
    ),
    "green": regions.BoundingBox(
        iymin=10, iymax=100, ixmin=200, ixmax=300
    ),
    "cyan": regions.BoundingBox(
        iymin=170, iymax=210, ixmin=90, ixmax=120
    ),
}

# -

wavranges = [
    (6000, 6100),
    (6150, 6195), (6205, 6250), (6330, 6345), (6400, 6450),
    (6480, 6490), (6520, 6530),
]

# +
fig, ax = plt.subplots(figsize=(12, 12))
for box in boxes.values():
    yslice, xslice = box.slices
    spec = oi6300cube[:, yslice, xslice].mean(axis=(1, 2))
    spec.plot()
    
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
sns.despine();
# -

yslice, xslice = boxes["SSN 152+168"].slices
subcube = oi6300cube[:, yslice, xslice]
contcube = extract.fit_continuum(
    subcube, wav_ranges=wavranges, deg=5, median=False,
)

subcube.mean(axis=(1, 2)).plot()
contcube.mean(axis=(1, 2)).plot()

contcube = extract.fit_continuum(
    oi6300cube, wav_ranges=wavranges, deg=5, median=False,
)

mcontcube = extract.fit_continuum(
    moi6300cube, wav_ranges=wavranges, deg=5, median=False,
)

# +
fig, ax = plt.subplots(figsize=(12, 12))
for box in boxes.values():
    yslice, xslice = box.slices
    spec = oi6300cube[:, yslice, xslice].mean(axis=(1, 2))
    cspec = contcube[:, yslice, xslice].mean(axis=(1, 2))
    spec.plot()
    cspec.plot(color="k")
    
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
sns.despine();
# -

# These are defined with respect to mcube
mboxes = {
#    "sw filament": regions.BoundingBox(
#        iymin=30, iymax=50, ixmin=300, ixmax=330,
#    ),
    "bow shock": regions.BoundingBox(
        iymin=165, iymax=205, ixmin=240, ixmax=290,
    ),
    "w filament": regions.BoundingBox(
        iymin=100, iymax=130, ixmin=25, ixmax=55,
    ),
    "c filament": regions.BoundingBox(
        iymin=195, iymax=210, ixmin=155, ixmax=195,
    ),
    "blue": regions.BoundingBox(
        iymin=75, iymax=140, ixmin=15, ixmax=40
    ),
    "red": regions.BoundingBox(
        iymin=200, iymax=250, ixmin=210, ixmax=300
    ),
    "magenta": regions.BoundingBox(
        iymin=10, iymax=50, ixmin=100, ixmax=150
    ),
    "green": regions.BoundingBox(
        iymin=10, iymax=100, ixmin=200, ixmax=300
    ),
    "cyan": regions.BoundingBox(
        iymin=170, iymax=210, ixmin=90, ixmax=120
    ),
}

# +
fig, ax = plt.subplots(figsize=(12, 12))
for box in mboxes.values():
    yslice, xslice = box.slices
    spec = moi6300cube[:, yslice, xslice].mean(axis=(1, 2))
    cspec = mcontcube[:, yslice, xslice].mean(axis=(1, 2))
    spec.plot()
    cspec.plot(color="k")
    
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.set(
#    ylim=[0, 750],
)
sns.despine();
# -

wavmin, wavmax = wavranges[0][0], wavranges[-1][1]
wavmin, wavmax

csubcube = (oi6300cube - contcube)
cdivcube = (oi6300cube / contcube)
csubcube.write(
    f"{prefix}-contsub.fits",
    savemask="nan",
    )
cdivcube.write(
    f"{prefix}-contdiv.fits",
    savemask="nan",
    )
contcube.write(
    f"{prefix}-cont.fits",
    savemask="nan",
    )

prefix = f"../big-data/ngc346-{wavmin:d}-{wavmax:d}-cube"
mcsubcube = (moi6300cube - mcontcube)
mcdivcube = (moi6300cube / mcontcube)
mcsubcube.write(
    f"{prefix}-contsub.fits",
    savemask="nan",
    )
mcdivcube.write(
    f"{prefix}-contdiv.fits",
    savemask="nan",
    )
mcontcube.write(
    f"{prefix}-cont.fits",
    savemask="nan",
    )

for label, box in mboxes.items():
    yslice, xslice = box.slices
    spec = csubcube[:, yslice, xslice].select_lambda(6280, 6330).mean(axis=(1, 2))
    try:
        spec.plot(color=label[0])
    except:
        pass

w1, w2

mom6300 = moments.find_moments(
    csubcube.select_lambda(6300, 6308)
)
mom6300m = moments.find_moments(
    mcsubcube.select_lambda(6300, 6308)
)

mom_pars_6300 = dict(
    restwav=6300.30,
    irange=[-150, 4.0e4],
    vrange=[85, 255],
    srange=[10, 120],    
)

plot_pars_6300=dict(
    ilabel="O I",
    label="6300",
    flabel="ngc346-oi",
    **mom_pars_6300,
)
g = moments.moments_corner_plot(
    mom6300, rebin=1, **plot_pars_6300
)

plot_pars_6300=dict(
    ilabel="O I",
    label="6300",
    flabel="ngc346-oi",
    **mom_pars_6300,
)
g = moments.moments_corner_plot(
    mom6300, rebin=4, **plot_pars_6300
)

g = moments.moments_corner_plot(
    mom6300m, rebin=4, **plot_pars_6300
)

# It looks like the possible symptoms of oversubtracted sky.

# +
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

(3e5*(mom6300[1] / 6300.30 - 1.0)).rebin(1).plot(
    vmin=100, vmax=220, 
    cmap="seismic", 
    ax=axes[0],
)
(3e5*(mom6300m[1] / 6300.30 - 1.0)).rebin(1).plot(
    vmin=100, vmax=220, 
    cmap="seismic", 
    ax=axes[1],
);

# +
fig, axes = plt.subplots(
    2, 2, 
    figsize=(10, 10),
    sharex=True, sharey=True,
)

imap = mom6300[0].copy()
vmap = 3e5*(mom6300[1] / 6300.30 - 1.0)
smap = 3e5*(mom6300[2] / 6300.30)

m = imap.data > 10.

vmap.mask = vmap.mask | (~m)
smap.mask = smap.mask | (~m)

imap.rebin(1).plot(
    vmin=-500, vmax=3e4, 
    cmap="turbo", 
    ax=axes[0, 0],
)

vmap.rebin(1).plot(
    vmin=100, vmax=220, 
    cmap="seismic", 
    ax=axes[0, 1],
)

smap.rebin(1).plot(
    vmin=0, vmax=120, 
    cmap="magma", 
    ax=axes[1, 0],
)

imap.rebin(1).plot(
    vmin=-350, vmax=-200, 
    cmap="viridis", 
    ax=axes[1, 1],
)
bg_6300 = contcube.select_lambda(w1, w2).mean(axis=0)
axes[1, 1].contour(
    bg_6300.data,
    levels=[300],
    colors="r",
)

fig.tight_layout();
# -

msky = (imap.data < -100) & (imap.data > -200) & (bg_6300.data < 300)
msky[310:, :] = False
msky[:, 310:] = False
msky[:10, :] = False
msky[:, :10] = False

msky.sum(), np.where(msky)

core_6300 = csubcube.select_lambda(6295, 6310)
sky_6300 = core_6300.copy() 
sky_6300.mask = sky_6300.mask | ~msky[None, : :]

sky_6300.mean(axis=(1, 2)).plot()

corr_6300 = core_6300 - sky_6300.mean(axis=(1, 2))

testpixels = [
    [250, 160], [150, 150], [160, 220],
    [70, 250], [75, 200], [310, 225],
    [25, 140], [250, 250], [140, 110], #[180, 290],
]
fig, axes = plt.subplots(
    3, 3, 
    figsize=(10, 8), 
    sharex=True,
    sharey="row",
)
for (j, i), ax in zip(testpixels, axes.flat):
    core_6300[:, j, i].plot(ax=ax)
    corr_6300[:, j, i].plot(ax=ax) 
    ax.set(xlabel="", ylabel="")
    ax.set_title(f"[{j}, {i}]")
fig.suptitle(
    "Before/after sky correction for faint/moderate/bright pixels"
)
sns.despine()
fig.tight_layout();

mom6300c = moments.find_moments(corr_6300.select_lambda(6300, 6308))

g = moments.moments_corner_plot(
    mom6300c, rebin=1, **plot_pars_6300
)

g = moments.moments_corner_plot(
    mom6300c, rebin=4, **plot_pars_6300
)

# +
fig, axes = plt.subplots(
    2, 2, 
    figsize=(10, 10),
    sharex=True, sharey=True,
)

imap = mom6300c[0].copy()
vmap = 3e5*(mom6300c[1] / 6300.30 - 1.0)
smap = 3e5*(mom6300c[2] / 6300.30)

m = imap.data > 0.

vmap.mask = vmap.mask | (~m)
smap.mask = smap.mask | (~m)

vmap_old = 3e5*(mom6300[1] / 6300.30 - 1.0)
vmap_old.mask = vmap_old.mask | (~m)


imap.rebin(1).plot(
    vmin=0, vmax=1e4, 
    cmap="turbo", 
    ax=axes[0, 0],
)

vmap.rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    ax=axes[0, 1],
)

smap.rebin(1).plot(
    vmin=40, vmax=120, 
    cmap="magma", 
    ax=axes[1, 0],
)

vmap_old.rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    ax=axes[1, 1],
)


fig.tight_layout();
# -

vmap.rebin(4).plot(
    vmin=130, vmax=180, 
    cmap="seismic", 
);

# There is clearly a problem still with the fainter pixels.  Witness the trend of sig with intensity.  But theis is not so important for studying the brighter parts.

mom_pars_6300 = dict(
    restwav=6300.30,
    irange=[-150, 4.0e4],
    vrange=[85, 255],
    srange=[30, 150],    
)
moments.save_moments_to_fits(
    mom6300c,
    label="6300",
    flabel="ngc346-sharp-oi",
    **mom_pars_6300,
)

# Now repeat for the co-added cube

# +
fig, axes = plt.subplots(
    2, 2, 
    figsize=(10, 10),
    sharex=True, sharey=True,
)

imap = mom6300m[0].copy()
vmap = 3e5*(mom6300m[1] / 6300.30 - 1.0)
smap = 3e5*(mom6300m[2] / 6300.30)

m = imap.data > 10.

vmap.mask = vmap.mask | (~m)
smap.mask = smap.mask | (~m)

imap.rebin(1).plot(
    vmin=-500, vmax=3e4, 
    cmap="turbo", 
    ax=axes[0, 0],
)

vmap.rebin(1).plot(
    vmin=100, vmax=220, 
    cmap="seismic", 
    ax=axes[0, 1],
)

smap.rebin(1).plot(
    vmin=0, vmax=120, 
    cmap="magma", 
    ax=axes[1, 0],
)

imap.rebin(1).plot(
    vmin=-200, vmax=-0, 
    cmap="viridis", 
    ax=axes[1, 1],
)
mbg_6300 = mcontcube.select_lambda(w1, w2).mean(axis=0)
axes[1, 1].contour(
    mbg_6300.data,
    levels=[300],
    colors="r",
)

fig.tight_layout();
# -

msky = (imap.data < -150) & (imap.data > -200) & (mbg_6300.data < 300)
msky[310:, :] = False
msky[:, 310:] = False
msky[:10, :] = False
msky[:, :10] = False

msky.sum(), np.where(msky)

core_6300 = mcsubcube.select_lambda(6295, 6310)
sky_6300 = core_6300.copy() 
sky_6300.mask = sky_6300.mask | ~msky[None, : :]

sky_6300.mean(axis=(1, 2)).plot()

corr_6300 = core_6300 - sky_6300.mean(axis=(1, 2))

testpixels = [
    [250, 160], [150, 150], [160, 220],
    [70, 250], [75, 200], [310, 225],
    [25, 140], [250, 250], [140, 110], #[180, 290],
]
fig, axes = plt.subplots(
    3, 3, 
    figsize=(10, 8), 
    sharex=True,
    sharey="row",
)
for (j, i), ax in zip(testpixels, axes.flat):
    core_6300[:, j, i].plot(ax=ax)
    corr_6300[:, j, i].plot(ax=ax) 
    ax.set(xlabel="", ylabel="")
    ax.set_title(f"[{j}, {i}]")
fig.suptitle(
    "Before/after sky correction for faint/moderate/bright pixels"
)
sns.despine()
fig.tight_layout();

mom6300mc = moments.find_moments(corr_6300.select_lambda(6300, 6308))

g = moments.moments_corner_plot(
    mom6300mc, rebin=1, **plot_pars_6300
)

g = moments.moments_corner_plot(
    mom6300mc, rebin=4, **plot_pars_6300
)

# +
fig, axes = plt.subplots(
    2, 2, 
    figsize=(10, 10),
    sharex=True, sharey=True,
)

imap = mom6300mc[0].copy()
vmap = 3e5*(mom6300mc[1] / 6300.30 - 1.0)
smap = 3e5*(mom6300mc[2] / 6300.30)

m = imap.data > 0.

vmap.mask = vmap.mask | (~m)
smap.mask = smap.mask | (~m)

vmap_old = 3e5*(mom6300m[1] / 6300.30 - 1.0)
vmap_old.mask = vmap_old.mask | (~m)


imap.rebin(1).plot(
    vmin=0, vmax=1e4, 
    cmap="turbo", 
    ax=axes[0, 0],
)

vmap.rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    ax=axes[0, 1],
)

smap.rebin(1).plot(
    vmin=40, vmax=120, 
    cmap="magma", 
    ax=axes[1, 0],
)

vmap_old.rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    ax=axes[1, 1],
)


fig.tight_layout();
# -

vmap.rebin(4).plot(
    vmin=130, vmax=180, 
    cmap="seismic", 
);

# This is way better than the sharp version

moments.save_moments_to_fits(
    mom6300mc,
    label="6300",
    flabel="ngc346-oi",
    restwav=6300.30,
)

# ### The other [O I] line: 6363 
#
# This time around, we will just work with the co-added cube

core_6363 = mcsubcube.select_lambda(6362, 6372)
mom6363 = moments.find_moments(core_6363)

mom_pars_6363 = dict(
    restwav=6363.78,
    irange=[-50, 1.3e4],
    vrange=[85, 255],
    srange=[10, 120],    
)

plot_pars_6363=dict(
    ilabel="O I",
    label="6363",
    flabel="ngc346-oi",
    **mom_pars_6363,
)
g = moments.moments_corner_plot(
    mom6363, rebin=1, **plot_pars_6363,
)

g = moments.moments_corner_plot(
    mom6363, rebin=4, **plot_pars_6363,
)

(3e5*(mom6363[1] / 6363.78 - 1.0)).rebin(1).plot(
    vmin=100, vmax=220, 
    cmap="seismic", 
    colorbar="v",
)

# +
fig, axes = plt.subplots(
    2, 2, 
    figsize=(10, 10),
    sharex=True, sharey=True,
)

imap = mom6363[0].copy()
vmap = 3e5*(mom6363[1] / 6363.78 - 1.0)
smap = 3e5*(mom6363[2] / 6363.78)

m = imap.data > 10.

vmap.mask = vmap.mask | (~m)
smap.mask = smap.mask | (~m)

imap.rebin(1).plot(
    vmin=-500/3, vmax=3e4/3, 
    cmap="turbo", 
    ax=axes[0, 0],
)

vmap.rebin(1).plot(
    vmin=100, vmax=220, 
    cmap="seismic", 
    ax=axes[0, 1],
)

smap.rebin(1).plot(
    vmin=0, vmax=120, 
    cmap="magma", 
    ax=axes[1, 0],
)

imap.rebin(1).plot(
    vmin=-250/3, vmax=-100/3, 
    cmap="viridis", 
    ax=axes[1, 1],
)
bg_6363 = contcube.select_lambda(6362, 6372).mean(axis=0)
axes[1, 1].contour(
    bg_6363.data,
    levels=[300],
    colors="r",
)

fig.tight_layout();
# -
# Use the same mask as we used for 6300


sky_6363 = core_6363.copy() 
sky_6363.mask = sky_6363.mask | ~msky[None, : :]

sky_6363.mean(axis=(1, 2)).plot()

corr_6363 = core_6363 - sky_6363.mean(axis=(1, 2))

testpixels = [
    [250, 160], [150, 150], [160, 220],
    [70, 250], [75, 200], [310, 225],
    [25, 140], [250, 250], [140, 110], #[180, 290],
]
fig, axes = plt.subplots(
    3, 3, 
    figsize=(10, 8), 
    sharex=True,
    sharey="row",
)
for (j, i), ax in zip(testpixels, axes.flat):
    core_6363[:, j, i].plot(ax=ax)
    corr_6363[:, j, i].plot(ax=ax) 
    ax.set(xlabel="", ylabel="")
    ax.set_title(f"[{j}, {i}]")
fig.suptitle(
    "Before/after sky correction for faint/moderate/bright pixels"
)
sns.despine()
fig.tight_layout();

mom6363c = moments.find_moments(corr_6363)

g = moments.moments_corner_plot(
    mom6363c, rebin=1, **plot_pars_6363
)

g = moments.moments_corner_plot(
    mom6363c, rebin=4, **plot_pars_6363
)

# +
fig, axes = plt.subplots(
    2, 2, 
    figsize=(10, 10),
    sharex=True, sharey=True,
)

imap = mom6363c[0].copy()
vmap = 3e5*(mom6363c[1] / 6363.78 - 1.0)
smap = 3e5*(mom6363c[2] / 6363.78)

m = imap.data > 0.

vmap.mask = vmap.mask | (~m)
smap.mask = smap.mask | (~m)

vmap_old = 3e5*(mom6363[1] / 6363.78 - 1.0)
vmap_old.mask = vmap_old.mask | (~m)


imap.rebin(1).plot(
    vmin=0, vmax=1e4/3, 
    cmap="turbo", 
    ax=axes[0, 0],
)

vmap.rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    ax=axes[0, 1],
)

smap.rebin(1).plot(
    vmin=40, vmax=120, 
    cmap="magma", 
    ax=axes[1, 0],
)

vmap_old.rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    ax=axes[1, 1],
)


fig.tight_layout();
# -

mom_pars_6363 = dict(
    restwav=6363.78,
    irange=[-150, 1.3e4],
    vrange=[85, 255],
    srange=[30, 150],    
)
moments.save_moments_to_fits(
    mom6363c,
    label="6363",
    flabel="ngc346-oi",
    **mom_pars_6363,
)

# ## Other lines in the wav range: [S III] and Si II

core_6312 = (moi6300cube - mcontcube).select_lambda(6311, 6321)
mom6312 = moments.find_moments(core_6312)



mom_pars_6312 = dict(
    restwav=6312.06,
    irange=[-150, 4e4],
    vrange=[85, 255],
    srange=[20, 120],    
)

plot_pars_6312=dict(
    ilabel="S III",
    label="6312",
    flabel="ngc346-siii",
    **mom_pars_6312,
)
g = moments.moments_corner_plot(
    mom6312, rebin=1, **plot_pars_6312,
)

g = moments.moments_corner_plot(
    mom6312, rebin=4, **plot_pars_6312,
)

(3e5*(mom6312[1] / 6312.06 - 1.0)).rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    colorbar="v",
);

moments.save_moments_to_fits(
    mom6312,
    label="6312",
    flabel="ngc346-siii",
    **mom_pars_6312,
)

cont_6312 = mcontcube.select_lambda(6311, 6321).mean(axis=0)

cont_6312.write("../data/ngc346-cont-6312-mean.fits")



# ### The Si II lines
#
# I will not re-run these for now because there is nothing to be gained for doing the co-added cube, since there is no strong diffuse component to any of them.

# Si II 6347

core_6347 = (oi6300cube - contcube).select_lambda(6345, 6355)
mom6347 = moments.find_moments(core_6347)

mom_pars_6347 = dict(
    restwav=6347.11,
    irange=[-150, 4e4],
    vrange=[85, 255],
    srange=[30, 150],    
)

plot_pars_6347=dict(
    ilabel="Si II",
    label="6347",
    flabel="ngc346-si-ii",
    **mom_pars_6347,
)
g = moments.moments_corner_plot(
    mom6347, rebin=1, **plot_pars_6347,
)

g = moments.moments_corner_plot(
    mom6347, rebin=4, **plot_pars_6347,
)

(3e5*(mom6347[1] / 6347.11 - 1.0)).rebin(4).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    colorbar="v",
)

moments.save_moments_to_fits(
    mom6347,
    label="6347",
    flabel="ngc346-sharp-si-ii",
    **mom_pars_6347,
)

core_6371 = (oi6300cube - contcube).select_lambda(6371, 6381)
mom6371 = moments.find_moments(core_6371)

mom_pars_6371 = dict(
    restwav=6371.36,
    irange=[-150, 4e4],
    vrange=[45, 305],
    srange=[30, 200],    
)

plot_pars_6371=dict(
    ilabel="Si II",
    label="6371",
    flabel="ngc346-si-ii",
    **mom_pars_6371,
)
g = moments.moments_corner_plot(
    mom6371, rebin=1, **plot_pars_6371,
)

g = moments.moments_corner_plot(
    mom6371, rebin=4, **plot_pars_6371,
)

(3e5*(mom6371[1] / 6371.36 - 1.0)).rebin(1).plot(
    vmin=120, vmax=180, 
    cmap="seismic", 
    colorbar="v",
)

moments.save_moments_to_fits(
    mom6371,
    label="6371",
    flabel="ngc346-sharp-si-ii",
    **mom_pars_6371,
)

# ~~Surprisingly, Si II 6371 has a diffuse component, but Si II 6347 does not~~  No, that was just because we were contaminated by [O I] 6363.  I have fixed that now.


