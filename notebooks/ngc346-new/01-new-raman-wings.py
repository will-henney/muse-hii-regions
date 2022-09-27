# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# # Looking for the H alpha Raman wings in SMC-N66 (NGC 346) with new cube from Peter Zeidler

# Start off the same as the other notebooks - import the libraries and load the data

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
sns.set_context("talk")

datapath = Path("/Users/will/Work/Muse-Hii-Data/SMC-NGC-346/")
fitsfilepath = datapath / "PeterZeidler" / "DATACUBE_FINAL_fwhm_cor.fits"
cube = Cube(str(fitsfilepath))
hacube = cube.select_lambda(6200.0, 6800.0)
cube_mask_orig = cube.mask.copy()
hacube_mask_orig = hacube.mask.copy()

# + pycharm={"name": "#%%\n"}
savepath = Path.cwd().parent.parent / "data-peter"
savepath
# -

# ## No longer need the sky correction
#
# Unlike with the standard pipeline cube, there is no pbvious problem with sky oversubtraction near H alpha in this cube.

cont6600 = hacube.select_lambda(6600.0, 6620.0).mean(axis=0)
hacore = (hacube.select_lambda(6560.0, 6572.0) - cont6600).sum(axis=0)
fig = plt.figure(figsize=(10, 10))
hacore.plot(scale="sqrt", vmin=5e4, vmax=2e5, cmap='viridis', colorbar='v')

# So this compares favorably with the previous version. Strangely, all the fluxes are higher. I am not sure why.

# ## Inspect the data cube
#
# The full data cube is in `cube`, while `hacube` is a 600 Å window around Hα.
#
# Looking at the data cube in DS9 I found a region that looks promising fro the Raman wings, so we will look at that first:

# +
# hacube.mask = (my_mask_3d) | hacube_mask_orig
fig, ax = plt.subplots(figsize=(10, 4))
spec_ha = hacube[:, 240:250, 230:250].mean(axis=(1, 2))
spec2 = spec_ha.copy().subspec(6400, 6650)
spec2.mask_region(6500, 6600)
cont = spec2.poly_spec(2)
spec_ha.plot()
cont.plot(color="m")
ax.axvline(6578, c="r", lw=0.3)
ax.axvline(6583, c="r", lw=0.3)
ax.axvline(6591, c="r", lw=0.3)
ax.axvline(6600, c="r", lw=0.3)
ax.axvline(6620, c="r", lw=0.3)

ax.axvline(6540, c="c", lw=0.3)
ax.axvline(6549, c="c", lw=0.3)



ax.set(ylim=[0, 500])
# -

# ## Images of different wavelength ranges
#
# So the red wing and blue wing look like they are definitely there. The red wing can be seen on both sides of the [N II] λ6583 line. In Orion, we only see it clearly for $\lambda > 6600$.  Perhaps this is because the C II λ6578 line is weaker in the LMC.

#

# First, make an image of the outer red wing:

cont6600 = hacube.select_lambda(6600.0, 6620.0).mean(axis=0)
wing = (hacube.select_lambda(6591.0, 6600.0) - cont6600).sum(axis=0)
fig = plt.figure(figsize=(10, 10))
wing.mask = wing.mask | (cont6600.data > 2000.0)
wing.plot(
    use_wcs=True,
    vmin=0,
    vmax=100,
    cmap="viridis",
    scale="linear",
    colorbar="v",
)
fig.axes[0].set_title(
    "Hα far red wing: full resolution",
    fontsize="large",
    pad=25,
);

# We see filamentary structure that is similar to the other emission lines (see 11-line-profiles notebook).  Note that I am masking out the stars by using a condition on `cont6600`.  It is a delicate balancing act between cutting out the PSF wings of the srars without losing too much of the nebular emission.
#
#
# The map is very noisy but we can try and improve things by rebinning to 8x8.  I also take the opportunity to trim off a 10-pixel margin all around, since there are some very noisy pixels there.

fig = plt.figure(figsize=(10, 10))
margin = 10
wing.mask[:margin, :] = True
wing.mask[-margin:, :] = True
wing.mask[:, :margin] = True
wing.mask[:, -margin:] = True
wing.rebin(4).plot(
    use_wcs=True,
    vmin=0,
    vmax=100,
    cmap="viridis",
    scale="linear",
    colorbar="v",
)
fig.axes[0].set_title(
    r"Hα far red wing: $4\times4$ binning",
    fontsize="large",
    pad=25,
);

# This looks a lot better.  We can see the central filament a lot more clearly.

inwing = (hacube.select_lambda(6578.0, 6583.0) - cont6600).sum(axis=0)
fig = plt.figure(figsize=(10, 10))
inwing.mask = inwing.mask | (cont6600.data > 2000.0)
margin = 10
inwing.mask[:margin, :] = True
inwing.mask[-margin:, :] = True
inwing.mask[:, :margin] = True
inwing.mask[:, -margin:] = True
inwing.rebin(4).plot(
    use_wcs=True,
    vmin=0,
    vmax=100,
    cmap="viridis",
    scale="linear",
    colorbar="v",
)
fig.axes[0].set_title(
    r"Hα near red wing: $4\times4$ binning",
    fontsize="large",
    pad=25,
);

bluewing = (hacube.select_lambda(6540.0, 6549.0) - cont6600).sum(axis=0)
fig = plt.figure(figsize=(10, 10))
bluewing.mask = bluewing.mask | (cont6600.data > 2000.0)
bluewing.mask[:margin, :] = True
bluewing.mask[-margin:, :] = True
bluewing.mask[:, :margin] = True
bluewing.mask[:, -margin:] = True
bluewing.rebin(4).plot(
    use_wcs=True,
    vmin=0,
    vmax=150,
    cmap="viridis",
    scale="linear",
    colorbar="v",
)
fig.axes[0].set_title(
    r"Hα blue wing: $4\times4$ binning",
    fontsize="large",
    pad=25,
);

hacore = (hacube.select_lambda(6555.0, 6578.0) - cont6600).sum(axis=0)
fig = plt.figure(figsize=(10, 10))
hacore.mask = hacore.mask | (cont6600.data > 10000.0)
hacore.rebin(1).plot(
    use_wcs=True,
    vmin=0,
    vmax=2e5,
    cmap="magma",
    scale="linear",
    colorbar="v",
)
fig.axes[0].set_title(
    "Hα line core: full resolution",
    fontsize="large",
    pad=25,
);

# ## Extract spectrum for rectangular regions
#
# Make some box regions for extracting the spectra.  I use the astropy affiliated `regions` package (see [docs](https://astropy-regions.readthedocs.io/en/latest/))

import regions

# I make four boxes. The first two are in regions where the Raman wing is strong, while the second two are where it is weak.

boxes = [
    # regions.BoundingBox(iymin=100, iymax=140, ixmin=15, ixmax=40),
    regions.BoundingBox(iymin=75, iymax=140, ixmin=15, ixmax=40),
    # regions.BoundingBox(iymin=180, iymax=240, ixmin=240, ixmax=300),
    regions.BoundingBox(iymin=200, iymax=250, ixmin=210, ixmax=300),
    regions.BoundingBox(iymin=10, iymax=50, ixmin=100, ixmax=150),
    regions.BoundingBox(iymin=10, iymax=100, ixmin=200, ixmax=300),
    regions.BoundingBox(iymin=170, iymax=210, ixmin=90, ixmax=120),
]

# Plot an image of the entire bandbass in pixel coordinates and plot the boxes on top of it:

fig, ax = plt.subplots(figsize=(10, 10))
wing.rebin(1).plot(
    vmin=0,
    vmax=50,
    scale="linear",
)
for box, c in zip(boxes, "brmgc"):
    box.plot(
        ax=ax,
        lw=3,
        edgecolor=c,
        linestyle="dashed",
        facecolor=(1.0, 1.0, 1.0, 0.4),
        fill=True,
    );

# We can get the pixel slices from each box like this:

boxes[0].slices

# So we can extract the spectrum for each box.

hacube.mask = hacube.mask
hacube.sum(axis=0).plot(scale='log')

# I plot the spectrum for each box and also fit the continuum, so we can easily see if there are any Raman wings present.

# +
fig, ax = plt.subplots(figsize=(10, 6))
offset = 0.0
for box, c in zip(boxes, "brmgc"):
    yslice, xslice = box.slices
    spec = hacube[:, yslice, xslice].mean(axis=(1, 2))
    spec /= spec[-1]
    spec += offset
    spec.plot(c=c, linewidth=1.5, label=f"box {c}")
    spec2 = spec.copy()
    #spec2.mask_region(6200, 6320)
    spec2.mask_region(6270, 6320)
    spec2.mask_region(6350, 6390)
    spec2.mask_region(6450, 6695)
    spec2.mask_region(6712, 6740)
    cont = spec2.poly_spec(4)
    cont.plot(c="k", linewidth=0.5)
    spec2.plot(c="y", linewidth=10, alpha=0.7, zorder=-100)
    offset += 0.1

ax.axvline(6633.0 * (1.0 + 160/3e5), color="k", lw=0.5)
ax.axvline(6664.0 * (1.0 + 160/3e5), color="k", lw=0.5)

ax.axvspan(6594.20, 6611.20, color="r", alpha=0.3, zorder=-100)
ax.axvspan(6612.05, 6628.20, color="r", alpha=0.2, zorder=-100)
ax.axvspan(6638.40, 6660.5, color="r", alpha=0.1, zorder=-100)

ax.legend()
ax.set(ylim=[0.9,2.0], ylabel=r"$F_\lambda / F_{6800}$ + offset")
sns.despine()
# -

# The two boxes that were selected to cover the central filament (blue and red) show clear wings.  Better seen on the red side, but also there on the blue.
#
# The other two boxes, which are off the filament (magenta and green) show no Raman wings at all.
#
# The thick yellow lines show the wavelength sections that are used for fitting the continuum (3rd-order polynomial).
#
# I have put vertical lines at the wavelengths of the 6633 and 6664 features, assuming redshift of 160 km/s.  The 6633 feature is in the middle of two emission lines. What are they? **Could be night sky airglow lines**

# ## Choose suitable bands to measure the Raman wings
#
# In the Orion paper, we have the 3 closest bands in the red wing: R040, R058, R087.  These are marked by pink boxes in the previous figure.
#
# In order to stand a chance of getting good maps, we need to get a better estimate of the continuum than the cont6600 that we used above, since this includes part of the Raman wing that we want to measure.
#
# On the other hand, fitting a 3rd of 4th order polynomial, like we just did for the rectangular boxes is not practical since the individual pixels are too noisy.
#
# A compromise would be to fit a linear trend between the continuum around 6400 and the continuum around 6700.
#
# The above figure shows that the continuum is approximately linear over this range. *We should check if there are any absorption edges that might cause discontinuities in the continuum*

# +
spec4fit = hacube.copy().select_lambda(6390, 6760)

testpixels = [
    [25, 250],
    [50, 80],
    [310, 225],
    [70, 250],
    [75, 200],
    [220, 250],
    [200, 100],
    [20, 230],
    [110, 26],
]
fig, axes = plt.subplots(
    3,
    3,
    figsize=(10, 8),
    sharex=True,
    sharey="row",
)
for (j, i), ax in zip(testpixels, axes.flat):
    spec1d = spec4fit.copy()[:, j, i]
    spec1d.mask_region(6445, 6690)
    spec1d.mask_region(6710, 6745)
    cont = spec1d.poly_spec(1)
    spec4fit[:, j, i].plot(ax=ax)
    cont.plot(ax=ax, linewidth=0.5, c="k")
    spec1d.plot(ax=ax, c="y", linewidth=10, alpha=0.7, zorder=-100)
    ax.set(xlabel="", ylabel="")
    ax.set_title(f"[{j}, {i}]")
axes[0, 0].set(ylim=[-10, 50])
axes[1, 0].set(ylim=[-10, 50])
axes[2, 0].set(ylim=[-10, 50])
fig.suptitle("Test of continuum fitting for selected pixels")
sns.despine()
fig.tight_layout()
# -

# That seems to have worked OK.  Now we need to do this for every pixel.

# ## Fit the continuum pixel-by-pixel to the whole cube

# +
from mpdaf.obj import iter_spe

spec4fit = hacube.copy().select_lambda(6390, 6760)
spec4fit.mask = spec4fit.mask | spec1d.mask[:, None, None]
cont_cube = spec4fit.clone(data_init=np.empty, var_init=np.zeros)

# +
# TOO SLOW - DO NOT RUN!
#
#for sp, co in zip(iter_spe(spec4fit), iter_spe(cont_cube)):
#    if not np.alltrue(sp.mask):
#        co[:] = sp.poly_spec(1)
# -

# This is ridiculously slow.  I need to rewrite it to do my own continuum fitting.  It turns out that the problem is that iterating over the individual `Spectrum` objects in a `Cube` is very inefficient.  So I have to make sure to iterate over just the `.data` components, which are numpy arrays.

from numpy.polynomial import Chebyshev as T
import itertools

nv, ny, nx = spec4fit.shape
wavs = spec4fit.wave.coord()

# +
spec1d = spec4fit.mean(axis=(1, 2))
spec1d.mask_region(6445, 6690)
spec1d.mask_region(6710, 6745)

m = ~spec1d.mask
m.sum(), m.shape
# -

for j, i in itertools.product(range(ny), range(nx)):
    spec1d = spec4fit.data[:, j, i]
    if i == 100 and j % 10 == 0:
        # Give an indication of progress every 10th row
        print(j, m.sum())
    try:
        p = T.fit(wavs[m], spec1d[m], deg=1)
        cont_cube.data[:, j, i] = p(wavs)
    except:
        pass

spec4fit.unmask()

spec4fit[:, 100, 100].mask

cont_cube.sum(axis=0).plot(scale="log")

# + pycharm={"name": "#%%\n"}
(spec4fit - cont_cube).sum(axis=0).plot(vmin=0, vmax=1e6, scale="sqrt")

# + pycharm={"name": "#%%\n"}

fig, ax = plt.subplots(figsize=(10, 6))
offset = 0.0
for box, c in zip(boxes, "brmgc"):
    yslice, xslice = box.slices
    spec = spec4fit[:, yslice, xslice].mean(axis=(1, 2))
    spec /= spec[-1]
    spec += offset
    cont = cont_cube[:, yslice, xslice].mean(axis=(1, 2))
    cont /= cont[-1]
    cont += offset
    spec.plot(c=c, linewidth=1.5, label=f"box {c}")
    cont.plot(c="k", linewidth=0.5)
    offset += 0.1

ax.axvline(6633.0 * (1.0 + 160/3e5), color="k", lw=0.5)
ax.axvline(6664.0 * (1.0 + 160/3e5), color="k", lw=0.5)

ax.axvspan(6594.20, 6611.20, color="r", alpha=0.3, zorder=-100)
ax.axvspan(6612.05, 6628.20, color="r", alpha=0.2, zorder=-100)
ax.axvspan(6638.40, 6660.5, color="r", alpha=0.1, zorder=-100)

ax.legend()
ax.set(ylim=[0.9,2.0], ylabel=r"$F_\lambda / F_{6800}$ + offset")
sns.despine()
# -

# Copy the limits of the Raman bands from the Orion project:

bands = {
    "B133": [6414.85, 6445.45],
    "B080": [6469.25, 6496.45],
    "B054": [6499.85, 6517.7],
    "B033": [6518.55, 6540.65],
    "R040": [6594.2, 6611.2],
    "R058": [6612.05, 6628.2],
    "R087": [6638.4, 6660.5],
    "R136": [6688.55, 6708.95],
}


bandcubes = {}
for band in bands:
    lam1, lam2 = bands[band]
    bandcubes[band] = (spec4fit - cont_cube).select_lambda(lam1, lam2)

# ## Continuum-subtracted images of the red bands

fig, axes = plt.subplots(
    2, 2,
    figsize=(10, 10),
    sharex=True, sharey=True,
)
rbands = list(bands)[4:]
for band, ax in zip(rbands, axes.flat):
    bandcubes[band].sum(axis=0).plot(
        ax=ax, vmin=0, vmax=50,
    )
    ax.set_title(band)
fig.tight_layout()

# So, we can se *something* in all the bands, but the first two are best.
#
# I am in two minds whether to mask out the point sources or not

# ## Continuum-subtracted images of the blue bands

fig, axes = plt.subplots(
    2, 2,
    figsize=(10, 10),
    sharex=True, sharey=True,
)
bbands = reversed(list(bands)[:4])
for band, ax in zip(bbands, axes.flat):
    bandcubes[band].sum(axis=0).plot(
        ax=ax, vmin=0, vmax=50,
    )
    ax.set_title(band)
fig.tight_layout();

# The last band does not show anything, but the other three do.

# + pycharm={"name": "#%%\n"}
for band in bandcubes:
    savefile = savepath / f"ngc346-raman-{band}-peter.fits"
    bandcubes[band].sum(axis=0).write(
        savefile,
        savemask="nan",
        checksum=True,
    )
# -

# ## Ratio of wing to core

hacore = (spec4fit - cont_cube).select_lambda(6560.0, 6572.0)

redsum = (bandcubes["R040"].sum(axis=0) +  bandcubes["R058"].sum(axis=0))
r = redsum / hacore.sum(axis=0)
starmask = cont_cube.sum(axis=0).data > 10*hacore.sum(axis=0).data
faintmask = hacore.sum(axis=0).data < 5e3
r.mask = r.mask | (r.data < 0.0) | starmask | (redsum.data < 10.0) | faintmask

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 10))

r.rebin(2).plot(vmin=0, vmax=0.002, cmap="gray_r", scale="linear", cb="v")
ax.contour(
    hacore.sum(axis=0).rebin(2).data,
    levels=[8e4, 12e4, 16e4, 32e4],
    cmap="YlOrRd_r",
)

# + pycharm={"name": "#%%\n"}
hacore.sum(axis=0).plot(vmin=0, vmax=2e5)
# -

# This shows a clear difference in distribution between the wings and the core.
#
# - [ ] **TODO** I should compare with the [S II]/Ha ratio

# ## Look at the very brightest Raman pixels
#
# These are all close to stars, but they are not actually stellar emission.

# +
cube_select = spec4fit - cont_cube
brightest_wing = redsum.data > 200.0
cube_select.mask = cube_select.mask | ~brightest_wing

# Trim off noisy bands close to edges
cube_select.mask[:, :, 310:] = True
#cube_select.mask[:, :10, :] = True
#cube_select.mask[:, :, :15] = True
cube_select.mask[:, -15:, :] = True

cont_cube_select = cont_cube.copy()
cont_cube_select.mask = cube_select.mask

lam1, lam2 = bands["R040"][0], bands["R058"][1]

fig, ax = plt.subplots(figsize=(10, 10))
(cube_select
 .select_lambda(lam1, lam2)
 .sum(axis=0)
 .plot(vmin=0, vmax=10000, scale="sqrt")
);

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 5))
cube_select.sum(axis=(1, 2)).plot(
    label="continuum-subtracted",
)
(0.03*cont_cube_select.sum(axis=(1, 2))).plot(
    label="0.03 x continuum",
)
ax.legend()
ax.axhline(0.0, c="k", linewidth=1)
ax.axhline(5e4, c="k", linewidth=0.5, linestyle="dashed")
ax.axvline(6633.0 * (1.0 + 160/3e5), color="k", lw=0.5)
ax.axvline(6664.0 * (1.0 + 160/3e5), color="k", lw=0.5)
ax.set(ylim=[-1e5, 5e5]);
# -

# ***We see the Raman absorption lines!!***
#
# Note that this spectrum is dominated by the bright knot seen in the emission line maps, which seems to be associated with a young high-mass star.
#
# There are two other bright knots too.
#
# The total brightness of the Raman wings from these knots is higher than that of the diffuse emission from the entire rest of the nebula (see below).

# *Note that multiplying a `Cube` by a constant destroys the mask, so we need to do the sum before multiplying*

# Another thing – some of this flux comes from some very bright pixels towards the edge of the map.  I had originally masked them all out (see the comment `# Trim off noisy bands close to edges` above), but I have reinstated the bottom and left edges, since they made quite a difference.
#
# - [ ] **TODO:** investigate this further

# ## Look at the diffuse but still bright Raman pixels

# + pycharm={"name": "#%%\n"}
cube_select2 = spec4fit - cont_cube

good_wing = (redsum.data > 80.0)
cube_select2.mask = cube_select2.mask | brightest_wing | ~good_wing | starmask

# Trim off noisy bands close to edges
cube_select2.mask[:, :, 310:] = True
cube_select2.mask[:, :10, :] = True
cube_select2.mask[:, :, :15] = True
cube_select2.mask[:, -15:, :] = True

cont_cube_select2 = cont_cube.copy()
cont_cube_select2.mask = cube_select2.mask

lam1, lam2 = bands["R040"][0], bands["R058"][1]

fig, ax = plt.subplots(figsize=(10, 10))
(cube_select2
 .select_lambda(lam1, lam2)
 .sum(axis=0).rebin(1)
 .plot(vmin=0, vmax=150, cmap="magma_r")
);

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 5))
cube_select2.sum(axis=(1, 2)).plot(
    label="diffuse Raman wings",
)
(0.03*cont_cube_select2.sum(axis=(1, 2))).plot(
    label="0.03 x continuum (diffuse)",
)
cube_select.sum(axis=(1, 2)).plot(
    label="compact Raman wings",
    linewidth=1, alpha=0.6,
    color="r",
)
ax.legend()
ax.axhline(0.0, c="k", linewidth=1)
ax.axvline(6633.0 * (1.0 + 160/3e5), color="k", lw=0.5)
ax.axvline(6664.0 * (1.0 + 160/3e5), color="k", lw=0.5)
ax.set(ylim=[-1e5, 5e5]);
# -

# ## Sources of interest
#
# 1. The bright knot (j, i, = 147, 122)
# 2. The cometary globule (j, i = 180, 109)
# 3. The corner globule (j, i = 40, 40) – has a faint embedded Hα point source with high linewidth, but the Ramn scattering is more extended.
# 4. The blue ball (j, i = 99, 59).  Much brighter in the innermost bands.  Perhaps it is not Raman scattering at all, but is a fast outflow?
#
# Actually, I will come back to this later. It would be easier to take notes in the org file.

# + pycharm={"name": "#%%\n"}
cube_save = spec4fit - cont_cube
cube_save.write(
    savepath / "ngc346-ha-plus-wings-cube-peter.fits",
    savemask="nan",
    checksum=True,
)
# -


