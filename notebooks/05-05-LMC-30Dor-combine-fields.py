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

# # Test of mosaicking together the 30 Dor fields: A, B, C, & D

# +
from pathlib import Path
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns

from mpdaf.obj import Image
from astropy.io import fits
from astropy.wcs import WCS
import reproject

sns.set_context("talk")
sns.set_color_codes()
# -

# ## Define the target WCS for mosaicking the 4 fields together
#
# Use the image from Castro to define the WCS for reprojection.

mosaic_ha_hdulist = fits.open("../data/MUSE_R136toWill/GAUS_Ha6562.8_060_Will.fits")
mosaic_ha_hdulist.info()

mosaic_hdr = mosaic_ha_hdulist[1].header

# Note that the header still refers to a 3rd dimension, even though this is just an image.  We need to remove all the keywords that mention dimension 3 so that the reprojection will work.

del mosaic_hdr["*3"]
del mosaic_hdr["CD3_*"]

# ## Test with the Hβ image.
#
# First we make an HDU of each field for the full mosaic.

pieces = {}
for field in "ABCD":
    infile = f"../data/lmc-30dor-{field}-hi-4861-bin01-sum.fits"
    hdu = fits.open(infile)["DATA"]
    del hdu.header["*3"]
    del hdu.header["CD3_*"]
    newdata, footprint = reproject.reproject_interp(
        hdu,
        mosaic_hdr,
    )
    pieces[field] = newdata

# Plot each piece on top of one another:

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
vmax = 5e5
ax.imshow(pieces["A"], vmin=0, vmax=vmax)
ax.imshow(pieces["B"], vmin=0, vmax=vmax)
ax.imshow(pieces["C"], vmin=0, vmax=vmax)
ax.imshow(pieces["D"], vmin=0, vmax=vmax)
ax.set_title("Hβ overlapping pieces")

# That looks fine - the edges of the pieces show slight artifacts, but they are not very noticeable.
#
# Now try combining the images with median:

combo = np.nanmedian(
    np.stack(list(pieces.values())),
    axis=0,
)
combo.shape

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
vmax = 5e5
ax.imshow(combo, vmin=0, vmax=vmax)
ax.set_title("Hβ median mosaic")

# This time the artifacts are even less visible.
#
# We do seem to cut off a small amount of the lower left field at the bottom.
#
# ## Compare with the Hα image from Castro:

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
haflux = mosaic_ha_hdulist[1].data * mosaic_ha_hdulist[3].data
ax.imshow(haflux, vmin=0, vmax=10e5)
ax.set_title("Hα from Castro")

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
haflux = mosaic_ha_hdulist[1].data * mosaic_ha_hdulist[3].data
decrement = haflux / combo
ax.imshow(decrement, vmin=1.5, vmax=6.0)
ax.contour(
    haflux,
    levels=[2.5e5, 5e5, 10e5],
    linewidths=[0.5, 1.0, 1.5],
    colors="w",
)
ax.contour(2 * combo, levels=[2.5e5, 5e5, 10e5], linewidths=[0.5, 1.0, 1.5], colors="r")
ax.set_title("Hα / Hβ")

# It looks like there is an astrometric offset here
#
# The Castro paper says that they based the astrometric calibration on catalog of Selman (1999).
#
# There is a line-emitting star near the center, which we could use to adjust the astrometry.

# ## Try and fix the astrometry
#
# I looked at the same star on the Castro images and on our images (field A) and got the coordinates from DS9:

# +
from astropy.coordinates import SkyCoord

c_castro = SkyCoord(ra=84.6602735, dec=-69.1035065, unit="deg")
c_eso = SkyCoord(ra=84.6599803, dec=-69.1036748, unit="deg")

Shift_1 = c_castro.ra - c_eso.ra
Shift_2 = c_castro.dec - c_eso.dec
Shift_1.deg, Shift_2.deg
# -

# It looks like field C has a slightly different shift.  I am using a different star for this:

# +
c_castro_C = SkyCoord(84.6689140, -69.0992165, unit="deg")
c_eso_C = SkyCoord(84.6684931, -69.0995408, unit="deg")

Shift_1C = c_castro_C.ra - c_eso_C.ra
Shift_2C = c_castro_C.dec - c_eso_C.dec
Shift_1C.deg, Shift_2C.deg
# -

# We need to add these to the `CRVAL` of our images in order to align with Castro:

pieces = []
for field in "ABCD":
    infile = f"../data/lmc-30dor-{field}-hi-4861-bin01-sum.fits"
    hdu = fits.open(infile)["DATA"]
    del hdu.header["*3"]
    del hdu.header["CD3_*"]
    if field == "C":
        hdu.header["CRVAL1"] += Shift_1C.deg
        hdu.header["CRVAL2"] += Shift_2C.deg
    else:
        hdu.header["CRVAL1"] += Shift_1.deg
        hdu.header["CRVAL2"] += Shift_2.deg

    newdata, footprint = reproject.reproject_interp(
        hdu,
        mosaic_hdr,
    )
    pieces.append(newdata)
combo = np.nanmedian(
    np.stack(pieces),
    axis=0,
)
combo.shape

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
haflux = mosaic_ha_hdulist[1].data * mosaic_ha_hdulist[3].data
decrement = haflux / combo
ax.imshow(decrement, vmin=1.5, vmax=6.0)
ax.contour(
    haflux,
    levels=[2.5e5, 5e5, 10e5],
    linewidths=[0.5, 1.0, 1.5],
    colors="w",
)
ax.contour(2 * combo, levels=[2.5e5, 5e5, 10e5], linewidths=[0.5, 1.0, 1.5], colors="r")
ax.set_title("Hα / Hβ")

# That looks a lot better.  Look at it again with emphasis on the lower extinction

fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
haflux = mosaic_ha_hdulist[1].data * mosaic_ha_hdulist[3].data
decrement = haflux / combo
ax.imshow(decrement, vmin=1.5, vmax=3.0)
ax.contour(
    haflux,
    levels=[2.5e5, 5e5, 10e5],
    linewidths=[0.5, 1.0, 1.5],
    colors="w",
)
ax.contour(2 * combo, levels=[2.5e5, 5e5, 10e5], linewidths=[0.5, 1.0, 1.5], colors="r")
ax.set_title("Hα / Hβ")

# In the left-hand 2 fields, the alignment is not perfect.  We could maybe fix this by looking at the continuum images, but I will leave it for now.

# ## Now process all the emission lines and save them to files

afiles = sorted(Path("../data").glob("lmc-30dor-A-*bin01-*.fits"))

afiles

for afile in afiles:
    pieces = []
    for field in "ABCD":
        infile = str(afile).replace("-A-", f"-{field}-")
        hdu = fits.open(infile)["DATA"]
        del hdu.header["*3"]
        del hdu.header["CD3_*"]
        if field == "C":
            hdu.header["CRVAL1"] += Shift_1C.deg
            hdu.header["CRVAL2"] += Shift_2C.deg
        else:
            hdu.header["CRVAL1"] += Shift_1.deg
            hdu.header["CRVAL2"] += Shift_2.deg

        newdata, footprint = reproject.reproject_interp(
            hdu,
            mosaic_hdr,
        )
        pieces.append(newdata)
    combo = np.nanmedian(
        np.stack(pieces),
        axis=0,
    )
    outfile = str(afile).replace("-A-", "-ABCD-")
    fits.PrimaryHDU(header=mosaic_hdr, data=combo).writeto(outfile, overwrite=True)

afile.stat().st_mtime

Path(outfile).stat().st_mtime


