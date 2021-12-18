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

# # Multi-resolution maps of CO in 30 Dor
#
# The idea is to combine the ALMA maps with the APEX maps to fill in the short baselines and get the extended emission.
#
# We already have the Alma map, so we need to get the APEX one

# ## APEX 12CO map from Okada et al (2019)

# We have a cube, so we will sum it in velocity
#
# ### Write out a summed image file

# +
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy import constants
from reproject import reproject_interp

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context("talk")
# -

DATAPATH = Path.cwd().parent / "big-data" / "30-Dor-Radio"

infile = "30Dor_CO21_map_30.fits"
hdulist = fits.open(DATAPATH / infile)
hdulist.info()

# Fix the velocity axis of the header before making the WCS

hdu = hdulist[0]
hdu.header["CRVAL3"] = hdu.header["VELO-LSR"]
hdu.header["CD3_3"] = hdu.header["CDELT3"] = hdu.header["DELTAV"]
w = WCS(hdu.header)
w

imsum = hdu.data[150:181, ...].sum(axis=0)
fits.PrimaryHDU(
    header=w.celestial.to_header(),
    data=imsum,
).writeto(
    DATAPATH / infile.replace(".fits", "-sum.fits"), 
    overwrite=True,
)

# ### Compare with the ALMA map

# We want this section to be independent of the previous one, so we will read the data in from the file again.

from astropy.convolution import convolve_fft, Gaussian2DKernel

infile = "30Dor_CO21_map_30-sum.fits"
hdu, = fits.open(DATAPATH / infile)
w = WCS(hdu.header)

DATAPATH2 = Path.cwd().parent / "data"
infile2 = "lmc-30dor-ABCD-12co-21-reproject-sum.fits"
hdu2, = fits.open(DATAPATH2 / infile2)
w2 = WCS(hdu2.header)

# Make a smooth version of alma image:

w2

# Find the equivalent RMS smoothing width of APEX map in Alma pixels

apex_fwhm = 30.0 * u.arcsec
alma_plate_scale = w2.wcs.cd[1, 1] * u.deg
apex_sigma_pixels = float(apex_fwhm / alma_plate_scale) / np.sqrt(8.0 * np.log(2))
apex_sigma_pixels

kernel = Gaussian2DKernel(apex_sigma_pixels)
im_smooth = convolve_fft(hdu2.data, kernel)

# +
fig, ax = plt.subplots(
    figsize=(10, 10),
    
    subplot_kw=dict(projection=w),
)
im = ax.imshow(hdu.data)
xlim, ylim = ax.get_xlim(), ax.get_ylim()
ax.contour(
    hdu.data, 
    levels=[1.5, 3, 6, 12, 25, 50, 100, 200, 400],
    colors="k",
)


if True:
    ax.contour(
        im_smooth, 
        transform=ax.get_transform(w2),
        levels=[1.5, 3, 6, 12, 25, 50, 100, 200, 400],
        colors="r",
    )
else:
    ax.imshow(
        im_smooth,
        transform=ax.get_transform(w2),
    )
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
fig.colorbar(im, ax=ax)
...;
# -

# So that looks good. The contours of the smoothed ALMA map look very similar to the contours of the APEX map.
#
# The biggest differences are close to the edges of the MUSE field because we haven't got any ALMA data outside of that. *Although we could get some!*

# The strange thing is that there is not a lot of room for any large-scale emission that would be left over from subtracting the smoothed ALMA map. 
#
# Could it be that the APEX map is missing those same small baselines that ALMA is? No - not possible since APEX is single-dish

# ## Make a bigger ALMA map
#
# Currently the ALMA map is restricted to the MUSE field.  We wil lmake another one that extends another 60 arcsec on each side. This will eliminate any edge effects in the MUSE field. 

# +
uids = ["218a", "2192", "219a"]
line_id = "12CO21"
label = "sum"
PREFIX = "Alma-2019.1.00843.S-30_doradus"

data_12co = {
    uid: fits.open(DATAPATH / f"{PREFIX}_{line_id}-{uid}-{label}.fits")["PRIMARY"]
    for uid in uids
}
# -

data_12co

# Make the output projection

hdr = hdu2.header.copy()
margin_arcsec = 60.0
pix_scale = 0.2
margin_pix = int(margin_arcsec / pix_scale)
hdr["NAXIS1"] += 2 * margin_pix
hdr["NAXIS2"] += 2 * margin_pix
hdr["CRPIX1"] += margin_pix
hdr["CRPIX2"] += margin_pix
hdr

images = [
    reproject_interp(hdu, hdr, return_footprint=False)
    for hdu in data_12co.values()
]

bigim = np.nanmedian(
    np.stack(images),
    axis=0,
)
#im[~np.isfinite(im)] = 0.0
fig, ax = plt.subplots(
    figsize=(10, 10),
    subplot_kw=dict(projection=WCS(hdr)),
)
ax.imshow(
    bigim,
     vmin=0.0, 
     vmax=1e2,
    cmap="inferno",
)
ax.contour(
    convolve_fft(bigim, kernel),
    levels=[1.5, 3, 6, 12, 25, 50, 100, 200, 400],
    colors="white",
)
ax.set_aspect("equal")

# This shows the original ALMA image in color map, and the image smoothed to 30 arcsec fwhm in contours. 

# ### Reproject smooth ALMA map to APEX pixels

imr = reproject_interp(
    (convolve_fft(bigim, kernel), hdr),
    hdu.header,
    return_footprint=False,
)

fig, [ax, axx] = plt.subplots(
    1, 2,
    figsize=(15, 8),
    sharex=True,
    sharey=True,
    subplot_kw=dict(projection=w),
)
levels = [1.5, 3, 6, 12, 25, 50, 100, 200, 400]
apex_map = hdu.data - np.nanmin(hdu.data)
im = ax.imshow(apex_map)
ax.contour(
    apex_map, 
    levels=levels,
    colors="k",
)
ax.contour(
    imr, 
    levels=levels,
    colors="r",
)
diff = apex_map -  0.0 * imr
axx.imshow(diff, vmin=im.norm.vmin, vmax=im.norm.vmax)
axx.contour(diff, levels=levels, colors="r")
axx.contour(diff, levels=[0], colors="w")
fig.colorbar(im, ax=[ax, axx])
...;

# So now we can reproject the difference map back onto the MUSE grid:

diff_big = reproject_interp(
    (diff, hdu.header),
    hdr,
    return_footprint=False,
)
diff_muse = reproject_interp(
    (diff, hdu.header),
    hdu2.header,
    return_footprint=False,
)

np.nanmax(diff_muse)

fig, ax = plt.subplots(
    figsize=(15, 15),
    subplot_kw=dict(projection=WCS(hdr)),
)
ax.imshow(
    diff_big + bigim,
#    bigim,
#    diff_muse,
    vmin=0.0,
    vmax=10.0,
    cmap="inferno",
)

fig, ax = plt.subplots(
    figsize=(15, 15),
    subplot_kw=dict(projection=w2),
)
ax.imshow(
    diff_muse + hdu2.data,
#    bigim,
#    diff_muse,
    vmin=0.0,
    vmax=10.0,
    cmap="inferno",
)
ax.contour(
    diff_muse + hdu2.data,
    levels=[0],
    colors="c",
)

fig, ax = plt.subplots(
    figsize=(15, 15),
    subplot_kw=dict(projection=w2),
)
ax.imshow(
    hdu2.data,
    vmin=0.0,
    vmax=10.0,
    cmap="inferno",
)
ax.contour(
    hdu2.data,
    levels=[0],
    colors="c",
)


