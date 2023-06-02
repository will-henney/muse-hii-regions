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

# + [markdown] pycharm={"name": "#%% md\n"}
# # Gaia astrometry of stars in NGC 346
#
# I want to get good positions for all the stars so I can align the HST and the MUSE images

# + pycharm={"name": "#%%\n"}
import pandas as pd
from pathlib import Path

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Load the Gaia sources

# + pycharm={"name": "#%%\n"}
datapath = Path("../data")

# + pycharm={"name": "#%%\n"}
df = pd.read_csv(datapath / "1621565655827O-result.csv")

# + pycharm={"name": "#%%\n"}
df

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Plot all the Gaia sources

# + pycharm={"name": "#%%\n"}
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()
sns.set_context("talk")

# + [markdown] pycharm={"name": "#%% md\n"}
# I set the color scale to be the magnitude

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 8))
scat = ax.scatter(
    x="ra", y="dec", c="phot_g_mean_mag", 
    data=df,
    s=20,
)
cb = fig.colorbar(scat, ax=ax)

# + [markdown] pycharm={"name": "#%% md\n"}
# So you can see the mass segregation straight away - the darker ones (brighter) tend to be more concentrated towards the center.  It looks like Gaia must have missed quite a lot of them though.
#
# *Note that the RA axis is the wrong way round* 

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Load an HST image
#
# We will start off with the UV 2200 Å broad-band image.  
#
# Originally I had downloaded the data in a strange format, where 3 different filter images were stacked in a cube.  But I have now gone back and grabbed the individual filters as separate files, since this is easier to work with. 
#
# All the HST data files are listed at this DOI [10.17909/t9-vtw1-c475](https://doi.org/10.17909/t9-vtw1-c475)

# + pycharm={"name": "#%%\n"}
import numpy as np

# + pycharm={"name": "#%%\n"}
from astropy.io import fits
from astropy.wcs import WCS

# + pycharm={"name": "#%%\n"}
bigdatapath = Path("../big-data")

# + pycharm={"name": "#%%\n"}
dataset = "hst_10248_03_acs_hrc_f220w"
hdulist = fits.open(bigdatapath / f"HST-NGC346/{dataset}/{dataset}_drz.fits")

# + pycharm={"name": "#%%\n"}
hdulist.info()

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Plot the image in celestial coordinates with Gais sources overlaid

# + pycharm={"name": "#%%\n"}
w = WCS(hdulist["SCI"].header)

# + pycharm={"name": "#%%\n"}
imdata = hdulist["SCI"].data

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(imdata, vmin=0.0, vmax=5.0, cmap="gray_r")
scat = ax.scatter(
    x="ra", y="dec", c="phot_g_mean_mag", 
    data=df,
    s=20, 
    alpha=0.4,
    transform=ax.get_transform('world'),
)

# + [markdown] pycharm={"name": "#%% md\n"}
# This shows that the HST coordinates are not quite right – there is an offset from the Gaia coordinates.
#
# We will zoom in on the central cluster to have a closer look.

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(imdata, vmin=0.0, vmax=5.0, cmap="gray_r")
scat = ax.scatter(
    x="ra", y="dec", edgecolors="r", 
    data=df,
    s=150, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.set(
    xlim=[1000, 1600],
    ylim=[600, 1200],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Adjust the WCS to align HST with Gaia

# + [markdown] pycharm={"name": "#%% md\n"}
# Now, try giving an offset to fix this.  We could either change `CRPIX` or the `CRVAL`, but the first seems simpler to reason about.  We want to move the HST stars to the right and down, so do the opposite to the reference pixel:

# + pycharm={"name": "#%%\n"}
ww = WCS(hdulist["SCI"].header)
ww.wcs.crpix

# + pycharm={"name": "#%%\n"}
ww.wcs.crpix -= np.array([19, -7])

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ww)
ax.imshow(imdata, vmin=0.0, vmax=5.0, cmap="gray_r")
scat = ax.scatter(
    x="ra", y="dec", edgecolors="r", 
    data=df,
    s=150, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.set(
    xlim=[1000, 1600],
    ylim=[600, 1200],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# That looks pretty good. 

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Repeat for the Hα image

# + pycharm={"name": "#%%\n"}
dataset_ha = "hst_10248_a3_acs_wfc_f658n"
hdulist_ha = fits.open(bigdatapath / f"HST-NGC346/{dataset_ha}/{dataset_ha}_drz.fits")

# + pycharm={"name": "#%%\n"}
hdulist_ha.info()

# + pycharm={"name": "#%%\n"}
w_ha = WCS(hdulist_ha["SCI"].header)
imdata_ha = hdulist_ha["SCI"].data

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w_ha)
ax.imshow(imdata_ha, vmin=0.0, vmax=1.0, cmap="gray_r")
scat = ax.scatter(
    x="ra", y="dec", edgecolors="r", 
    data=df,
    s=50, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.set(
    xlim=[2000, 3200],
    ylim=[1600, 2800],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# Zoom in some more and switch to a logarithmic brightness scaling:

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w_ha)
ax.imshow(
    np.log10(imdata_ha), 
    vmin=-1.0, vmax=3.0, 
    cmap="gray_r",
)
scat = ax.scatter(
    x="ra", y="dec", edgecolors="r", 
    data=df,
    s=50, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.set(
    xlim=[2400, 2800],
    ylim=[2000, 2400],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# This time, the whift is up and to the right. Let's try and fix it:

# + pycharm={"name": "#%%\n"}
w_ha_fix = WCS(hdulist_ha["SCI"].header)
w_ha_fix.wcs.crpix += np.array([6, 7.5])
w_ha_fix.wcs.crpix

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w_ha_fix)
ax.imshow(
    np.log10(imdata_ha), 
    vmin=-0.7, vmax=0.7, 
    cmap="gray_r",
)
scat = ax.scatter(
    x="ra", y="dec", edgecolors="r", 
    data=df,
    s=50, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.set(
    xlim=[2400, 2800],
    ylim=[2000, 2400],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Save the corrected versions as FITS files
#
# First the Hα image:

# + pycharm={"name": "#%%\n"}
for hdu in hdulist_ha:
    if hdu.is_image:
        hdu.header.update(w_ha_fix.to_header())

# + pycharm={"name": "#%%\n"}
hdulist_ha.writeto(
    bigdatapath / "ngc346-hst-acs-f658n-wcsgaia.fits",
    overwrite=True,
)

# + [markdown] pycharm={"name": "#%% md\n"}
# Now the UV continuum image

# + pycharm={"name": "#%%\n"}
for hdu in hdulist:
    if hdu.is_image:
        hdu.header.update(ww.to_header())

# + pycharm={"name": "#%%\n"}
hdulist.writeto(
    bigdatapath / "ngc346-hst-acs-f220w-wcsgaia.fits",
    overwrite=True,
)

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Load the DAOPHOT source list
#
# The HST data came with some data tables that seem to be the results of point source extraction.
#
# First look at the Hα ones. 

# + pycharm={"name": "#%%\n"}
from astropy.io import ascii

# + pycharm={"name": "#%%\n"}
stars = ascii.read(
    str(bigdatapath / f"HST-NGC346/{dataset_ha}/{dataset_ha}_daophot_trm.cat"),
    format="commented_header",
    header_start=56,
)

# + pycharm={"name": "#%%\n"}
stars

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w_ha_fix)
ax.imshow(
    np.log10(imdata_ha), 
    vmin=-1.0, vmax=3.0, 
    cmap="gray_r",
)
scat = ax.scatter(
    x="c1", y="c2", edgecolors="r", 
    data=stars,
    s=50, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('pixel'),
)
ax.set(
    xlim=[2400, 2800],
    ylim=[2000, 2400],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# Hmm, lots of stars are missing.  It seems to be avoiding regions with diffuse Hα emission.  We will zoom out to check.

# + pycharm={"name": "#%%\n"}
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=w_ha_fix)
ax.imshow(
    np.log10(imdata_ha), 
    vmin=-0.7, vmax=0.7, 
    cmap="gray_r",
)
scat = ax.scatter(
    x="c1", y="c2", edgecolors="r", 
    data=stars,
    s=50, 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('pixel'),
)
ax.set(
    xlim=[1600, 3600],
    ylim=[1200, 3200],
#    xlim=[2400, 2800],
#    ylim=[2000, 2400],
);

# + [markdown] pycharm={"name": "#%% md\n"}
# Yes, this confirms that the source extraction avoids regions where the nebula is bright.
#
# We will look at the results from the broad-band filters instead, starting with the UV. *Cancel that – I can't make sense of the source list file*

# + pycharm={"name": "#%%\n"}
#uvstars = ascii.read(
#    str(bigdatapath / f"HST-NGC346/{dataset}/{dataset}_daophot_trm.cat"),
#    format="fixed_width_no_header",
#    data_start=409,
#)        

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Try getting stars from Hubble Source Catalog instead
#
# This is now in a separate notebook: `04-01-ngc-346-hsc-sources.ipynb`
