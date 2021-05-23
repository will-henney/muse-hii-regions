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

# # NGC 346 stars from Sabbi et al (2007)
#
# In my never-ending quest to get a catalog of all the stars in NGC 346, I am now going to try using the astroquery Vizier interface to get the tables from: 
#
# ```bibtex
# @article{Sabbi:2007h,
# 	author = {{Sabbi}, E. and {Sirianni}, M. and {Nota}, A. and {Tosi}, M. and {Gallagher}, J. and {Meixner}, M. and {Oey}, M.~S. and {Walterbos}, R. and {Pasquali}, A. and {Smith}, L.~J. and {Angeretti}, L.},
# 	journal = {\aj},
# 	month = jan,
# 	number = {1},
# 	pages = {44-57},
# 	title = {{Past and Present Star Formation in the SMC: NGC 346 and its Neighborhood}},
# 	volume = {133},
# 	year = 2007}
#
# ```

# ## Find which catalog we want
#
# I had to install the `astroquery` library from conda-forge.

from astroquery.vizier import Vizier

# I am following a similar method to the examples in the [astroquery docs](https://astroquery.readthedocs.io/en/latest/vizier/vizier.html):

catalog_list = Vizier.find_catalogs('Sabbi 2007')

# We now have a list of matching catalogs:

for k,v in catalog_list.items():
    print(k)
    print(v.description, end="\n\n")

# This gave more than I was expecting, but the one we want is `J/AJ/133/44`:

catalog_id = "J/AJ/133/44"

# ## Download the tables
#
# We have to set `Vizier.ROW_LIMIT = -1` to make sure we get the complete table.  Otherwise, the default is to return only 50 rows.

Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs(catalog_id)
catalogs

# This gives two tables.  The first is a list of all 79960 stars that they have found in NGC 346.  The second is a list of the 16 sub-clusters that they have identified:

stars, subclusters = catalogs
subclusters

# ## Look at the star catalog

stars

# ## Put the coordinate columns in a more useful form
#
# We want to convert the coordinate columns into degrees for ease of plotting, so we need some astropy libraries:

from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude
import astropy.units as u

# I had originally thought of using SkyCoord, but this is not necessary

# +
# SkyCoord(stars["RAJ2000"], stars["DEJ2000"], unit=(u.hourangle, u.deg))
# -

# Easier to just convert the RA and DEC columns separately (note, these cells take a while to run – 10s of seconds on my laptop):

ra = Longitude(stars["RAJ2000"], unit=u.hourangle)

dec = Latitude(stars["DEJ2000"], unit=u.deg)

# Check that the values in degrees look reasonable:

ra.to(u.deg), dec

# Add them in as new columns to the table and set display precision to about 10 mas. 

stars["ra"] = ra.to(u.deg)
stars["dec"] = dec
stars["ra"].format = "{:.6f}"
stars["dec"].format = "{:.6f}"
stars

# ## Add $V - I$ color to the table

# For some reason the units get lost, so I add them explicitly and set the display precision to be the same as for the original photometry columns.  For good measure, I also do the error propagation. 

import numpy as np

# +
stars["V-I"] = stars["F555W"] - stars["F814W"]
stars["V-I"].format = "{:.3f}"
stars["V-I"].unit = u.mag

stars["e_V-I"] = np.hypot(stars["e_F555W"], stars["e_F814W"])
stars["e_V-I"].format = "{:.3f}"
stars["e_V-I"].unit = u.mag

stars
# -

# The stars seem to be listed in order of $V$-band brightness.  The brightest have $V < 12$ and have relatively blue colors, with very small errors on the color.  The faintest have $V \approx 27$ with uniformly red colors, $V - I \approx +1.5$ with color uncertainty $\approx 0.1$ mag. 
#
# Some of the bright stars have $V - I < 0$, indicating that they are OB stars.  The bright ones with $V - I > 0$ could be red giants/AGB/supergiants or could be reddened OB stars.  The foreground reddening is very low, but there is patchy dense molecular gas.

# Comparing with the O stars in Orion.  The brightest stars there have $V \approx 5$.  The distance to Orion is 0.41 pc, versus 61.7 kpc for the SMC, so the equivalent magnitude in NGC 346 would be:

5 + 5*np.log10(61.7/0.41)

# So anything with $V < 16$ is brighter than the Orion stars.  We have lots of those, and the brightest is nearly 100 times brighter, which seems unlikely.  Perhaps it is a foreground contaminant, or maybe a small cluster.  Also, in Orion there is about 1 mag of foreground extinction.

# ## Plot the distribution of stars

# +
from pathlib import Path

import pandas as pd
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.table import Table

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()
sns.set_context("talk")
# -

# Use the same nominal coordinates of NGC 346 that we used in the previous notebook for consistency:

ra0, dec0 = 14.77101, -72.1771

# + active=""
# Find the aspect ratio to keep the pixels square:
# -

aspect = 1.0/np.cos(np.deg2rad(dec0))
print("Angle distorsion factor:", aspect)

# +
fig, ax = plt.subplots(figsize=(20, 20))

points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    s=0.1,
    c=stars["V-I"].data,
    vmin=-1.0,
    vmax=1.0,
    cmap="icefire",
    alpha=1.0,
    label=f'{len(stars)} stars'
)
fig.colorbar(points, ax=ax, label="$V - I$")
ax.invert_xaxis()
ax.set_aspect(aspect)
ax.set(
    xlabel='RA [deg]',
    ylabel='Dec [deg]',
)
ax.legend()
fig.tight_layout()
# -

# We can see an evolved cluster at the top center, which has red stars only.  The blue stars are concentrated in the core of NGC 346, which is in the center of the image. 
#
# Now, we zoom in on the core:

# +
dec_range = 0.5/60.0
ra_range = dec_range * aspect

fig, ax = plt.subplots(figsize=(20, 20))

points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=200.0/(stars["F555W"] - 11.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
)
fig.colorbar(points, ax=ax, label="$V - I$")
ax.plot(ra0, dec0, 'kx', markersize=20)
ax.invert_xaxis()
ax.set_aspect(aspect)
ax.set(
    xlim=[ra0 - ra_range, ra0 + ra_range],
    ylim=[dec0 - dec_range, dec0 + dec_range],
    xlabel='RA [deg]',
    ylabel='Dec [deg]',
)
# -

# Comparing with the HSC results from the previous notebook, there seems to be an offset in the color values, which is odd.  The Sabbi colors are about +0.3 redder than the HSC ones.  I have adjusted the `vmin, vmax` values for the color bar so that the colors look similar to in the previous notebook. 
#
# Apart from that, the results are comparable.  It looks like we have more sources here, which is good.

# ## Combine with Hα map

bigdatapath = Path("../big-data")
hdu = fits.open(bigdatapath / "ngc346-hst-acs-f658n-wcsgaia.fits")["SCI"]
w = WCS(hdu.header)
imha = hdu.data

# I have loaded the ACS H alpha image, which has been aligned to the Gaia coordinate frame. I also load the table of Gaia source positions. So, we will be able to see if there is any problem with the Sabbi coordinates:

datapath = Path("../data")
dfgaia = pd.read_csv(datapath / "1621565655827O-result.csv")

# Now plot everything together.  Just to make misalignments super obvious, I will zoom in even more than I did before:

x0, y0 = 2600, 2200
#dx, dy = 200, 200
dx, dy = 100, 100

# +
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.7, vmax=0.7, 
    cmap="gray_r",
)
points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=400.0/(stars["F555W"] - 11.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
#fig.colorbar(points, ax=ax, label="$V - I$")
gpoints = ax.scatter(
    x="ra", y="dec", 
    edgecolors="w", 
    data=dfgaia,
    s=100/(dfgaia["phot_g_mean_mag"] - 15), 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
);
# -

# It is clear there is an offset in RA with the Sabbi position. It looks to be about 2 ACS pixels, or 0.1 arcsec. So, we will fix that:

rafix = ra + aspect * 0.08 * u.arcsec
rafix

stars["ra"] = rafix.to(u.deg)
stars["ra"].format = "{:.6f}"
stars[:2]

# So that looks OK,  The RA column has changed in the last three digits.

# +
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.7, vmax=0.7, 
    cmap="gray_r",
)
points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=400.0/(stars["F555W"] - 13.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)
gpoints = ax.scatter(
    x="ra", y="dec", 
    edgecolors="w", 
    data=dfgaia,
    s=1200/(dfgaia["phot_g_mean_mag"] - 13), 
    alpha=1.0,
    linewidth=3,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
);
# -

# That's looking a lot better now. 
#
# Now, take a zoomed out view as well:

dx, dy = 1000, 1000

# +
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-1.0, vmax=0.3, 
    cmap="gray_r",
)
points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=40.0/(stars["F555W"] - 13.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=0.5,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)
gpoints = ax.scatter(
    x="ra", y="dec", 
    edgecolors="w", 
    data=dfgaia,
    s=120/(dfgaia["phot_g_mean_mag"] - 13), 
    alpha=1.0,
    linewidth=1.5,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
);
# -

# ## Add in the sub-cluster positions
#
# Put the coordinates in degrees. And also convert the cluster radius to degrees:

ra_sc = Longitude(subclusters["RAJ2000"], unit=u.hourangle)
dec_sc = Latitude(subclusters["DEJ2000"], unit=u.deg)
subclusters["ra"] = ra_sc.to(u.deg)
subclusters["dec"] = dec_sc
subclusters["ra"].format = "{:.6f}"
subclusters["dec"].format = "{:.6f}"

D_smc = 61.7*u.kpc
R_deg = (subclusters["Rad"].to(u.au).value / (D_smc.to(u.pc)).value) * u.arcsec#.to(u.deg)
R_deg

subclusters["R"] = R_deg
subclusters

from astropy.visualization.wcsaxes import SphericalCircle

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-1.2, vmax=0.3, 
    cmap="gray_r",
)
for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=4,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
    ax.text(
        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
        ha="center", va="center", color="w",
        fontweight="black",
        transform=ax.get_transform("world"),
    )


# +
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-1.0, vmax=0.3, 
    cmap="gray_r",
)
points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=1.0/(stars["F555W"] - 12.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
#fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)
#gpoints = ax.scatter(
#    x="ra", y="dec", 
#    edgecolors="w", 
#    data=dfgaia,
#    s=120/(dfgaia["phot_g_mean_mag"] - 13), 
#    alpha=1.0,
#    linewidth=1.5,
#    facecolor="none",
#    transform=ax.get_transform('world'),
#)
#ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=4,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
    ax.text(
        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
        ha="center", va="center", color="w",
        fontweight="black",
        transform=ax.get_transform("world"),
    )


#ax.set(
#    xlim=[x0 - dx, x0 + dx],
#    ylim=[y0 - dy, y0 + dy],
#);

# +
dx, dy = 1000, 1000

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-1.0, vmax=0.3, 
    cmap="gray_r",
)
points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=3.0/(stars["F555W"] - 12.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)
#gpoints = ax.scatter(
#    x="ra", y="dec", 
#    edgecolors="w", 
#    data=dfgaia,
#    s=120/(dfgaia["phot_g_mean_mag"] - 13), 
#    alpha=1.0,
#    linewidth=1.5,
#    facecolor="none",
#    transform=ax.get_transform('world'),
#)
#ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=4,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
    ax.text(
        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
        ha="center", va="center", color="w",
        fontweight="black",
        transform=ax.get_transform("world"),
        clip_on=True,
    )


ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.tight_layout();
# -

dx, dy = 400, 400

# +
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.5, vmax=0.8, 
    cmap="gray_r",
)
points = ax.scatter(
    stars['ra'].data, 
    stars['dec'].data, 
    c=stars["V-I"].data,
    s=100.0/(stars["F555W"] - 12.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)
#gpoints = ax.scatter(
#    x="ra", y="dec", 
#    edgecolors="w", 
#    data=dfgaia,
#    s=120/(dfgaia["phot_g_mean_mag"] - 13), 
#    alpha=1.0,
#    linewidth=1.5,
#    facecolor="none",
#    transform=ax.get_transform('world'),
#)
#ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

for data in subclusters[0:3]:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='k', facecolor='none', 
        linestyle="dashed",
        linewidth=4,
        clip_on=True,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
    ax.text(
        data["ra"], data["dec"], data["__SSN2007_"],
        ha="center", va="center", color="k",
        fontweight="black", fontsize="x-large",
        transform=ax.get_transform("world"),
    )


ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.tight_layout();
# -

# ## Try and apply CMD cuts to isolate YSOs
#
# Annoyingly, the catalog does not include the estimated age of the stars.  
#
# But it looks like $V < 20 + 4 (V - I)$ would make a good cut.  This would select the upper main sequence plus the young stars above the lower main sequence (plus evolved stars, but there shouldn't be many of those).

mcut = stars["F555W"] < 20.0 + 4*stars["V-I"]

mcut.sum(), (~mcut).sum()

ystars = stars[mcut][::-1]
ystars

lmsstars = stars[~mcut]

# +
dx, dy = 2000, 2000
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-1.0, vmax=0.3, 
    cmap="gray_r",
)
points = ax.scatter(
    ystars['ra'].data, 
    ystars['dec'].data, 
    c=ystars["V-I"].data,
    s=10.0/(ystars["F555W"] - 12.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
#fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)
#gpoints = ax.scatter(
#    x="ra", y="dec", 
#    edgecolors="w", 
#    data=dfgaia,
#    s=120/(dfgaia["phot_g_mean_mag"] - 13), 
#    alpha=1.0,
#    linewidth=1.5,
#    facecolor="none",
#    transform=ax.get_transform('world'),
#)
#ax.plot(ra0, dec0, 'kx', markersize=20, transform=ax.get_transform('world'),)

for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=4,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
    ax.text(
        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
        ha="center", va="center", color="w",
        fontweight="black",
        clip_on=True,
        transform=ax.get_transform("world"),
    )



ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
#fig.tight_layout();
fig.savefig("../figs/ngc-346-star-map-zoom-A.pdf");

# +
dx, dy = 500, 500
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    imha, 
    vmin=0.2, vmax=1.0, 
    cmap="gray_r",
)
points = ax.scatter(
    ystars['ra'].data, 
    ystars['dec'].data, 
    c=ystars["V-I"].data,
    s=100.0/(ystars["F555W"] - 12.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)

for data in subclusters[[0, 1, 2, 6, 7]]:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='k', facecolor='none', 
        linestyle="dashed",
        linewidth=4,
        clip_on=True,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
    ax.text(
        data["ra"], data["dec"], data["__SSN2007_"],
        ha="center", va="center", color="k",
        fontweight="black", fontsize="x-large",
        transform=ax.get_transform("world"),
    )



ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.savefig("../figs/ngc-346-star-map-zoom-B.pdf")

# +
dx, dy = 100, 100
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.5, vmax=1.0, 
    cmap="gray_r",
)
points = ax.scatter(
    ystars['ra'].data, 
    ystars['dec'].data, 
    c=ystars["V-I"].data,
    s=400.0/(ystars["F555W"] - 13.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)

points2 = ax.scatter(
    lmsstars['ra'].data, 
    lmsstars['dec'].data, 
    edgecolors="w",
    facecolor="none",
    s=400.0/(lmsstars["F555W"] - 13.0),
    alpha=1.0,
    linewidth=3,
    transform=ax.get_transform('world'),    
)


ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.savefig("../figs/ngc-346-star-map-zoom-C.pdf");

# +
dx, dy = 30, 30
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
#    np.log10(imha), 
#    vmin=-0.5, vmax=2.0, 
    np.sqrt(imha), 
    vmin=0.0, vmax=8.0, 
    cmap="gray_r",
)
points = ax.scatter(
    ystars['ra'].data, 
    ystars['dec'].data, 
    c=ystars["V-I"].data,
    s=4*400.0/(ystars["F555W"] - 13.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=0.8)

points2 = ax.scatter(
    lmsstars['ra'].data, 
    lmsstars['dec'].data, 
    edgecolors="y",
    facecolor="none",
    s=4*400.0/(lmsstars["F555W"] - 13.0),
    linewidth=3,
    alpha=1.0,
    transform=ax.get_transform('world'),    
)


ax.set(
    xlim=[x0 + 50 - dx, x0 + 50 + dx],
    ylim=[y0 - 40 - dy, y0 - 40 + dy],
)
fig.tight_layout();
# -


