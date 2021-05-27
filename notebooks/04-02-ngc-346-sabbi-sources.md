---
jupyter:
  jupytext:
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# NGC 346 stars from Sabbi et al (2007)

In my never-ending quest to get a catalog of all the stars in NGC 346, I am now going to try using the astroquery Vizier interface to get the tables from: 

```bibtex
@article{Sabbi:2007h,
	author = {{Sabbi}, E. and {Sirianni}, M. and {Nota}, A. and {Tosi}, M. and {Gallagher}, J. and {Meixner}, M. and {Oey}, M.~S. and {Walterbos}, R. and {Pasquali}, A. and {Smith}, L.~J. and {Angeretti}, L.},
	journal = {\aj},
	month = jan,
	number = {1},
	pages = {44-57},
	title = {{Past and Present Star Formation in the SMC: NGC 346 and its Neighborhood}},
	volume = {133},
	year = 2007}

```


## Find which catalog we want

I had to install the `astroquery` library from conda-forge.

```python
from astroquery.vizier import Vizier
```

I am following a similar method to the examples in the [astroquery docs](https://astroquery.readthedocs.io/en/latest/vizier/vizier.html):

```python
catalog_list = Vizier.find_catalogs('Sabbi 2007')
```

We now have a list of matching catalogs:

```python
for k,v in catalog_list.items():
    print(k)
    print(v.description, end="\n\n")
```

This gave more than I was expecting, but the one we want is `J/AJ/133/44`:

```python
catalog_id = "J/AJ/133/44"
```

## Download the tables

We have to set `Vizier.ROW_LIMIT = -1` to make sure we get the complete table.  Otherwise, the default is to return only 50 rows.

```python
Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs(catalog_id)
catalogs
```

This gives two tables.  The first is a list of all 79960 stars that they have found in NGC 346.  The second is a list of the 16 sub-clusters that they have identified:

```python
stars, subclusters = catalogs
subclusters
```

## Look at the star catalog

```python
stars
```

## Put the coordinate columns in a more useful form

We want to convert the coordinate columns into degrees for ease of plotting, so we need some astropy libraries:

```python
from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude
import astropy.units as u
```

I had originally thought of using SkyCoord, but this is not necessary

```python
# SkyCoord(stars["RAJ2000"], stars["DEJ2000"], unit=(u.hourangle, u.deg))
```

Easier to just convert the RA and DEC columns separately (note, these cells take a while to run – 10s of seconds on my laptop):

```python
ra = Longitude(stars["RAJ2000"], unit=u.hourangle)
```

```python
dec = Latitude(stars["DEJ2000"], unit=u.deg)
```

Check that the values in degrees look reasonable:

```python
ra.to(u.deg), dec
```

Add them in as new columns to the table and set display precision to about 10 mas. 

```python
stars["ra"] = ra.to(u.deg)
stars["dec"] = dec
stars["ra"].format = "{:.6f}"
stars["dec"].format = "{:.6f}"
stars
```

## Add $V - I$ color to the table


For some reason the units get lost, so I add them explicitly and set the display precision to be the same as for the original photometry columns.  For good measure, I also do the error propagation. 

```python
import numpy as np
```

```python
stars["V-I"] = stars["F555W"] - stars["F814W"]
stars["V-I"].format = "{:.3f}"
stars["V-I"].unit = u.mag

stars["e_V-I"] = np.hypot(stars["e_F555W"], stars["e_F814W"])
stars["e_V-I"].format = "{:.3f}"
stars["e_V-I"].unit = u.mag

stars
```

The stars seem to be listed in order of $V$-band brightness.  The brightest have $V < 12$ and have relatively blue colors, with very small errors on the color.  The faintest have $V \approx 27$ with uniformly red colors, $V - I \approx +1.5$ with color uncertainty $\approx 0.1$ mag. 

Some of the bright stars have $V - I < 0$, indicating that they are OB stars.  The bright ones with $V - I > 0$ could be red giants/AGB/supergiants or could be reddened OB stars.  The foreground reddening is very low, but there is patchy dense molecular gas.


Comparing with the O stars in Orion.  The brightest stars there have $V \approx 5$.  The distance to Orion is 0.41 pc, versus 61.7 kpc for the SMC, so the equivalent magnitude in NGC 346 would be:

```python
5 + 5*np.log10(61.7/0.41)
```

So anything with $V < 16$ is brighter than the Orion stars.  We have lots of those, and the brightest is nearly 100 times brighter, which seems unlikely.  Perhaps it is a foreground contaminant, or maybe a small cluster.  Also, in Orion there is about 1 mag of foreground extinction.


## Plot the distribution of stars

```python
from pathlib import Path

import pandas as pd
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.table import Table

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()
sns.set_context("talk")
```

Use the same nominal coordinates of NGC 346 that we used in the previous notebook for consistency:

```python
ra0, dec0 = 14.77101, -72.1771
```

<!-- #raw -->
Find the aspect ratio to keep the pixels square:
<!-- #endraw -->

```python
aspect = 1.0/np.cos(np.deg2rad(dec0))
print("Angle distorsion factor:", aspect)
```

```python
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
```

We can see an evolved cluster at the top center, which has red stars only.  The blue stars are concentrated in the core of NGC 346, which is in the center of the image. 

Now, we zoom in on the core:

```python
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
```

Comparing with the HSC results from the previous notebook, there seems to be an offset in the color values, which is odd.  The Sabbi colors are about +0.3 redder than the HSC ones.  I have adjusted the `vmin, vmax` values for the color bar so that the colors look similar to in the previous notebook. 

Apart from that, the results are comparable.  It looks like we have more sources here, which is good.


## Combine with Hα map

```python
bigdatapath = Path("../big-data")
hdu = fits.open(bigdatapath / "ngc346-hst-acs-f658n-wcsgaia.fits")["SCI"]
w = WCS(hdu.header)
imha = hdu.data
```

I have loaded the ACS H alpha image, which has been aligned to the Gaia coordinate frame. I also load the table of Gaia source positions. So, we will be able to see if there is any problem with the Sabbi coordinates:

```python
datapath = Path("../data")
dfgaia = pd.read_csv(datapath / "1621565655827O-result.csv")
```

Now plot everything together.  Just to make misalignments super obvious, I will zoom in even more than I did before:

```python
x0, y0 = 2600, 2200
#dx, dy = 200, 200
dx, dy = 100, 100
```

```python
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
```

It is clear there is an offset in RA with the Sabbi position. It looks to be about 2 ACS pixels, or 0.1 arcsec. So, we will fix that:

```python
rafix = ra + aspect * 0.08 * u.arcsec
rafix
```

```python
stars["ra"] = rafix.to(u.deg)
stars["ra"].format = "{:.6f}"
stars[:2]
```

So that looks OK,  The RA column has changed in the last three digits.

```python
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
```

That's looking a lot better now. 

Now, take a zoomed out view as well:

```python
dx, dy = 1000, 1000
```

```python
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
```

## Add in the sub-cluster positions

Put the coordinates in degrees. And also convert the cluster radius to degrees:

```python
ra_sc = Longitude(subclusters["RAJ2000"], unit=u.hourangle)
dec_sc = Latitude(subclusters["DEJ2000"], unit=u.deg)
subclusters["ra"] = ra_sc.to(u.deg)
subclusters["dec"] = dec_sc
subclusters["ra"].format = "{:.6f}"
subclusters["dec"].format = "{:.6f}"
```

```python
D_smc = 61.7*u.kpc
R_deg = (subclusters["Rad"].to(u.au).value / (D_smc.to(u.pc)).value) * u.arcsec#.to(u.deg)
R_deg
```

```python
subclusters["R"] = R_deg
subclusters
```

```python
from astropy.visualization.wcsaxes import SphericalCircle
```

```python
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

```

```python
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
```

```python
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
```

```python
dx, dy = 400, 400
```

```python
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
```

## Try and apply CMD cuts to isolate YSOs

Annoyingly, the catalog does not include the estimated age of the stars.  

But it looks like $V < 20 + 4 (V - I)$ would make a good cut.  This would select the upper main sequence plus the young stars above the lower main sequence (plus evolved stars, but there shouldn't be many of those).

```python
mcut = stars["F555W"] < 20.0 + 4*stars["V-I"]
```

```python
mcut.sum(), (~mcut).sum()
```

```python
ystars = stars[mcut][::-1]
ystars
```

```python
lmsstars = stars[~mcut]
```

```python
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
```

```python
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
```

```python
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
```

```python
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
```

## Better cuts in the CMD

Actually, it turns out that Gouliermis:2006h already have the regions defined for UMS, LMS, PMS, and RGB.  So we should use those. 

First, let's plot the CMD to look at it.

```python
mag0 = 22.1
mUMS = (stars["F555W"] < mag0) & (stars["F555W"] < 29 - 15 * stars["V-I"])
mRGB = (stars["F555W"] < mag0) & (stars["F555W"] >= 29 - 15 * stars["V-I"])
mLMS = (stars["F555W"] >= mag0) & (stars["F555W"] > 18.5 + 5.5 * stars["V-I"])
mPMS = (stars["F555W"] >= mag0) & (stars["F555W"] <= 18.5 + 5.5 * stars["V-I"])


stars['reg'] = 0.0
stars['reg'][mUMS] = 1.0
stars['reg'][mRGB] = 2.0
stars['reg'][mLMS] = 3.0
stars['reg'][mPMS] = 4.0


fig, ax = plt.subplots(figsize=(15, 15))
ax.scatter(
    stars["V-I"],
    stars["F555W"],
    marker=".",
    s=1,
    c=stars['reg'],
    alpha=0.15,
    cmap="hsv",
)
ax.set(
    xlim=[-1.5, 3.0],
    ylim=[27.0, 13.0],
)
```

```python
ystars = stars[mPMS | mUMS][::-1]
ystars
```

```python
lmsstars = stars[mLMS]
rgbstars = stars[mRGB]
```

```python

```

```python
dx, dy = 2000, 2000
fig = plt.figure(figsize=(8, 8))
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
    s=10.0/(ystars["F555W"] - 13.0),
    edgecolors="k",
    lw=0.05,
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
        linewidth=2,
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
fig.tight_layout();
fig.savefig("../figs/ngc-346-star-map-zoom-A.pdf");
```

Also indicate the coordinates of Source C from Table 1 of Rubio:2018f, which come from near-IR imaging.

```python
c0_fk5 = SkyCoord(
    "00 59 05.43",
    "-72 10 35.5",
    unit=(u.hourangle, u.deg),
    frame="fk5",
)
```

```python
c0_fk5
```

```python
c0 = c0_fk5.icrs
c0
```

```python
c0.ra.to_string(u.hourangle), c0.dec.to_string(u.deg)
```

Note that I have converted the coordinates from J2000 (as stated in Rubio:2018f) to ICRS for consistency with the Gaia frame.  But this makes a difference of only 0.01 arcsec, whereas the precision of the data is only 0.1 arcsec (based on significant figures given, who knows wjat the accuracy is).

```python
dx, dy = 500, 500
fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    imha, 
    vmin=0.2, vmax=1.0, 
    cmap="gray_r",
)

irpoint = ax.scatter(
    c0.ra.deg,
    c0.dec.deg,
    marker="+",
    s=200,
    c="w",
    linewidth=3,
    transform=ax.get_transform('world'), 
)

points = ax.scatter(
    ystars['ra'].data, 
    ystars['dec'].data, 
    c=ystars["V-I"].data,
    s=30.0/(ystars["F555W"] - 12.0),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    edgecolors="k",
    linewidth=0.1,
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=1.0)

for data in subclusters[[0, 1, 2, 6, 7]]:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='k', facecolor='none', 
        linestyle="dashed",
        linewidth=2,
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
fig.tight_layout()
fig.savefig("../figs/ngc-346-star-map-zoom-B.pdf")
```

```python

```

```python
dx, dy = 100, 100
fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.5, vmax=1.0, 
    cmap="gray_r",
)


irpoint = ax.scatter(
    c0.ra.deg,
    c0.dec.deg,
    marker="+",
    s=1000,
    c="w",
    linewidth=3,
    transform=ax.get_transform('world'), 
)

points = ax.scatter(
    ystars['ra'].data, 
    ystars['dec'].data, 
    c=ystars["V-I"].data,
    s=200.0/(ystars["F555W"] - 13.5),
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    edgecolors="k",
    linewidth=0.3,
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
fig.colorbar(points, ax=ax, label="$V - I$", shrink=1.0)

points2 = ax.scatter(
    lmsstars['ra'].data, 
    lmsstars['dec'].data, 
    edgecolors="w",
    facecolor="none",
    s=200.0/(lmsstars["F555W"] - 13.5),
    alpha=1.0,
    linewidth=1,
    transform=ax.get_transform('world'),    
)

points3 = ax.scatter(
    rgbstars['ra'].data, 
    rgbstars['dec'].data, 
    edgecolors="r",
    facecolor="none",
    s=200.0/(rgbstars["F555W"] - 13.5),
    linewidth=1,
    alpha=1.0,
    transform=ax.get_transform('world'),    
)

ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.tight_layout()
fig.savefig("../figs/ngc-346-star-map-zoom-C.pdf");
```

```python
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

irpoint = ax.scatter(
    c0.ra.deg,
    c0.dec.deg,
    marker="+",
    s=4*1000,
    c="r",
    linewidth=10,
    transform=ax.get_transform('world'), 
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
points3 = ax.scatter(
    rgbstars['ra'].data, 
    rgbstars['dec'].data, 
    edgecolors="r",
    facecolor="none",
    s=4*400.0/(rgbstars["F555W"] - 13.0),
    linewidth=3,
    alpha=1.0,
    transform=ax.get_transform('world'),    
)



ax.set(
    xlim=[x0 + 50 - dx, x0 + 50 + dx],
    ylim=[y0 - 40 - dy, y0 - 40 + dy],
)
fig.tight_layout();
```

So, it doesn't make much difference really, but it is good to do the cuts properly.  We do see quite a lot of the RGB-area sources, but these might not be RGB, but instead could be intermediate-mass young stars


## Look at CMD of the central cluster

We can take a cut in RA and Dec. 

```python
x0, y0
```

```python
w.pixel_to_world_values(x0, y0)
```

```python
dx, dy = 100, 100
ra1, dec1 = w.pixel_to_world_values(x0 + dx, y0 - dy)
ra2, dec2 = w.pixel_to_world_values(x0 - dx, y0 + dy)

mbox = (
    (stars["ra"] >= ra1)
    & (stars["ra"] <= ra2)
    & (stars["dec"] >= dec1)
    & (stars["dec"] <= dec2)
)
boxstars = stars[mbox]
boxstars
```

```python
fig, ax = plt.subplots(figsize=(15, 15))
ax.scatter(
    stars["V-I"],
    stars["F555W"],
    marker=".",
    s=1,
    c="k",
    alpha=0.15,
    cmap="hsv",
)

ax.scatter(
    boxstars["V-I"],
    boxstars["F555W"],
    marker="o",
    s=20,
    c="r",
    alpha=0.5,
)


ax.set(
    xlim=[-1.5, 3.0],
    ylim=[27.0, 13.0],
)
```

```python
mYSO = (
    (stars["F555W"] < 17) 
    & (stars["V-I"] > 0.5)
    & (stars["V-I"] < 1.0)
)
stars[mYSO]
```

Looks like ID152, which is MPG454 is our star!


And find the separation of each other star from infrared source C, which we calculated earlier.

```python
boxstars["sep"] = SkyCoord(
    boxstars["ra"], boxstars["dec"], unit=u.deg
).separation(c0).to(u.arcsec)
boxstars["sep"].format = "{:.2f}"
```

Now put them in order of that separation:

```python
boxstars.sort("sep")
boxstars
```

This shows that the "other" star is ID168, which doesn't have an MPG number.  Although MPG did not have the resolution to distinguish between the two sources, which are separated by 0.3 arcsec. 

Note that the identification of these two stars (152 and 168) was already mentioned in Sec 2.3, para 1 of Rubio:2018f

The NE star (we could call it MPG454-NE) is within 0.02 arcsec of the IR source, which is coincident within the astrometric precision.

Now, we can redo the CMD, but showing the separation with symbol color.

```python
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(
    stars["V-I"],
    stars["F555W"],
    marker=".",
    edgecolors="none",
    s=1,
    c="k",
    alpha=0.15,
    cmap="hsv",
)

points = ax.scatter(
    boxstars["V-I"],
    boxstars["F555W"],
    marker="o",
    s=50/(0.15 + boxstars["sep"]),
    c=boxstars["sep"],
    cmap="plasma",
    vmin=0.0,
    vmax=8.0,
    alpha=1.0,
)
fig.colorbar(
    points, ax=ax,
    orientation="horizontal",
    label="Separation from IR source, arcsec"
)

ax.set(
    xlim=[-1.5, 3.0],
    ylim=[27.0, 13.0],
    xlabel="$V - I$",
    ylabel="$V$",
)
fig.tight_layout()
fig.savefig("../figs/ngc-346-cmd-zoom.pdf");
```

<!-- #region -->
It is clear that the the closer stars to the IR source tend to be younger – they are further from the main sequence line.

Note that the Sun would have $V = 24$ at the distance of the SMC.

The reddest stars on MS have $V - I = 1.4$, which corresponds to $T_\mathrm{eff} \approx 4800$ K (see Fig 8 of Da-Rio:2010a). 

The reddest PMS stars in this area have $V - I = 2.2$, which corresponds to $T_\mathrm{eff} \approx 3800$ K.  



One of the young-looking stars has $V = 24$, $V - I = 2$. So, a bit more than $1\,L_\odot$ (because BC is bigger for cooler stars), and probably $T_\mathrm{eff} \approx 4100$ K. This could be a 0.6 Msun star that is less than 1 Myr old. 
<!-- #endregion -->

## Look at the brightest stars

At the separation column in to the full table:

```python
stars["sep"] = SkyCoord(
    stars["ra"], stars["dec"], unit=u.deg
).separation(c0).to(u.arcsec)
stars["sep"].format = "{:.2f}"
```

```python
mbright = stars["F555W"] < 17.0
mbright = mbright & mUMS
mbright.sum(), (~mbright).sum()
```

```python
mblue = stars["V-I"] < 0.0
(mblue & mbright).sum()
```

```python
obstars = stars[mbright & mblue]
obstars
```

```python
obstars[obstars["sep"] < 40.0]
```

## Look at stars that are similar to SSN 14



```python
stars[stars["ID"] == 14]
```

Look at difference in V and V-I

```python
stars["diff SSN14"] = np.hypot(
    stars["F555W"] - 14.09,
    stars["V-I"] - (-0.228),
)
m14 = stars["diff SSN14"] <= 0.5
```

```python
stars[m14]
```

Kind of unsurprising that they are the consecurtive ones on the list!


## Plot the positions of the 150 brightest stars



```python
dx, dy = 3000, 3000
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-1.0, vmax=0.3, 
    cmap="gray_r",
)
bbstars = stars[:150]
points = ax.scatter(
    bbstars['ra'].data, 
    bbstars['dec'].data, 
    c=bbstars["V-I"].data,
    s=10**(4.5-0.2*bbstars["F555W"]),
    edgecolors="k",
    lw=0.2,
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)

for i, star in enumerate(bbstars):
    xpos = ["left", "right"][i % 2]
    ypos = ["top", "center", "bottom"][i % 3]
    tcol = "cmbrg"[i % 5]
    ax.text(
        star["ra"], star["dec"],
        f"{star['ID']:03d}", 
        fontsize="xx-small",
        ha=xpos,
        va=ypos,
        color=tcol,
        fontweight="black",
        clip_on=False,
        transform=ax.get_transform('world'),   
    )

for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=2,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
#    ax.text(
#        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
#        ha="center", va="center", color="w",
#        fontweight="black",
#        clip_on=True,
#        transform=ax.get_transform("world"),
#    )



ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.tight_layout();
fig.savefig("../figs/ngc-346-star-map-bright-A.pdf");
```

```python
bbstars
```

```python
dx, dy = 500, 500
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.7, vmax=0.6, 
    cmap="gray_r",
)
bbstars = stars[:150]
points = ax.scatter(
    bbstars['ra'].data, 
    bbstars['dec'].data, 
    c=bbstars["V-I"].data,
    s=10**(5.0-0.2*bbstars["F555W"]),
    edgecolors="k",
    lw=0.2,
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)

for i, star in enumerate(bbstars):
    xpos = ["left", "right"][i % 2]
    ypos = ["top", "center", "bottom"][i % 3]
    tcol = "cmbrg"[i % 5]
    ax.text(
        star["ra"], star["dec"],
        f"{star['ID']:03d}", 
        fontsize="small",
        ha=xpos,
        va=ypos,
        color=tcol,
        fontweight="black",
        clip_on=True,
        transform=ax.get_transform('world'),   
    )

for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=2,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
#    ax.text(
#        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
#        ha="center", va="center", color="w",
#        fontweight="black",
#        clip_on=True,
#        transform=ax.get_transform("world"),
#    )



ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.tight_layout();
fig.savefig("../figs/ngc-346-star-map-bright-B.pdf");
```

```python
dx, dy = 100, 100
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.4, vmax=1.0, 
    cmap="gray_r",
)
bbstars = stars[:300]
points = ax.scatter(
    bbstars['ra'].data, 
    bbstars['dec'].data, 
    c=bbstars["V-I"].data,
    s=10**(5.5-0.2*bbstars["F555W"]),
    edgecolors="k",
    lw=0.2,
    vmin=-0.7 + 0.3,
    vmax=1.0 + 0.3,
    cmap="rainbow",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)

for i, star in enumerate(bbstars):
    xpos = ["left", "right"][i % 2]
    ypos = ["top", "center", "bottom"][i % 3]
    tcol = "cmbrg"[i % 5]
    ax.text(
        star["ra"], star["dec"],
        f"{star['ID']:03d}", 
        fontsize="small",
        ha=xpos,
        va=ypos,
        color=tcol,
        fontweight="black",
        clip_on=True,
        transform=ax.get_transform('world'),   
    )

for data in subclusters:
    circ = SphericalCircle(
        (data["ra"]*u.deg, data["dec"]*u.deg), 
        data["R"]*u.arcsec,
        edgecolor='w', facecolor='none', 
        linestyle="dotted",
        linewidth=2,
        transform=ax.get_transform("world"),
    )
    ax.add_patch(circ)
#    ax.text(
#        data["ra"], data["dec"], data["__SSN2007_"].split()[-1],
#        ha="center", va="center", color="w",
#        fontweight="black",
#        clip_on=True,
#        transform=ax.get_transform("world"),
#    )



ax.set(
    xlim=[x0 - dx, x0 + dx],
    ylim=[y0 - dy, y0 + dy],
)
fig.tight_layout();
fig.savefig("../figs/ngc-346-star-map-bright-C.pdf");
```

```python

```
