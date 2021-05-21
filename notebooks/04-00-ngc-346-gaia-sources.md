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

# Gaia astrometry of stars in NGC 346

I want to get good positions for all the stars so I can align the HST and the MUSE images

```python
import pandas as pd
from pathlib import Path
```

## Load the Gaia sources

```python
datapath = Path("../data")
```

```python
df = pd.read_csv(datapath / "1621565655827O-result.csv")
```

```python
df
```

## Plot all the Gaia sources

```python
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()
sns.set_context("talk")
```

I set the color scale to be the magnitude

```python
fig, ax = plt.subplots(figsize=(10, 8))
scat = ax.scatter(
    x="ra", y="dec", c="phot_g_mean_mag", 
    data=df,
    s=20,
)
cb = fig.colorbar(scat, ax=ax)
```

So you can see the mass segregation straight away - the darker ones (brighter) tend to be more concentrated towards the center.  It looks like Gaia must have missed quite a lot of them though.

*Note that the RA axis is the wrong way round* 


## Load an HST image

We will start off with the UV one.

```python
import numpy as np
```

```python
from astropy.io import fits
from astropy.wcs import WCS
```

```python
bigdatapath = Path("../big-data")
```

```python
hdulist = fits.open(
    bigdatapath / 
    "HST-NGC346/hst_10248_03_acs_hrc_f330w_f220w/hst_10248_03_acs_hrc_f330w_f220w_sci.fits")
```

```python
hdulist.info()
```

## Plot the image in celestial coordinates with Gais sources overlaid

```python
w = WCS(hdulist[0].header).celestial
```

```python
imdata = hdulist[0].data[0, :, :]
```

```python
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
```

This shows that the HST coordinates are cbot quite right – there is an offset from the Gaia coordinates.

We will zoom in on the central cluster to have a closer look.

```python
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
```

## Adjust the WCS to align HST with Gaia


Now, try giving an offset to fix this.  We could either change `CRPIX` or the `CRVAL`, but the first seems simpler to reason about.  We want to move the HST stars to the right and down, so do the opposite to the reference pixel:

```python
ww = WCS(hdulist[0].header).celestial
ww.wcs.crpix
```

```python
ww.wcs.crpix -= np.array([19, -7])
```

```python
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
```

That looks pretty good. 


## Repeat for the Hα image

```python
acspath = bigdatapath / "HST-NGC346/hst_10248_a3_acs_wfc_f658n"
hdulist_acs = fits.open(
    acspath / "hst_10248_a3_acs_wfc_f658n_drz.fits"
)
```

```python
hdulist_acs.info()
```

```python
wacs = WCS(hdulist_acs["SCI"].header)
imdata_acs = hdulist_acs["SCI"].data.astype("float")
```

```python
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=wacs)
ax.imshow(imdata_acs, vmin=0.0, vmax=1.0, cmap="gray_r")
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
```

Zoom in some more and switch to a logarithmic brightness scaling:

```python
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=wacs)
ax.imshow(
    np.log10(imdata_acs), 
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
```

This time, the whift is up and to the right. Let's try and fix it:

```python
wacs_fix = WCS(hdulist_acs["SCI"].header)
wacs_fix.wcs.crpix += np.array([6, 7.5])
```

```python
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=wacs_fix)
ax.imshow(
    np.log10(imdata_acs), 
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
```

```python

```

```python

```

### Load the DAOPHOT source list

The HST data came with some data tables that seem to

```python
from astropy.io import ascii
```

```python
stars = ascii.read(
    str(acspath / "hst_10248_a3_acs_wfc_f658n_daophot_trm.cat"),
    format="commented_header",
    header_start=56,
)
```

```python
stars
```

```python
ascii.read?
```

```python

```
