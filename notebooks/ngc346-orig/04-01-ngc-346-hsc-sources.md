---
jupyter:
  jupytext:
    encoding: '# -*- coding: utf-8 -*-'
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

# NGC 346 stars from Hubble Source Catalog

This is based off the notebook example at https://archive.stsci.edu/hst/hsc/help/api/hscv3_smc_api.html

But I have taken the HSC API routines and put them in a file.  I had to do a `conda install requests`. 

```python
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy.table import Table

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()
sns.set_context("talk")

sys.path.append("../lib")
import hubble_source_catalog as hsc
```

Set up the target:

```python
target = 'NGC 346'
ra, dec = hsc.resolve(target)
print(target, ra, dec)
```


And do the actual search:

```python
import time
```

```python
columns = """MatchID,MatchRA,MatchDec,CI,A_F555W,A_F814W""".split(",")
columns = [x.strip() for x in columns]
columns = [x for x in columns if x and not x.startswith('#')]

# select objects with at least one ACS F555W and ACS F814W measurement
# and with concentration index 0.9 < CI < 1.6, consistent with point sources
# search a small 3x3 arcmin box in RA and Dec centered on the SMC
ddec = 3.0 / 60.0
dra = ddec / np.cos(np.radians(dec))
constraints = {'A_F555W_N.gte': 1, 'A_F814W_N.gte': 1, 'CI.gt':0.5, 'CI.lt':1.6,
               'MatchDec.gt': dec - ddec, 'MatchDec.lt': dec + ddec,
               'MatchRA.gt': ra - dra, 'MatchRA.lt': ra + dra}

# do a search with a large number of rows allowed
t0 = time.time()
tab = ascii.read(
    hsc.hscsearch(
        table="summary",
        release='v3',
        columns=columns,
        verbose=True,
        pagesize=2000000,
        **constraints
    )
)
print("{:.1f} s: retrieved data and converted to {}-row astropy table".format(time.time()-t0, len(tab)))

# compute color column and select for objects in more limited color range
tab['V-I'] = tab['A_F555W'] - tab['A_F814W']
tab = tab[(tab['V-I'] < 2.5) & (tab['V-I'] > -2.5)]
print("{:.1f} s: selected {} objects with -2.5 < V-I < 2.5".format(time.time()-t0, len(tab)))

# clean up the output format
tab['A_F555W'].format = "{:.3f}"
tab['A_F814W'].format = "{:.3f}"
tab['V-I'].format = "{:.3f}"
tab['CI'].format = "{:.3f}"
tab['MatchRA'].format = "{:.6f}"
tab['MatchDec'].format = "{:.6f}"

tab
```

```python
fig, ax = plt.subplots(figsize=(20, 20))

stars = ax.scatter(
    tab['MatchRA'], 
    tab['MatchDec'], 
    s=0.1,
    c=tab["V-I"].data,
    cmap="icefire_r",
    alpha=1.0,
    label=f'{len(tab)} HSC measurements'
)
ax.plot(ra, dec, 'rx', label=target, markersize=10)
fig.colorbar(stars, ax=ax)
ax.invert_xaxis()
ax.set_aspect(1.0/np.cos(np.radians(dec)))
ax.set(
    xlabel='RA [deg]',
    ylabel='Dec [deg]',
)
ax.legend()
```

Zoom in on cluster

```python
fig, ax = plt.subplots(figsize=(20, 20))

stars = ax.scatter(
    tab['MatchRA'], 
    tab['MatchDec'], 
    s=100.0/(tab["A_F555W"] - 15.1),
    c=tab["V-I"].data,
    vmin=-0.7,
    vmax=1.0,
    cmap="rainbow",
    alpha=1.0,
    label=f'{len(tab)} HSC measurements'
)
ax.plot(ra, dec, 'rx', label=target, markersize=10)
fig.colorbar(stars, ax=ax)
ax.invert_xaxis()
ax.set_aspect(1.0/np.cos(np.radians(dec)))
dec_range = 0.5/60.0
ra_range = dec_range / np.cos(np.radians(dec))
ax.set(
    xlim=[ra - ra_range, ra + ra_range],
    ylim=[dec - dec_range, dec + dec_range],
    xlabel='RA [deg]',
    ylabel='Dec [deg]',
)
ax.legend()
```

I am using hue to represent the $V - I$ color and symbol size to represent the $V$ magnitude.


## Combine with Hα map

```python
bigdatapath = Path("../big-data")
```

```python
hdu = fits.open(bigdatapath / "ngc346-hst-acs-f658n-wcsgaia.fits")["SCI"]
w = WCS(hdu.header)
imha = hdu.data
```

```python
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.7, vmax=0.7, 
    cmap="gray_r",
)
stars = ax.scatter(
    tab['MatchRA'], 
    tab['MatchDec'], 
    s=30.0/(tab["A_F555W"] - 15.1),
    c=np.array(tab["V-I"].data),
    #edgecolors="w",
    vmin=-0.7,
    vmax=1.0,
    cmap="rainbow",
    facecolor="none",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
ax.set(
    xlim=[2400, 2800],
    ylim=[2000, 2400],
);
```

## Add the Gaia sources too

```python
datapath = Path("../data")
dfgaia = pd.read_csv(datapath / "1621565655827O-result.csv")
```

```python
fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(1, 1, 1, projection=w)
ax.imshow(
    np.log10(imha), 
    vmin=-0.7, vmax=0.7, 
    cmap="gray_r",
)
stars = ax.scatter(
    tab['MatchRA'], 
    tab['MatchDec'], 
    s=30.0/(tab["A_F555W"] - 15.1),
    c=np.array(tab["V-I"].data),
    #edgecolors="w",
    vmin=-0.7,
    vmax=1.0,
    cmap="rainbow",
    facecolor="none",
    alpha=1.0,
    transform=ax.get_transform('world'),    
)
gstars = ax.scatter(
    x="ra", y="dec", 
    edgecolors="w", 
    data=dfgaia,
    s=100/(dfgaia["phot_g_mean_mag"] - 15), 
    alpha=1.0,
    facecolor="none",
    transform=ax.get_transform('world'),
)
ax.set(
    xlim=[2400, 2800],
    ylim=[2000, 2400],
);
```

So we can see that neither catalog is complete.

What about the bespoke studies, like Nota+ (2006) and Sabbi+ (2007)?

- [ ] ***Yes, there is a Vizier table for Sabbi, which I should get***

```python

```
