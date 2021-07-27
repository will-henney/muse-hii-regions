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

```python
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
```

Use the mosaiced image from Castro to define the WCS for reprojection.

```python
mosaic_ha_hdulist = fits.open("../data/MUSE_R136toWill/GAUS_Ha6562.8_060_Will.fits")
mosaic_ha_hdulist.info()
```

```python
mosaic_hdr = mosaic_ha_hdulist[1].header
```

Note that the header still refers to a 3rd dimension, even though this is just an image.  We need to remove all the keywords that mention dimension 3 so that the reprojection will work.

```python
del mosaic_hdr["*3"]
del mosaic_hdr["CD3_*"]
```

Test with the Hβ image.  First we make an HDU of each field for the full mosaic.

```python
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
```

```python
pieces
```

```python
fig, ax = plt.subplots(
    figsize=(12, 12),
    subplot_kw=dict(projection=WCS(mosaic_hdr)),
)
ax.imshow(pieces["A"], vmin=0, vmax=1e6)
ax.imshow(pieces["B"], vmin=0, vmax=1e6)
ax.imshow(pieces["C"], vmin=0, vmax=1e6)
ax.imshow(pieces["D"], vmin=0, vmax=1e6)
```

```python
Image??
```

```python
field = "A"
infile = f"../data/lmc-30dor-{field}-hi-4861-bin01-sum.fits"
hdu = fits.open(infile)["DATA"]
del hdu.header["*3"]
del hdu.header["CD3*"]
WCS(hdu.header)

```

```python
WCS(mosaic_hdr)
```

```python
infiles = Path("../data").glob("lmc-30dor-[ABCD]-*bin01-sum.fits")
```

```python
list(infiles)
```

```python

```
