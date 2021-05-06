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

# Extracting subregions from a MUSE cube

I want to learn how to extract regions along the spatial and wavelength axes and manipulate them to do things continuum subtraction.

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
```

```python
datapath = Path("/Users/will/Work/Muse-Hii-Data/SMC-NGC-346/")
fitsfilepath = datapath / "ADP.2017-10-16T11_04_19.247.fits"
cube = Cube(str(fitsfilepath))
cube.info()
```

## Extracting wavelength slices


We make a summed spectrum of the entire field again, just like in the previous notebook.

```python
sp0 = cube.sum(axis=(1, 2))
```

Now, we take a look at the region around the H alpha line. The `Spectrum.subspec` method allows us to use Å instead of pixels.

```python
fig, ax = plt.subplots(figsize=(10, 4))
sp0.subspec(6520.0, 6620.0).plot()
ax.set(yscale="log")
```

This shows that a good range for extracting Ha would be 6560 to 6575 Å.  So, we just do a simple-minded sum of that range. For some reason, the equivalent of `Spectrum.subspec` for a `Cube` has a different name, `Cube.select_lambda`:

```python
im_ha = cube.select_lambda(6560.0, 6575.0).sum(axis=0)
```

And we have a look at it:

```python
fig, ax = plt.subplots(figsize=(10, 10))
im_ha.plot(use_wcs=True, cmap="gray_r", scale="log", colorbar="v")
```

So, that looks good, but it includes the continuum.  A very simple-minded way of removing the continuum would just be to take an average value in line-free wavelength ranges.

We can use 6555-60 and 6575-80:

```python
im_blue_cont = cube.select_lambda(6555.0, 6560.0).mean(axis=0)
im_red_cont = cube.select_lambda(6575.0, 6580.0).mean(axis=0)
im_mean_cont = (im_blue_cont + im_red_cont) / 2.0
```

So we have a map of the average continuum from the red and blue sides of the line, which we will now visualize:

```python
fig, ax = plt.subplots(figsize=(10, 10))
im_mean_cont.plot(use_wcs=True, cmap="gray_r", scale="log", colorbar="v")
```

Now I subtract the average continuum level before summing again to find the BG-subtracted line map. We have to do it this way because `im_mean_cont` still has the per-Å units since it is a mean, so we can't just subtract it from `im_ha`.

```python
im_ha_bgsub = (cube.select_lambda(6560.0, 6575.0) - im_mean_cont).sum(axis=0)
```

```python
fig, ax = plt.subplots(figsize=(10, 10))
im_ha_bgsub.plot(
    use_wcs=True,
    vmin=-2.0e4,
    vmax=1.0e5,
    cmap="gray_r",
    scale="linear",
    colorbar="v",
)
```

So. now we can see beautiful details of the Hα emission. There are chains of globules and elephant-trunk structures that cross the nebula. These are physically quite large – about 2 arcsec, which is about 0.5 parsec.  The entire field is about 15 pc square.

The white dots are presumably from the underlying photospheric absorption from some of the stars.  Other stars are seen in emission.


Strangely, though, the brightnesses go negative

```python
im_ha.data.min()
```

```python
im_mean_cont.data.min()
```

```python
im_ha_ew = (im_ha_bgsub - im_ha.data.min()) / (im_mean_cont - im_mean_cont.data.min())
```

```python
fig, ax = plt.subplots(figsize=(10, 10))
im_ha_ew.plot(
    vmin=-10.0,
    vmax=10.0,
    use_wcs=True,
    cmap="gray_r",
    scale="linear",
    colorbar="v",
)
```

 ## Extracting spatial regions


We will go back to the ha+continuum image and try and plot it without the stars

```python
savemask = im_ha.mask
savemask
```

We mask out the bright continuum sources using `mask_selection`

```python
im_ha.mask_selection(np.where(im_mean_cont.data > 1000.0))
```

```python
fig, ax = plt.subplots(figsize=(10, 10))
im_ha.plot(use_wcs=True, cmap="gray_r", scale="log", colorbar="v")
```

What this uses is an array of pixels, such as provided by `np.where`:

```python
bright_cont_selection = np.where(im_mean_cont.data > 1000.0)
bright_cont_selection
```

```python
savecubemask = cube.mask
```

```python
cube.mask_selection(bright_cont_selection)
```

```python

```
