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


*Note that the `Cube.get_image()` function can do the continuum subtraction semi-automatically, which may be more convenient in some instances.*


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


### Masking of 2D images

We will go back to the ha+continuum image and try and plot it without the stars

```python
savemask = im_ha.mask.copy()
savemask
```

We mask out the bright continuum sources using `mask_selection`

```python
im_ha.mask_selection(np.where(im_mean_cont.data > 1000.0))
```

```python
fig, ax = plt.subplots(figsize=(10, 10))
im_ha.plot(use_wcs=True, cmap="viridis", scale="log", colorbar="v")
```

What this uses is an array of pixels, such as provided by `np.where`:

```python
bright_cont_selection = np.where(im_mean_cont.data > 1000.0)
bright_cont_selection
```

### Masking of a 3D cube (small subcube version)

My first attempts at doing this did not work as I expected, so I am going to extract a small sub-cube that will be easy to look at.  


There are at least three ways that look like they might be used to make a sub-cube:

1. `Cube.subcube(center=..., size=..., lbda=..., ...)` cuts out a rectangle.  It returns a view so we would have to do `.copy()`
2. `Cube.truncate(coord=[w1, x1, y1, w2, x2, y2], ...)` does return a copy, and with a slightly different interface
3. `Cube[k1:k2, j1:j2, i1:i2]` is the simplest approach using pixel slice notation, and should return a copy

We will try the third way first, but just for the celestial coordinates, then use `.select_lambda()` to extract a small section of the wavelength axis:

```python
cc = cube[:, 5:15, 130:140].select_lambda(6554.0, 6581.0)
```

```python
cc.info()
```

We have a 10x10 pixel postage stamp, with 22 pixels in wavelength. Some of the pixels are from the image border and are masked, as we can see here:

```python
cc.sum(axis=0).plot(use_wcs=True)
```

I chose this region because it contains a star that we can try and mask out.  Here is the spectrum averaged over all pixels and plotted on a linear scale:

```python
fig, ax = plt.subplots(figsize=(10, 4))
cc.mean(axis=(1, 2)).plot()
ax.set(ylim=[0, None]);
```

So, we have roughly equal parts continuum and line emission.  

Now, we look at the mask (steps by multiple pixels so we can see it better): 

```python
cc.mask[::5, ::2, ::2]
```

We see a stack of wavelength images. Each has y-axis inverted, so the `True` values in each top row are the bottom line of pixels in each image.  


I am going to try and mask the pixels that have relatively strong continuum, and then recalculate the average spectrum.  With luck, that will make the line stand out more.  First, save the original mask (I use copy here because otherwise they are just different views of the same object - we should also use copy to put the original mask back if we ever want to). 

```python
cc_mask_orig = cc.mask
```

```python
id(cc_mask_orig), id(cc.mask)
```

```python
cc_mask_orig = cc.mask.copy()
```

```python
id(cc_mask_orig), id(cc.mask)
```

#### Masking via simple slicing on the celestial pixel axes

First, I try doing it by hand with pixel slicing of the mask.  I will attempt to mask out the bottom right of the cube (8x8 pixels), so we just leave a 2 pixel L-shaped strip:

```python
cc.mask[:, :8, -8:] = True
```

```python
cc.sum(axis=0).plot(use_wcs=True)
```

So that looks OK. Now I check if it has had any effect on the average spectrum:

```python
fig, ax = plt.subplots(figsize=(10, 4))
cc.mean(axis=(1, 2)).plot()
ax.set(ylim=[0, None]);
```

Yes, the continuum level is now much lower, although this is at a cost of removing most of the pixels. 

#### Masking via logical operation to make a boolean array

Now, we will see if we can use a logical condition to set the mask.  But first, restore the original

```python
cc.mask = cc_mask_orig.copy()
```

We can make a 2D mask, which is true for all pixels where the continuum intensity is more than half the average intensity in the passband:

```python
my_mask = cc[:6, :, :].mean(axis=0).data > 0.5*cc.mean(axis=0).data
```

Note that we had to extract the `.data` before comparing, otherwise we get something very strange.  Even so, the mask we get back is itself a masked array – so it is a masked mask!

```python
my_mask[::2, ::2]
```

The `.data` part has a logical array where our condition is true.  The `.mask` part is another logical array where the cube's mask is true.  We can simplify life by just ORing them together to make a simple array:

```python
my_mask = my_mask.mask | my_mask.data
my_mask[::2, ::2]
```

In order to apply this mask to the cube, we need to make it 3-dimensional:

```python
my_mask_3d = np.repeat(my_mask[None, :, :], cc.shape[0], axis=0)
```

```python
my_mask_3d[::5, ::2, ::2]
```

Now directly set the cube's mask to be this boolean array.

```python
cc.mask = my_mask_3d.copy()
```

Have a look at the result. First collapse the spectral axis to make an image:

```python
cc.sum(axis=0).plot(use_wcs=True)
```

Looks good.  We see that the original border is still masked, plus a roughly circular patch around the star.

Now collapse the celestial axes to get a spectrum, which should be of the unmasked pixels only:

```python
fig, ax = plt.subplots(figsize=(10, 4))
cc.mean(axis=(1, 2)).plot()
ax.set(ylim=[0, None]);
```

Yep. that looks fine.


### Masking of a 3D cube (big cube version)

Now that I have a method that seems to work, I will go back to the full cube and try things there.  Although, I will restrict it to the an extended wavelength range around H alpha. We will use the range 6200 to 6800 Å, which is the same as Fig 3 of my Raman paper.

```python
hacube = cube.select_lambda(6200.0, 6800.0)
hacube.info()
```

Save a copy of the original mask again

```python
hacube_mask_orig = hacube.mask.copy()
```

```python
fig, ax = plt.subplots(figsize=(10, 10))
hacube.sum(axis=0).plot(use_wcs=True, cmap="magma", scale="log", colorbar="v")
```

```python
fig, ax = plt.subplots(figsize=(10, 4))
hacube.mean(axis=(1, 2)).plot()
ax.set(ylim=[0, 400]);
```

There really is a lot of continuum here.  It would be good to get rid of some of it so we can see the lines better!

Make a 3D mask for `hacube` where each image plane is the mask of the strong continuum sources from above.

```python
my_mask_3d = np.repeat(
    im_ha.mask[None, :, :], 
    hacube.shape[0], 
    axis=0,
)
```

Have a look at a sparse sample of the cube mask:

```python
my_mask_3d[:2, ::20, ::50]
```

Looks OK. Now we apply it to the cube.  We need to make sure we are also masking out the voxels that were originally masked, since some of them have invalid data, which can cause problems with the spectra.

```python
hacube.mask = my_mask_3d.copy() | hacube_mask_orig
```

And redo the summed image and average spectrum:

```python
fig, ax = plt.subplots(figsize=(10, 10))
hacube.sum(axis=0).plot(use_wcs=True, cmap="magma", scale="log", colorbar="v")
```

```python
fig, ax = plt.subplots(figsize=(10, 4))
hacube.mean(axis=(1, 2)).plot()
ax.set(ylim=[0, 400]);
```

```python
#hacube.mask = (my_mask_3d) | hacube_mask_orig
hacube.mask = hacube_mask_orig.copy()
fig, ax = plt.subplots(figsize=(10, 4))
spec_ha = hacube[:, 240:250, 230:250].mean(axis=(1, 2))
spec2 = spec_ha.subspec(6400, 6650)
spec2.mask_region(6500, 6600)
cont = spec2.poly_spec(2)
spec_ha.unmask()
spec_ha.plot()
cont.plot(color="m")
ax.axvline(6578,c="r", lw=0.3)
ax.axvline(6583, c="r", lw=0.3)
ax.axvline(6591,c="r", lw=0.3)
ax.axvline(6600, c="r", lw=0.3)
ax.axvline(6620, c="r", lw=0.3)

ax.axvline(6540, c="c", lw=0.3)
ax.axvline(6549, c="c", lw=0.3)



ax.set(ylim=[0, 100]);
```

```python
cont6600 = cube.select_lambda(6600.0, 6620.0).mean(axis=0)
wing = (cube.select_lambda(6591.0, 6600.0) -  cont6600).sum(axis=0)
fig, ax = plt.subplots(figsize=(10, 10))
wing.mask = wing.mask | (cont6600.data > 70.0)
wing.plot(
    use_wcs=True, 
    vmin=-10, vmax=50,
    cmap="viridis", 
    scale="linear", 
    colorbar="v",
)
```

```python
inwing = (cube.select_lambda(6578.0, 6583.0) -  cont6600).sum(axis=0)
fig, ax = plt.subplots(figsize=(10, 10))
inwing.mask = inwing.mask | (cont6600.data > 70.0)
inwing.plot(
    use_wcs=True, 
    vmin=-10, vmax=50,
    cmap="viridis", 
    scale="linear", 
    colorbar="v",
)
```

```python
hacube.mask = hacube_mask_orig.copy() | wing.mask
fig, ax = plt.subplots(figsize=(10, 10))
hacube.sum(axis=0).plot()
```

```python
fig, ax = plt.subplots(figsize=(10, 4))
hacube[:, 100:140, 15:40].mean(axis=(1, 2)).plot()
hacube[:, 180:240, 240:300].mean(axis=(1, 2)).plot(c="r")
hacube[:, 10:50, 100:150].mean(axis=(1, 2)).plot(c="m")
hacube[:, 10:100, 200:300].mean(axis=(1, 2)).plot(c="g")

ax.set(ylim=[0.0, 50])
```

```python
hacube.mask = hacube_mask_orig.copy()
fig, ax = plt.subplots(figsize=(10, 4))
hacube[:, -30:, 50:70].mean(axis=(1, 2)).plot()
ax.set(ylim=[-100, 20])
```

```python
fig, ax = plt.subplots(figsize=(10, 4))
cube[:, -30:, 50:70].mean(axis=(1, 2)).plot()
ax.set(ylim=[-1000, 100])
```

```python

```
