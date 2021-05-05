# -*- coding: utf-8 -*-
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

# # Extracting subregions from a MUSE cube
#
# I want to learn how to extract regions along the spatial and wavelength axes and manipulate them to do things continuum subtraction.

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube

datapath = Path("/Users/will/Work/Muse-Hii-Data/SMC-NGC-346/")
fitsfilepath = datapath / "ADP.2017-10-16T11_04_19.247.fits"
cube = Cube(str(fitsfilepath))
cube.info()

# ## Extracting wavelength slices

# We make a summed spectrum of the entire field again, just like in the previous notebook.

sp0 = cube.sum(axis=(1, 2))

# Now, we take a look at the region around the H alpha line. The `Spectrum.subspec` method allows us to use Å instead of pixels. 

fig, ax = plt.subplots(figsize=(10, 4))
sp0.subspec(6520.0, 6620.0).plot()
ax.set(yscale="log");

# This shows that a good range for extracting Ha would be 6560 to 6575 Å.  So, we just do a simple-minded sum of that range. For some reason, the equivalent of `Spectrum.subspec` for a `Cube` has a different name, `Cube.select_lambda`: 

im_ha = cube.select_lambda(6560.0, 6575.0).sum(axis=0)

# And we have a look at it:

fig, ax = plt.subplots(figsize=(10, 10))
im_ha.plot(use_wcs=True, cmap="gray_r", scale="log", colorbar="v");

# So, that looks good, but it includes the continuum.  A very simple-minded way of removing the continuum would just be to take an average value in line-free wavelength ranges.  
#
# We can use 6555-60 and 6575-80:

im_blue_cont = cube.select_lambda(6555.0, 6560.0).mean(axis=0)
im_red_cont = cube.select_lambda(6575.0, 6580.0).mean(axis=0)
im_mean_cont = (im_blue_cont + im_red_cont) / 2.0

# So we have a map of the average continuum from the red and blue sides of the line, which we will now visualize:

fig, ax = plt.subplots(figsize=(10, 10))
im_mean_cont.plot(use_wcs=True, cmap="gray_r", scale="log", colorbar="v");

# Now I subtract the average continuum level before summing again to find the BG-subtracted line map. We have to do it this way because `im_mean_cont` still has the per-Å units since it is a mean, so we can't just subtract it from `im_ha`. 

im_ha_bgsub = (
    cube.select_lambda(6560.0, 6575.0) - im_mean_cont
).sum(axis=0)

fig, ax = plt.subplots(figsize=(10, 10))
im_ha_bgsub.plot(
    use_wcs=True, 
    vmin=-2.0e4,
    vmax=1.0e5,
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);

# So. now we can see beautiful details of the Hα emission. There are chains of globules and elephant-trunk structures that cross the nebula. These are physically quite large – about 2 arcsec, which is about 0.5 parsec.  The entire field is about 15 pc square.
#
# The white dots are presumably from the underlying photospheric absorption from some of the stars.  Other stars are seen in emission.

# Strangely, though, the brightnesses go negative 

im_ha.data.min()

im_mean_cont.data.min()

im_ha_ew = (im_ha_bgsub - im_ha.data.min()) / (im_mean_cont - im_mean_cont.data.min())

fig, ax = plt.subplots(figsize=(10, 10))
im_ha_ew.plot(
    vmin=-10.0, vmax=10.0,
    use_wcs=True, 
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);

fig, ax = plt.subplots(figsize=(10, 4))
sp0.subspec(6700.0, 6750.0).plot()
ax.set(yscale="log");

im_c1 = cube.select_lambda(6710.0, 6715.0).mean(axis=0)
im_c2 = cube.select_lambda(6725.0, 6730.0).mean(axis=0)
im_c6716 = (im_c1 + im_c2) / 2.0
im_c3 = cube.select_lambda(6740.0, 6745.0).mean(axis=0)
im_c6731 = (im_c3 + im_c2) / 2.0

im_sii16_bgsub = (
    cube.select_lambda(6715.0, 6725.0) - im_c6716
).sum(axis=0)
im_sii31_bgsub = (
    cube.select_lambda(6730.0, 6740.0) - im_c6731
).sum(axis=0)

im_sii31_bgsub -= im_sii31_bgsub.data.min()
im_sii16_bgsub -= im_sii16_bgsub.data.min()

fig, ax = plt.subplots(figsize=(10, 10))
(im_sii16_bgsub + im_sii31_bgsub).plot(
    vmin=0.0,
    vmax=20000.0,
    use_wcs=True, 
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);

fig, ax = plt.subplots(figsize=(10, 10))
(im_sii31_bgsub / im_sii16_bgsub).plot(
    vmin=0.3,
    vmax=1.4,
    use_wcs=True, 
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);

fig, ax = plt.subplots(figsize=(10, 4))
sp0.subspec(5700.0, 5800.0).plot()
ax.set(yscale="log");

fig, ax = plt.subplots(figsize=(10, 4))
sp0.subspec(6250.0, 6350.0).plot()
ax.set(yscale="log");

fig, ax = plt.subplots(figsize=(10, 4))
sp0.subspec(9050, 9100).plot()
ax.set(yscale="log");

im_c6312 = cube.select_lambda(6320.0, 6340.0).mean(axis=0)
im_c2 = cube.select_lambda(9065.0, 9070.0).mean(axis=0)
im_c3 = cube.select_lambda(9080.0, 9085.0).mean(axis=0)
im_c9069 = (im_c3 + im_c2) / 2.0

im_siii6312_bgsub = (
    cube.select_lambda(6310.0, 6320.0) - im_c6312
).sum(axis=0)
im_siii9069_bgsub = (
    cube.select_lambda(9070.0, 9080.0) - im_c9069
).sum(axis=0)

im_siii6312_bgsub += 150.0
im_siii9069_bgsub -= im_siii9069_bgsub.data.min()

fig, ax = plt.subplots(figsize=(10, 10))
im_siii9069_bgsub.plot(
    vmin=0.0,
    vmax=10000.0,
    use_wcs=True, 
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);

fig, ax = plt.subplots(figsize=(10, 10))
(im_siii6312_bgsub / im_siii9069_bgsub).plot(
    vmin=0.0,
    vmax=0.1,
    use_wcs=True, 
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);

fig, ax = plt.subplots(figsize=(10, 10))
(im_siii9069_bgsub / im_sii31_bgsub).plot(
    vmin=0.0,
    vmax=12.0,
    use_wcs=True, 
    cmap="gray_r", 
    scale="linear", 
    colorbar="v",
);


