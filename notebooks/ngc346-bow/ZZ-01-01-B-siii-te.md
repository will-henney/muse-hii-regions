---
jupyter:
  jupytext:
    encoding: '# -*- coding: utf-8 -*-'
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.18.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Notebook B of line ratios from saved maps from the PZ cube

This is from the splitting up of ZZ-01-01-more-line-ratios.ipynb. It includes calculating the S III temperature and various ionization diagnostics with forbidden lines

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
import regions
import sys
import pandas as pd
import cmasher as cmr
import pyneb as pn

sns.set_context("talk")
sns.set_color_codes()
```
## Path to the root of this repo

```python
ROOT = Path.cwd().parent.parent
```

```python
mapspath = ROOT / "data/n346-bow-lines/refined"
```

```python
def get_image(lineid, combo="P-007", ext="SUM"):
    im = Image(str(mapspath / f"map-{lineid}-{combo}.fits"), ext=ext)
    im.data = im.data.astype(float)
    return im


def get_image_raw(lineid, variant="csub", channel="ABC"):
    p = ROOT / "data" / "n346-bow-lines" / f"maps-n346-PZ-2pass-{variant}-007"
    candidates = list(p.glob(f"**/*-{lineid}-{channel}.fits"))
    assert candidates
    imfile = candidates[0]
    im = Image(str(imfile))
    im.data = im.data.astype(float)
    return im
```


## Reddening law for SMC

```python
rc = pn.RedCorr()
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
rc.law = "F99"
```

## Calculate the [S III] temperature


```python


im6312 = get_image_raw("s-iii-6312-06")
im9069 = get_image_raw("s-iii-9068-90")
cont6312 = get_image_raw("s-iii-6312-06", variant="cont")

# The raw ratio:

fig, ax = plt.subplots(figsize=(12, 12))
(im6312 / im9069).plot(vmin=0.08, vmax=0.15, cmap="gray", colorbar="v")

# That looks more or less like the last time

imax = 15000
slope = 0.1
m = im9069.data < imax
m = m & (im9069.data > -100)
m = m & (im6312.data < slope * imax)
m = m & ~im9069.mask & ~im6312.mask
df = pd.DataFrame(
    {
        "9069": im9069.data[m].astype(float),
        "6312": im6312.data[m].astype(float),
    }
)
```


```python
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)

g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope * imax], "--", color="r")
g.fig.suptitle("Correlation between [S III] 9069 and 6312 brightness")
```

Nothing needs to get added to anything.

```python
imax = 5000
slope = 0.1
x = im9069.data
y = im6312.data
m = x < imax
m = m & (x > -100) & (y > -100)
m = m & (y < 1.2 * slope * imax)
m = m & ~im9069.mask & ~im6312.mask
df = pd.DataFrame(
    {
        "9069": x.data[m],
        "6312": y.data[m],
    }
)
```

```python
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)

g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope * imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, 1.2 * slope * imax], "--", color="r")
g.fig.suptitle("ZOOMED Correlation between [S III] 9069 and 6312 brightness")
```

```python
df.describe()
```

### Now we need to correct for reddening.

```python
imEBV = Image(str(ROOT / "data/ngc346-ZZ-reddening-E_BV.fits"))
```

```python
A9069 = rc.X(9069) * imEBV
im9069c = im9069.copy()
im9069c.data = im9069.data * 10 ** (0.4 * A9069.data)
```

```python
A6312 = rc.X(6312) * imEBV
im6312c = im6312.copy()
im6312c.data = (im6312.data) * 10 ** (0.4 * A6312.data)
```

```python
n = 1
fig, ax = plt.subplots(figsize=(12, 12))
(im6312c.rebin(n) / im9069c.rebin(n)).plot(
    vmin=0.1, vmax=0.15, cmap="magma", colorbar="v"
)
```

```python
n = 4
x = np.log10(im9069c.rebin(n).data)
y = np.log10(im6312c.rebin(n).data / im9069c.rebin(n).data)
m = (x > 3.0) & (x < 4.5)
m = m & (y > -1.05) & (y < -0.75)
m = m & ~im9069c.rebin(n).mask & ~im6312c.rebin(n).mask
df = pd.DataFrame(
    {
        "log10 9069": x[m],
        "log 10 6312/9069": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
```

This is totally different from what we got from the ESO cube.  This version is much more reliable, since with the other one we had no idea where the zeropoint was


Now, make a mask of EW(6312).  But first, we need to correct the zero point of the continuum.

```python
im6312_zero = 0.0
cont6312_zero = 0.0
imax = 600
x = im6312.data - im6312_zero
y = cont6312.data - cont6312_zero
m = x < imax
m = m & (x > 100)
m = m & (y < 3000) & (y > 100)
m = m & ~cont6312.mask & ~im6312.mask
df = pd.DataFrame(
    {
        "6312": x[m],
        "cont": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)

g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
```

```python
fig, ax = plt.subplots(figsize=(10, 10))
ew6312 = 3 * 1.25 * (im6312 - im6312_zero) / (cont6312 - cont6312_zero)
ew6312.plot(vmin=1.0, vmax=5.0, scale="sqrt")
ax.contour(ew6312.data, levels=[0.5], colors="r", linewidths=3)
ax.contour(im9069.data, levels=[6000.0], colors="k", linewidths=1)
```

```python
fixmask = (ew6312.data < 1.0) | (im9069.data < 3000.0)
fixmask = fixmask & (im6312c.data < 0.1 * im9069c.data)
fixmask = fixmask & (im6312c.data > 0.2 * im9069c.data)
iborder = 12
fixmask[:iborder, :] = True
fixmask[-iborder:, :] = True
fixmask[:, :iborder] = True
fixmask[:, -iborder:] = True
```

```python
im6312c.mask = im6312c.mask | fixmask
```

```python
n = 4
fig, ax = plt.subplots(figsize=(12, 12))
(im6312c.rebin(n) / im9069c.rebin(n)).plot(
    vmin=0.1, vmax=0.15, cmap="magma", colorbar="v"
)
```

Honestly, this looks identical to the last one. So I do not know what I am even doing here.

```python
n = 4
x = np.log10(im9069c.rebin(n).data)
y = np.log10(im6312c.rebin(n).data / im9069c.rebin(n).data)
z = im9069c.rebin(n).data
m = (x > 2.5) & (x < 5.5)
m = m & (y > -1.1) & (y < -0.6)
m = m & ~im9069c.rebin(n).mask & ~im6312c.rebin(n).mask
df = pd.DataFrame(
    {
        "log10 9069": x[m],
        "log 10 6312/9069": y[m],
    }
)
kws = dict(weights=z[m], bins=30)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=kws,
    diag_kws=kws,
)
```

### Convert to actual temperatures with pyneb

```python
s3 = pn.Atom("S", 3)
```

```python
s3.getTemDen([0.1, 0.2], den=100.0, wave1=6300, wave2=9069)
```

```python
r_s3_grid = np.linspace(0.05, 0.25, 201)
T_s3_grid = s3.getTemDen(r_s3_grid, den=100.0, wave1=6300, wave2=9069)
```

```python
imT_siii = im6312c.clone(data_init=np.empty)
imT_siii.data[~fixmask] = np.interp(
    im6312c.data[~fixmask] / im9069c.data[~fixmask],
    r_s3_grid,
    T_s3_grid,
    left=np.nan,
    right=np.nan,
)
# imT_siii.mask = imT_siii.mask | fixmask
# imT_siii.data[imT_siii.mask] = np.nan
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
imT_siii.rebin(2).plot(colorbar="v", cmap="hot", vmin=12000, vmax=16000)
```

```python
badpix = ~np.isfinite(imT_siii.data)
imT_siii.mask = imT_siii.mask | badpix
```

```python
imT_siii.write(str(ROOT / "data/ngc346-ZZ-T-siii.fits"), savemask="nan")
```

The rather disappointing conclusion of this is that the [S III] temperatures do vary from about 13 to 16 kK, but they don't show anything special at the bow shock, being about 13.7 +/- 0.4 kK there.

If anything, the T is lower in the bow shock.

Average over whole FOV is 14.2 +/- 0.8 kK after smoothing to eliminate the noise contribution.  This implies $t^2 = 0.003$ in plane of sky, which is small.


## Calculate [O III]/[S III]

```python
im5007 = get_image_raw("o-iii-5006-84")
```

Correct for extinction:

```python
A5007 = rc.X(5007) * imEBV
im5007c = im5007.copy()
im5007c.data = im5007.data * 10 ** (0.4 * A5007.data)
```

```python
np.nanmedian(imEBV[150:250, 200:300].data.data)
```

```python
median_EBV = np.nanmedian(imEBV[150:250, 200:300].data.data)
median_EBV
```

```python
im5007cc = im5007.copy()
im5007cc.data = im5007.data * 10 ** (0.4 * rc.X(5007) * median_EBV)
im9069cc = im9069.copy()
im9069cc.data = im9069.data * 10 ** (0.4 * rc.X(9069) * median_EBV)
```

Quick look:

```python
fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
((im5007cc) / im9069cc).plot(
    vmin=10,
    vmax=55,
    colorbar="v",
    scale="linear",
    cmap="Purples",
    ax=axes[0],
)
((im5007c) / im9069c).plot(
    vmin=10,
    vmax=55,
    colorbar="v",
    scale="linear",
    cmap="Purples",
    ax=axes[1],
)
```

The left panel is with a median reddening. The right panel is with a pixel-by-pixel reddening but that fails at the position of stars.


Check zero points:

```python
imax = 15000
imin = 2000
slope = 45
slope2 = 35
x = im9069cc.data
y = im5007cc.data
m = x < imax
m = m & (x > imin)
m = m & (y < slope * imax)
m = m & (y > slope * imin)
m = m & ~im9069c.mask & ~im5007c.mask
df = pd.DataFrame(
    {
        "9069": x[m],
        "5007": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
# g.axes[1, 0].axvline(0.0, color="r")
# g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope * imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, slope2 * imax], "--", color="r")
g.fig.suptitle("Correlation between [S III] 9069 and [O III] 5007 brightness")
```

```python
slope = 45
slope2 = 35
x = im9069cc.data
y = im5007cc.data
m = x > imin
m = m & (y > 10 * x) & (y < 60 * x)
m = m & ~im9069cc.mask & ~im5007cc.mask
df = pd.DataFrame(
    {
        "log 9069": np.log10(x[m]),
        "log 5007/9069": np.log10(y[m] / x[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.fig.suptitle("Correlation between [S III] 9069 and [O III] / [S III] ratio")
```

```python
imR_oiii_siii = (im5007cc) / im9069cc
imR_oiii_siii.write(
    str(ROOT / "data/ngc346-ZZ-R-oiii-5007-siii-9069.fits"), savemask="nan"
)
```

## Calculate [O III] / Hβ

This might be better since at least it is not affected by reddening.

```python
imhb = get_image_raw("h-i-4861-32")
```

```python
imax = 60000
imin = 2000
slope = 5.0
slope2 = 6.0
x = imhb.data
y = im5007.data
m = x < imax
m = m & (x > imin)
m = m & (y < slope * imax)
m = m & (y > slope * imin)
m = m & ~imhb.mask & ~im5007.mask
df = pd.DataFrame(
    {
        "4861": x[m],
        "5007": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.axes[1, 0].plot([0, imax], [0, slope * imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, slope2 * imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and [O III] 5007 brightness")
```

```python
n = 1
x = imhb.rebin(n).data
y = im5007.rebin(n).data
m = (x > imin) & (y > x) & (y < 10 * x)
m = m & ~imhb.rebin(n).mask & ~im5007.rebin(n).mask
df = pd.DataFrame(
    {
        "log10(4861)": np.log10(x[m]),
        "log10(5007 / 4861)": np.log10(y[m] / x[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
)
g.fig.suptitle("Correlation between Hβ 4861 and [O III] / Hβ ratio")
```

```python
imR_oiii_hb = (im5007) / (imhb)
imR_oiii_hb.write(str(ROOT / "data/ngc346-ZZ-R-oiii-5007-hi-4861.fits"), savemask="nan")
```

```python
im4740 = get_image_raw("ar-iv-4740-17")
```

```python
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel
```

```python
im = convolve_fft(im4740.data, Gaussian2DKernel(5.0))
```

```python
fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
imR_oiii_siii.plot(
    vmin=10,
    vmax=55,
    colorbar="v",
    scale="linear",
    ax=axes[0],
    cmap="Purples",
)
imR_oiii_hb.plot(
    vmin=2,
    vmax=6.5,
    colorbar="v",
    scale="linear",
    ax=axes[1],
    cmap="Greens",
)
for ax in axes:
    ax.contour(im, levels=[100, 150, 200], linewidths=1, colors="y")
axes[0].set_title("[O III] / [S III]")
axes[1].set_title("[O III] / Hβ")
```

So, the bow shock shows up much better in [O III]/[S III] than it does in [O III]/Hb.  This implies that it must be a lack of [S III], rather than an excess of [O III] that characterises the bow shock.
