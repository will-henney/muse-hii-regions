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

# Notebook C of line ratios from saved maps from the PZ cube

This is from the splitting up of ZZ-01-01-more-line-ratios.ipynb. It includes analysis of the He I and He II lines

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

## Calculate He I / Hβ

Let us see if this has a hole in it where the He II is coming from.
+

```python
im5875 = get_image_raw("he-i-5875-62")
im4922 = get_image_raw("he-i-4921-93")
im5048 = get_image_raw("he-i-5047-74")
im5016 = get_image_raw("he-i-5015-68")
im6678 = get_image_raw("he-i-6678-15")
im7065 = get_image_raw("he-i-7065-28")
im7281 = get_image_raw("he-i-7281-35")
imcont = get_image_raw("he-i-5875-62", variant="cont")
ew5875 = 3 * 1.25 * im5875 / imcont
imhb = get_image_raw("h-i-4861-32")
imha = get_image_raw("h-i-656*")
```

```python
fig, axes = plt.subplots(3, 2, sharey=True, figsize=(12, 18))
im5875.plot(ax=axes[0, 0], vmin=0, vmax=7500)
axes[0, 0].contour(ew5875.data, levels=[3], colors="r", linewidths=1)
im4922.plot(ax=axes[0, 1], vmin=0, vmax=550)
(imhb).plot(ax=axes[1, 0], vmin=0, vmax=60000)
im5016.plot(ax=axes[1, 1], vmin=0, vmax=700)
im6678.plot(ax=axes[2, 0], vmin=0, vmax=2000)
im7065.plot(ax=axes[2, 1], vmin=0, vmax=1400)
axes[0, 0].set_title("He I 5875")
axes[0, 1].set_title("He I 4922")
axes[1, 0].set_title("H I 4861")
axes[1, 1].set_title("He I 5016")
axes[2, 0].set_title("He I 6678")
axes[2, 1].set_title("He I 7065")
```


So 5875 is 10 to 100 times brighter than the other two. And it is almost identical to Hβ!

```python
imax = 30000
imin = 1000
slope = 0.12
x = imhb.data
y = im5875.data
m = x < imax
m = m & (x > imin)
m = m & (y < slope * imax)
m = m & (y > imin * slope)
m = m & ~imhb.mask & ~im5875.mask
df = pd.DataFrame(
    {
        "4861": x[m],
        "5875": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(weights=x[m], bins=200),
    diag_kws=dict(weights=x[m], bins=200),
)
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope * imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and He I 5875 brightness")
```

```python
imax = 100000
imin = 1000
slope = 0.12
x = imhb.data
y = im5875.data
m = x < imax
m = m & (x > imin)
m = m & (y < slope * imax)
m = m & (y > imin * slope)
m = m & ~imhb.mask & ~im5875.mask
df = pd.DataFrame(
    {
        "4861": x[m],
        "5875": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(weights=x[m], bins=200),
    diag_kws=dict(weights=x[m], bins=200),
)
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope * imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and He I 5875 brightness")
```

```python
imEBV = Image(str(ROOT / "data/ngc346-ZZ-reddening-E_BV.fits"))
```

```python
imR_hei_hb = im5875 / (imhb)

fig, ax = plt.subplots(figsize=(12, 12))
imR_hei_hb.plot(colorbar="v", cmap="gray", vmin=0.11, vmax=0.14)
```

```python
red_R_hei_hb = imR_hei_hb.copy()
red_R_hei_hb.data = 10 ** (0.4 * imEBV.data * (rc.X(4861) - rc.X(5875)))

fig, ax = plt.subplots(figsize=(12, 12))
(imR_hei_hb / red_R_hei_hb).plot(colorbar="v", cmap="gray", vmin=0.1, vmax=0.115)
```

So if we correct it for reddening, then lots of spurious structure disappears.  But we are left with very little variation at all, except for at the mYSO and the top right corner, which both show low He I.

```python
fig, ax = plt.subplots(figsize=(12, 12))
(im6678 / imha).plot(colorbar="v", cmap="gray", vmin=0.009, vmax=0.011)
```

```python
1/3.5
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
n = 4
(im6678.rebin(n) / im5875.rebin(n)).plot(colorbar="v", cmap="gray", vmin=0.25, vmax=0.35)
```

### He I emissivities

```python
hei = pn.RecAtom("He", 1)

dens = [10, 50, 100, 200, 400, 800]
tems = [11000, 13000, 15000, 18000]
e4713 = hei.getEmissivity(tems, dens, wave=4713)
e4921 = hei.getEmissivity(tems, dens, wave=4922)
e5016 = hei.getEmissivity(tems, dens, wave=5016)
e5047 = hei.getEmissivity(tems, dens, wave=5047)
e5876 = hei.getEmissivity(tems, dens, wave=5876)
e6678 = hei.getEmissivity(tems, dens, wave=6678)
e7065 = hei.getEmissivity(tems, dens, wave=7065)
e7281 = hei.getEmissivity(tems, dens, wave=7281)
```

```python
np.round(e5016 / e5876, 3)
```

So 5016/5876 is a T indicator for low densities (becomes insensitive at high densities). 

At 50 pcc it is +15% with T

```python
np.round(e5876 / e6678, 3)
```

5876 / 6678 is flat with n (+0.3% at low-T, -1.4% at high-T) and flat with n (+0.5% at n=50)

```python
np.round((e7281 + e7065) / e6678, 3)
```

So (7281 + 7065) / 6678 has positive variation with T (30%) and n (20%)

```python
np.round(e7065 / e6678, 3)
```

Just 7065/6678 is positive with T (+37%) and n (+32%)

```python
np.round(e7065 / e7281, 3)
```

7065 / 7281 is flattish with T (+5%) and positive with density (+24%)

```python
np.round(e7281 / e6678, 3)
```

7281 / 6678 is flattish with n (+7%) and positive with temperature (+29%)

```python
100 * (253 - 196)/196
```

```python
np.round(e4713 / e5876, 3)
```

4713 / 5876 is too weak

```python
np.round(e4921/ e5876, 3)
```

4921/5876 is rather weak too and not very sensitive: +2% with T at 50pcc, -5% with n at 13,000 K


#### Conclusions on He I ratios

* 5016/5876 is a reasonable T indicator (for densities < 500 pcc) *but it is much lower than predicted – maybe radiative transfer effects?*
* 7065/7281 is a possible density indicator
* 7281/6678 is a possible T indicator

But we need to do reddening correction for all of them. 

```python
_ebv = imEBV.data[np.isfinite(imEBV.data)]

np.mean(_ebv), np.std(_ebv), np.median(_ebv), np.median(np.abs(_ebv - np.median(_ebv)))
```

```python
10 ** (0.4 * 0.156 * (rc.X(7065) - rc.X(7281))), 10 ** (0.4 * 0.156 * (rc.X(7281) - rc.X(6678)))
```

### He I n-sensitive ratio



```python
fig, ax = plt.subplots(figsize=(12, 12))
n = 4
(1.012 * im7065.rebin(4) / im7281.rebin(4)).plot(colorbar="v", cmap="magma", vmin=3, vmax=5)
```

```python
n = 4
x = im6678.rebin(n).data
y = 1.012 * im7065.rebin(n).data / im7281.rebin(n).data 
z = im6678.rebin(n).data
m = (x > 0) & (x < 2000)
m = m & (y > 3) & (y < 5)
m = m & ~im7065.rebin(n).mask & ~im7281.rebin(n).mask
df = pd.DataFrame(
    {
        "He I 6678": x[m],
        "7065/7281": y[m],
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

```python
np.round(e7065 / e7281, 3)
```

### He I T-sensitive ratio




The 7281/6678 ratio depends primarily on T

```python
fig, ax = plt.subplots(figsize=(12, 12))
n = 4
(0.96 * im7281.rebin(4) / im6678.rebin(4)).plot(colorbar="v", cmap="magma", vmin=0.17, vmax=0.22)
```

```python

```

```python
n = 4
x = im6678.rebin(n).data
y = 0.96 * im7281.rebin(n).data / im6678.rebin(n).data
z = im6678.rebin(n).data
m = (x > 0) & (x < 2000)
m = m & (y > 0.1) & (y < 0.25)
m = m & ~im6678.rebin(n).mask & ~im7281.rebin(n).mask
df = pd.DataFrame(
    {
        "He I 6678": x[m],
        "7281/6678": y[m],
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

```python
np.round(e7281 / e6678, 3)
```

```python
tems
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
(im5875 / im6678).plot(colorbar="v", cmap="viridis", vmin=3.4, vmax=3.8)
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
(im5016 / im5875).plot(colorbar="v", cmap="viridis", vmin=0.05, vmax=0.15)
```


## Calculate He II / Hβ



```python
im4686 = get_image_raw("he-ii-4685-68")
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
im4686.plot(colorbar="v", cmap="gray", vmin=-30.0, vmax=300)
```

```python
imR_heii_hb = im4686 / imhb

fig, ax = plt.subplots(figsize=(12, 12))
imR_heii_hb.plot(colorbar="v", cmap="gray", vmin=-0.003, vmax=0.01)
```

```python
n = 2
xslice, yslice = slice(200, 300), slice(100, 250)
x = imR_heii_hb[yslice, xslice].rebin(n).data
y = (
    imR_hei_hb[yslice, xslice].rebin(n).data
    / red_R_hei_hb[yslice, xslice].rebin(n).data
)
z = im4686[yslice, xslice].rebin(n).data
m = x < 0.03
m = m & (x > 0)
m = m & (y < 0.113)
m = m & (y > 0.103)
m = (
    m
    & ~imR_heii_hb[yslice, xslice].rebin(n).mask
    & ~imR_hei_hb[yslice, xslice].rebin(n).mask
)
df = pd.DataFrame(
    {
        "4686 / 4861": x[m],
        "5875 / 4861": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(weights=z[m], bins=30),
    diag_kws=dict(weights=z[m], bins=30),
)
g.fig.suptitle("Correlation between He II / Hβ and He I / Hβ ratios")
```

So there is a *tiny* change in 5875/4861 from 0.109 to 0.107 as 4686/4861 increases.

```python
df["high"] = df["4686 / 4861"] > 0.005
df
```

```python
sns.histplot(
    data=df,
    x="5875 / 4861",
    hue="high",
    multiple="dodge",
    shrink=1.0,
    stat="probability",
    common_norm=False,
    alpha=0.8,
    bins=6,
)
```

```python
df["high"] = df["4686 / 4861"] > 0.005
sns.histplot(
    data=df,
    x="5875 / 4861",
    hue="high",
    # multiple="dodge",
    shrink=1.0,
    stat="probability",
    cumulative=True,
    common_norm=False,
    alpha=0.8,
    bins=1000,
    element="step",
).axhline(0.5, linestyle="dotted", color="k")
```

Try to remove the pattern noise

```python
imR_heii_hb_fake = imR_heii_hb.copy()
m = imR_heii_hb_fake.data > 0.01
imR_heii_hb_fake.data[m] = np.nan
vv = np.nanmedian(imR_heii_hb_fake.data, axis=0, keepdims=True)
hh = np.nanmedian(imR_heii_hb_fake.data, axis=1, keepdims=True)
imR_heii_hb_fake.data = vv + hh
imR_heii_hb_fake.data -= np.nanmedian(vv + hh)
imR_heii_hb_c = imR_heii_hb.copy()
imR_heii_hb_c.data -= imR_heii_hb_fake.data
```

```python
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
imR_heii_hb_fake.plot(ax=axes[0, 0], colorbar="v", cmap="gray", vmin=-0.01, vmax=0.01)
(imR_hei_hb / red_R_hei_hb).plot(
    ax=axes[0, 1], colorbar="v", cmap="gray", vmin=0.1, vmax=0.115
)
imR_heii_hb.plot(ax=axes[1, 0], colorbar="v", cmap="gray", vmin=-0.01, vmax=0.01)
imR_heii_hb_c.plot(ax=axes[1, 1], colorbar="v", cmap="gray", vmin=-0.01, vmax=0.01)
```

```python
n = 2
xslice, yslice = slice(200, 300), slice(100, 250)
# xslice, yslice = slice(230, 300), slice(130, 220)
x = imR_heii_hb_c[yslice, xslice].rebin(n).data
y = (
    imR_hei_hb[yslice, xslice].rebin(n).data
    / red_R_hei_hb[yslice, xslice].rebin(n).data
)
z = im4686[yslice, xslice].rebin(n).data
m = x < 0.03
m = m & (x > 0)
m = m & (y < 0.113)
m = m & (y > 0.103)
m = (
    m
    & ~imR_heii_hb[yslice, xslice].rebin(n).mask
    & ~imR_hei_hb[yslice, xslice].rebin(n).mask
)
df = pd.DataFrame(
    {
        "4686 / 4861": x[m],
        "5875 / 4861": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(weights=z[m], bins=10),
    diag_kws=dict(weights=z[m], bins=10),
)
g.fig.suptitle("Correlation between He II / Hβ and He I / Hβ ratios")
```

```python
df["high"] = df["4686 / 4861"] > 0.003
sns.histplot(
    data=df,
    x="5875 / 4861",
    hue="high",
    multiple="dodge",
    shrink=1.0,
    stat="probability",
    common_norm=False,
    alpha=0.8,
    bins=6,
)
```

```python
df["high"] = df["4686 / 4861"] > 0.005
sns.histplot(
    data=df,
    x="5875 / 4861",
    hue="high",
    # multiple="dodge",
    shrink=1.0,
    stat="probability",
    cumulative=True,
    common_norm=False,
    alpha=0.8,
    bins=1000,
    element="step",
).axhline(0.5, linestyle="dotted", color="k")
```

So, the depatterning improved the map but did not help with the statistics.



```python

```
