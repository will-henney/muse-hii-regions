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

<!-- #region -->
# Spatial profiles along slit from FORS1 spectra

I would like to achieve the following:

1. [X] Identify W 3 and the bow shock in the spectra
2. [ ] Trace the full western extent of the bow shock in [Ar IV]
    - the MUSE field extends only 22 arcsec from the star
    - the FORS1 spectrum shows a roughly linear ramp that extends about 34 arcsec (see image below)
3. [ ] Calculate the [O III] 4363/5007 temperature profile of the bow shock
    - See if there is any evidence for temperatures as high as the [Ar IV] temperature
4. [ ] Measure the [Ne III] 3869 profile.  Ne$^+$ has the same ionization potential as Ar$^{+2}$ (40 eV), but [Ne III] seems to be much more broadly distributed than [Ar IV]



![FORS1 profile of Ar IV](assets/ngc346-screenshot-fors1-ariv-profile.png)


<!-- #endregion -->

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import regions
from astropy.io import fits
from astropy.wcs import WCS
import sys

sys.path.append("../lib")
import extract

sns.set_context("talk")
sns.set_color_codes()
```

Load line and continuum images for the two blue spectral ranges.  We may use the continuum for masking out regions that have excessive stellar contamination.

```python
hdu1 = fits.open("../data/ngc346-fors1-A-3500-4600-contsub.fits")[0]
hdu2 = fits.open("../data/ngc346-fors1-A-4400-6000-contsub.fits")[0]
hdu1c = fits.open("../data/ngc346-fors1-A-3500-4600-cont.fits")[0]
hdu2c = fits.open("../data/ngc346-fors1-A-4400-6000-cont.fits")[0]
```

```python
hdu1B = fits.open("../data/ngc346-fors1-B-3562-4600-contsub.fits")[0]
hdu2B = fits.open("../data/ngc346-fors1-B-4400-6000-contsub.fits")[0]
hdu1cB = fits.open("../data/ngc346-fors1-B-3562-4600-cont.fits")[0]
hdu2cB = fits.open("../data/ngc346-fors1-B-4400-6000-cont.fits")[0]
```

Set up the world coords for the two ranges

```python
wcs1 = WCS(hdu1)
wcs2 = WCS(hdu2)
wcs1B = WCS(hdu1B)
wcs2B = WCS(hdu2B)
for w in wcs1, wcs2, wcs1B, wcs2B:
    w.wcs.cdelt[1] *= -1.0
```

Set up the lines that we want with wavelength limits:

```python
pv1A = extract.PositionVelocityImage(hdu1.data, wcs1)
pv2A = extract.PositionVelocityImage(hdu2.data, wcs2)
pvc1A = extract.PositionVelocityImage(hdu1c.data, wcs1)
pvc2A = extract.PositionVelocityImage(hdu2c.data, wcs2)
pv1B = extract.PositionVelocityImage(hdu1B.data, wcs1B)
pv2B = extract.PositionVelocityImage(hdu2B.data, wcs2B)
pvc1B = extract.PositionVelocityImage(hdu1cB.data, wcs1B)
pvc2B = extract.PositionVelocityImage(hdu2cB.data, wcs2B)
```

```python
restwav = {
    "He II 4686": 4685.68,
    "[Ar IV] 4740": 4740.17,
    "[O III] 5007": 5006.84,
    "[O III] 4363": 4363.209,
    "[Ne III] 3869": 3869.07,
    "H I 4861": 4861.32,
    "H I 4340": 4340.463,
    "H I 4102": 4101.735,
    "He I 5876": 5875.62,
    "[Cl III] 5518": 5517.71,
    "[Cl III] 5538": 5537.88,
    "[N II] 5755": 5755.08,
    "[Fe III] 4658": 4658.10,
    "[Fe III] 4987": 4987.20,
    "[Fe III] 5270": 5270.4,
    "O II 4650": 4650.00,
}
```

```python
emlines = {k: extract.EmissionLine(k, v) for k, v in restwav.items()}
```

```python
emlines
```

```python
for em in emlines.values():
    if em.wav0 > 4400:
        em.pvA = pv2A
        em.pvB = pv2B
        em.pvcA = pvc2A
        em.pvcB = pvc2B
    else:
        em.pvA = pv1A
        em.pvB = pv1B
        em.pvcA = pvc1A
        em.pvcB = pvc1B
```

```python
for em in emlines.values():
    em.A = em.pvA.slit_profile(em)
    em.ewA = em.pvA.slit_ew_profile(em, em.pvcA)
    em.B = em.pvB.slit_profile(em)
    em.ewB = em.pvB.slit_ew_profile(em, em.pvcB)
```

```python
fig, axes = plt.subplots(
    2,
    1,
    figsize=(12, 12),
    sharex=True,
)
e = emlines["[Ar IV] 4740"]
axes[0].plot(e.A.position, e.A.data, linewidth=1.0, alpha=1.0)
axes[1].plot(e.B.position, e.B.data, linewidth=1.0, alpha=1.0)

e = emlines["He II 4686"]
axes[0].plot(e.A.position, e.A.data, linewidth=2.0, alpha=1.0)
axes[1].plot(e.B.position, e.B.data, linewidth=2.0, alpha=1.0)

axes[0].axhline(0.0, linestyle="dashed", color="k")
axes[0].axvline(0.0, linestyle="dotted", color="k")

axes[1].axhline(0.0, linestyle="dashed", color="k")
axes[1].axvline(0.0, linestyle="dotted", color="k")

for ax in axes:
    ax.set(
        xlim=[-25, 75],
        ylim=[-30, 60],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 12))
e = emlines["[O III] 5007"]
axes[0].plot(e.A.position, 1e-3 * e.A.data, linewidth=1.0, alpha=1.0)
axes[1].plot(e.B.position, 1e-3 * e.B.data, linewidth=1.0, alpha=1.0)

e = emlines["[O III] 4363"]
axes[0].plot(e.A.position, 0.075 * e.A.data, linewidth=1.0, alpha=1.0)
axes[1].plot(e.B.position, 0.075 * e.B.data, linewidth=1.0, alpha=1.0)

e = emlines["He II 4686"]
axes[0].plot(e.A.position, e.A.data, linewidth=2.0, alpha=1.0)
axes[1].plot(e.B.position, e.B.data, linewidth=2.0, alpha=1.0)

for ax in axes:
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.set(
        xlim=[-25, 75],
        ylim=[-30, 60],
    )
```

```python
Amask = emlines["[O III] 5007"].ewA.data > 200.0
Bmask = emlines["[O III] 5007"].ewB.data > 200.0

for e in emlines.values():
    e.multiA = e.A.multibin(mask=Amask)
    e.multiB = e.B.multibin(mask=Bmask)
```


```python
e.multiB[16].data
```

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
e1 = emlines["[O III] 4363"]
e2 = emlines["[O III] 5007"]
xmin, xmax = -10, 50

n = 1
ratio = e1.multiA[n].data / e2.multiA[n].data
pos = e1.multiA[n].position
ax.plot(pos, ratio, linewidth=1.0, alpha=1.0, drawstyle="steps-mid")

n = 8
ratio = e1.multiA[n].data / e2.multiA[n].data
pos = e1.multiA[n].position
ax.plot(pos, ratio, linewidth=4.0, alpha=1.0, drawstyle="steps-mid")

ax.plot(e2.A.position, 3e-7*e2.A.data)

ax.axhline(0.0, linestyle="dashed", color="k")
m = (pos > xmin) & (pos < xmax)
ax.axhline(
    np.nanmedian(ratio[m]),
    linestyle="dotted", color="k",
)
ax.axvline(0.0, linestyle="dotted", color="k")
ax.set(
    xlim=[xmin, xmax],
    ylim=[-0.003, 0.02],
);
```

```python
fig, axes = plt.subplots(
    2,
    1,
    figsize=(12, 12),
    sharex=True,
)
n = 8
e = emlines["[Ar IV] 4740"]
axes[0].plot(
    e.multiA[n].position, e.multiA[n].data, 
    linewidth=1.0, alpha=1.0, drawstyle="steps-mid",
)
axes[1].plot(
    e.multiB[n].position, e.multiB[n].data, 
    linewidth=1.0, alpha=1.0, drawstyle="steps-mid",
)

e = emlines["He II 4686"]
axes[0].plot(
    e.multiA[n].position, e.multiA[n].data, 
    linewidth=2.0, alpha=1.0, drawstyle="steps-mid",
)
axes[1].plot(
    e.multiB[n].position, e.multiB[n].data, 
    linewidth=2.0, alpha=1.0, drawstyle="steps-mid",
)

axes[0].axhline(0.0, linestyle="dashed", color="k")
axes[0].axvline(0.0, linestyle="dotted", color="k")

axes[1].axhline(0.0, linestyle="dashed", color="k")
axes[1].axvline(0.0, linestyle="dotted", color="k")

for ax in axes:
    ax.set(
        xlim=[-25, 75],
        ylim=[-30, 60],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, s in zip(axes, [posA, posB], [sA, sB]):
    ax.plot(
        pos,
        s["[O III] 4363"]
        / s["[O III] 5007"]
        / np.median(s["[O III] 4363"] / s["[O III] 5007"]),
        linewidth=0.3,
        alpha=1.0,
        label="[O III] 4363 / 5007",
    )
    ax.plot(
        pos,
        s["[Ne III] 3869"]
        / s["[O III] 5007"]
        / np.median(s["[Ne III] 3869"] / s["[O III] 5007"]),
        linewidth=0.3,
        alpha=1.0,
        label="[Ne III] 3869 / [O III] 5007",
    )
    ax.plot(
        pos,
        0.1 * s["[O III] 5007"] / np.median(s["[O III] 5007"]),
        #    0.1* s4341 / np.median(s4341),
        linewidth=1.0,
        alpha=1.0,
        label="[O III] 5007",
    )
    ax.plot(
        pos,
        s["H I 4340"] / s["H I 4861"] / np.median(s["H I 4340"] / s["H I 4861"]),
        linewidth=0.5,
        alpha=1.0,
        label="Hγ / Hβ",
    )
    # ax.plot(positions, 0.0003 * s4686, linewidth=2.0, alpha=1.0)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")

    ax.legend(fontsize="x-small")
    ax.set(
        # xlim=[-10, 50],
        xlim=[-230, 200],
        ylim=[-0.01, 1.5],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line in "[O III] 5007", "H I 4861":
        ax.plot(pos, ew[line], linewidth=3.0, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        # xlim=[-10, 50],
        # ylim=[-0.003, 0.02],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line in "He I 5876", "[Ne III] 3869":
        ax.plot(pos, ew[line], linewidth=3.0, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        # xlim=[-10, 50],
        # ylim=[-0.003, 0.02],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line in "[O III] 4363", "[N II] 5755":
        ax.plot(pos, ew[line], linewidth=3.0, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        # xlim=[-10, 50],
        # ylim=[-0.003, 0.02],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    m = ew["[O III] 5007"] > 250.0
    for line in "[Ar IV] 4740", "He II 4686":
        y = ew[line][:]
        y[~m] = np.nan
        ax.plot(pos, y, linewidth=1.5, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        xlim=[-10, 50],
        ylim=[-2, 3],
    )
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, s, ew in zip(axes, [posA, posB], [sA, sB], [ewA, ewB]):
    m = ew["[O III] 5007"] > 250.0
    for line in "[Fe III] 4658", "O II 4650":
        y = s[line][:]
        y[~m] = np.nan
        ax.plot(pos, y, linewidth=0.5, alpha=1.0, label=line)
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend()
    ax.set(
        xlim=[-100, 50],
        ylim=[-20, 40],
    )
```

```python
a = np.arange(4)
np.kron(a, np.ones((4,)))
```

```python
sys.path.append("../multibin-maps")
import rebin_utils
```

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
for ax, pos, ew in zip(axes, [posA, posB], [ewA, ewB]):
    for line, color in zip(["[Ar IV] 4740", "He II 4686"], "cr"):
        m = ew["[O III] 5007"] > 250.0
        x, y = pos[:], ew[line][:]
        w = np.ones_like(y)
        lw = 0.5
        for n in [1, 2, 4, 8]:
            ax.plot(
                x,
                y,
                linewidth=lw,
                alpha=1.0,
                label=line,
                color=color,
                drawstyle="steps-mid",
            )
            [x, y], m, w = rebin_utils.downsample1d([x, y], m, weights=w)
            y[~m] = np.nan
            lw += 0.5
            line = "_nolabel"
    ax.axhline(0.0, linestyle="dashed", color="k")
    ax.axvline(0.0, linestyle="dotted", color="k")
    ax.legend(fontsize="x-small")
    ax.set(
        xlim=[-20, 50],
        # xlim=[-100, 50],
        ylim=[-2, 3],
    )
```
