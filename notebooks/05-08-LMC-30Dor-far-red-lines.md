---
jupyter:
  jupytext:
    encoding: '# -*- coding: utf-8 -*-'
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# 30 Dor far red lines: 7800–8700 Å, 8600–9500 Å

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import regions
import sys

sys.path.append("../lib")
import moments
import extract

sns.set_context("talk")
sns.set_color_codes()
```

## Start with 7800–8700

```python tags=[]
cubeA = Cube("../big-data/lmc-30dor-A-subcube-78-87.fits")
cubeB = Cube("../big-data/lmc-30dor-B-subcube-78-87.fits")
cubeC = Cube("../big-data/lmc-30dor-C-subcube-78-87.fits")
cubeD = Cube("../big-data/lmc-30dor-D-subcube-78-87.fits")
```

### Inspect the average spectrum for each field

```python
fig, ax = plt.subplots(figsize=(12, 6))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cube.sum(axis=(1, 2)).plot(label=label)
ax.legend()
ax.set(
    ylim=[0.1e8, 0.9e8],
)
sns.despine()
```

We can see the optical ghost around 8600 Å.  But the main problem here is going to be the Paschen jump.  Let us hope that a 5th order polynomial can fit it OK.


### Define continuum wavelength ranges

We can try the same wav ranges as we used for NGC 346

```python
wavranges = [
    (7800, 7805), (7825, 7835), (7900, 7910),
    (7935, 7950), (8000, 8010), (8025, 8045), (8067, 8082),
    (8107, 8125), (8160, 8180),
    (8200, 8220), (8250, 8270), (8354, 8364), (8368, 8378), (8385, 8395),
    (8405, 8415), (8432, 8442), (8480, 8490),
    (8525, 8535), (8560, 8570), (8590, 8600),
    (8620, 8640), (8680, 8700),
]
```

```python
nv, ny, nx = cubeA.data.shape
ny, nx
```

Split each field up into 4 roughly equal tiles:

```python
mm = 160  # middle of each image
fig, ax = plt.subplots(figsize=(12, 8))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cu11 = cube[:, :mm, :mm]
    cu12 = cube[:, :mm, mm:]
    cu21 = cube[:, mm:, :mm]
    cu22 = cube[:, mm:, mm:]
    cu11.mean(axis=(1, 2)).plot(label=f"{label}11")
    cu12.mean(axis=(1, 2)).plot(label=f"{label}12")
    cu21.mean(axis=(1, 2)).plot(label=f"{label}21")
    cu22.mean(axis=(1, 2)).plot(label=f"{label}22")

for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.legend(ncol=4, fontsize="x-small")
ax.set(
    yscale="log",
    ylim=[50, 5e3],
)
sns.despine()
```

### Test the polynomial fitting with field A

This takes about a minute for each field.

```python
contA = extract.fit_continuum(
    cubeA,
    wav_ranges=wavranges,
    deg=8,
    median=False,
)
```

#### Inspect the results for different portions of the field.

```python
fig, ax = plt.subplots(figsize=(12, 8))
cubeA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
contA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.set(ylim=[0, 800])
```

I had to go back and forth a few times adjusting the wav ranges.  It is difficult to get a good fit on the blue side because of the WR features.

As can be seen here, the final version is not perfect – it slightly overpredicts the continuum around 4700 to 4800.  This might affect some of the weak [Fe III] lines, but the [Ar IV] 4740 does not seem to be much affected.


### Now do the other fields

```python
contB = extract.fit_continuum(
    cubeB,
    wav_ranges=wavranges,
    deg=8,
    median=False,
)
```

```python
contC = extract.fit_continuum(
    cubeC,
    wav_ranges=wavranges,
    deg=8,
    median=False,
)
```

```python
contD = extract.fit_continuum(
    cubeD,
    wav_ranges=wavranges,
    deg=8,
    median=False,
)
```

### Save the continuum-subtracted cubes

```python
csub = {}
cdiv = {}
for cube, cont, label in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    "DCBA",
):
    prefix = f"../big-data/lmc-30dor-{label}-subcube-78-87"
    csub[label] = cube - cont
    cdiv[label] = cube / cont
    csub[label].write(
        f"{prefix}-contsub.fits",
        savemask="nan",
    )
    cdiv[label].write(
        f"{prefix}-contdiv.fits",
        savemask="nan",
    )
    cont.write(
        f"{prefix}-cont.fits",
        savemask="nan",
    )
```


## Repeat for the next range: 8600 to 9500

```python tags=[]
cubeA = Cube("../big-data/lmc-30dor-A-subcube-86-95.fits")
cubeB = Cube("../big-data/lmc-30dor-B-subcube-86-95.fits")
cubeC = Cube("../big-data/lmc-30dor-C-subcube-86-95.fits")
cubeD = Cube("../big-data/lmc-30dor-D-subcube-86-95.fits")
```

### 86–95 Inspect the average spectrum for each field

```python
fig, ax = plt.subplots(figsize=(12, 6))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cube.sum(axis=(1, 2)).plot(label=label)
ax.legend()
ax.set(
    ylim=[0.3e8, 0.7e8],
)
sns.despine()
```

Again, we get WR features (I think)  around 7120, and some over compensation for the atmospheric asbsorption.


### 86–95 Define continuum wavelength ranges

No large ranges of clear continuum here

```python
wavranges = [
    (8658, 8668), (8675, 8685), (8696, 8705),
    (8720, 8730), (8745, 8755), (8800, 8820), 
    (8890, 8900), (8930, 8940), (8970, 8980),
    (9140, 9150),
    (9200, 9210), (9250, 9270),
]
```

```python
nv, ny, nx = cubeA.data.shape
ny, nx
```

Split each field up into 4 roughly equal tiles:

```python
mm = 160  # middle of each image
fig, ax = plt.subplots(figsize=(12, 8))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cu11 = cube[:, :mm, :mm]
    cu12 = cube[:, :mm, mm:]
    cu21 = cube[:, mm:, :mm]
    cu22 = cube[:, mm:, mm:]
    cu11.mean(axis=(1, 2)).plot(label=f"{label}11")
    cu12.mean(axis=(1, 2)).plot(label=f"{label}12")
    cu21.mean(axis=(1, 2)).plot(label=f"{label}21")
    cu22.mean(axis=(1, 2)).plot(label=f"{label}22")

for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.legend(ncol=4, fontsize="x-small")
ax.set(
    yscale="log",
    ylim=[40, 8e3],
)
sns.despine()
```

### 86-94 Test the polynomial fitting with field A

This takes about a minute for each field.

```python
contA = extract.fit_continuum(
    cubeA,
    wav_ranges=wavranges,
    deg=2,
    median=False,
)
```

#### 86–94 Inspect the results for different portions of the field.

```python
fig, ax = plt.subplots(figsize=(12, 8))
cubeA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
contA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.set(ylim=[0, 800])
```

### 86–94 Now do the other fields

```python
contB = extract.fit_continuum(
    cubeB,
    wav_ranges=wavranges,
    deg=2,
    median=False,
)
```

```python
contC = extract.fit_continuum(
    cubeC,
    wav_ranges=wavranges,
    deg=2,
    median=False,
)
```

```python
contD = extract.fit_continuum(
    cubeD,
    wav_ranges=wavranges,
    deg=2,
    median=False,
)
```

### 86-94 Save the continuum-subtracted cubes

```python
csub = {}
cdiv = {}
for cube, cont, label in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    "DCBA",
):
    prefix = f"../big-data/lmc-30dor-{label}-subcube-86-94"
    csub[label] = cube - cont
    cdiv[label] = cube / cont
    csub[label].write(
        f"{prefix}-contsub.fits",
        savemask="nan",
    )
    cdiv[label].write(
        f"{prefix}-contdiv.fits",
        savemask="nan",
    )
    cont.write(
        f"{prefix}-cont.fits",
        savemask="nan",
    )
```

```python

```

