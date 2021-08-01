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

# 30 Dor red lines: 6200–7100 Å, 7000–7900, etc

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

## Start with 6200–7100

```python tags=[]
cubeA = Cube("../big-data/lmc-30dor-A-subcube-62-71.fits")
cubeB = Cube("../big-data/lmc-30dor-B-subcube-62-71.fits")
cubeC = Cube("../big-data/lmc-30dor-C-subcube-62-71.fits")
cubeD = Cube("../big-data/lmc-30dor-D-subcube-62-71.fits")
```

### Inspect the average spectrum for each field

```python
fig, ax = plt.subplots(figsize=(12, 6))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cube.sum(axis=(1, 2)).plot(label=label)
ax.legend()
ax.set(
    ylim=[0.5e8, 1.5e8],
)
sns.despine()
```

Again, we get WR features (I think)  around 6600, and of course the Raman wings


### Define continuum wavelength ranges

We can try the same wav ranges as we used for NGC 346

```python
wavranges = [
    (6220, 6280),
    (6400, 6450),
    (6700, 6710),
    (6760, 6810), (7015, 7050),
    (7090, 7100), 
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
    ylim=[120, 9e3],
)
sns.despine()
```

### Test the polynomial fitting with field A

This takes about a minute for each field.

```python
contA = extract.fit_continuum(
    cubeA,
    wav_ranges=wavranges,
    deg=5,
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


#### Look at some line images

```python
fig, axes = plt.subplots(3, 2, figsize=(12, 18), sharex=True, sharey=True)
(cubeA - contA).select_lambda(6300, 6310).sum(axis=0).plot(
    ax=axes[0, 0],
    vmin=-10,
    vmax=6000,
)
axes[0, 0].set_title("[O I] 6300")

(cubeA - contA).select_lambda(6312, 6322).sum(axis=0).plot(
    ax=axes[0, 1],
    vmin=-10,
    vmax=10000,
)
axes[0, 1].set_title("[S III] 6312")

(cubeA - contA).select_lambda(7002, 7012).sum(axis=0).plot(
    ax=axes[1, 0],
    vmin=-3,
    vmax=250,
)
axes[1, 0].set_title("O I 7002")

(cubeA - contA).select_lambda(6462, 6472).sum(axis=0).plot(
    ax=axes[1, 1],
    vmin=-3,
    vmax=150,
)
axes[1, 1].set_title("C II 6462")

(cubeA - contA).select_lambda(6731, 6741).sum(axis=0).plot(
    ax=axes[2, 0],
    vmin=-10,
    vmax=30000,
)
axes[2, 0].set_title("[S II] 6731")

(cubeA - contA).select_lambda(7065, 7075).sum(axis=0).plot(
    ax=axes[2, 1],
    vmin=-10,
    vmax=20000,
)
axes[2, 1].set_title("He I 7065")
```

### Now do the other fields

```python
contB = extract.fit_continuum(
    cubeB,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)
```

```python
contC = extract.fit_continuum(
    cubeC,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)
```

```python
contD = extract.fit_continuum(
    cubeD,
    wav_ranges=wavranges,
    deg=5,
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
    prefix = f"../big-data/lmc-30dor-{label}-subcube-62-71"
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


## Repeat for the next range: 7000 to 7900

```python tags=[]
cubeA = Cube("../big-data/lmc-30dor-A-subcube-70-79.fits")
cubeB = Cube("../big-data/lmc-30dor-B-subcube-70-79.fits")
cubeC = Cube("../big-data/lmc-30dor-C-subcube-70-79.fits")
cubeD = Cube("../big-data/lmc-30dor-D-subcube-70-79.fits")
```

### 70–79 Inspect the average spectrum for each field

```python
fig, ax = plt.subplots(figsize=(12, 6))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cube.sum(axis=(1, 2)).plot(label=label)
ax.legend()
ax.set(
    ylim=[0.4e8, 1.0e8],
)
sns.despine()
```

Again, we get WR features (I think)  around 7120, and some over compensation for the atmospheric asbsorption.


### 70–79 Define continuum wavelength ranges

No large ranges of clear continuum here

```python
wavranges = [
    (7015, 7040), (7090, 7100), (7150, 7160),
    (7210, 7225), (7260, 7275), 
    (7380, 7390),
    (7410, 7420), (7450, 7460),
    (7510, 7520),
    (7720, 7740),
    (7780, 7790), (7830, 7840),
    (7890, 7900),
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

### 70–79 Test the polynomial fitting with field A

This takes about a minute for each field.

```python
contA = extract.fit_continuum(
    cubeA,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)
```

#### 70–79 Inspect the results for different portions of the field.

```python
fig, ax = plt.subplots(figsize=(12, 8))
cubeA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
contA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.set(ylim=[0, 800])
```

### 70–79 Now do the other fields

```python
contB = extract.fit_continuum(
    cubeB,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)
```

```python
contC = extract.fit_continuum(
    cubeC,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)
```

```python
contD = extract.fit_continuum(
    cubeD,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)
```

### 70–79 Save the continuum-subtracted cubes

```python
csub = {}
cdiv = {}
for cube, cont, label in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    "DCBA",
):
    prefix = f"../big-data/lmc-30dor-{label}-subcube-70-79"
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
