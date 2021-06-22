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

# NGC 346 stellar velocities from He II 5041.5 Å absorption line

Castro+ (2018) used this line for Tarantula stars.  It has the advantage that it is strong absorption in O stars and there is no nebular contamination.

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
import masktools

sns.set_context("talk")
sns.set_color_codes()
```

```python
cube = Cube("../big-data/ngc346-5300-6100-cube-contsub.fits")
cont = Cube("../big-data/ngc346-5300-6100-cube-cont.fits")
```

```python
win5411 = cube.select_lambda(5400, 5430)
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
win5411.sum(axis=0).plot(
    vmin=-10000, 
    vmax=100, 
    cmap="gray",
    colorbar="v"
);
```

```python
mom5411 = moments.find_moments(win5411)
```

```python
cont5411 = cont.select_lambda(5400, 5430).mean(axis=0)
dlam, _, _ = cube.get_step()
ew5411 = -dlam * mom5411[0] / cont5411
```

Mask out region where there are no stars

```python
masktools.trim_edges(ew5411, 10)
masktools.trim_edges(cont5411, 10)
ew5411.mask = ew5411.mask | (cont5411.data < 3000)
```

```python
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
ew5411.plot(
    ax=axes[0],
    vmin=0, 
    vmax=1.2, 
    cmap="viridis",
    colorbar="v"
)
cont5411.plot(
    ax=axes[1],
    vmin=0, 
    vmax=1e5, 
    cmap="gray",
    scale="log",
    colorbar="v"
)
axes[1].contour(cont5411.data, levels=[3000], colors="r")
fig.tight_layout();
```

```python

```

```python
wav0 = 5411.52 # from Atomic Line List
vel = 3e5 * (mom5411[1] - wav0) / wav0
vel.mask = vel.mask | (cont5411.data < 3000) | (ew5411.data < 0.2)
masktools.trim_edges(vel, 10)
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
vel.rebin(1).plot(
    vmin=100, 
    vmax=240, 
    cmap="seismic",
    colorbar="v",
    use_wcs=True,
);
ax.contour(cont5411.data, levels=[3000], colors="g");
```

```python
seg = cont5411.segment(minsize=4, background=5000)
```

```python
bseg = [_ for _ in seg if _.data.max() > 10000]

len(bseg)
```

```python
from astropy.wcs import WCS
```

```python
ny, nx = 5, 4 
fig = plt.figure(figsize=(12, 12))
for i, subim in enumerate(bseg):
    w = WCS(subim.get_data_hdu().header).celestial
    ax = fig.add_subplot(ny, nx, i + 1, projection=w)
    subim.plot(ax=ax, vmin=0, vmax=20000, use_wcs=True)
fig.tight_layout()
```

```python
import pandas as pd
```

```python
df = pd.DataFrame([_.peak() for _ in bseg])
df
```

```python
fig = plt.figure(figsize=(12, 12))
ax = ax = fig.add_subplot(
    1, 1, 1, 
    projection=vel.wcs.wcs,
)
vel.rebin(1).plot(
    ax=ax,
    vmin=100, 
    vmax=240, 
    cmap="seismic",
    colorbar="v",
)
ax.scatter(
    "x", "y", 
    data=df, 
    marker="o", color="k", s=500, facecolor="none",
    transform=ax.get_transform("icrs"),
);
```

```python
source = bseg[0]
```

```python
source.get_range()
```

```python
x, y = df["x"][10], df["y"][10]
vel.wcs.sky
```

```python
x, y
```

```python
type(vel.wcs.wcs)
```

```python
vel.wcs.wcs.world_to_pixel_values(x, y)
```

```python

```
