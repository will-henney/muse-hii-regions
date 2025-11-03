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

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
# Notebook I of more line ratios from saved maps from the PZ cubes

This is adapted from 01-02-yet-more-line-ratios. It calculates the Ar IV / Ar III ratio

<!-- #endregion -->


```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
from zz_utils import get_image_raw, ROOT


# ## Ar ionization balance
#
# We will calculate the conditions at the very peak of the [Ar IV] emission.
#
# We need a very careful slection of the background (BG) and bow shock (BS) samples, since we want to make sure we are in the little triangle window where the intermediate ionization lines are not too contaminated by globule i-fronts and unrelated filaments.
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
i0, j0, w, h = 234, 200, 12, 8
bgbox = regions.RegionBoundingBox(
    iymin=j0 - h // 2,
    iymax=j0 + h // 2,
    ixmin=i0 - w // 2,
    ixmax=i0 + w // 2,
)
i0, j0, w, h = 250, 193, 8, 8
bsbox = regions.RegionBoundingBox(
    iymin=j0 - h // 2,
    iymax=j0 + h // 2,
    ixmin=i0 - w // 2,
    ixmax=i0 + w // 2,
)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
imariv.plot(ax=axes[0], vmin=0, vmax=700, colorbar="v")
imariii.plot(ax=axes[1], vmin=0, vmax=2500, colorbar="v")
(imariv / imariii).plot(ax=axes[2], vmin=0, vmax=0.45, cmap="magma_r", colorbar="v")
for ax in axes:
    bsbox.plot(ax=ax, color="w")
    bgbox.plot(ax=ax, color="w")
    ax.set(
        xlim=[200, 300],
        ylim=[100, 250],
    )
fig.tight_layout()
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
bs_slices, _ = bsbox.get_overlap_slices(imariv.shape)
bg_slices, _ = bgbox.get_overlap_slices(imariv.shape)

bs_ariv = imariv[bs_slices].data.mean()
bg_ariv = imariv[bg_slices].data.mean()
sbs_ariv = imariv[bs_slices].data.std()
sbg_ariv = imariv[bg_slices].data.std()

bs_slices, _ = bsbox.get_overlap_slices(imariii.shape)
bg_slices, _ = bgbox.get_overlap_slices(imariii.shape)
bs_ariii = imariii[bs_slices].data.mean()
bg_ariii = imariii[bg_slices].data.mean()
sbs_ariii = imariii[bs_slices].data.std()
sbg_ariii = imariii[bg_slices].data.std()
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
f"[Ar IV]: BS = {bs_ariv:.2f} +/- {sbs_ariv:.2f}, BG = {bg_ariv:.2f} +/- {sbg_ariv:.2f}"
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
f"[Ar III]: BS = {bs_ariii:.2f} +/- {sbs_ariii:.2f}, BG = {bg_ariii:.2f} +/- {sbg_ariii:.2f}"
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
BS_ariv = bs_ariv - bg_ariv
sBS_ariv = np.hypot(sbs_ariv, sbg_ariv)

BS_ariii = bs_ariii - bg_ariii
sBS_ariii = np.hypot(sbs_ariii, sbg_ariii)

BS_ar_iv_iii = BS_ariv / BS_ariii
sBS_ar_iv_iii = BS_ar_iv_iii * np.hypot(sBS_ariii / BS_ariii, sBS_ariv / BS_ariv)

f"BG-subtracted BS: [Ar IV] / [Ar III] = {BS_ar_iv_iii:.3f} +/- {sBS_ar_iv_iii:.3f}"
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
ar4 = pn.Atom("Ar", 4)
ar3 = pn.Atom("Ar", 3)
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
First T is for BG nebula. Second and third are lower and upper limits for bow shock.  See analysis of [Ar IV] temperature in `10-01` notebook
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
Ts = [12500, 14000, 16000]
e4711 = ar4.getEmissivity(tem=Ts, den=10.0, wave=4711)
e4740 = ar4.getEmissivity(tem=Ts, den=10.0, wave=4740)
e7136 = ar3.getEmissivity(tem=Ts, den=10.0, wave=7136)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e4711 + e4740
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e7136
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
So the emissivity is very similar, strangely.

Anyway, we should have:
$$
\frac {n(\mathrm{Ar^{3+}})} {n(\mathrm{Ar^{2+}})}
=
\frac{I(4711 + 4740)}{I(7136)}
\, \frac{e(7136)}{e(4711 + 4740)}
$$
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e7136 / (e4711 + e4740)
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
So the T uncertainty of +/- 1000 K would give +/- 10% uncertainty in the emissivity ratio. For the time being we take the middle value:
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
Ar3p_over_Ar2p = BS_ar_iv_iii * (e7136 / (e4711 + e4740))[1]
Ar3p_over_Ar2p
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
Or, close enough to unity.  So, at the inner edge of the bow shock we have 50% Ar++ and 50% Ar+++. We can use this to constrain the stellar spectrum if we run some Cloudy models

Now do the same, but for the BG nebula
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
Ar3p_over_Ar2p_BG = (bg_ariv / bg_ariii) * (e7136 / (e4711 + e4740))[0]
Ar3p_over_Ar2p_BG
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
So this implies 20% Ar+++ in the BG region.

However, these are both still integrals along the line of sight, so the actual variation in ionization might be larger.

For instance, the Ar+++ fraction might be higher than 50% at the inner edge.  We could do a LOS integration on the Cloudy emissivities to investigate this.

Likewise, the Ar+++ fraction in the BG might be lower than 0.2 since that could be due to contamination by the bow shock wing emission.
<!-- #endregion -->
