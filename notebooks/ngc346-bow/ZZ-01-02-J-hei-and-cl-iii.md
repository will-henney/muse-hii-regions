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
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
## Can we get a He I density?
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
dgrid = [1.0, 10.0, 100.0, 1000.0]
T0 = [11000, 13000, 18000]
he1.getEmissivity(tem=T0, den=dgrid, wave=5876) / he1.getEmissivity(
    tem=T0, den=dgrid, wave=6678
)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
he1.getEmissivity(tem=T0, den=dgrid, wave=4922) / he1.getEmissivity(
    tem=T0, den=dgrid, wave=5876
)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
he1.getEmissivity(tem=T0, den=dgrid, wave=5048) / he1.getEmissivity(
    tem=T0, den=dgrid, wave=5876
)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
im5048 = Image(str(datadir / "ngc346-hei-5048-bin01-sum.fits"))
im5876 = Image(str(datadir / "ngc346-hei-5875-bin01-sum.fits"))
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
n = 1
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
im5048.rebin(n).plot(vmin=-10, vmax=60, ax=axes[0, 0], colorbar="v")
im5876.rebin(n).plot(vmin=-10, vmax=3000, ax=axes[0, 1], colorbar="v")
imcont.rebin(n).plot(vmin=0, vmax=1e4, ax=axes[1, 0], colorbar="v")
(im5048.rebin(n) / im5876.rebin(n)).plot(
    ax=axes[1, 1],
    vmin=0.0,
    vmax=0.03,
    cmap="magma",
    colorbar="v",
)
for ax in axes.flat:
    bsbox.plot(ax=ax, color="w")
    bgbox.plot(ax=ax, color="w")
    ax.set(
        xlim=[200, 300],
        ylim=[100, 250],
    )
fig.tight_layout()
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
yyslice = slice(164, 204)  # original
# yyslice = slice(160, 210) # broader
# yyslice = slice(170, 200) # narrower
# yyslice = slice(180, 200) # top half ultra narrow
hei_5048_profile = make_profile(im5048)
hei_5876_profile = make_profile(im5876)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(hei_profile)
pos = (np.arange(nx) - ix0) * 0.2

ax.plot(pos, 0.01 * hei_5876_profile / np.median(hei_5876_profile), label="He I", lw=4)
ax.plot(pos, hei_5048_profile / hei_5876_profile, label="5048 / 5875", lw=3)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=3, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Surface brightness",
    xlim=[-12, 22],
    ylim=[0, 0.03],
)
sns.despine()
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
bs_5048 = im5048[bs_slices].data.mean()
bs_5876 = im5876[bs_slices].data.mean()
bs_5048 / bs_5876
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
bg_5048 = im5048[bg_slices].data.mean()
bg_5876 = im5876[bg_slices].data.mean()
bg_5048 / bg_5876
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
These are all way lower than the theoretical values for reasonable temperatures.  Maybe the 5048 line is affected by underlying stellar absorption.
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
# yyslice = slice(164, 204) # original
# yyslice = slice(160, 210) # broader
# yyslice = slice(170, 200) # narrower
yyslice = slice(180, 200)  # top half ultra narrow
n_sii_profile = make_profile(im_n_sii)
T_siii_profile = make_profile(im_T_siii)
sii_profile = make_profile(im6731)
im9069 = Image(str(datadir / "ngc346-siii-9069-bin01-sum.fits"))
siii_profile = make_profile(im9069)
```

```python pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 4))
ix0 = 227.5
nx = len(sii_profile)
pos = (np.arange(nx) - ix0) * 0.2
pos2 = (np.arange(len(siii_profile)) - ix0) * 0.2

ax.plot(
    pos, 0.01 * n_sii_profile, ds="steps-mid", label="$n$([S II]) / 100 cm$^{-3}$", lw=2
)
ax.plot(
    pos2, 0.0001 * T_siii_profile, ds="steps-mid", label="$T$([S III]) / 10,000 K", lw=2
)
ax.plot(
    pos, 1.0 * sii_profile / np.median(sii_profile), label="[S II] brightness", lw=3
)
ax.plot(
    pos2, 1.8 * siii_profile / np.median(siii_profile), label="[S III] brightness", lw=3
)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=2, loc="upper left", fontsize="x-small")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[0, 3.95],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-sii-siii-ne-te.pdf")
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} tags=["temperature"] -->
So the [S III] temperature is significantly larger than the [O III] temperature from the nebula. It is around 14000 K and has a drop towards W 3 coming from the east side, and then a step and a very constant region.  But the step occurs before the bow shock, so is probably unrelated.
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(sii_profile)
pos = (np.arange(nx) - ix0) * 0.2
pos2 = (np.arange(len(siii_profile)) - ix0) * 0.2

ax.plot(pos, 0.01 * n_sii_profile, label="$n$([S II]) / 100 cm$^{-3}$", lw=4)
ax.plot(pos2, 0.0001 * T_siii_profile, label="$T$([S III]) / 10,000 K", lw=3)
ax.plot(pos2, 1.0 * sii_profile / siii_profile, label="[S II] / [S III]", lw=3)
ax.plot(pos2, 4.0 * oi_profile / sii_profile, label="[O I] / [S II]", lw=3)


ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=2, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[0, 2.1],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-sii-siii-ratio-ne-te.pdf")
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
imhei_c = Image(str(datadir / "ngc346-hei-5875-correct.fits"))
imhi_c = Image(str(datadir / "ngc346-hi-4861-correct.fits"))
imheii_c = Image(str(datadir / "ngc346-heii-4686-correct.fits"))
imcont2 = Image(str(datadir / "ngc346-cont-4686-mean.fits"))
```


```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
m = imcont2.data > 500.0
for im in imhei_c, imheii_c, imhi_c:
    im.mask = im.mask | m
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
# yyslice = slice(164, 204) # original
yyslice = slice(160, 210)  # broader
# yyslice = slice(170, 200) # narrower
# yyslice = slice(180, 200) # top half ultra narrow
hei_c_profile = make_profile(imhei_c)
hi_c_profile = make_profile(imhi_c)
heii_c_profile = make_profile(imheii_c)
cont_profile = make_profile(imcont2)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(10, 4))
ix0 = 227.5
nx = len(hei_profile)
pos = (np.arange(nx) - ix0) * 0.2

ax.plot(
    pos,
    heii_c_profile / hi_c_profile,
    ds="steps-mid",
    label="He II λ4686 / H I λ4861",
    lw=2,
)
ax.plot(
    pos,
    hei_c_profile / hi_c_profile - 0.10,
    ds="steps-mid",
    label="(He I λ5875 / H I λ4861) – 0.1",
    lw=2,
)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=1, loc="upper right", fontsize="x-small")
ax.set_yticks([0.0, 0.005, 0.010])
ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[-0.005, 0.015],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-he-ratios.pdf")
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(hei_profile)
pos = (np.arange(nx) - ix0) * 0.2

ax.plot(
    pos,
    heii_c_profile / hi_c_profile,
    ds="steps-mid",
    label="He II λ4686 / H I λ4861",
    lw=2,
)
ax.plot(
    pos,
    hei_c_profile / hi_c_profile - 0.10,
    ds="steps-mid",
    label="(He I λ5875 / H I λ4861) – 0.1",
    lw=2,
)
ax.plot(
    pos,
    0.012 * cont_profile / np.max(cont_profile[np.abs(pos) < 10]),
    ds="steps-mid",
    label="continuum",
    lw=2,
    color="k",
)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=1, loc="upper right")
ax.set_yticks([0.0, 0.005, 0.010])
ax.set(
    xlabel="Offset west from W 3, arcsec",
    xlim=[-12, 22],
    ylim=[-0.005, 0.015],
)
sns.despine()
```

These show the continuum as well, so we can see the PSF width of about 5 pixels, which is less thn the shell thickness of about 10 pixels.

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
im5518 = Image(str(datadir / "ngc346-cliii-5518-bin01-sum.fits"))
im5538 = Image(str(datadir / "ngc346-cliii-5538-bin01-sum.fits"))
imha = Image(str(datadir / "ngc346-hi-6563-bin01-sum.fits"))
imcont = Image(str(datadir / "ngc346-cont-4686-mean.fits"))
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
cl3 = pn.Atom("Cl", 3)
Rlo = cl3.getLowDensRatio(wave1=5518, wave2=5538)
Rhi = cl3.getHighDensRatio(wave1=5518, wave2=5538)
Rlo, Rhi
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
m = imcont.data > 3e2
im5518.mask = im5518.mask | m
im5538.mask = im5538.mask | m
trim_edges(im5518, 20)
trim_edges(im5538, 20)
trim_edges(imcont, 20)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
shift5538 = 15.0
shift5518 = 23.0
im5538.data += shift5538
im5518.data += shift5518
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
n = 16
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
im5538.rebin(n).plot(vmin=-10, vmax=120, ax=axes[0, 0], colorbar="v")
im5518.rebin(n).plot(vmin=-10, vmax=120, ax=axes[0, 1], colorbar="v")
imcont.rebin(n).plot(vmin=0, vmax=1e4, ax=axes[1, 0], colorbar="v")
(im5518.rebin(n) / im5538.rebin(n)).plot(
    ax=axes[1, 1],
    vmin=0.0,
    vmax=2.0,
    cmap="gray",
    colorbar="v",
)
fig.tight_layout()
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
n = 8

imx = im5538.rebin(n)
imy = (im5518 / Rlo).rebin(n)

imin, imax = -10, 100
m = (imx.data < imax) & (imy.data < imax)
m = m & (imx.data > imin) & (imy.data > imin)
m = m & (imy.data / imx.data > -2) & (imy.data / imx.data < 4.0)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "5538": imx.data[m],
        "5518": imy.data[m],
        "5518 / 5538": imy.data[m] / imx.data[m],
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
g.axes[1, 0].plot([imin, imax], [imin, imax], "--", color="r")
g.axes[1, 0].plot([imin, imax], [imin, imax], "--", color="r")
g.fig.suptitle("Correlation between [Cl III] 5538 and 5518 brightness")
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
n = 32

imx = imha.rebin(n)
imy = im5518.rebin(n) / im5538.rebin(n)

imin, imax = -10, 100
m = (imy.data > 0.5) & (imy.data < 2.0)
m = m & (imx.data < 1.5e5)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": imx.data[m],
        "5518 / 5538": imy.data[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=50,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=50,
    ),
)
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(1.5, color="r")
g.fig.suptitle("Correlation between [Cl III] 5538, 5518 sum and ratio")
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
xslice, yslice = slice(230, 300), slice(144, 245)

n = 2

imx = imha[yslice, xslice].rebin(n)
imy = im5518[yslice, xslice].rebin(n) / im5538[yslice, xslice].rebin(n)

imin, imax = -10, 100
m = (imy.data > 0.15) & (imy.data < 15.0)
m = m & (imx.data < 1.5e5)
m = m & ~imx.mask & ~imy.mask
df = pd.DataFrame(
    {
        "6563": np.log10(imx.data[m]),
        "5518 / 5538": np.log10(imy.data[m]),
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=imx.data[m],
        bins=30,
    ),
    diag_kws=dict(
        weights=imx.data[m],
        bins=30,
    ),
)
g.axes[1, 0].axhline(np.log10(1.5), color="r")
g.axes[1, 1].axvline(np.log10(1.5), color="r")


g.fig.suptitle("Correlation between [Cl III] 5538 / 5518 ratio and Ha brightness")
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
m = ~im5518[yslice, xslice].mask
npix = m.sum()
y = im5518[yslice, xslice].data[m]
x = im5538[yslice, xslice].data[m]
w = imha[yslice, xslice].data[m]

# unweighted
ym = y.mean()
xm = x.mean()
R1 = ym / xm
dR1 = np.sqrt((y.var() / ym**2 + x.var() / xm**2) / (npix - 1))

# weighted
ymw = np.average(y, weights=w)
xmw = np.average(x, weights=w)
R2 = ymw / xmw
dR2 = np.sqrt(
    (
        np.average(((y - ymw) / ymw) ** 2, weights=w)
        + np.average(((x - ymw) / xmw) ** 2, weights=w)
    )
    / (npix - 1)
)
f"Unweighted R = {R1:.4f} +/- {dR1:.4f}; Weighted R = {R2:.4f} +/- {dR2:.4f}"
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
cl3.getTemDen(1.44, tem=12000, wave1=5518, wave2=5538)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e5518 = cl3.getEmissivity(tem=12000, den=[1, 10, 100, 1000], wave=5518)
e5538 = cl3.getEmissivity(tem=12000, den=[1, 10, 100, 1000], wave=5538)
e5518 / e5538
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e5518 = cl3.getEmissivity(tem=15000, den=[1, 10, 100, 1000], wave=5518)
e5538 = cl3.getEmissivity(tem=15000, den=[1, 10, 100, 1000], wave=5538)
e5518 / e5538
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
e5518 = cl3.getEmissivity(tem=8000, den=[1, 10, 100, 1000], wave=5518)
e5538 = cl3.getEmissivity(tem=8000, den=[1, 10, 100, 1000], wave=5538)
e5518 / e5538
```

<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->
Check the 2 sigma lower limit
<!-- #endregion -->

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
cl3.getTemDen(R2 - 5 * dR2, tem=15000, wave1=5518, wave2=5538)
```

```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
cl3.getTemDen(R2 - 5 * dR2, tem=1000, wave1=5518, wave2=5538)
```
