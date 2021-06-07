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

# Line ratios from saved moment images

In this notebook, I will work with only pre-calculated line maps, which have already been extracted from the cube and have been corrected for the sky problems (at least in theory).  The extraction process is carried out in the 03 series of notebooks.

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
import regions
import sys
import pandas as pd
import pyneb as pn

sns.set_context("talk")
sns.set_color_codes()
```

## Calculate reddening from Balmer decrement


Load the Hα and Hβ maps

```python
imha = Image("../data/ngc346-hi-6563-bin01-sum.fits")
```

```python
imhb = Image("../data/ngc346-hi-4861-bin01-sum.fits")
```

### Look at the raw Hα/Hβ ratio:

```python
fig, ax = plt.subplots(figsize=(12, 12))
(imha/imhb).plot(vmin=1.0, vmax=5.0, cmap="gray", colorbar="v")
```

So we see a lot of structure there.  In priniciple, lighter means more extinction.  This seems to be real at the bottom of the image, where we see clear signs of the foreground filament. 

But in other parts, the ratio is suspiciously well correlated with the brightness.  So we need to fix that.


### PyNeb calculation of intrinsic Balmer decrement

```python
hi = pn.RecAtom("H", 1)
```

Calculate the theoretical Balmer decrement from PyNeb. Density and temperature from Valerdi:2019a

```python
tem, den = 12500, 100
R0 = hi.getEmissivity(tem, den, wave=6563) / hi.getEmissivity(tem, den, wave=4861)
R0
```

### Look at correlation between Hα and Hβ in the faint limit

To make thinks easier, I multiply the Hb values by R0 so we have a square plot.  I zoom in on the faint parts:

```python
imax = 30000
m = imha.data < imax
m = m & (R0*imhb.data < imax)
m = m & ~imha.mask & ~imhb.mask
df = pd.DataFrame(
    {
        "ha": imha.data[m],
        "hb": R0*imhb.data[m],
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
g.axes[1, 0].plot([0, imax], [0, imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, 0.9*imax], "--", color="r")
g.fig.suptitle("Correlation between Ha and Hb brightness");
```

So, the slope is not unity, meaning the extinction is not zero.  But the intercept is not zero either. Clearly, we must fix that


### Correct zero point of Hβ map

```python
x = imha.data
y = R0*imhb.data - 5500
imax = 30000
m = (x < imax) & (y < imax)
m = m & ~imha.mask & ~imhb.mask
df = pd.DataFrame(
    {
        "ha": x[m],
        "hb": y[m],
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
g.axes[1, 0].plot([0, imax], [0, imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, 0.88*imax], "--", color="r")
g.fig.suptitle("Correlation between Ha and corrected Hb brightness");
```

That looks way better.  Expand out to brighter pixels:

```python
imax = 150000
m = (x < imax) & (y < imax)
m = m & ~imha.mask & ~imhb.mask
df = pd.DataFrame(
    {
        "ha": x[m],
        "hb": y[m],
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
g.axes[1, 0].plot([0, imax], [0, imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, 0.88*imax], "--", color="r")
g.fig.suptitle("Correlation between Ha and corrected Hb brightness");
```

### Final corrected Hα/Hβ ratio image

Now fix the Hb image and try again:

```python
hbfix = 5500 / R0
fig, ax = plt.subplots(figsize=(12, 12))
(imha/(imhb - hbfix)).plot(vmin=1.0, vmax=5.0, cmap="gray", colorbar="v")
```

That is on the same brightness map as before, and it now completely eliminates the spurious structure associate with surface brightness – hurray!

Now define some regions to take averages

```python
boxes = {
    "sw filament": regions.BoundingBox(
        iymin=30, iymax=50, ixmin=300, ixmax=330,
    ),
    "bow shock": regions.BoundingBox(
        iymin=165, iymax=205, ixmin=240, ixmax=290,
    ),
    "w filament": regions.BoundingBox(
        iymin=100, iymax=130, ixmin=25, ixmax=55,
    ),
    "c filament": regions.BoundingBox(
        iymin=195, iymax=210, ixmin=155, ixmax=195,
    ),
}
```

Plot on a better scale and show the regions:

```python
fig, ax = plt.subplots(figsize=(12, 12))
(imha/(imhb - hbfix)).plot(
    vmin=R0, vmax=1.5*R0, 
    scale="linear",
    cmap="gray_r", colorbar="v",
)

for box, c in zip(boxes.values(), "yrmgc"):
    box.plot( 
        ax=ax, 
        lw=3, 
        edgecolor=c, 
        linestyle="dashed",
#        facecolor=(1.0, 1.0, 1.0, 0.4), 
        fill=False,
    );
```

We can see some very high extinction at in the S filaments.  And some small increase in extinction in the main diagonal filament.  This is probably limited having foreground emission to some extent.

Look at average values in the sample boxes

```python
for label, box in boxes.items():
    yslice, xslice = box.slices
    ha = np.median(imha[yslice, xslice].data.data)
    hb = np.median(imhb[yslice, xslice].data.data - hbfix)
    print(f"{label}: {ha/hb:.3f}")
```

I tried mean and median, and it made very little difference.  Lowest in the bow shock region; slightly higher in the west and central filaments.  Much higher in the southwest filament. 





### The reddening law

```python
pn.RedCorr().getLaws()
```

PyNeb does not seem to have anything specifically tailored to the SMC.  The average SMC extinction law is supposedly simply $1/\lambda$.  

But, it is possible to get a SMC curve by using the "F99-like" option, which uses the curve of Fitzpatrick & Massa 1990, ApJS, 72, 163. This depends on $R_V$ and 6 other parameters (!!!).  Most of the parameters only affect the UV part of the curve, which does not concern us. 

Then, we can use the average values of $R_V$ and the other parameters, which were fit by Gordon:2003l to SMC stars. This is $R_V = 2.74 \pm 0.13$.

So here I compare that SMC curve with $1/\lambda$ and with the Clayton curve for Milky Way (but also adjusted to $R_V = 2.74$):

```python
def A_lam(wave):
    return 4861.32/wave

def my_X(wave, params=[]):
    """A_lam / E(B - V) ~ lam^-1"""
    return A_lam(wave) / (A_lam(4400) - A_lam(5500))

rc = pn.RedCorr()
rc.UserFunction = my_X
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
f, ax = plt.subplots(figsize=(10,10))
rc.plot(laws = ["user", "F99", "CCM89"], ax=ax)

ax.set(
    xlim=[4000, 9300],
    ylim=[-2.5, 1.0],
#    xlim=[4000, 7000],
#    ylim=[-1, 1],

)
```

So the Gordon curve is flatter in the blue, steeper in green, and flatter in red, as compared to $1/\lambda$. 

```python
rc = pn.RedCorr()
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
rc.law = "F99"
```

Test it out for the bow shock region:

```python
rc.setCorr(obs_over_theo=3.123 / R0, wave1=6563., wave2=4861.)
rc.E_BV, rc.cHbeta
```

And for the highest extinction region

```python
rc.setCorr(obs_over_theo=4.165 / R0, wave1=6563., wave2=4861.)
rc.E_BV, rc.cHbeta
```

So $E(B - V)$ varies from about 0.1 to about 0.35. This is similar to what is found for the stars. 


### The reddening map

We can now make a map of $E(B - V)$

```python
R = imha/(imhb - hbfix)
rc.setCorr(
    obs_over_theo=R.data / R0, 
    wave1=6563., 
    wave2=4861.
)
imEBV = R.copy()
imEBV.data = rc.E_BV
imEBV.mask = imha.mask | imhb.mask
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
imEBV.plot(
    vmin=0.0, vmax=1.0, 
    scale="linear",
    cmap="magma_r", colorbar="v",
);
```

Looks like I would expect. Check values in the boxes:

```python
for label, box in boxes.items():
    yslice, xslice = box.slices
    ebv = np.median(imEBV[yslice, xslice].data.data)
    print(f"{label}: {ebv:.3f}")
```

These seem the same as before.  But we want to eliminate extreme values. 

```python
badpix = (imEBV.data > 1.0) | (imEBV.data < 0.0)
imEBV.mask = imEBV.mask | badpix
```

Save it to a file:

```python
imEBV.write("../data/ngc346-reddening-E_BV.fits", savemask="nan")
```

Lots of regions are affected by the stellar absorption.  There are apparent increases in reddening at the position of each star.  This is not real, but is due to the photospheric absorption having more of an effect on Hb (mainly because the emission line is weaker). 

At some point, I am going to have to deal with that. But it is not an issue for the bow shock emission, since this is in an area free of stars.  We should just use the median bow shock reddening of $E(B-V) = 0.087$ so that we don't introduce any extra noise.

```python
rc.E_BV = 0.087
wavs = np.arange(4600, 9300)
Alam = rc.E_BV*rc.X(wavs)

fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(wavs, Alam)
ax.set(
    xlabel="Wavelength, $\lambda$, Å",
    ylabel="Extinction, $A_\lambda$, mag",
    title=f"Extinction curve for $E(B - V) = {rc.E_BV:.3f}$",
)
sns.despine();
```

## Calculate the [S III] temperature

```python
im6312 = Image("../data/ngc346-siii-6312-bin01-sum.fits")
im9069 = Image("../data/ngc346-siii-9069-bin01-sum.fits")
cont6312 = Image("../data/ngc346-cont-6312-mean.fits")
```

The raw ratio:

```python
fig, ax = plt.subplots(figsize=(12, 12))
(im6312/im9069).plot(vmin=-0.2, vmax=0.2, cmap="gray", colorbar="v")
```

Oh dear, that is completely dominated by the zero-point error.  But we can fix it!

```python
imax = 5000
slope = 0.1
m = im9069.data < imax
m = m & (im9069.data > -100)
m = m & (im6312.data < slope*imax)
m = m & ~im9069.mask & ~im6312.mask
df = pd.DataFrame(
    {
        "9069": im9069.data[m],
        "6312": im6312.data[m],
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
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.fig.suptitle("Correlation between [S III] 9069 and 6312 brightness");
```

So it is obvious that we need to add about 100 to the 6312 brightness.  But we can already see that the slope of the relation gets shallower for 9069 > 1000 – so fainter is hotter!

```python
imax = 2000
slope = 0.1
x = im9069.data
y = im6312.data + 105.0
m = x < imax
m = m & (x > -100)
m = m & (y < 1.2*slope*imax)
m = m & ~im9069.mask & ~im6312.mask
df = pd.DataFrame(
    {
        "9069": x[m],
        "6312": y[m],
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
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, 1.2*slope*imax], "--", color="r")
g.fig.suptitle("CORRECTED Correlation between [S III] 9069 and 6312 brightness");
```

Now we need to correct for reddening.

```python
A9069 = rc.X(9069) * imEBV
im9069c = im9069.copy()
im9069c.data = im9069.data * 10**(0.4 * A9069.data)
```

```python
A6312 = rc.X(6312) * imEBV
im6312c = im6312.copy()
im6312c.data = (im6312.data + 105) * 10**(0.4 * A6312.data)
```

```python
n = 1
fig, ax = plt.subplots(figsize=(12, 12))
(im6312c.rebin(n) / im9069c.rebin(n)).plot(vmin=0.1, vmax=0.2, cmap="magma", colorbar="v")
```

```python
n = 4
x = np.log10(im9069c.rebin(n).data)
y = np.log10(im6312c.rebin(n).data / im9069c.rebin(n).data)
m = (x > 2.0) & (x < 5.0)
m = m & (y > -1.1) & (y < -0.6)
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
);
```

Now, make a mask of EW(6312).  But first, we need to correct the zero point of the continuum.

```python
im6312_zero = -105.0
cont6312_zero = -15.0
imax = 200
x = im6312.data - im6312_zero
y = cont6312.data - cont6312_zero
m = x < imax
m = m & (x > -10)
m = m & (y < 50) & (y > -10)
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
ew6312 = 1.25*(im6312 - im6312_zero)/ (cont6312 - cont6312_zero)
ew6312.plot(vmin=0.0, vmax=30.0, scale="sqrt")
ax.contour(ew6312.data, levels=[0.5], colors="r", linewidths=3)
ax.contour(im9069.data, levels=[500.0], colors="k", linewidths=1)
```

```python
fixmask = (ew6312.data < 1.0) | (im9069.data < 500.0)
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
n = 1
fig, ax = plt.subplots(figsize=(12, 12))
(im6312c.rebin(n) / im9069c.rebin(n)).plot(vmin=0.1, vmax=0.2, cmap="magma", colorbar="v")
```

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
    plot_kws=kws, diag_kws=kws,
);
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
    r_s3_grid, T_s3_grid,
    left=np.nan, right=np.nan,
)
#imT_siii.mask = imT_siii.mask | fixmask
#imT_siii.data[imT_siii.mask] = np.nan
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
imT_siii.plot(colorbar="v", cmap="hot", vmin=11000, vmax=22000); 
```

```python
badpix = ~np.isfinite(imT_siii.data)
imT_siii.mask = imT_siii.mask | badpix
```

```python
imT_siii.write("../data/ngc346-T-siii.fits", savemask="nan")
```

The rather disappointing conclusion of this is that the [S III] temperatures do vary from about 13 to 16 kK, but they don't show anything special at the bow shock, being about 13.7 +/- 0.4 kK there. 

Average over whole FOV is 14.2 +/- 0.8 kK after smoothing to eliminate the noise contribution.  This implies $t^2 = 0.003$ in plane of sky, which is small.


## Calculate [O III]/[S III]

```python
im5007 = Image("../data/ngc346-oiii-5007-bin01-sum.fits")
```

Correct for extinction:

```python
A5007 = rc.X(5007) * imEBV
im5007c = im5007.copy()
im5007c.data = im5007.data * 10**(0.4 * A5007.data)
```

```python
median_EBV = np.median(imEBV[150:250, 200:300].data)
median_EBV
```

```python
im5007cc = im5007.copy()
im5007cc.data = im5007.data * 10**(0.4 * rc.X(5007) * median_EBV)
im9069cc = im9069.copy()
im9069cc.data = im9069.data * 10**(0.4 * rc.X(9069) * median_EBV)
```

Quick look:

```python
fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
((im5007cc - 55000)/im9069cc).plot(
    vmin=20, vmax=100, 
    colorbar="v", scale="linear",
    ax=axes[0],
);
((im5007c - 55000)/im9069c).plot(
    vmin=20, vmax=100, 
    colorbar="v", scale="linear",
    ax=axes[1],
);
```

Fix zero points:

```python
imax = 6000
slope = 45
slope2 = 35
x = im9069c.data
y = im5007c.data - 55000
m = x < imax
m = m & (x > -100)
m = m & (y < slope*imax)
m = m & (y > -100*slope)
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
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, slope2*imax], "--", color="r")
g.fig.suptitle("Correlation between [S III] 9069 and [O III] 5007 brightness");
```

```python
imax = 6000
slope = 45
slopw = 35
x = im9069cc.data
y = im5007cc.data - 55000
m = x < imax
m = m & (x > -100)
m = m & (y < slope*imax)
m = m & (y > -100*slope)
m = m & ~im9069cc.mask & ~im5007cc.mask
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
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.axes[1, 0].plot([0, imax], [0, slope2*imax], "--", color="r")
g.fig.suptitle("Correlation between [S III] 9069 and [O III] 5007 brightness");
```

```python
slope = 45
slope2 = 35
x = im9069cc.data
y = im5007cc.data - 55000
m = (x > 100)
m = m & (y > x) & (y < 100*x)
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
g.fig.suptitle("Correlation between [S III] 9069 and [O III] / [S III] ratio");
```

```python
imR_oiii_siii = (im5007cc - 55000)/im9069cc
imR_oiii_siii.write("../data/ngc346-R-oiii-5007-siii-9069.fits", savemask="nan")
```

## Calculate [O III] / Hβ

This might be better since at least it is not affected by reddening. 

```python
imax = 10000
slope = 5.
x = imhb.data - hbfix
y = im5007.data - 33000
m = x < imax
m = m & (x > -100)
m = m & (y < slope*imax)
m = m & (y > -100*slope)
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
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and [O III] 5007 brightness");
```

```python
imax = 50000
slope = 5.
x = imhb.data - hbfix
y = im5007.data - 33000
m = x < imax
m = m & (x > -100)
m = m & (y < slope*imax)
m = m & (y > -100*slope)
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
g.axes[1, 0].axvline(0.0, color="r")
g.axes[1, 0].axhline(0.0, color="r")
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and [O III] 5007 brightness");
```

```python
n = 1
x = imhb.rebin(n).data - hbfix
y = im5007.rebin(n).data - 33000
m = (x > 2000) & (y > x) & (y < 10 * x)
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
g.fig.suptitle("Correlation between Hβ 4861 and [O III] / Hβ ratio");
```

```python
imR_oiii_hb = (im5007 - 33000)/(imhb - hbfix)
imR_oiii_hb.write("../data/ngc346-R-oiii-5007-hi-4861.fits", savemask="nan")
```

```python
fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 6))
imR_oiii_siii.plot(
    vmin=20, vmax=100, 
    colorbar="v", scale="linear",
    ax=axes[0],
);
imR_oiii_hb.plot(
    vmin=3, vmax=7, 
    colorbar="v", scale="linear",
    ax=axes[1],
)
axes[0].set_title("[O III] / [S III]")
axes[1].set_title("[O III] / Hβ")
```

## Calculate He I / Hβ

Let us see if this has a hole in it where the He II is coming from.

```python
im5875 = Image("../data/ngc346-hei-5875-bin01-sum.fits")
im4922 = Image("../data/ngc346-hei-4922-bin01-sum.fits")
im5048 = Image("../data/ngc346-hei-5048-bin01-sum.fits")
```

```python
fig, axes = plt.subplots(2, 2, sharey=True, figsize=(12, 12))
im5875.plot(ax=axes[0, 0], vmin=0, vmax=5000)
im4922.plot(ax=axes[0, 1], vmin=0, vmax=300)
(imhb - hbfix).plot(ax=axes[1, 0], vmin=0, vmax=40000)
im5048.plot(ax=axes[1, 1], vmin=0, vmax=50)
axes[0, 0].set_title("He I 5875")
axes[0, 1].set_title("He I 4922")
axes[1, 0].set_title("H I 4861")
axes[1, 1].set_title("He I 5048")
```

So 5875 is 10 to 100 times brighter than the other two. And it is almost identical to Hβ! 

```python
imax = 10000
slope = 0.12
x = imhb.data - hbfix
y = im5875.data
m = x < imax
m = m & (x > -100)
m = m & (y < slope*imax)
m = m & (y > -100*slope)
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
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and He I 5875 brightness");
```

```python
imax = 100000
slope = 0.12
x = imhb.data - hbfix
y = im5875.data
m = x < imax
m = m & (x > -100)
m = m & (y < slope*imax)
m = m & (y > -100*slope)
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
g.axes[1, 0].plot([0, imax], [0, slope*imax], "--", color="r")
g.fig.suptitle("Correlation between Hβ 4861 and He I 5875 brightness");
```

```python
imR_hei_hb = im5875 / (imhb - hbfix)
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
imR_hei_hb.plot(colorbar="v", cmap="gray", vmin=0.11, vmax=0.14); 
```

```python
red_R_hei_hb = imR_hei_hb.copy()
red_R_hei_hb.data = 10**(0.4*imEBV.data*(rc.X(4861) - rc.X(5875)))
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
(imR_5875_4861 / red_R_hei_hb).plot(colorbar="v", cmap="gray", vmin=0.1, vmax=0.115); 
```

So if we correct it for reddening, then lots of spurious structure disappears.  But we are left with very little variation at all, except for at the mYSO and the top right corner, which both show low He I.


## Calculate He II / Hβ



```python
im4686 = Image("../data/ngc346-heii-4686-bin01-sum.fits")
```

```python
fig, ax = plt.subplots(figsize=(12, 12))
im4686.plot(colorbar="v", cmap="gray", vmin=0.0, vmax=300);
```

```python
imR_heii_hb = im4686 / (imhb - hbfix)

fig, ax = plt.subplots(figsize=(12, 12))
imR_heii_hb.plot(colorbar="v", cmap="gray", vmin=0.0, vmax=0.02);
```

```python
n = 2
xslice, yslice = slice(200, 300), slice(100, 250)
x = imR_heii_hb[yslice, xslice].rebin(n).data
y = imR_hei_hb[yslice, xslice].rebin(n).data / red_R_hei_hb[yslice, xslice].rebin(n).data
z = im4686[yslice, xslice].rebin(n).data
m = x < 0.03
m = m & (x > 0)
m = m & (y < 0.113)
m = m & (y > 0.103)
m = m & ~imR_heii_hb[yslice, xslice].rebin(n).mask & ~imR_hei_hb[yslice, xslice].rebin(n).mask
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
g.fig.suptitle("Correlation between He II / Hβ and He I / Hβ ratios");
```

So there is a *tiny* change in 5875/4861 from 0.109 to 0.107 as 4686/4861 increases.

```python
df["high"] = df["4686 / 4861"] > 0.003
df
```

```python
sns.histplot(
    data=df, 
    x="5875 / 4861", 
    hue="high",
    multiple="stack", 
    shrink=1.0,
    stat="probability",
    common_norm=False,
    bins=10,
)
```

## Ratio of [Ar IV] / [Ar III]

```python
im4711 = Image("../data/ngc346-ariv-4711-bin01-sum.fits")
im4740 = Image("../data/ngc346-ariv-4740-bin01-sum.fits")
im7171 = Image("../data/ngc346-ariv-7171-bin01-sum.fits")
im7263 = Image("../data/ngc346-ariv-7263-bin01-sum.fits")
im7136 = Image("../data/ngc346-ariii-7136-bin01-sum.fits")
```

```python
fig, axes = plt.subplots(2, 2, sharey=True, figsize=(12, 12))
im4711.plot(ax=axes[0, 0], vmin=0, vmax=400)
im4740.plot(ax=axes[0, 1], vmin=0, vmax=250)
(im7171 + im7263).plot(ax=axes[1, 0], vmin=0, vmax=30)
im7136.plot(ax=axes[1, 1], vmin=0, vmax=3500)
axes[0, 0].set_title("[Ar IV] 4711 + He I 4713")
axes[0, 1].set_title("[Ar IV] 4740")
axes[1, 0].set_title("[Ar IV] 7171 + 7263")
axes[1, 1].set_title("[Ar III] 7136")
```

```python
n = 8
fig, axes = plt.subplots(2, 2, sharey=True, figsize=(12, 12))
(im4711.rebin(n) / im4740.rebin(n)).plot(ax=axes[0, 0], vmin=0, vmax=4)
(im4740.rebin(n) / im7136.rebin(n)).plot(ax=axes[0, 1], vmin=0, vmax=0.15)
((im7171.rebin(n) + im7263.rebin(n)) / (im4740.rebin(n) + im4711.rebin(n))).plot(ax=axes[1, 0], vmin=0, vmax=0.08)
im7136.rebin(n).plot(ax=axes[1, 1], vmin=0, vmax=3500)
axes[0, 0].set_title("([Ar IV] 4711 + He I 4713) / [Ar IV] 4740")
axes[0, 1].set_title("[Ar IV] 4740 / [Ar III] 7136")
axes[1, 0].set_title("[Ar IV] (7171 + 7263) / (4711 + 4740)")
axes[1, 1].set_title("[Ar III] 7136")
```

Now we must subtract the He I line!

```python
hei = pn.RecAtom("He", 1)
```

```python
dens = [50, 100, 200]
tems = [13000, 18000]
e4713 = hei.getEmissivity(tems, dens, wave=4713)
e5876 = hei.getEmissivity(tems, dens, wave=5876)
e4713 / e5876
```

There is a slight temperature dependence, but almost no density dependence if we use the 5876 line.  This is probably the best because it has good signal to noise.  

We can assume that the He I temperature is the same as the [S III] temperature.

But we will check the other lines as well. 

```python
e4922 = hei.getEmissivity(tems, dens, wave=4922)
e5048 = hei.getEmissivity(tems, dens, wave=5048)
e4713 / e4922, e4713 / e5048
```

The 4922 has the same T-dependence as 5876, just 10 times weaker. The 5048 has a constant ratio, but it is so weak that we cannot use it.  So 5876 it is ...


### Average values of T and reddening to use in the corrections

We would introduce too much noise by using the pixel-by-pixel values of $T$ and $E(B - V)$, so we will construct an average value by using the He I brightenss as a weight, but masking out the mYSO


```python
def trim_edges(im, m):
    """Trim m pixels of each edge of image in place by setting mask"""
    im.mask[:m, :] = True
    im.mask[-m:, :] = True
    im.mask[:, :m] = True
    im.mask[:, -m:] = True
    return None
```

```python
im_hei_weight = im5875.copy()
im_hei_weight.mask[140:157, 110:141] = True
im_hei_weight.mask[94:104, 55:65] = True
trim_edges(im_hei_weight, 10)
im_hei_weight.data.mask = im_hei_weight.mask
fig, ax = plt.subplots(figsize=(10, 10))
im_hei_weight.plot(cmap="gray_r")
ax.set_title("Weight mask for He I emission");
```

That looks OK. Now calculate some averages:

```python
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
im1 = im_hei_weight.copy()
im1.data = im_hei_weight.data * imEBV.data
im1.mask = im_hei_weight.mask | imEBV.mask
im1.plot(ax=axes[0])
im2 = im_hei_weight.copy()
im2.data = im_hei_weight.data * imT_siii.data
im2.mask = im_hei_weight.mask | imT_siii.mask
im2.plot(ax=axes[1])
```

```python
avHe_EBV = np.average(imEBV.data, weights=im_hei_weight.data)
avHe_Tsiii = np.average(imT_siii.data, weights=im_hei_weight.data)
f"He I brightness-weighted averages: E(B-V) = {avHe_EBV:.2f}, T = {avHe_Tsiii/1000:.2f} kK"
```

```python
median_EBV
```

The `avHe_Tsiii` looks good. But to be honest, I am a bit suspicious of the `avHe_EBV` reddening, since there are lots of anomalous spots of high $E(B-V)$ that correspond to stars (presumably underlying stellar absorption affecting the Balmer decrement).  

I couls use the median instead, but I have ended up using the weighted one after all, since we seem to be oversubtracting if anything.

```python
avHe_reddening_4713_5876 = 10**(0.4 * avHe_EBV * (rc.X(4713) - rc.X(5876)))
median_reddening_4713_5876 = 10**(0.4 * median_EBV * (rc.X(4713) - rc.X(5876)))
avHe_reddening_4713_5876, median_reddening_4713_5876
```

```python
dens = 100.0
avHe_e4713_5876 = (
    hei.getEmissivity(avHe_Tsiii, dens, wave=4713) 
    / hei.getEmissivity(avHe_Tsiii, dens, wave=5876)
)
avHe_e4713_5876
```

Now do the correction by faking the 4713 line and subtracting it:

```python
#im_fake_4713 = (avHe_e4713_5876 / avHe_reddening_4713_5876) * im5875
im_fake_4713 = 0.0255 * im5875
im4711c = im4711 - im_fake_4713
```

Make a common minimal mask to use for all the [Ar IV] lines, which we will then combine with a brightness-based mask for the weaker lines and ratios:

```python
cont4686 = Image("../data/ngc346-cont-4686-mean.fits")
```

I need to decide how bright a star needs to be before I mask out that bit of the image. 5000 in the `cont4686` image seems a reasonable value. 

```python
fig, ax = plt.subplots(figsize=(10, 10))
cont4686.plot(colorbar="v", vmin=0, vmax=1e5, scale="sqrt")
ax.contour(cont4686.data, levels=[5e3], colors="r")
```

```python
im_ariv_sum = im4711c + im4740
trim_edges(im_ariv_sum, 12)
im_ariv_sum.mask[78:88, 190:199] = True
im_ariv_sum.mask[234:236, 266:271] = True
im_ariv_sum.mask[81:84, 52:55] = True
im_ariv_sum.mask = im_ariv_sum.mask | (cont4686.data > 1.5e3)
im_ariv_sum.mask = im_ariv_sum.mask | (im_ariv_sum.data > 650)
im_ariv_sum.mask = im_ariv_sum.mask | (im_ariv_sum.data < -250)
fig, ax = plt.subplots(figsize=(10, 10))
im_ariv_sum.rebin(1).plot(
    colorbar="v", 
    vmin=-10, vmax=600, 
    cmap="gray_r", 
    scale="sqrt"
);
```

That is looking good.  Apply the mask to all the other images

```python
for im in im4711c, im4740, im7171, im7263:
    im.mask = im.mask | im_ariv_sum.mask
```

```python
ariv_R1 = im4711c / im4740
ariv_R1.mask = ariv_R1.mask | (im_ariv_sum.data < 250)
ariv_R3_plus_R4 = (im7171 + im7263) / (im4740 + im4711c)
ariv_R3_plus_R4.mask = ariv_R3_plus_R4.mask | (im_ariv_sum.data < 300)
```

```python
fig, axes = plt.subplots(2, 2, figsize=(12, 12), sharex="row", sharey="row")
im4711c.plot(ax=axes[0, 0], vmin=-20, vmax=300, colorbar="v")
im4740.plot(ax=axes[0, 1], vmin=-20, vmax=300/1.4, colorbar="v")
n = 8
#(im4711c.rebin(n) / im4740.rebin(n)).plot(ax=axes[1, 0], vmin=0, vmax=2, colorbar="v")
#(
#    (im7171.rebin(n) + im7263.rebin(n)) 
#     / (im4740.rebin(n) + im4711c.rebin(n))
#).plot(ax=axes[1, 1], vmin=0, vmax=0.08, colorbar="v")
ariv_R1.rebin(n).plot(ax=axes[1, 0], vmin=0, vmax=2, cmap="mako_r", colorbar="v")
ariv_R3_plus_R4.rebin(n).plot(ax=axes[1, 1], vmin=0, vmax=0.08, cmap="inferno", colorbar="v")
```

```python
n = 4
xslice, yslice = slice(200, 300), slice(100, 250)
ratio = 1.4
xmax = 250
ymax = ratio*xmax
x = im4740[yslice, xslice].rebin(n).data
y = im4711c[yslice, xslice].rebin(n).data 
z = im_ariv_sum[yslice, xslice].rebin(n).data
m = (x > -100.0) & (y > -100.0) & (x < xmax) & (y < ymax)
#m = m & (x > 0)
#m = m & (y < 0.113)
#m = m & (y > 0.103)
m = m & ~im4740[yslice, xslice].rebin(n).mask & ~im4711[yslice, xslice].rebin(n).mask
df = pd.DataFrame(
    {
        "4740": x[m],
        "4711": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        #weights=z[m], 
        bins=64//n),
    diag_kws=dict(
        #weights=z[m], 
        bins=128//n),
)
g.axes[1, 0].plot([0, xmax], [0, ymax])
g.fig.suptitle("Correlation between [Ar IV] 4711 and 4740");
```

```python
n = 4
xslice, yslice = slice(200, 300), slice(100, 250)
ratio = 1.35
xmax = 300
ymax = ratio * xmax
x = im4740[yslice, xslice].rebin(n).data
y = im4711c[yslice, xslice].rebin(n).data 
z = im_ariv_sum[yslice, xslice].rebin(n).data
m = (x > -100.0) & (y > -100.0) & (x < xmax) & (y < ymax)
y = y / x
m = m & (y < 2.5*ratio) & (y > 0.0)
#m = m & (x > 0)
#m = m & (y < 0.113)
#m = m & (y > 0.103)
m = m & ~im4740[yslice, xslice].rebin(n).mask & ~im4711[yslice, xslice].rebin(n).mask
df = pd.DataFrame(
    {
        "4740": x[m],
        "4711 / 4740": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(weights=z[m], bins=200//n),
    diag_kws=dict(weights=z[m], bins=100//int(np.sqrt(n))),
)
g.axes[1, 0].axhline(ratio, linestyle="dashed", color="k", linewidth=2)
g.axes[1, 1].axvline(ratio, linestyle="dashed", color="k", linewidth=2)
g.fig.suptitle("Correlation between [Ar IV] 4740 and 4711 / 4740");
```

<!-- #raw -->
So the density-sensitive ratio is $R_1 = 1.35 \pm 0.1$ 
<!-- #endraw -->

```python
ariv = pn.Atom("Ar", 4)
```

```python
ariv.getTemDen([1.30, 1.35, 1.40], tem=17500, wave1=4711, wave2=4740)
```

```python
ariv.getSources()
```

```python
n = 4
xslice, yslice = slice(200, 300), slice(100, 250)
ratio = 0.025
xmax = 600
ymax = ratio*xmax
x = im_ariv_sum[yslice, xslice].rebin(n).data
y = im7171[yslice, xslice].rebin(n).data 
z = im_ariv_sum[yslice, xslice].rebin(n).data
m = (x > -0.5*xmax) & (y > -2*ymax) & (x < xmax) & (y < 3*ymax)
#m = m & (x > 0)
#m = m & (y < 0.113)
#m = m & (y > 0.103)
m = m & ~im_ariv_sum[yslice, xslice].rebin(n).mask & ~im7171[yslice, xslice].rebin(n).mask
df = pd.DataFrame(
    {
        "4711 + 4740": x[m],
        "7171": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=z[m], 
        bins=128//n),
    diag_kws=dict(
        weights=200 + z[m], 
        bins=128//n),
)
g.axes[1, 0].plot([0, xmax], [0, ymax])
g.fig.suptitle("Correlation between [Ar IV] 4711+40 and 7171");
```

```python
n = 4
xslice, yslice = slice(200, 300), slice(100, 250)
ratio = 0.025
xmax = 600
ymax = ratio*xmax
x = im_ariv_sum[yslice, xslice].rebin(n).data
y = im7263[yslice, xslice].rebin(n).data 
z = im_ariv_sum[yslice, xslice].rebin(n).data
m = (x > -0.5*xmax) & (y > -2*ymax) & (x < xmax) & (y < 3*ymax)
#m = m & (x > 0)
#m = m & (y < 0.113)
#m = m & (y > 0.103)
m = m & ~im_ariv_sum[yslice, xslice].rebin(n).mask & ~im7171[yslice, xslice].rebin(n).mask
df = pd.DataFrame(
    {
        "4711 + 4740": x[m],
        "7263": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=z[m], 
        bins=128//n),
    diag_kws=dict(
        weights=200 + z[m], 
        bins=128//n),
)
g.axes[1, 0].plot([0, xmax], [0, ymax])
g.fig.suptitle("Correlation between [Ar IV] 4711+40 and 7263");
```

```python
n = 4
xslice, yslice = slice(200, 300), slice(100, 250)
ratio = 0.025
xmax = 600
ymax = ratio*xmax
x = im_ariv_sum[yslice, xslice].rebin(n).data
y = im7171[yslice, xslice].rebin(n).data 
z = im_ariv_sum[yslice, xslice].rebin(n).data
m = (x > -0.5*xmax) & (y > -2*ymax) & (x < xmax) & (y < 3*ymax)
y = y / x
m = m & (y > -ratio) & (y < 3*ratio)
#m = m & (x > 0)
#m = m & (y < 0.113)
#m = m & (y > 0.103)
m = m & ~im_ariv_sum[yslice, xslice].rebin(n).mask & ~im7171[yslice, xslice].rebin(n).mask
df = pd.DataFrame(
    {
        "4711 + 4740": x[m],
        "7171 / (4711 + 4740)": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=z[m], 
        bins=128//n),
    diag_kws=dict(
        weights=200 + z[m], 
        bins=128//n),
)
g.fig.suptitle("Correlation between [Ar IV] 4711+40 and 7171 / (4711 + 4740)");
```

```python
n = 4
xslice, yslice = slice(200, 300), slice(100, 250)
ratio = 0.025
xmax = 600
ymax = ratio*xmax
x = im_ariv_sum[yslice, xslice].rebin(n).data
y = im7263[yslice, xslice].rebin(n).data 
z = im_ariv_sum[yslice, xslice].rebin(n).data
m = (x > -0.5*xmax) & (y > -2*ymax) & (x < xmax) & (y < 3*ymax)
y = y / x
m = m & (y > -ratio) & (y < 3*ratio)
#m = m & (x > 0)
#m = m & (y < 0.113)
#m = m & (y > 0.103)
m = m & ~im_ariv_sum[yslice, xslice].rebin(n).mask & ~im7263[yslice, xslice].rebin(n).mask
df = pd.DataFrame(
    {
        "4711 + 4740": x[m],
        "7263 / (4711 + 4740)": y[m],
    }
)
g = sns.pairplot(
    df,
    kind="hist",
    height=4,
    corner=True,
    plot_kws=dict(
        weights=z[m], 
        bins=128//n),
    diag_kws=dict(
        weights=200 + z[m], 
        bins=128//n),
)
g.fig.suptitle("Correlation between [Ar IV] 4711+40 and 7263 / (4711 + 4740)");
```

```python

```
