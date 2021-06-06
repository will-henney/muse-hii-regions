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

These seem the same as before.  Save it to a file:

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
n = 1
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
);
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

```python

```
