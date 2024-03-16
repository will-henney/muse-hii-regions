---
jupyter:
  jupytext:
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Cloudy models of NGC 346 bow shock around Walborn 3

```python
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import cmasher as cmr
import astropy.units as u
from pathlib import Path
import sys
from cloudytab import cloudytab
```

```python
cloudytab?
```

```python
sns.set_context("talk")
sns.set_color_codes()
```

```python
ROOT = Path.cwd().parent.parent
m1 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n010")
m2 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n010-p")
m3 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n030-p")
```

```python
m1.data.keys()
```

```python

```

## Optical lines

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12))


# colnames = m.data["emis"].colnames[1:]

embands = [
 'He 2 4685.70A',
 'Ar 4 4740.12A',
 'Ne 3 3868.76A',
 'O  3 5006.84A',
 'Blnd 5875.66A',
 'Ar 3 7135.79A',
 'H  1 4861.33A',
 'Ca B 6562.82A',
 'O  2 7319.99A',
]

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.neon', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m1, m2, m3], axes):
    radius = m.data["rad"]["radius"] * u.cm.to(u.pc) 
    hb = m.data["emis"]['H  1 4861.33A'] 
    for emband, color in zip(embands, colors):
        em = m.data["emis"][emband] 
        ax.plot(radius, em / hb.max(), label=emband, color=color)
    ax.set(
        yscale="log",
        ylim=[0.001, 10.1],
        xlabel="Radius, pc",
        ylabel="Emissivity",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant density, n = 10")
axes[1].set_title("Constant pressure, n = 10")
axes[2].set_title("Constant pressure, n = 30")
sns.despine()
fig.tight_layout();
```

## Infrared lines and bands

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12))


# colnames = m.data["emis"].colnames[1:]

embands = [
 "Ne 3 15.5509m",
 "Ne 2 12.8101m",
 "S  4 10.5076m",
 "S  3 18.7078m",
 "S  3 33.4704m",
 "Si 2 34.8046m",
]

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.neon', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m1, m2, m3], axes):
    radius = m.data["rad"]["radius"] * u.cm.to(u.pc) 
    hb = m.data["emis"]['H  1 4861.33A'] 
    for emband, color in zip(embands, colors):
        em = m.data["emis"][emband] 
        ax.plot(radius, em / hb.max(), label=emband, color=color)
    ax.set(
        yscale="log",
        ylim=[0.001, 10.1],
        xlabel="Radius, pc",
        ylabel="Emissivity",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant density, n = 10")
axes[1].set_title("Constant pressure, n = 10")
axes[2].set_title("Constant pressure, n = 30")
sns.despine()
fig.tight_layout();
```

## Physical variables


Class to make dataframe of all cloudy files

```python
class C:
    """Dictionary of pandas dataframes"""
    def __init__(self, d):
        for k, v in d.items():
            setattr(self, k, v.to_pandas())
```

```python
m1.p = C(m1.data)
m2.p = C(m2.data)
m3.p = C(m3.data)
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12))
for m, ax in zip([m1, m2, m3], axes):
    m.radius = m.p.rad.radius * u.cm.to(u.pc) 
    ax.plot(m.radius, m.p.ovr.eden, label="eden")
    ax.plot(m.radius, m.p.ovr.hden * m.p.ovr.HII, label="H II")
    ax.plot(m.radius, m.p.ovr.hden * m.p.ovr.HeII, label="He II")
    ax.plot(m.radius, m.p.ovr.hden * m.p.ovr.HeIII, label="He III")
    ax.plot(m.radius, m.p.ovr.hden * m.p.Ar["Ar+3"], label="Ar IV")
    ax.plot(m.radius, m.p.ovr.hden * m.p.Ne["Ne+2"], label="Ne III")
    ax.plot(m.radius, m.p.ovr.hden * m.p.O["O+2"], label="O III")
    ax.plot(m.radius, m.p.ovr.hden * m.p.S["S+2"], label="S III")
    ax.plot(m.radius, m.p.ovr.hden * m.p.S["S+3"], label="S IV")
    ax.plot(m.radius, 0.001 * m.p.ovr.Te, label="Te, kK")
axes[0].legend(ncol=3)
axes[0].set_title("Constant density, n = 10")
axes[1].set_title("Constant pressure, n = 10")
axes[2].set_title("Constant pressure, n = 30")
axes[-1].set(
    xlabel="Radius, pc",
)
sns.despine()
fig.tight_layout();
```

```python
m.p.Ar
```

```python
m.p.cont
```

```python
import astropy.constants as const
```

```python
wavnorm = (const.h * const.c / u.rydberg).to(u.micron)
freqnorm = (u.rydberg / const.h).to(u.Hz)
sednorm = (u.erg / u.s) / const.L_sun.cgs
#m.p.cont["Cont  nu"] * wavnorm 
wavs = wavnorm / m1.data["cont"]["Cont  nu"]
freqs = freqnorm * m1.data["cont"]["Cont  nu"]
sednorm
```

```python
fig, ax = plt.subplots(figsize=(15, 10))
ax.plot(wavs, sednorm * m3.data["cont"]["incident"])
ax.plot(wavs, sednorm * m3.data["cont"]["DiffOut"])
ax.plot(wavs, sednorm * m3.data["cont"]["trans"])

ax.axvline(24.0, lw=5, color="k", alpha=0.3)
for ip in 1.0, 1.8, 4.0:
    ax.axvspan(0, 0.0912/ip, lw=0, color="r", alpha=0.1)
ax.set(
    xscale="log",
    yscale="log",
    xlim=[1e-2, 1e3],
    ylim=[300, 1e6],
    xlabel="Wavelength, micron",
    ylabel=r"$\nu L_\nu$, L$_\odot$",
)
sns.despine()
fig.tight_layout();
```
```python
m4 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n010-p-r08")
m5 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n005-p-r08")
m6 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n100-p-r08")
m7 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n050-p-r08")
m8 = cloudytab.CloudyModel(ROOT / "cloudy/models/w3-n010-d01-r08")
```


```python
m4.data.keys()
```

```python
fig, ax = plt.subplots(figsize=(15, 10))
ax.plot(wavs, sednorm * m4.data["cont"]["incident"])
ax.plot(wavs, sednorm * m4.data["cont"]["DiffOut"], label="n010-p-r08")
ax.plot(wavs, sednorm * m1.data["cont"]["DiffOut"], label="n010")
ax.plot(wavs, sednorm * m8.data["cont"]["DiffOut"], label="n010-d01-r08")
ax.plot(wavs, sednorm * m7.data["cont"]["DiffOut"], label="n050-p-r08")

ax.axvline(24.0, lw=5, color="k", alpha=0.3)
ax.axvspan(0, 0.0912/4, lw=0, color="r", alpha=0.3)
ax.set(
    xscale="log",
    yscale="log",
    xlim=[1e-2, 1e3],
    ylim=[1.0, 1.0e7],
    xlabel="Wavelength, micron",
    ylabel=r"$\nu L_\nu$, L$_\odot$",
)
ax.legend()
sns.despine()
fig.tight_layout();
```

```python
for m in m4, m5, m6, m7, m8:
    m.p = C(m.data)
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)
for m, ax in zip([m4, m7, m8], axes):
    m.radius = m.p.rad.radius * u.cm.to(u.pc) 
    ax.plot(m.radius, m.p.ovr.eden, label="eden")
    ax.plot(m.radius, m.p.ovr.hden * m.p.ovr.HII, label="H II")
    ax.plot(m.radius, m.p.ovr.hden * m.p.ovr.HeII, label="He II")
    ax.plot(m.radius, m.p.ovr.hden * m.p.ovr.HeIII, label="He III")
    ax.plot(m.radius, m.p.ovr.hden * m.p.Ar["Ar+3"], label="Ar IV")
    ax.plot(m.radius, m.p.ovr.hden * m.p.Ne["Ne+2"], label="Ne III")
    ax.plot(m.radius, m.p.ovr.hden * m.p.O["O+2"], label="O III")
    ax.plot(m.radius, 0.001 * m.p.ovr.Te, label="Te, kK")
axes[-1].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
axes[-1].set(
    xlabel="Radius, pc",
)
sns.despine()
fig.tight_layout();
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)


# colnames = m.data["emis"].colnames[1:]

embands = [
 'He 2 4685.70A',
 'Ar 4 4740.12A',
 'Ne 3 3868.76A',
 'O  3 5006.84A',
 'Blnd 5875.66A',
 'Ar 3 7135.79A',
 'H  1 4861.33A',
 'H  1 6562.82A',
 'O  2 7319.99A',
]

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m4, m7, m8], axes):
    radius = m.data["rad"]["radius"] * u.cm.to(u.pc) 
    hb = m.data["emis"]['H  1 4861.33A'] 
    for emband, color in zip(embands, colors):
        em = m.data["emis"][emband] 
        ax.plot(radius, em / hb.max(), label=emband, color=color)
    ax.set(
        yscale="log",
        ylim=[0.001, 10.1],
        xlabel="Radius, pc",
        ylabel="Emissivity",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
sns.despine()
fig.tight_layout();
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)


# colnames = m.data["emis"].colnames[1:]

embands = [
 'He 2 4685.70A',
 'Ar 4 4740.12A',
 'Ne 3 3868.76A',
 'O  3 5006.84A',
 'Blnd 5875.66A',
 'Ar 3 7135.79A',
# 'H  1 4861.33A',
 'H  1 6562.82A',
 'O  2 7319.99A',
]

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m4, m7, m8], axes):
    radius = m.data["rad"]["radius"] * u.cm.to(u.pc) 
    hb = m.data["emis"]['H  1 4861.33A'] 
    for emband, color in zip(embands, colors):
        em = m.data["emis"][emband] 
        ax.plot(radius, em / em.max(), label=emband, color=color)
    ax.set(
        yscale="linear",
        ylim=[0.00, 1.1],
        xlabel="Radius, pc",
        ylabel="Emissivity",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
sns.despine()
fig.tight_layout();
```

## Brightness versus projected radius in spherical symmetry

First approximation to bow shock shape is that it is a hemisphere

Therefore at each projected radius $b$, the brightness is given by:
$$
S(b) = \int_{-\infty}^\infty j(r) \, dz
$$
where
$$
r^2 = b^2 + z^2 
\quad \Rightarrow \quad
2 r\, dr = 2 z\, dz
\quad \Rightarrow \quad
dz = \frac{r}{z}\, dr
$$
Therefore
$$
S(b) = 2 \int_b^\infty j(r) \, \frac{r}{(r^2 - b^2)^{1/2}} \, dr
$$


### Naive implementation of surface brightness intergral

This version makes a regular grid of impact parameter and a regular grid of dummy radii for the integration

```python
nb = 200
def brightness_regrid(r, dr, e, nb, verbose=False):
    b = np.linspace(0.0, r.max(), nb)
    drmin = np.min(dr)
    drmin = np.percentile(dr, 0.5)
    _r = np.arange(0.0, r.max(), step=drmin)
    nr = len(_r)
    if verbose:
        print(f"{nr=}, {drmin=:.2e}")
#    _r = np.linspace(0.0, r.max(), 30 * nb + 5)
    _e = np.interp(_r, r, e, left=0.0, right=0.0)
    _dr = [r.max() / nr] * nr
    bgrid = np.stack([b] * nr, axis=0)
    rgrid = np.stack([_r] * nb, axis=1)
    egrid = np.stack([_e] * nb, axis=1)
    drgrid = np.stack([_dr] * nb, axis=1)
    rgrid[rgrid <= bgrid + 0.5 * drgrid] = np.nan
    sb = 2 * np.nansum(egrid * rgrid * drgrid / (drgrid + np.sqrt(rgrid**2 - bgrid**2)), axis=0)
    return b, sb
```

### Alternatively, we can do the integral in z

This is better, since we have no singularity.  Since we have discrete points, we can do it without having to do any interpolation, using either trapezium or simpson. 

```python
from scipy import integrate
def brightness_discrete(r, dr, e, n_inner=50, verbose=False, integrator=np.trapz):
    """Perform integral of surface brightness along line of sight

    Suitable values for `integrator` are numpy.trapz or scipy.integrate.simpson
    """

    # Use the Cloudy radial points with additional uniform grid from origin to inner boundary
    b_inner = np.linspace(0.0, r.min(), num=n_inner, endpoint=False)
    b = np.concatenate((b_inner, r))
    
    sb = np.zeros_like(b)
    # For each impact parameter
    for i, _b in enumerate(b):
        # Select all radii greater than impact parameter
        m = r >= _b
        # Array of LOS positions for each of these radii
        z = np.sqrt(r[m]**2 - _b**2)
        _e = e[m]
        # Integrate along z to find brightness
        sb[i] = 2 * integrator(_e, z)
    return b, sb
```

```python
brightness = brightness_discrete
```

```python
m = m8
nb = None
r = m.data["rad"]["radius"]
dr = m.data["rad"]["dr"]
e = m.data["emis"]["He 2 4685.70A"]
b, s = brightness(r, dr, e)
len(s)
```

```python
m7.data["emis"].colnames[]
```

### Pre-calculate the surface brightness profiles for all models


```python
from astropy.table import Table

def sb_table(model):
    r = model.data["rad"]["radius"]
    dr = model.data["rad"]["dr"]    
    elabels = model.data["emis"].colnames[1:] # skip the depth column
    sbdict = {}
    for elabel in elabels:
        sbdict["b"], sbdict[elabel] = brightness(r, dr, model.data["emis"][elabel])
    return Table(sbdict)
```

```python
%%timeit -n 1 -r 1
for m in m1, m2, m3, m4, m5, m6, m7, m8:
    m.sb = sb_table(m)
```

### Optical line surface brightness

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)


# colnsmes = m.data["emis"].colnames[1:]

embands = [
 'He 2 4685.70A',
 'Ar 4 4740.12A',
 'Ne 3 3868.76A',
 'O  3 5006.84A',
 'Blnd 5875.66A',
 'Ar 3 7135.79A',
 'H  1 6562.82A',
 'O  2 7319.99A',
]

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m4, m7, m8], axes):
    radius = m.sb["b"] * u.cm.to(u.pc) 
    for emband, color in zip(embands, colors):
        sb = m.sb[emband]
        ax.plot(radius, sb / sb.max(), "-", label=emband, color=color)
    ax.set(
        yscale="linear",
        ylim=[0.00, 1.1],
        xlabel="Radius, pc",
        ylabel="Surface brightness",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
sns.despine()
fig.tight_layout();
```

### IR line and continuum surface brightness

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)


# colnames = m.data["emis"].colnames[1:]

embands = [
 "S  4 10.5076m",
 "Ne 3 15.5509m",
 "S  3 18.7078m",
 "S  3 33.4704m",
 "Ne 2 12.8101m",
 "Si 2 34.8046m",
]
normband = "S  3 18.7078m"
sbnorm = m.sb[normband].mean()

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m4, m7, m8], axes):
    radius = m.sb["b"] * u.cm.to(u.pc) 
    for emband, color in zip(embands, colors):
        sb = m.sb[emband]
        ax.plot(radius, sb, label=emband, color=color)
    ax.set(
        yscale="linear",
        ylim=[0.00, None],
        xlabel="Radius, pc",
        ylabel="Surface brightness",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
sns.despine()
fig.tight_layout();
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)


# colnames = m.data["emis"].colnames[1:]

embands = [
 "PAHC 10.9000m",
 "nFnu 15.6901m",
 "nFnu 19.6199m",
 "nFnu 24.7829m",
 "nFnu 30.8695m",
 "nFnu 41.2152m",
 "nFnu 60.8322m",
]
normband = "nFnu 15.6901m"
sbnorm = m.sb[normband].mean()

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(embands), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m4, m7, m8], axes):
    radius = m.sb["b"] * u.cm.to(u.pc)
    for emband, color in zip(embands, colors):
        sb = m.sb[emband]
        ax.plot(radius, sb / sbnorm, label=emband, color=color)
    ax.set(
        yscale="linear",
        ylim=[0.00, None],
        xlabel="Radius, pc",
        ylabel="Surface brightness",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
sns.despine()
fig.tight_layout();
```

## Line ratio diagnostics

I can look at the same infrared line ratios as in the Spitzer data.  And also at optical ratios

```python
ratio_dict = {
    "c15-11-mir": ("nFnu 15.6901m", "PAHC 10.9000m"),
    "c25-15-mir": ("nFnu 24.7829m", "nFnu 15.6901m"),
    "s43-mir": ("S  4 10.5076m", "S  3 18.7078m"),
    "ne3s3-mir": ("Ne 3 15.5509m", "S  3 18.7078m"),
    "ar43-opt": ("Ar 4 4740.12A", "Ar 3 7135.79A"),
    "o3hb-opt": ("O  3 5006.84A", "H  1 4861.33A"),
    "He21-opt": ("He 2 4685.70A", "Blnd 5875.66A"),
}
```

```python
def sb_ratio_table(model, ratios):
    result = {"b": model.sb["b"]}
    for key, value in ratios.items():
        result[key] = model.sb[value[0]] / model.sb[value[1]]
    return Table(result)
```

```python
%%timeit -n 1 -r 1
for m in m1, m2, m3, m4, m5, m6, m7, m8:
    m.sbratios = sb_ratio_table(m, ratio_dict)
```

```python
m4.sbratios.to_pandas().describe()
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(ratio_dict), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m4, m7, m8], axes):
    radius = m.sbratios["b"] * u.cm.to(u.pc)
    rlabels = m.sbratios.colnames[1:]
    for rlabel, color in zip(rlabels, colors):
        ax.plot(radius, m.sbratios[rlabel], label=rlabel, color=color)
    ax.set(
        yscale="log",
        ylim=[0.01, 500],
        xlabel="Radius, pc",
        ylabel="Ratio",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant pressure, n = 10, Rmax = 8 pc")
axes[1].set_title("Constant pressure, n = 50, Rmax = 8 pc")
axes[2].set_title("Density law $r^{-1}$, n = 10, Rmax = 8 pc")
sns.despine()
fig.tight_layout();
```

```python
fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)

# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.chroma_r', 
    len(ratio_dict), 
    cmap_range=(0.15, 0.85), 
    return_fmt='hex'
)

for m, ax in zip([m1, m2, m3], axes):
    radius = m.sbratios["b"] * u.cm.to(u.pc)
    rlabels = m.sbratios.colnames[1:]
    for rlabel, color in zip(rlabels, colors):
        ax.plot(radius, m.sbratios[rlabel], label=rlabel, color=color)
    ax.set(
        yscale="log",
        ylim=[0.01, 100],
        xlabel="Radius, pc",
        ylabel="Ratio",
    )
axes[0].legend(ncol=3)
axes[0].set_title("Constant density, n = 10")
axes[1].set_title("Constant pressure, n = 10")
axes[2].set_title("Constant pressure, n = 30")
sns.despine()
fig.tight_layout();
```

So these are looking pretty different from the observed ratios, with the exception of the 25/15 micron continuum ratio, which is around 10 in both the observations and the model. 

```python

```
