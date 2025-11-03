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
# Notebook G of more line ratios from saved maps from the PZ cubes

This is adapted from 01-02-yet-more-line-ratios. It includes the Ha density from emission measure

<!-- #endregion -->


```python jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Image
from zz_utils import get_image_raw, ROOT


# ## Alternative route to density from Balmer lines
#
# We can estimate the bow shock shell emission measure from the jump in the surface brightness of H beta, say. This has the advantage that we know that the hydrogen is fully ionized.

A_hb = 0.33
jump_hb = 31.7 * 300 * muse_bright_unit * 10 ** (0.4 * A_hb)
f"{jump_hb=}"

# The 300 is my estimate from the above figure of the Hb jump at the bow shock inner edge, while the 31.7 is what we divided Hb by before plotting it.
#
# Thhen we also correct for the dust extinction `A_hb`.
#
# So this should be the intrinsic surface brightness of Hb, which should be EM times emission coefficient / 4 pi
#

EM = jump_hb * (4 * np.pi * u.sr).to(u.arcsec**2) / (e4861 * pn_e_units)
EM

# There are three answers, corrending to assumed temperature of 12,500, 13,800, or 15,500 K

np.sqrt(EM / depth_heii).to(u.cm**-3)

# So that gives a bigger density than I was getting before. On the other hand, it might be a slight over-estimate since the H+ path length may be largerby a factor of 2 (see next section). There is also a small correction for the ionized helium contribution to the electron density, but that is much smaller than the uncertainty in the path length.
#
# So, we would finally get
#

np.sqrt(EM / (2 * depth_heii)).to(u.cm**-3)

# So, taking into account the uncertainty in the path length would give us a hydrogen density of 35 +/- 10 pcc

# ### Sanity check from He II / H I
#
# In principle, we can come to the same conclusion starting from the 4686/4861 ratio, so long as we account for the fraction of the H beta emission that comes from the bow shock.
```

```python
fig, ax = plt.subplots(figsize=(15, 6))
ix0 = 227.5
nx = len(heii_profile)
pos = (np.arange(nx) - ix0) * 0.2
imargin = 10
ratio = heii_profile / hb_profile
ax.plot(pos[imargin:], ratio[imargin:], label="He II / H beta", lw=4)

ax.axhline(0, color="k")

ax.axvline(0, color="k", lw=1, ls="dashed")
ax.axvspan(2.0, 9.0, 0.4, 0.8, color="k", alpha=0.1, linewidth=0, zorder=-100)
ax.legend(ncol=3, loc="upper left")

ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Line ratio",
    #    xlim=[-12, 22],
    ylim=[-0.005, 0.018],
)
sns.despine()
fig.savefig(figdir / "ngc346-bow-shock-heii-hbeta-cuts.pdf")
```

So the peak ratio is about 0.015. Note that we already did this once above in the *He II emission measure* section. We got exactly the same value.

We now have to divide this by the fraction of H beta that comes frm the bow shock.

```python
hb_bow_frac = 0.3
heii_hb_bow_max = np.max(ratio[imargin:]) / hb_bow_frac
heii_hb_bow_max
```

Then we should have
$$
\texttt{heii\_hb\_bow\_max} =
\frac{I(4686)}{I(4681)} =
\frac{
\int e(4686)\, n(\mathrm{He^{++}})\, n_e\, dz
}{
\int e(4861)\, n(\mathrm{H^+})\, n_e\, dz
}
\approx
\frac{e(4686)}{e(4861)}\,
y(\mathrm{He})\,
x(\mathrm{He}^{++})
$$

```python
e4686, e4861, e4686 / e4861
```

So the He II emission coefficient is about 12 times larger than H beta, which almost exactly cancels out the He abundance factor. And the temperature dependence is very slight.

```python
yhe = 0.08325295
xheiii = (heii_hb_bow_max * e4861) / (yhe * e4686)
xheiii
```

Note that this assumes that the line-of-sight path length is the same for He++ and for H+. In reality it may be longer for H+ because the thickness of the emitting shell is greater. In the thin shell approximation, the ratio of path lengths is the sqrt of the ratio of thickness.

The ratio of thicknesses is between 2 and 5, depending on just how we measure it. So this would give a ratio of path lengths of about (2 +/- 1), which is also what you would have guessed by looking at the chord lengths on the plane of the sky (for [Ar IV] and He II).

This would give a final value of $x(\mathrm{He^{++}}) \approx 0.1$


<!-- #region jupyter={"outputs_hidden": false} pycharm={"name": "#%% md\n"} -->

<!-- #endregion -->
