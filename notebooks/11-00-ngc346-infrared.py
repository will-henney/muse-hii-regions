# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # NGC 346 infrared profiles
#
# Broad slit profiles across the bow shock

from pathlib import Path
import numpy as np
from mpdaf.obj import Image
from matplotlib import pyplot as plt
import seaborn as sns
import regions
import cmasher as cmr
sns.set_context("talk")
sns.set_color_codes()

irfiles = Path("../data").glob("ngc346-ir-*.fits")
irfiles = sorted(irfiles)
irfiles


class IRim:
    def __init__(self, path):
        self.path = path
        self.im = Image(str(self.path))
        s = str(path).replace("../data/ngc346-ir-", "").replace(".fits", "")
        self.wav = np.round(0.1 * float(s[:4]), 1)
        self.label = s[5:]
        self.profile = self.im.data[710:790, :].mean(axis=0)
        nx = len(profile)
        self.pos = 0.2 * (np.arange(nx) - nx/2)
        self.bg = (self.profile[:50].mean() + self.profile[-50:].mean()) / 2
        self.profile -= self.bg
        self.im.data -= self.bg


        
    def __repr__(self):
        return f"IRim({self.label}, {self.wav})"


irdata = {}
for path in sorted(irfiles):
    irim = IRim(path)
    irdata[irim.label] = irim
irdata

irdata["WISE3"].im.plot(vmin=0.0, vmax=35.0, use_wcs=True)

fig, ax = plt.subplots(figsize=(12, 6))
for irim in irdata.values():
    ax.plot(irim.pos, irim.profile, label=irim.label)
ax.legend(ncol=6, fontsize="x-small")
ax.axvline(0.0)
ax.set(
    xlabel="Offset, arcsec",
    ylabel="Surface brightness, MJy / sr",
    yscale="log",
    ylim=[0.1, 9000],
)
sns.despine();

goodbands = "WISE3", "MSX-E", "MSX-C", "WISE4", "MIPS1"

fig, ax = plt.subplots(figsize=(12, 6))
# Take N colors from named colormap in [0.15, 0.85] range in HEX
colors = cmr.take_cmap_colors(
    'cmr.apple_r', 
    len(irdata), 
    cmap_range=(0.25, 0.75), 
    return_fmt='hex'
)
offset = 0.0
for irim, color in zip(irdata.values(), colors):
    lw = 3.0 if irim.label in goodbands else 1.0
    alpha = 0.5 if "MSX" in irim.label else 1.0
    norm = np.mean(irim.profile.data)
    ax.plot(irim.pos, offset + irim.profile / norm, color=color, lw=lw, alpha=alpha)
    ax.text(155.0, offset, irim.label, va="center", fontsize="x-small", color=color)
    ax.text(-155.0, offset, f"{irim.wav:.1f}", ha="right", va="center", fontsize="x-small", color=color)
    offset += 1.0
#ax.legend(ncol=6, fontsize="x-small")
ax.axvline(0.0)
ax.set(
    xlabel="Offset west from W 3, arcsec",
    ylabel="Relative brightness + offset",
    yscale="linear",
    ylim=[-2.0, 3.0 + offset],
    xlim=[-180, 180],
)
sns.despine()
fig.savefig("../figs/ngc346-infrared-profiles.pdf");




