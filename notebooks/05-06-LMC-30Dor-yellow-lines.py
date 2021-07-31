# -*- coding: utf-8 -*-
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

# +
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
# -

cubeA = Cube("../big-data/lmc-30dor-A-subcube-54-63.fits")
cubeB = Cube("../big-data/lmc-30dor-B-subcube-54-63.fits")
cubeC = Cube("../big-data/lmc-30dor-C-subcube-54-63.fits")
cubeD = Cube("../big-data/lmc-30dor-D-subcube-54-63.fits")

# ## Inspect the average spectrum for each field

fig, ax = plt.subplots(figsize=(12, 6))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cube.sum(axis=(1, 2)).plot(label=label)
ax.legend()
ax.set(
    ylim=[0.5e8, 1.5e8],
)
sns.despine()

# Again, we get WR features (I think)  around 5800.  We can also see the DIB 5781 absorption. 

# ## Define continuum wavelength ranges
#
# We can try the same wav ranges as we used for NGC 346

wavranges = [
    (5470, 5510),
    (5550, 5570),
    (5600, 5700),
    (6000, 6030), (6070, 6100),
    (6180, 6220), (6260, 6280),
]

nv, ny, nx = cubeA.data.shape
ny, nx

# Split each field up into 4 roughly equal tiles:

# +
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
# -

# ## Test the polynomial fitting with field A
#
# This takes about a minute for each field.

contA = extract.fit_continuum(
    cubeA,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)

# ### Inspect the results for different portions of the field.

fig, ax = plt.subplots(figsize=(12, 8))
cubeA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
contA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.set(ylim=[0, 800])

# I had to go back and forth a few times adjusting the wav ranges.  It is difficult to get a good fit on the blue side because of the WR features.
#
# As can be seen here, the final version is not perfect – it slightly overpredicts the continuum around 4700 to 4800.  This might affect some of the weak [Fe III] lines, but the [Ar IV] 4740 does not seem to be much affected.

# ### Look at some line images

# +
fig, axes = plt.subplots(3, 2, figsize=(12, 18), sharex=True, sharey=True)
(cubeA - contA).select_lambda(5517, 5527).sum(axis=0).plot(
    ax=axes[0, 0],
    vmin=-10,
    vmax=2000,
)
axes[0, 0].set_title("[Cl III] 5517")
(cubeA - contA).select_lambda(5538, 5548).sum(axis=0).plot(
    ax=axes[0, 1],
    vmin=-10,
    vmax=2000,
)
axes[0, 1].set_title("[Cl III] 5538")
(cubeA - contA).select_lambda(5958, 5968).sum(axis=0).plot(
    ax=axes[1, 0],
    vmin=-10,
    vmax=200,
)
axes[1, 0].set_title("Si II 5958")
(cubeA - contA).select_lambda(5755, 5765).sum(axis=0).plot(
    ax=axes[1, 1],
    vmin=-10,
    vmax=600,
)
axes[1, 1].set_title("[N II] 5755")
(cubeA - contA).select_lambda(5876, 5886).sum(axis=0).plot(
    ax=axes[2, 0],
    vmin=-10,
    vmax=90000,
)
axes[2, 0].set_title("[He I 5876")
(cubeA - contA).select_lambda(5979, 5989).sum(axis=0).plot(
    ax=axes[2, 1],
    vmin=-10,
    vmax=300,
)
axes[2, 1].set_title("Si II 5979")


# -

# ## Now do the other fields

contB = extract.fit_continuum(
    cubeB,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)

contC = extract.fit_continuum(
    cubeC,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)

contD = extract.fit_continuum(
    cubeD,
    wav_ranges=wavranges,
    deg=5,
    median=False,
)

# ## Save the continuum-subtracted cubes

csub = {}
cdiv = {}
for cube, cont, label in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    "DCBA",
):
    prefix = f"../big-data/lmc-30dor-{label}-subcube-54-63"
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

