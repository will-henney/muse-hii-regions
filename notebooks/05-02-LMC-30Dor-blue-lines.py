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

cubeA = Cube("../big-data/lmc-30dor-A-subcube-46-55.fits")
cubeB = Cube("../big-data/lmc-30dor-B-subcube-46-55.fits")
cubeC = Cube("../big-data/lmc-30dor-C-subcube-46-55.fits")
cubeD = Cube("../big-data/lmc-30dor-D-subcube-46-55.fits")

fig, ax = plt.subplots(figsize=(12, 6))
for cube, label in zip([cubeA, cubeB, cubeC, cubeD], "ABCD"):
    cube.sum(axis=(1, 2)).plot(label=label)
ax.legend()
ax.set(
    ylim=[0.8e8, 2.2e8],
)
sns.despine();

# We can see the blue-bump from WR stars arounf 4700 Å.  This includes broad He II plus C III/IV, N III, and N V.  Hopefull, this will be restricted to particular regions in each field.
#
# We can try the same wav ranges as we used for NGC 346

wavranges = [
    (4600, 4635), (4725, 4738), (4780, 4815), 
    #(4890, 4900), (4937, 4947), (4974, 4984), (5028, 5038), 
    (5070, 5100), (5120, 5145), (5220, 5260), (5300, 5340), (5470, 5500),
]

nv, ny, nx = cubeA.data.shape
ny, nx

# Split each field up into 4 roughly equal tiles:

# +
mm = 160 # middle of each image
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
    ylim=[None, 2e4],
)
sns.despine();
# -

contA = extract.fit_continuum(
    cubeA, wav_ranges=wavranges, deg=5, median=False,
)

fig, ax = plt.subplots(figsize=(12, 8))
cubeA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
contA[:, 275:300, 120:150].mean(axis=(1, 2)).plot()
for wavrange in wavranges:
    ax.axvspan(*wavrange, alpha=0.3)
ax.set(ylim=[0, 800])

fig, axes = plt.subplots(
    3, 2, 
    figsize=(12, 18), 
    sharex=True, sharey=True
)
(cubeA - contA).select_lambda(4740, 4750).sum(axis=0).plot(
    ax=axes[0, 0], vmin=-10, vmax=1000,
)
(cubeA - contA).select_lambda(4685, 4692).sum(axis=0).plot(
    ax=axes[0, 1], vmin=-10, vmax=500,
)
(cubeA - contA).select_lambda(4630, 4650).sum(axis=0).plot(
    ax=axes[1, 0], vmin=-10, vmax=700,
)
(cubeA - contA).select_lambda(5200, 5210).sum(axis=0).plot(
    ax=axes[1, 1], vmin=-10, vmax=1000,
)
(cubeA - contA).select_lambda(4958, 4964).sum(axis=0).plot(
    ax=axes[2, 0], vmin=-10, vmax=900000,
)
(cubeA - contA).select_lambda(4820, 4850).sum(axis=0).plot(
    ax=axes[2, 1], vmin=-10, vmax=300,
)

contB = extract.fit_continuum(
    cubeB, wav_ranges=wavranges, deg=5, median=False,
)

contC = extract.fit_continuum(
    cubeC, wav_ranges=wavranges, deg=5, median=False,
)

contD = extract.fit_continuum(
    cubeD, wav_ranges=wavranges, deg=5, median=False,
)

fig, axes = plt.subplots(
    2, 2, 
    figsize=(12, 12), 
    sharex=True, sharey=True
)
for cube, cont, ax in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    axes.flat,
):
    (cube - cont).select_lambda(4740, 4750).sum(axis=0).plot(
        ax=ax, 
        vmin=-10, 
        vmax=1500,
    )
fig.suptitle("[Ar IV] 4740")
fig.tight_layout();

fig, axes = plt.subplots(
    2, 2, 
    figsize=(12, 12), 
    sharex=True, sharey=True
)
for cube, cont, ax in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    axes.flat,
):
    (cube - cont).select_lambda(4685, 4692).sum(axis=0).plot(
        ax=ax, 
        vmin=-10, 
        vmax=1500,
    )
fig.suptitle("He II 4686")
fig.tight_layout();

fig, axes = plt.subplots(
    2, 2, 
    figsize=(12, 12), 
    sharex=True, sharey=True
)
for cube, cont, ax in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    axes.flat,
):
    (cube - cont).select_lambda(5200, 5210).sum(axis=0).plot(
        ax=ax, 
        vmin=-10, 
        vmax=1500,
    )
fig.suptitle("[N I] 5200")
fig.tight_layout();

fig, axes = plt.subplots(
    2, 2, 
    figsize=(12, 12), 
    sharex=True, sharey=True
)
for cube, cont, ax in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    axes.flat,
):
    (cube - cont).select_lambda(4820, 4855).sum(axis=0).rebin(2).plot(
        ax=ax, 
        vmin=-10, 
        vmax=1000,
    )
fig.suptitle("Hβ blue wing")
fig.tight_layout();

fig, axes = plt.subplots(
    2, 2, 
    figsize=(12, 12), 
    sharex=True, sharey=True
)
for cube, cont, ax in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    axes.flat,
):
    (cube - cont).select_lambda(4855, 4870).sum(axis=0).rebin(1).plot(
        ax=ax, 
        vmin=-10, 
        vmax=500000,
    )
fig.suptitle("Hβ core")
fig.tight_layout();

fig, axes = plt.subplots(
    2, 2, 
    figsize=(12, 12), 
    sharex=True, sharey=True
)
for cube, cont, ax in zip(
    [cubeD, cubeC, cubeB, cubeA],
    [contD, contC, contB, contA],
    axes.flat,
):
    (cube - cont).select_lambda(4871, 4876).sum(axis=0).rebin(1).plot(
        ax=ax, 
        vmin=-10, 
        vmax=3000,
    )
fig.suptitle("Hβ red wing")
fig.tight_layout();


