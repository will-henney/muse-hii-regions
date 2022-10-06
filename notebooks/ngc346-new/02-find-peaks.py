# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# + [markdown] pycharm={"name": "#%% md\n"}
# # Automatically find peaks in the spectrum

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Imports
#
# Standard library imports

# + pycharm={"name": "#%%\n"}
from pathlib import Path
from typing import Union

# + [markdown] pycharm={"name": "#%% md\n"}
# Third-party  imports

# + pycharm={"name": "#%%\n"}
import numpy as np
import scipy.signal as si
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
import seaborn as sns
import cmasher as cmr
from mpdaf.obj import Cube, Image, Spectrum
import regions as rg

sns.set_context("talk")

# + [markdown] pycharm={"name": "#%% md\n"}
# Will's libraries

# + pycharm={"name": "#%%\n"}
import wcsfile

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Paths for data files
#
# Files in `small_data_path` should not be too large to be checked into git

# + pycharm={"name": "#%%\n"}
data_path = Path.cwd().parent.parent / "big-data" / "ngc346new"
small_data_path = Path.cwd().parent.parent / "data"

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Read in data files

# + pycharm={"name": "#%%\n"}
cube = Cube(str(data_path / "n346-muse-csub-101.fits"))

# + [markdown] pycharm={"name": "#%% md\n"}
# Switch to reading the fits files directly

# + pycharm={"name": "#%%\n"}
star_mask = fits.open(small_data_path / "n346-mask-stars.fits")[1].data
yso_mask = fits.open(small_data_path / "n346-mask-yso.fits")[1].data

# + pycharm={"name": "#%%\n"}
Image(data=star_mask - yso_mask).plot()

# + pycharm={"name": "#%%\n"}
mask = np.where(yso_mask > 0.0, 1.0, np.nan)
spec_yso = np.nanmean(cube.data * mask[None, :, :], axis=(1, 2))

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_yso, wave=cube.wave).plot()
ax.set_ylim(-50, 200)
ax.set_xlim(8100, 8600)
...;

# + pycharm={"name": "#%%\n"}
mask = np.where(yso_mask + star_mask == 0.0, 1.0, np.nan)
spec_all = np.nanmean(cube.data * mask[None, :, :], axis=(1, 2))

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_all, wave=cube.wave).plot()
ax.set_ylim(-5, 20)
ax.set_xlim(6200, 6800)
...;

# + pycharm={"name": "#%%\n"}
peaks, props = si.find_peaks(spec_yso, prominence=5.0,   distance=4)

# + pycharm={"name": "#%%\n"}
peaks

# + pycharm={"name": "#%%\n"}
props['prominences']

# + pycharm={"name": "#%%\n"}
wave = cube.wave

# + pycharm={"name": "#%%\n"}
wave.coord()[peaks]

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_yso, wave=cube.wave).plot()
ymax=200
ax.scatter(wave.coord()[peaks], np.minimum(ymax, props["prominences"]), marker=".", color="r")
ax.set_ylim(-50, ymax * 1.1)
ax.set_xlim(8100, 8600)
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_yso, wave=cube.wave).plot()
ymax=200
ax.scatter(wave.coord()[peaks], np.minimum(ymax, props["prominences"]), marker=".", color="r")
ax.set_ylim(-50, ymax * 1.1)
ax.set_xlim(6200, 6800)
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_yso, wave=cube.wave).plot()
ymax=200
ax.scatter(wave.coord()[peaks], np.minimum(ymax, props["prominences"]), marker=".", color="r")
ax.set_ylim(-50, ymax * 1.1)
ax.set_xlim(4600, 5400)
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_yso, wave=cube.wave).plot()
ymax=200
ax.scatter(wave.coord()[peaks], np.minimum(ymax, props["prominences"]), marker=".", color="r")
ax.set_ylim(-0.1 * ymax, ymax * 1.1)
ax.set_xlim(6950, 7300)
...;

# + pycharm={"name": "#%%\n"}
peaks_all, props_all = si.find_peaks(spec_all, prominence=0.3, distance=4)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_all, wave=cube.wave).plot()
ymax=5
ax.scatter(wave.coord()[peaks_all], np.minimum(ymax, props_all["prominences"]), marker="x", color="m")
ax.set_ylim(-0.1 * ymax, ymax * 1.1)
ax.set_xlim(6950, 7300)
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_all, wave=cube.wave).plot()
ymax=10
ax.scatter(wave.coord()[peaks_all], np.minimum(ymax, props_all["prominences"]), marker="x", color="m")
ax.set_ylim(-0.1 * ymax, ymax * 1.1)
ax.set_xlim(4600, 5400)
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_all, wave=cube.wave).plot()
ymax=20
ax.scatter(wave.coord()[peaks_all], np.minimum(ymax, props_all["prominences"]), marker="x", color="m")
ax.set_ylim(-0.1 * ymax, ymax * 1.1)
ax.set_xlim(8100, 8600)
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
Spectrum(data=spec_all, wave=cube.wave).plot()
ymax=20
ax.scatter(wave.coord()[peaks_all], np.minimum(ymax, props_all["prominences"]), marker="x", color="m")
ax.set_ylim(-0.1 * ymax, ymax * 1.1)
ax.set_xlim(8600, 9300)
...;

# + pycharm={"name": "#%%\n"}
props_all['prominences']


# + pycharm={"name": "#%%\n"}
import pandas as pd

# + pycharm={"name": "#%%\n"}
df = pd.DataFrame(props_all, index=wave.coord()[peaks_all])
df

# + pycharm={"name": "#%%\n"}

