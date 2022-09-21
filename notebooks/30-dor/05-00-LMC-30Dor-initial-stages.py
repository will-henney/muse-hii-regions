# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Initial stages of data reduction for 30 Dor
#
# The plan is to split the big cubes into more manageable chunks, so that Mabel can more easily work with them.

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube

sns.set_context("talk")

# I have already looked at field A, so this time I will look at field C first:

# +
datapath = Path("/Users/will/Work/Muse-Hii-Data/30Dor/2021-06/")
fileA = "ADP.2016-07-14T14:17:17.826.fits"
fileB = "ADP.2016-07-19T06:16:41.316.fits"
fileC = "ADP.2016-07-19T11:27:30.987.fits"
fileD = "ADP.2016-07-15T13:22:09.432.fits"

# cubeA = Cube(str(datapath / fileA))
# cubeB = Cube(str(datapath / fileB))
cubeC = Cube(str(datapath / fileC))
# cubeD = Cube(str(datapath / fileD))
# -

# We will divide the wavelength range into chunks of length 900 Å, but with overlap region of of 100 Å between each.

wav_min, _, _, wav_max, _, _ = cubeC.get_range()
wav_min, wav_max, wav_max - wav_min

wav_window = 900.0
wav_overlap = 100.0
nwin = 1 + np.ceil((wav_max - wav_min) / (wav_window + wav_overlap))
wav_win_min = wav_min + np.arange(nwin) * (wav_window - wav_overlap)
wav_win_max = wav_win_min + wav_window
np.round(wav_win_min), np.round(wav_win_max)


def divide_in_subcubes(cube, width=900.0, overlap=100.0):
    """
    Divide MPDAF `cube` into subcubes of `width` Angstrom,
    overlapped by `overlap` Angstrom.

    Returns dict of the cubes with labels of the wavelength range
    in units of 100 Angstrom, e.g., '46-55' for 4600 to 5500 Angstrom
    """
    wav_min, _, _, wav_max, _, _ = cube.get_range()
    nwin = 1 + np.ceil((wav_max - wav_min) / (width + overlap))
    wav_win_min = wav_min + np.arange(nwin) * (width - overlap)
    wav_win_max = wav_win_min + width
    subcubes = {}
    for wav1, wav2 in zip(wav_win_min, wav_win_max):
        label = f"{int(round(wav1/100)):d}-{int(round(wav2/100)):d}"
        subcubes[label] = cube.select_lambda(wav1, wav2)
    return subcubes


subcubes = divide_in_subcubes(cubeC)

subcubes


def write_subcubes(
    subcubes,
    prefix,
    outdir="../big-data",
):
    """
    Write each subcube to a separate FITS file

    Returns a list of all the files
    """
    file_paths = []
    for label, subcube in subcubes.items():
        file_path = f"{outdir}/{prefix}-subcube-{label}.fits"
        subcube.write(file_path)
        file_paths.append(file_path)
    return file_paths


write_subcubes(subcubes, "lmc-30dor-C")

# Now try to do all the steps at once for field A

# %%timeit -n1 -r1
_files = write_subcubes(
    divide_in_subcubes(
        Cube(str(datapath / fileA))
    ),
    "lmc-30dor-A",
)

# That seemed to work, and it took 42.6 seconds on my machine, which includes reading in the big cube, splitting it up and then writing it out again.  Hopefully it doesn't leak any memory.
#
# It seems that variables set inside the `%%timeit` block are not available outside it, so I will try printing the file list instead.
#
# Now do the remaining two cubes:

# %%timeit -n1 -r1
_files = write_subcubes(
    divide_in_subcubes(
        Cube(str(datapath / fileB))
    ),
    "lmc-30dor-B",
)
print(_files)


# That was even faster: 30 seconds.
#
# Looks like the D cube did not download fully. I am trying to get it again, but is going to take 30 minutes.

_files = write_subcubes(
    divide_in_subcubes(
        Cube(str(datapath / fileD))
    ),
    "lmc-30dor-D",
)
print(_files)

# Hurray, that worked finally.  So this notebooks work is now finished. 
