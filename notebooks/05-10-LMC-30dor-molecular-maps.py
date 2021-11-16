# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Maps of molecular gas around 30 Dor

# We want to compare the distribution of molecular gas with that of the Raman wings. 
#
# Mabel is going to be doing most of that, but here I am going to create some integrated line maps from the original data cubes obtained from observatory archives.

from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

# ## Alma maps of 30 Dor-10 GMC
#
# These maps are described in [Indebetouw et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...888...56I).
#
#

DATAPATH = Path.cwd().parent / "big-data" / "30-Dor-Radio"

DATAPATH

PATH_2013 = DATAPATH / "2013.1.00346.S/science_goal.uid___A001_X12a_X279/group.uid___A001_X12a_X27a/member.uid___A001_X12a_X27b/product"

hdulist = fits.open(PATH_2013 / "member.uid___A001_X12a_X27b.30_doradus_13CO21.image.fits")

hdulist.info()

w = WCS(hdulist["PRIMARY"].header)

w.celestial

# The original data has 4 axes (4th is degenerate axis of Stokes parameters).  So, we have to sum over 2 axes to get down to the pure image.

image = np.nansum(
    hdulist["PRIMARY"].data,
    axis=(0, 1),
)
image.shape

# Now write out the image with a new header:

fits.PrimaryHDU(
    header=w.celestial.to_header(),
    data=image,
).writeto(
    DATAPATH / "Alma-2013.1.00346.S-30_doradus_13CO21-sum.fits",
    overwrite=True,
)

image = np.nanmax(
    hdulist["PRIMARY"].data,
    axis=(0, 1),
)
fits.PrimaryHDU(
    header=w.celestial.to_header(),
    data=image,
).writeto(
    DATAPATH / "Alma-2013.1.00346.S-30_doradus_13CO21-peak.fits",
    overwrite=True,
)

# ## The 2019 Alma observations from the `2019.1.00843.S` program

RAW_DATAPATH = Path.home() / "Work"/ "Alma-Data" / "LMC-30-Dor"


def fitspath(uid: str, spw: int=25):
    """Find an ALMA spectral cube file with given `uid`"""
    filename = f"member.uid___A001_X1465_X{uid}.30_Doradus_sci.spw{spw}.cube.I.pbcor.fits"
    matches = list(RAW_DATAPATH.rglob(filename))
    assert len(matches) == 1
    return matches[0]


# Check that we can find a cube:

fitspath("219a")

paths_12co = {uid: fitspath(uid) for uid in ("218a", "2192", "219a")}

for p in paths_12co.values():
    fits.open(p).info()

for uid, p in paths_12co.items():
    hdu = fits.open(p)["PRIMARY"]
    w = WCS(hdu.header)
    image = np.nansum(hdu.data, axis=(0, 1))
    savepath = DATAPATH / f"Alma-2019.1.00843.S-30_doradus_12CO21-{uid}-sum.fits"
    fits.PrimaryHDU(header=w.celestial.to_header(), data=image).writeto(savepath, overwrite=True)
    image = np.nanmax(hdu.data, axis=(0, 1))
    savepath = DATAPATH / f"Alma-2019.1.00843.S-30_doradus_12CO21-{uid}-peak.fits"
    fits.PrimaryHDU(header=w.celestial.to_header(), data=image).writeto(savepath, overwrite=True)
    


