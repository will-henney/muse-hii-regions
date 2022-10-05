"""Determine spectral continuum by median filtering

These are based on earlier routines that I wrote for the
00-check-peter-cube notebook. The difference is that this time I will
not use the mpdaf library, so it should be more general and more
efficient

Will Henney: 2022-10-04

"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
import scipy.ndimage as ndi

DATA_PATH = Path.home() / "Work/Muse-Hii-Data/SMC-NGC-346"

cubename = "ADP.2017-10-16T11_04_19.247.fits"
outprefix = "n346-muse"


def get_median_continuum(data: np.ndarray, window_size=11):
    """Take windowed median along first axis of data array"""
    ndim = len(data.shape)
    # Conform the window size parameter to the shape of the data
    size = (window_size,) + (1,) * (ndim - 1)
    return ndi.median_filter(data, size=size)


if __name__ == "__main__":

    try:
        window_size = int(sys.argv[1])
    except (IndexError, ValueError):
        sys.exit(f"Usage: {sys.argv[0]} WINDOW_SIZE")

    hdulist = fits.open(DATA_PATH / cubename)
    cont_data = get_median_continuum(hdulist[1].data, window_size)

    fits.PrimaryHDU(header=hdulist[1].header, data=cont_data,).writeto(
        f"{outprefix}-cont-{window_size:03d}.fits",
        overwrite=True,
    )
    fits.PrimaryHDU(
        header=hdulist[1].header,
        data=hdulist[1].data - cont_data,
    ).writeto(
        f"{outprefix}-csub-{window_size:03d}.fits",
        overwrite=True,
    )
