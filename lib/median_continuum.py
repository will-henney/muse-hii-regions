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
import typer


def get_median_continuum(data: np.ndarray, window_size=11):
    """Take windowed median along first axis of data array"""
    ndim = len(data.shape)
    # Conform the window size parameter to the shape of the data
    size = (window_size,) + (1,) * (ndim - 1)
    return ndi.median_filter(data, size=size)


def main(
    window_size: int,
    cube_name: str = "ADP.2017-10-16T11_04_19.247.fits",
    out_prefix: str = "n346-muse",
    data_path=Path.home() / "Work/Muse-Hii-Data/SMC-NGC-346",
    two_pass: bool = False,
    first_window_size: int = 11,
    shave_threshold: float = 1.0,
):
    """Find and remove continuum from cube by median filtering"""

    hdulist = fits.open(data_path / cube_name)

    if two_pass:
        # First filter pass
        cont_data = get_median_continuum(hdulist[1].data, first_window_size)
        # Save off the lines that go more than shave_threshold above continuum
        shaved_data = np.minimum(hdulist[1].data, cont_data + shave_threshold)
        # Second filter pass
        cont_data = get_median_continuum(shaved_data, window_size)
    else:
        cont_data = get_median_continuum(hdulist[1].data, window_size)

    # Write out the new cubes
    for data, label in [
        (cont_data, "cont"),  # Continuum
        (hdulist[1].data - cont_data, "csub"),  # Original minus continuum
        (hdulist[1].data / cont_data, "cdiv"),  # Original over continuum
    ]:
        fits.PrimaryHDU(header=hdulist[1].header, data=data).writeto(
            f"{out_prefix}-{label}-{window_size:03d}.fits",
            overwrite=True,
        )


if __name__ == "__main__":
    typer.run(main)
