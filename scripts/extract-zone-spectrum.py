from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import typer
import sys

def get_spectrum_from_cube(
        cube: np.ndarray,
        reduction_method: callable = np.nanmean,
) -> np.ndarray:
    """Extract 1D spectrum from cube"""
    assert len(cube.shape) == 3
    spec = reduction_method(
        cube,
        axis=(1, 2),
    )
    assert len(spec.shape) == 1 and len(spec) == cube.shape[0]
    return spec


def main(
        cube_file: str,
        output_id: str,
        output_dir: str="zone_spectra",
        method: str="mean",
        jmin: int=0,
):
    """Extract 1D spectra from cube for each region in file"""

    reduction_method_options = {
        "mean": np.nanmean,
        "median": np.nanmedian,
    }
    assert method in reduction_method_options

    # Read the spectral cube
    hdu = fits.open(cube_file)[1]
    yslice = slice(jmin, None)
    spec = get_spectrum_from_cube(
        hdu.data[:, yslice, :],
        reduction_method=reduction_method_options[method],
    )

    # Make sure the output folder exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    # And save spectrum as a FITS file
    fits.PrimaryHDU(
        header=WCS(hdu.header).spectral.to_header(),
        data=spec,
    ).writeto(f"{output_dir}/{output_id}-spec1d.fits", overwrite=True)


if __name__ == "__main__":
    typer.run(main)
