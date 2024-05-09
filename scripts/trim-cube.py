from astropy.io import fits
import numpy as np
import typer
from pathlib import Path


def main(
    hdu_index: int = 1,
    cube_name: str = "PeterZeidler/DATACUBE_FINAL_fwhm_cor.fits",
    data_path=Path.home() / "Work/Muse-Hii-Data/SMC-NGC-346",
    extra_suffix: str = "_trim",
    target_shape: tuple[int, int, int] = (3801, 326, 346),
):
    # Open the FITS file
    hdu = fits.open(data_path / cube_name)[hdu_index]
    # Trim off the rightmost column
    hdu.data = hdu.data[..., :-1]
    # Check shape of result
    assert hdu.data.shape == target_shape

    # Write the new FITS file
    hdu.writeto(data_path / cube_name.replace(".fits", f"{extra_suffix}.fits"), overwrite=True)


if __name__ == "__main__":
    typer.run(main)
