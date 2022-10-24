from pathlib import Path
import numpy as np
from astropy.io import fits
import typer

def combine_maps(
        pattern: str="*-uil-*",
        prefix: str="all-uil",
        min_strength: float=0.0,
):
    """Make sum and median images by combining several maps"""

    fits_paths = sorted(Path.cwd().glob(f"{pattern}.fits"))

    hdus = [fits.open(p)[0] for p in fits_paths]

    data_stack = np.stack(
        [
            hdu.data for hdu in hdus
            if float(hdu.header["STRENGTH"]) >= min_strength
        ],
        axis=0,
    )

    print(f"Combining {len(data_stack)} images")

    fits.PrimaryHDU(
        header=hdus[-1].header,
        data=np.nansum(data_stack, axis=0),
    ).writeto(
        f"{prefix}-sum.fits",
        overwrite=True,
    )
    fits.PrimaryHDU(
        header=hdus[-1].header,
        data=np.nanmedian(data_stack, axis=0),
    ).writeto(
        f"{prefix}-median.fits",
        overwrite=True,
    )


if __name__ == "__main__":
    typer.run(combine_maps)
