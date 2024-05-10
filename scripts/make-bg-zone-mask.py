from pathlib import Path
import numpy as np
from astropy.io import fits
import yaml
import typer

def trim_pixel_border_from_mask(mask: np.ndarray, border: int) -> np.ndarray:
    """Set to False all pixels within certain border of mask array"""
    assert mask.ndim == 2
    mask2 = np.zeros_like(mask)
    # Easier to set the pixels that we do not want to set to FALSE
    mask2[border:-border,  border:-border] = True
    # and just AND it with the original
    return mask & mask2

BG_THRESHOLDS = [("III", -0.1), ("S", 0.2)]


def main(
        trim_border: int=5,
):
    """Write a fits image file for the BG zone mask"""

    # Initialise the mask to all True
    mask = None
    # For each criterion, set to False any pixels that do not meet it
    for zlabel, threshold in BG_THRESHOLDS:
        zone_bright_map = fits.open(f"zone-{zlabel}-bright-map.fits")[0]
        if mask is None:
            mask = zone_bright_map.data < threshold
        else:
            mask = mask & (zone_bright_map.data < threshold)
    # And trim around the border to avoid noisy pixels
    mask = trim_pixel_border_from_mask(mask, trim_border)
    maskfilename = f"zone-BG-mask.fits"
    fits.PrimaryHDU(
        header=zone_bright_map.header,
        data=mask.astype(int),
    ).writeto(maskfilename, overwrite=True)
    print("Saved mask to", maskfilename)

if __name__ == "__main__":
    typer.run(main)
