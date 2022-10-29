from typing import Union
from pathlib import Path
import numpy as np
from mpdaf.obj import Cube, Spectrum, Image
from astropy.io import fits
import typer
import regions as rg
import slugify

def trim_pixel_border_from_mask(mask: np.ndarray, border: int) -> np.ndarray:
    """Set to False all pixels within certain border of mask array"""
    assert mask.ndim == 2
    mask2 = np.zeros_like(mask)
    # Easier to set the pixels that we do not want to set to FALSE
    mask2[border:-border,  border:-border] = True
    # and just AND it with the original
    return mask & mask2



def main(
        zone_indices_file: str,
        zone_list: str="0 I II III IV MYSO S",
        bright_threshold: float=0.3,
        trim_border: int=5,
):
    """Write a fits image file for each zone mask"""

    # Read in the map of zone indices
    zone_index_hdu = fits.open(zone_indices_file)[0]
    # What zones do we have?
    zone_labels = zone_list.split(" ")

    # For each zone, make a mask that selects only that zone's pixels,
    # then write it to a file
    for izone, label in enumerate(zone_labels):
        mask = (zone_index_hdu.data == izone) 
        zone_bright_map = fits.open(f"zone-{label}-bright-map.fits")[0].data
        # Additionally require brightness to exceed some threshold
        mask = mask & (zone_bright_map > bright_threshold)
        # And trim around the border to avoid noisy pixels
        mask = trim_pixel_border_from_mask(mask, trim_border)
        maskfilename = f"zone-{label}-mask.fits"
        fits.PrimaryHDU(
            header=zone_index_hdu.header,
            data=mask.astype(int),
        ).writeto(maskfilename, overwrite=True)
        print("Saved mask to", maskfilename)

if __name__ == "__main__":
    typer.run(main)
