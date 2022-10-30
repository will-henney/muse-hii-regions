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



def main(
        zone_indices_file: str,
        zone_file: str="zones.yaml",
        trim_border: int=5,
):
    """Write a fits image file for each zone mask"""

    # Read in the map of zone indices
    zone_index_hdu = fits.open(zone_indices_file)[0]
    # And erad in the zone metadata
    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    # For each zone, make a mask that selects only that zone's pixels,
    # then write it to a file
    for izone, zone in enumerate(zones):
        _, label = zone["label"].split("-")
        mask = (zone_index_hdu.data == izone)
        zone_bright_map = fits.open(f"zone-{label}-bright-map.fits")[0].data
        # Additionally require brightness to exceed some threshold
        mask = mask & (zone_bright_map > zone["min_bright"])
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
