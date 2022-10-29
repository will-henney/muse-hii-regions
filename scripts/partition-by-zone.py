import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
import yaml
import typer

def main(
        zone_file: str="zones.yaml",
        output_file: str="zone_indices.fits",
        smooth: float=1.5,
):
    """Divide map into different zones according to which line type predominates"""
    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    bright_maps = []
    for zone in zones:
        # Load typical brightness map for this zone
        hdu = fits.open(zone["fits_file"])[0]
        if hdu.data is None:
            hdu = fits.open(zone["fits_file"])[1]
        # Smooth it a bit
        hdu.data = convolve(hdu.data, Gaussian2DKernel(smooth))
        # Find brightness limits corresponding to per-zone specified percentiles
        vmin, vmax = np.nanpercentile(hdu.data, zone["percentiles"])
        # Linear rescaling of limits to range [0, 1]
        bright_map = (hdu.data - vmin) / (vmax - vmin)
        # Construct file name to save normalized brightness of this zone
        bright_file = f"{zone['label']}-bright-map.fits"
        # Save the normalized brightness map
        fits.PrimaryHDU(header=hdu.header, data=bright_map).writeto(bright_file, overwrite=True)
        # And add to the list
        bright_maps.append(bright_map)

    # Make three-dimensional stack of maps
    bright_map_stack = np.stack(bright_maps, axis=0)
    # Find which map in the stack is brightest for each pixel
    izone_map = np.argmax(bright_map_stack, axis=0).astype(float)
    # Make a combined map of all these maximum values
    bmax_map = np.max(bright_map_stack, axis=0)
    # And use it to eliminate pixels that are very faint in ALL maps
    izone_map[bmax_map < 0.0] = np.nan
    # Save the map of zone indices
    fits.PrimaryHDU(header=hdu.header, data=izone_map).writeto(output_file, overwrite=True)



if __name__ == "__main__":
    typer.run(main)
