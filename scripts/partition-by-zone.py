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
        hdu = fits.open(zone["fits_file"])[0]
        hdu.data = convolve(hdu.data, Gaussian2DKernel(smooth))
        vmin, vmax = np.nanpercentile(hdu.data, zone["percentiles"])
        bright_maps.append((hdu.data - vmin) / (vmax - vmin))
    bright_map_stack = np.stack(bright_maps, axis=0)
    izone_map = np.argmax(bright_map_stack, axis=0).astype(float)
    bmax_map = np.max(bright_map_stack, axis=0)
    izone_map[bmax_map < 0.1] = np.nan
    fits.PrimaryHDU(header=hdu.header, data=izone_map).writeto(output_file, overwrite=True)



if __name__ == "__main__":
    typer.run(main)
