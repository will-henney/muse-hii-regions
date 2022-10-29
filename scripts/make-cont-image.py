from typing import Union
from pathlib import Path
import numpy as np
from mpdaf.obj import Cube, Spectrum, Image
from astropy.io import fits
import typer
import slugify

def main(
        cube_file: str,
        output_id: str,
        wave_range: tuple[float, float]=(4600.0, 9300.0),
        percentiles: tuple[float, float]=(5.0, 95.0),
):
    """Write a continuum image scaled between percentiles"""

    # Read the spectral cube
    cube = Cube(cube_file)
    image = cube.get_image(wave_range)
    vmin, vmax = np.nanpercentile(image.data, percentiles)
    image = (image - vmin) / (vmax - vmin)
    image.write(f"cont-image-{output_id}.fits")

if __name__ == "__main__":
    typer.run(main)
