from typing import Union
from pathlib import Path
import numpy as np
from mpdaf.obj import Cube, Spectrum, Image
from astropy.io import fits
import typer
import slugify

def main(
        cube_file: str,
        mask_file: str,
        output_id: str,
):
    """Apply an image mask to a cube. Set voxels to NaN where mask is zero/False"""

    # Read the spectral cube
    cube = Cube(cube_file)

    # Read the spaxel mask
    mask_image = Image(mask_file)

    # Extend cube mask to include everywhere that the spaxel mask is
    # false
    cube.mask = cube.mask | ~mask_image.data.astype(bool)[:, ...]
    # Also, mask out strange values from cube
    cube.mask = cube.mask | (cube.data == 0.0) | (cube.data == -1.0)

    # Write the result to the same dir that the cube came from
    cube.write(cube_file.replace(".fits", f"-{output_id}.fits"), savemask="nan")

if __name__ == "__main__":
    typer.run(main)
