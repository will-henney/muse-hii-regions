from typing import Union
from pathlib import Path
import numpy as np
from mpdaf.obj import Cube, Spectrum
import typer
import regions as rg

def get_spectrum_from_region(
        cube: Cube,
        region: Union[rg.PixelRegion, rg.SkyRegion],
        reduction_method: callable = np.sum,
        debug: bool = False,
) -> Spectrum:
    try:
        region_mask = region.to_mask()
    except AttributeError:
        region_mask = region.to_pixel(cube.wcs.wcs).to_mask()
    nv, ny, nx = cube.shape
    # Slices into 2D arrays
    slices_large, slices_small = (
        region_mask
        .get_overlap_slices((ny, nx))
    )
    if debug:
        print('2D slice:', slices_large)
    slices_cube = (slice(None, None),) + slices_large
    image_mask_large = region_mask.to_image((ny, nx))
    image_mask_small = image_mask_large[slices_large]
    cube_cutout = cube.data[slices_cube]
    cube_cutout[cube.mask[slices_cube]] = 0.0
    spec = reduction_method(
        cube_cutout * image_mask_small[None, :, :],
        axis=(1, 2),
    )
    return Spectrum(wave=cube.wave, data=spec, unit=cube.unit)


def main(
        region_file: str,
        out_prefix: str="n346-muse",
        data_path: Path=(
            Path.home() / "Work/Muse-Hii-Data/SMC-NGC-346"
        ),
        cube_name: str="ADP.2017-10-16T11_04_19.247.fits",
        star_mask_file: Union[str, None]=None,
):
    """Extract 1D spectra from cube for each region in file"""

    sky_regions = rg.Regions.read(region_file)
    region_dict = {reg.meta["label"]: reg for reg in sky_regions}



if __name__ == "__main__":
    typer.run(main)
