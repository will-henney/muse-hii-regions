from typing import Union
from pathlib import Path
import numpy as np
from mpdaf.obj import Cube, Spectrum, Image
import typer
import regions as rg
import slugify

def get_spectrum_from_region(
        cube: Cube,
        region: Union[rg.PixelRegion, rg.SkyRegion],
        reduction_method: callable = np.nanmean,
        extra_image_mask: Union[np.ndarray, None] = None,
        debug: bool = False,
) -> Spectrum:
    """Extract 1D spectrum from region in an efficient way"""
    try:
        region_mask = region.to_mask()
    except AttributeError:
        region_mask = region.to_pixel(cube.wcs.wcs).to_mask()
    nv, ny, nx = cube.shape
    if extra_image_mask is not None:
        assert extra_image_mask.shape == ny, nx
    # Slices into 2D arrays
    slices_large, slices_small = (
        region_mask
        .get_overlap_slices((ny, nx))
    )
    if debug:
        print('2D slice:', slices_large)
    slices_cube = (slice(None, None),) + slices_large
    image_mask_large = region_mask.to_image((ny, nx))
    if extra_image_mask is not None:
        # extra_image_mask should be true for pixels that we want to
        # include
        image_mask_large[~extra_image_mask] = np.nan

    image_mask_small = image_mask_large[slices_large]
    cube_cutout = cube.data[slices_cube]
    cube_cutout[cube.mask[slices_cube]] = np.nan
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
        out_path: Path=Path.cwd(),
        exclude_mask_file: Union[str, None]=typer.Option(
            None,
            help="Optional file with image to mask out data. Only pixels with zero will be included."
        ),
):
    """Extract 1D spectra from cube for each region in file"""

    sky_regions = rg.Regions.read(region_file)
    region_dict = {reg.meta["label"]: reg for reg in sky_regions}

    cube = Cube(str(data_path / cube_name))

    # Set the extra image mask to be true where the exclude_mask_file image is zero
    if exclude_mask_file is not None:
        extra_image_mask = np.where(
            Image(exclude_mask_file).data == 0.0,
            True,
            False
        )
    else:
        extra_image_mask = None

    # Now do the work to get the spectra
    spec_dict = {
        label: get_spectrum_from_region(cube, reg, extra_image_mask=extra_image_mask)
        for label, reg in region_dict.items()
    }

    # And save them all
    for label, spec in spec_dict.items():
        spec.write(str(out_path / f"{out_prefix}-{slugify.slugify(label)}.fits"))


if __name__ == "__main__":
    typer.run(main)
