from typing import Union
from pathlib import Path
import numpy as np
from mpdaf.obj import Cube, Spectrum, Image
from astropy.io import fits
import typer
import regions as rg
import slugify

def get_spectrum_from_masked_cube(
        cube: Cube,
        spaxel_mask: np.ndarray,
        reduction_method: callable = np.nanmean,
        debug: bool = False,
) -> Spectrum:
    """Extract 1D spectrum from cube, but only spaxels where mask is True"""
    nv, ny, nx = cube.shape
    # Check that mask is 2D and the right size
    assert spaxel_mask.shape == (ny, nx)
    spec = reduction_method(
        # If we index by the mask over last two axes, then we get a 1d
        # array of spectra, which we can combine
        cube.data[:, spaxel_mask],
        axis=-1,
    )
    return Spectrum(wave=cube.wave, data=spec, unit=cube.unit)


def main(
        zone_indices_file: str,
        cube_file: str,
        exclude_mask_file: str,
        output_id: str,
        output_dir: str="zone_spectra",
        zone_list: str="0 I II III IV",
):
    """Extract 1D spectra from cube for each region in file"""

    # Read in the map of zone indices
    zone_index_map = fits.open(zone_indices_file)[0].data
    # What zones do we have?
    zone_labels = zone_list.split(" ")
    # For each zone, make a mask that selects only that zone's pixels
    zone_masks = [zone_index_map == izone for izone in range(len(zone_labels))]

    # Read the spectral cube
    cube = Cube(cube_file)
    # Nan-ify the data array, since we do not use the mask array again later
    cube.data[cube.mask] = np.nan

    # Set the extra image mask to be true where the exclude_mask_file
    # image is zero
    extra_image_mask = np.where(
        Image(exclude_mask_file).data == 0.0,
        True,
        False
    )

    # Now do the get the 1D spectrum for each zoen
    spec_dict = {
        label: get_spectrum_from_masked_cube(
            cube,
            zone_mask | extra_image_mask
        )
        for label, zone_mask in zip(zone_labels, zone_masks)
    }

    # Make sure the output folder exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    # And save each spectrum as a separate FITS file
    for label, spec in spec_dict.items():
        spec.write(f"{output_dir}/zone-{label}-{output_id}-spec1d.fits")


if __name__ == "__main__":
    typer.run(main)
