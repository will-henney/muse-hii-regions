import numpy as np
import sys
from pathlib import Path
import typer
import yaml
import slugify
from text_unidecode import unidecode
from astropy.io import fits
from astropy.wcs import WCS

unwanted_types = ["sky", "telluric", "noise", "nan", "stellar"]


def get_line_type(s):
    if s is None:
        return None
    ltype = slugify.slugify(str(s).rstrip("?"))
    if not ltype or ltype in unwanted_types:
        return None
    else:
        return ltype


def load_cube_hdu(
    cube_path: Path,
):
    return fits.open(cube_path)[0]


def get_id_string(data):
    s = f"{data['Index']:04d}-"
    s += slugify.slugify(data["ID"])
    if "UIL" in data["ID"]:
        s += "-" + slugify.slugify(f"{data['lambda_HM']:.2f}")
    return s


def main(
    cube_file: str = "n346-muse-2pass-csub-007.fits",
    project_root: Path = Path("../"),
    yaml_folder: str = "data/n346-lines/all-lines-orig",
    maps_folder: str = "data/n346-bow-lines",
):
    """Create ABC channel maps of all emission lines from data in YAML files"""


    # Get the spectral cube
    big_data_folder = project_root / "big-data" / "ngc346new"
    cube_path = big_data_folder / f"{cube_file}"
    cube = load_cube_hdu(cube_path)

    # The YAML files contain metadata for each line
    yaml_files = (project_root / yaml_folder).glob("*.yaml")
    # Loop over all the lines
    for yaml_file in yaml_files:
        with open(yaml_file) as f:
            metadata = yaml.safe_load(f)
        # Group all lines of same type into their own folder
        line_type = get_line_type(metadata["Type"])
        if line_type is None:
            # Skip unwanted types
            continue

        # Make folder for this line type if necessary
        save_path = project_root / maps_folder / f"maps-{cube_path.stem}" / f"type-{line_type}" 
        save_path.mkdir(exist_ok=True, parents=True)

        ipeak = metadata["Index"]
        # Save each of 3 channels as a separate FITS file
        images = {}
        for chan_label, ichan in zip(["A", "B", "C"], [ipeak - 1, ipeak, ipeak + 1]):
            images[chan_label] = cube.data[ichan, ...]
        # And also the moments
        images["ABC"] = images["A"] + images["B"] + images["C"]
        images["m1"] = (images["C"] - images["A"]) / images["ABC"]
        images["m2"] = (images["C"] + images["A"]) / images["ABC"]
        header = WCS(cube.header).celestial.to_header()
        # FITS headers allow only ASCII strings
        header.update({k: unidecode(str(v)) for k, v in metadata.items()})
        for label, image in images.items():
            fits_file = get_id_string(metadata) + f"-{label}.fits"
            fits.PrimaryHDU(header=header, data=image).writeto(
                save_path / fits_file, overwrite=True
            )
        print("Image saved to", save_path / f"{get_id_string(metadata)}-*.fits")


if __name__ == "__main__":
    typer.run(main)
