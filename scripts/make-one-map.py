import numpy as np
import sys
from pathlib import Path
import typer
import openpyxl
import yaml
import slugify
from text_unidecode import unidecode
from astropy.io import fits
from astropy.wcs import WCS

unwanted_types = ["sky",  "telluric", "noise", "nan"]

def get_line_type(s):
    ltype = slugify.slugify(str(s).rstrip("?"))
    if ltype in unwanted_types:
        return None
    else:
        return ltype


def load_cube_hdu(
        cwindow: int,
        prefix: str="n346-muse-csub",
        big_data_folder: Path=Path("../../big-data/ngc346new"),
):
    cube_path = big_data_folder / f"{prefix}-{cwindow:03d}.fits"
    return fits.open(cube_path)[0]

def get_id_string(data):
    s = f"{data['Index']:04d}-"
    s += slugify.slugify(data["ID"])
    if "UIL" in data["ID"]:
        s += "-" + slugify.slugify(f"{data['lambda_HM']:.2f}")
    return s

def choose_cont_window(data: dict) -> tuple[int, bool]:
    """Decide which type of continuum subtraction is preferred

    Returns tuple: width of window, and whether to subtract baseline
    """

    # Case of no preference given
    if not data["Cont_method"]:
        # Just use the wide window
        return 101, False

    try:
        # Case of only one method listed and it is an integer
        return int(data["Cont_method"]), False
    except ValueError:
        # Case of various methods listed, or one that contains letters. Take the first
        cont_methods = data["Cont_method"].split(",")
        first_cont_method = cont_methods[0]
        if first_cont_method.endswith("B"):
            # Case that we want to subtract the baseline
            return int(first_cont_method.rstrip("B")), True
        else:
            # Case that we do not
            return int(first_cont_method), False



def main(
        yaml_file : str,
):
    """Create map of a single emission line from data in YAML file
    """
    if not yaml_file.endswith(".yaml"):
        yaml_file = yaml_file + ".yaml"
    with open(yaml_file) as f:
        metadata = yaml.load(f)

    # Group all lines of same type into their own folder
    line_type = get_line_type(metadata["Type"])
    save_path = Path("type-" + line_type)
    save_path.mkdir(exist_ok=True)

    cwindow, yes_sub_base = choose_cont_window(metadata)
    cube = load_cube_hdu(cwindow)
    ipeak = metadata["Index"]
    # First try: just use 3 pixels along wave axis
    cube_window = cube.data[ipeak-1:ipeak+2, ...]
    if yes_sub_base:
        # This will fail if the line is broad
        base = 0.5 * (cube.data[ipeak-2, ...] + cube.data[ipeak+2, ...])
        cube_window -= base
    image = np.sum(cube_window, axis=0)
    header = WCS(cube.header).celestial.to_header()
    # FITS headers allow only ASCII strings
    header.update({k: unidecode(str(v)) for k, v in metadata.items()})

    fits_file = get_id_string(metadata) + ".fits"
    fits.PrimaryHDU(header=header, data=image).writeto(save_path / fits_file, overwrite=True)
    print("Image saved to", save_path / fits_file)

if __name__ == "__main__":
    typer.run(main)
