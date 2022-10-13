from mpdaf.obj import Spectrum
import numpy as np
import scipy.signal as si
from  astropy.table import Table
import typer

def fmt_float(x):
    return f"{int(10*x):04d}"


def main(
        spec_file: str,
        min_prominence: float=3.0,
        min_distance: float=3.0,
):
    """Find peaks in a 1D spectrum"""
    assert spec_file.endswith(".fits")
    spec = Spectrum(spec_file)
    # Find pixel positions of peaks
    peaks, props = si.find_peaks(
        spec.data,
        prominence=min_prominence,
        distance=min_distance,
        # FWHM must be between 2 and 10 pixels to allow blends
        width=(1.5, 10.0),
    )
    # Convert to wavelengths
    waves = spec.wave.coord()[peaks]
    # Make a table of the results
    tab = Table(
        {"Wavelength": waves, "Pixel": peaks, **props}
    )
    # Do not use too many decimal places for the float columns
    for col in [
            "Wavelength", "prominences", "widths",
            "width_heights", "left_ips", "right_ips",
    ]:
        tab[col] = np.round(tab[col], 4)
    # And save it as CSV format
    suffix = ("-peaks"
              f"-p{fmt_float(min_prominence)}"
              f"-d{fmt_float(min_distance)}.csv")
    tab_file = spec_file.replace(".fits", suffix)
    tab.write(tab_file, format="ascii.ecsv", overwrite=True)
    print(tab_file)



if __name__ == "__main__":
    typer.run(main)
