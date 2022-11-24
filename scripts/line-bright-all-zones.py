from mpdaf.obj import Spectrum
import numpy as np
from pathlib import Path
import sys
import yaml
import typer
import slugify
from astropy.modeling import models, fitting
from discrete_gaussian_model import DiscreteGaussianModel

INDEX_PATTERN = "[0-9]" * 4
FITTER = fitting.LevMarLSQFitter()

ZONES_TO_FIX_SKY = ["0", "I", "II", "III", "MYSO"]
SKY_SPECIES = ["OH", "O_2"]

def get_id_string(data):
    s = f"{data['Index']:04d}-{str(data['Type']).rstrip('?')}"
    if data["ID"]:
        s += "-" + data["ID"]
        if "UIL" in data["ID"]:
            s += "-" + f"{int(np.round(data['lambda_HM']))}"
    return slugify.slugify(s)

def sanitize(x):
    return float(np.round(x, 4))

def sanitize_array(x: np.ndarray) -> list:
    return [sanitize(_) for _ in x.tolist()]


def fit_gauss7(wave: np.ndarray, spec: np.ndarray, mask: np.ndarray):
    assert wave.shape == spec.shape == mask.shape == (7,)
    assert mask.dtype == bool
    # Make a copy so that we do not affect the global array
    _mask = mask.copy()
    # Unmask the central pixels for this line, so we can fit them
    _mask[2:5] = True
    # Initial guess for Gaussian
    g0 = DiscreteGaussianModel(
        amplitude=spec[3],
        mean=wave[3],
        stddev=1.0,
        bin_width=wave[1] - wave[0],
        fixed={"bin_width": True},
    )
    # Try and fit it
    return FITTER(g0, wave[_mask], spec[_mask])
    # return g0

def is_sky_blend(d: dict):
    """Check if this line is a UIL or similar that is blended with sky"""
    if d.get("ID", "") and d.get("ID", "").endswith("+"):
        # Case where we have a blended line: check that it is with
        # a sky line by looking at the note to the ID column
        try:
            # We need to look in first element of list inside dict inside dict
            is_blend = any(s in d["Notes"]["ID"][0] for s in SKY_SPECIES)
        except (KeyError, IndexError):
            # But either of the dicts may be absent
            is_blend = False
        if d["Type"].startswith(("Med", "High")):
            # Also, do not try this with Medium or higher
            # ionization lines, since they will have real emission
            # in Zone IV
            is_blend = False
    else:
        # Or, not even a blended line
        is_blend = False
    return is_blend


def main(
        spec_id_label: str="c007-chop-mean",
        orig_data_dir: str="all-lines-orig",
        zone_file: str="zones.yaml",
        zone_spectra_dir: str="zone_spectra",
        debug: bool=False,
):
    """Find line strength of all lines from all zones"""

    # First, get the zones
    with open(zone_file) as f:
        zones = yaml.safe_load(f)
    # And load the corresponding spectra
    nwave = None
    for zone in zones:
        specfile = f"{zone_spectra_dir}/{zone['label']}-{spec_id_label}-spec1d.fits"
        zone["spec"] = Spectrum(str(specfile))
        # Save the Zone IV spec when we find it, since we will use it later to fix the sky
        if zone["label"] == "zone-IV":
            spec_IV = zone["spec"]
        if nwave is None:
            nwave = len(zone["spec"].data)
        else:
            assert len(zone["spec"].data) == nwave

    # Next, get all the lines
    line_files = sorted(Path(orig_data_dir).glob(f"{INDEX_PATTERN}.yaml"))

    # And the list of indices to avoid when calculating the BG
    avoid_indices = np.loadtxt(
        Path(orig_data_dir) / "line-indices.txt",
        dtype=int,
    )
    # Make a 1D mask of acceptable BG pixels
    bg_mask = np.ones((nwave,), bool)
    bg_mask[avoid_indices] = False

    # And an array of the wave pixel indices
    indices = np.arange(nwave)

    # Now loop over all lines
    for line_file in line_files:
        with open(line_file) as f:
            metadata = yaml.safe_load(f)
        ipeak = metadata["Index"]
        # Check if we are blended with sky
        metadata["sky_blend"] = is_sky_blend(metadata)
        # Get the strength for each zone
        for zone in zones:
            spec = zone["spec"].copy()
            # Case where we have a sky blend
            if metadata["sky_blend"] and zone["label"].split("-")[-1] in ZONES_TO_FIX_SKY:
                # Use Zone IV to model the sky
                spec = spec - spec_IV

            # Wide window includes BG and line
            win7 = spec.data[ipeak-3:ipeak+4]
            # Corresponding wavelengths, converted from m to Angstrom
            wave7 = 1.0e10 * spec.wave.coord(indices[ipeak-3:ipeak+4])
            # Corresponding mask
            m7 = bg_mask[ipeak-3:ipeak+4]
            # Narrow window includes only line
            win3 = spec.data[ipeak-1:ipeak+2]
            bg_npix = np.sum(m7)
            if bg_npix > 0:
                bg_mean = np.mean(win7[m7])
                bg_sig = np.std(win7[m7])
            else:
                bg_mean = bg_sig = 0.0
            line_sum = np.sum(win3 - bg_mean)
            gfit = fit_gauss7(wave7, win7 - bg_mean, m7)
            metadata[zone["label"]] = {
                "Strength": sanitize(line_sum),
                "Sigma": sanitize(bg_sig),
                "BG npix": int(bg_npix),
                "BG mean": sanitize(bg_mean),
                "Pixel Wave": sanitize(wave7[3]),
            }
            metadata[zone["label"]]["Gauss Fit"] = {
                "Amplitude": sanitize(gfit.amplitude.value),
                "Mean Wave": sanitize(gfit.mean.value),
                "RMS Width": sanitize(gfit.stddev.value),
            }
            if debug:
                metadata[zone["label"]]["Window"] = {
                    "Wave": sanitize_array(wave7),
                    "Spectrum": sanitize_array(win7 - bg_mean),
                    "Fit":  sanitize_array(gfit(wave7)),
                    "BG Mask": m7.tolist(),
                }
                #metadata[zone["label"]]["Fit Info"] = str(FITTER.fit_info)

        out_dir = Path(f"all-lines-{spec_id_label}")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / (get_id_string(metadata) + ".yaml")
        with open(out_path, "w") as f:
            yaml.dump(metadata, f)

if __name__ == "__main__":
    typer.run(main)
