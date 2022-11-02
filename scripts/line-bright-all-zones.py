from mpdaf.obj import Spectrum
import numpy as np
from pathlib import Path
import yaml
import typer
import slugify

INDEX_PATTERN = "[0-9]" * 4

def get_id_string(data):
    s = f"{data['Index']:04d}-{str(data['Type']).rstrip('?')}"
    if data["ID"]:
        s += "-" + data["ID"]
        if "UIL" in data["ID"]:
            s += "-" + f"{int(np.round(data['lambda_HM']))}"
    return slugify.slugify(s)

def sanitize(x):
    return float(np.round(x, 4))

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

    # And the list of 
    # Now loop over all lines
    for line_file in line_files:
        with open(line_file) as f:
            metadata = yaml.safe_load(f)
        ipeak = metadata["Index"]
        # Get the strength for each zone
        for zone in zones:
            spec = zone["spec"]
            # Wide window includes BG and line
            win7 = spec.data[ipeak-3:ipeak+4]
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
            metadata[zone["label"]] = {
                "Strength": sanitize(line_sum),
                "Sigma": sanitize(bg_sig),
                "BG npix": int(bg_npix),
                "BG mean": sanitize(bg_mean),
            }

        out_dir = Path(f"all-lines-{spec_id_label}")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / (get_id_string(metadata) + ".yaml")
        if debug:
            print("Saving line data to", out_path)
        with open(out_path, "w") as f:
            yaml.dump(metadata, f)

if __name__ == "__main__":
    typer.run(main)
