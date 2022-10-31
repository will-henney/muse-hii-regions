from mpdaf.obj import Spectrum
import numpy as np
from pathlib import Path
import yaml
import typer
import slugify

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
        line_file: str,
        spec_id_label: str="c007-chop-mean",
        zone_file: str="zones.yaml",
        zone_spectra_dir: str="zone_spectra",
):
    """Find line strength of a single line from all zones"""

    with open(line_file) as f:
        metadata = yaml.safe_load(f)
    ipeak = metadata["Index"]
    metadata["Zones Strength"] = {}
    metadata["Zones BG"] = {}
    metadata["Zones BG sig"] = {}

    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    for zone in zones:
        specfile = f"{zone_spectra_dir}/{zone['label']}-{spec_id_label}-spec1d.fits"
        spec = Spectrum(str(specfile))
        line_window = spec.data[ipeak-1:ipeak+2]
        bg_window = spec.data[[ipeak-3, ipeak-2, ipeak+2, ipeak+3]]
        bg_mean = np.mean(bg_window)
        bg_sig = np.std(bg_window)
        line_sum = np.sum(line_window)
        label = zone["label"]
        metadata["Zones Strength"][label] = sanitize(line_sum)
        # Average background extended over the line window
        metadata["Zones BG"][label] = sanitize(bg_mean * len(line_window))
        # Standard error of the same
        metadata["Zones BG sig"][label] = sanitize(
            bg_sig * len(line_window) / np.sqrt(len(bg_window) - 1)
        )

    out_dir = Path(f"all-lines-{spec_id_label}")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (get_id_string(metadata) + ".yaml")
    print("Saving line data to", out_path)
    with open(out_path, "w") as f:
        yaml.dump(metadata, f)

if __name__ == "__main__":
    typer.run(main)
