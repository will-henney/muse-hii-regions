from mpdaf.obj import Spectrum
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from pathlib import Path
import yaml
import typer

def main(
        id_label: str,
        linthresh: float=20.0,
        star_scale: float=5.0,
        zone_file: str="../zones.yaml",
):
    """Plot of spectra from zones"""
    fig, ax = plt.subplots(figsize=(100, 3))
    figfile = f"spec1d-all-zones-{id_label}.pdf"

    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    nzones = len(zones)
    yspan = 1.5 * linthresh
    dy = yspan / (nzones - 1)
    offset = yspan / 2
    for zone in zones:
        specfile = f"{zone['label']}-{id_label}-spec1d.fits"
        spec = Spectrum(str(specfile))
        scale = star_scale if zone["label"].endswith(("S", "MYSO")) else 1.0
        color = sns.dark_palette(
            tuple(zone["husl"]),
            input="husl",
        )[-1]
        (spec / scale + offset).plot(label=zone["label"], linewidth=1, color=color)
        ax.axhline(offset, linewidth=0.5, color=color)
        offset -= dy

    ax.legend(ncol=nzones)
    ax.minorticks_on()
    ax.grid(which="major", linewidth=0.5)
    ax.grid(which="minor", linewidth=0.2)
    ax.set_yscale("symlog", linthresh=linthresh, linscale=2.0)
    ax.xaxis.set_major_locator(MaxNLocator(100))
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
