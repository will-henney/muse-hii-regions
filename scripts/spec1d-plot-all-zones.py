from mpdaf.obj import Spectrum
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path
import typer

def main(
        id_label: str,
        linthresh: float=20.0,
        star_scale: float=5.0,
):
    """Plot of spectra from zones"""
    fig, ax = plt.subplots(figsize=(50, 3))
    figfile = f"spec1d-all-zones-{id_label}.pdf"

    specfiles = sorted(Path.cwd().glob(f"zone-*-{id_label}-spec1d.fits"))
    nzones = len(specfiles)
    yspan = 1.5 * linthresh
    dy = yspan / (nzones - 1)
    offset = yspan / 2
    for specfile in specfiles:
        zone = specfile.stem.split("-")[1]
        spec = Spectrum(str(specfile))
        scale = star_scale if zone in ["S", "MYSO"] else 1.0
        (spec / scale + offset).plot(label=f"Zone {zone}", linewidth=1)
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
