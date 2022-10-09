from mpdaf.obj import Spectrum
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import typer

def main(
        region: str,
        peak_suffix: str="p0010-d0030",
        linthresh: float=20.0,
):
    """Plot of spectra with identified peaks"""
    fig, ax = plt.subplots(figsize=(50, 3))
    figfile = f"spec1d-peaks-{region}-{peak_suffix}.pdf"


    spec = Spectrum(f"n346-nostar-{region}.fits")
    spec_bg = Spectrum(f"n346-nostar-{region}-bg.fits")
    spec.plot(label=region, linewidth=1)
    spec_bg.plot(label=f"{region} BG", linewidth=0.5)

    tab = Table.read(
        f"n346-nostar-{region}-peaks-{peak_suffix}.csv",
        format="ascii.ecsv",
    )
    ax.scatter("Wavelength", "prominences", data=tab,
               marker="x", color="r", s=15)
    ax.legend(ncol=4)
    ax.minorticks_on()
    ax.grid(which="major", linewidth=0.5)
    ax.grid(which="minor", linewidth=0.2)
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.xaxis.set_major_locator(MaxNLocator(100))
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
