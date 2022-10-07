from mpdaf.obj import Spectrum
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import typer

def main(region: str, linthresh: float=20.0):
    """Comparison plot of with/without stars, plus BG"""
    fig, ax = plt.subplots(figsize=(50, 3))
    figfile = f"spec1d-compare-{region}.pdf"
    spec = Spectrum(f"n346-nostar-{region}.fits")
    spec_nomask = Spectrum(f"n346-all-{region}.fits")
    spec_bg = Spectrum(f"n346-nostar-{region}-bg.fits")
    spec_bgsub = spec - spec_bg
    spec.plot(label="star mask", linewidth=2)
    spec_nomask.plot(label="no mask", linewidth=1)
    spec_bg.plot(label="BG star mask", linewidth=0.5)
    spec_bgsub.plot(label="BG-subtracted", linewidth=0.5)
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
