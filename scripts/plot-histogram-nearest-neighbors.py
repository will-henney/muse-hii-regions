from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore


types = [
    "Identified",
    "Unidentified",
]
colors = ["r", "k"]
sizes = [50, 25]
lmin, lmax = 0.0, 3.0

def poisson(wn, wn0):
    return 1 - np.exp(-2 * wn / wn0)

def main(
        id_label: str,
        nbins: int=10,
):
    """Plot histograms of the nearest neighbor distribution"""

    sns.set_color_codes()
    fig, axes = plt.subplots(3,  1, figsize=(8, 6), sharex=True)
    figfile = f"hist-nearest-neighbors-{id_label}.pdf"
    bins = np.logspace(lmin, lmax, nbins)
    dl = (lmax - lmin) / nbins
    centers = (bins[:-1] + bins[1:]) / 2
    finegrid = np.logspace(0.0, 3.0, 200)

    for ax, _type, color, size in zip(axes[1:], types, colors, sizes):
        df = pd.read_csv(f"all-lines-{id_label}/nearest-neighbors-{_type}.csv").set_index("Index")
        df = df[1:-1]
        H, _, _ = ax.hist(df.dwave_nn, bins=bins, color=color, cumulative=True, label=_type)
        # mean distance between lines
        wave0 = np.mean(df.dwave_mean)
        #wn0 = np.median(df.dwn)
        ax.plot(
            finegrid,
            H[-1] * poisson(finegrid, wave0),
            color="y",
            linewidth=2,
            label=f"{wave0:.1f} Å$^{{-1}}$",
        )

        HH, _ = np.histogram(df.dwave_nn[df.mutual], bins=bins)
        P1 = HH.cumsum() / H
        axes[0].plot(centers, P1, c=color, label=_type)
        ax.legend()
        ax.set_ylabel("# of lines")

    axes[-1].set(
        xscale="log",
        xlabel=r"Nearest neighbor wavenumber difference, $\mathrm{Å}^{-1}$",
    )
    axes[0].axhline(2/3, color="b", linestyle="dashed")
    axes[0].legend()
    axes[0].set(
        ylim=[0.45, 1.1],
        ylabel="Degree of reciprocity",
    )
    sns.despine()
    fig.savefig(figfile)
    print(figfile, end="")





if __name__ == "__main__":
    typer.run(main)
