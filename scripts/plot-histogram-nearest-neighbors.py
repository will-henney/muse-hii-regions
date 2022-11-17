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
    "Cross",
]
colors = ["r", "k", "c"]
lmin, lmax = -0.1, 2.9
wave_floor = 3.0
clabels = ["Identified lines", "Unidentified lines", "Cross comparison"]
labels = ["_nolabel", "_nolabel", "Poisson process"]

def poisson(wn, wn0):
    return 1 - np.exp(-2 * wn / wn0)
def main(
        id_label: str,
        nbins: int=10,
):
    """Plot histograms of the nearest neighbor distribution"""

    sns.set_color_codes()
    fig, axes = plt.subplots(len(types),  1, figsize=(5, 5), sharex=True)
    figfile = f"hist-nearest-neighbors-{id_label}.pdf"
    bins = np.logspace(lmin, lmax, nbins)
    dl = (lmax - lmin) / nbins
    centers = (bins[:-1] + bins[1:]) / 2
    finegrid = np.logspace(lmin, lmax, 200)

    wave0 = {}
    for ax, _type, color, label, clabel in zip(axes, types, colors, labels, clabels):
        df = pd.read_csv(f"all-lines-{id_label}/nearest-neighbors-{_type}.csv").set_index("Index")
        df = df[1:-1]
        H, _, _ = ax.hist(
            df.dwave_nn, bins=bins, color=color, cumulative=True,
            # label=clabel,
        )
        # mean distance between lines
        if _type == "Cross":
            wave0[_type] = wave0["Identified"]
        else:
            wave0[_type] = np.mean(df.dwave_mean)
        #wn0 = np.median(df.dwn)
        ax.plot(
            finegrid,
            H[-1] * np.where(
                finegrid > wave_floor,
                (poisson(finegrid, wave0[_type]) - poisson(wave_floor, wave0[_type]))
                / (1.0 - poisson(wave_floor, wave0[_type])),
                np.nan
            ),
            color="b",
            linestyle="dashed",
            linewidth=3,
            label=label,
            # label=fr"Poisson: $\rho_\lambda = {1/wave0:.3f}$ Å$^{{-1}}$",
        )
        ax.text(
            150, 0.5 * H[-1], clabel,
            ha="center", va="center", fontweight="bold",
            color=color, bbox={"facecolor": "w", "edgecolor": "none", "pad": 2},
        )

        # HH, _ = np.histogram(df.dwave_nn[df.mutual], bins=bins)
        # P1 = HH.cumsum() / H
        # axes[0].plot(centers, P1, c=color, label=_type)
    axes[1].set_ylabel("Cumulative\n# of lines", rotation="horizontal", ha="right", va="center")

    axes[-1].legend(fontsize="small", handlelength=3.0, framealpha=1.0)

    axes[-1].set(
        xscale="log",
        xlabel=r"Nearest neighbor wavelength difference: $\Delta\lambda_{\mathrm{nn}}$, Å",
        xticks=[1, 10, 100],
        xticklabels=["1", "10", "100"],
    )
    # axes[0].axhline(2/3, color="b", linestyle="dashed")
    # axes[0].legend()
    # axes[0].set(
    #     ylim=[0.45, 1.1],
    #     ylabel="Degree of reciprocity",
    # )
    sns.despine()
    fig.savefig(figfile, bbox_inches="tight")
    print(figfile, end="")





if __name__ == "__main__":
    typer.run(main)
