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
#     "Cross",
]
colors = [
    "r",
    "k",
    # "c",
]

# 
wavgrid = np.arange(4600, 9400, dtype=float)

def finesse_rolling(waves, values, window_size):
    """Take rolling mean in window"""
    return values.rolling(window=window_size, center=True)

def finesse_reindex(waves, values, window_size):
    """Reindex, then interpolate, then rolling window"""
    return (values
            .reindex(range(values.index.min(),  values.index.max()))
            .interpolate(method="linear")
            .rolling(window=window_size, center=True)
            )

finesse = finesse_rolling
#finesse = finesse_reindex

def main(
        id_label: str,
        window_size: int=10,
        species_file: str="species.yaml",
):
    """Plot statistics of the nearest neighbor distances versus wavelength"""
    with open(species_file) as f:
        info = yaml.safe_load(f)
    # Drop the Deep type since it only has one line
    info["types"].pop("Deep")
    # But add in the combination of all Identified lines
    info["types"]["Identified"] = info["types"]["Medium"]

    # Three panel stack of plots
    sns.set_color_codes()
    fig, axes = plt.subplots(3,  1, figsize=(5, 5), sharex=True)
    figfile = f"nearest-neighbor-stats-{id_label}.pdf"
    for _type, color in zip(types, colors):
        df = pd.read_csv(
            f"all-lines-{id_label}/nearest-neighbors-{_type}.csv"
        ).set_index("Index").sort_values(by="wave")
        # Make window-averaged arrays to plot
        wave = finesse(df.wave, df.wave, window_size).mean()
        density = finesse(df.wave, 1 / df.dwave_mean, window_size)
        ratio_nn = finesse(df.wave, df.dwave_nn / df.dwave_mean, window_size)
        reciprocity = finesse(df.wave, df.mutual.astype(float), window_size)

        axes[0].plot(wave, density.mean(), ds="steps-mid", color=color, label=f"{_type} lines")
        sem = density.std() / np.sqrt(window_size - 1)
        axes[0].fill_between(
            wave,
            density.mean() - sem,
            density.mean() + sem,
            step="mid",
            color=color,
            alpha=0.2,
            linewidth=0.0,
        )

        axes[1].plot(wave, ratio_nn.mean(), ds="steps-mid", color=color)
        sem = ratio_nn.std() / np.sqrt(window_size - 1)
        axes[1].fill_between(
            wave,
            ratio_nn.mean() - sem,
            ratio_nn.mean() + sem,
            step="mid",
            color=color,
            alpha=0.2,
            linewidth=0.0,
        )

        axes[2].plot(wave, reciprocity.mean(), ds="steps-mid", color=color)
        sem = reciprocity.std() / np.sqrt(window_size - 1)
        axes[2].fill_between(
            wave,
            reciprocity.mean() - sem,
            reciprocity.mean() + sem,
            step="mid",
            color=color,
            alpha=0.2,
            linewidth=0.0,
        )



    axes[-1].set(
        xlabel=r"Wavelength, Å",
    )

    axes[0].set(
        ylim=[0.0, None],
        yscale="linear",
    )
    axes[0].set_ylabel(
        "Spectral\ndensity:\n" + r" $\rho_\lambda = 1 / \langle \Delta\lambda \rangle$, Å$^{-1}$",
        rotation="horizontal", ha="center", va="center",
        labelpad=50,
    )

    axes[1].axhline(1/2, color="b", linestyle="dashed")
    axes[1].set(
        ylim=[0.0, 1.05],
    )
    axes[1].set_ylabel(
        "Nearest-neighbor\n" + r"ratio: $\langle \Delta\lambda_{\mathrm{nn}} \rangle / \langle \Delta\lambda \rangle$",
        rotation="horizontal", ha="center", va="center",
        labelpad=50,
    )

    axes[2].axhline(2/3, color="b", linestyle="dashed")
    axes[2].set(
        ylim=[0.0, 1.05],
    )
    axes[2].set_ylabel(
        "Degree of\nreciprocity:\n" + r"$P_{\mathrm{nn-nn}}$",
        rotation="horizontal", ha="center", va="center",
        labelpad=50,
    )



    axes[0].legend()
    for ax in axes:
        ax.grid(axis="x", color="k", linestyle="dotted")
    sns.despine()
    fig.savefig(figfile, bbox_inches="tight")
    print(figfile, end="")





if __name__ == "__main__":
    typer.run(main)
