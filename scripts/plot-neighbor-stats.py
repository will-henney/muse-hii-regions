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

# 
wavgrid = np.arange(4600, 9400, dtype=float)

def finesse_regrid(waves, values, window_size=None):
    """Fix up a quantity as a function of wavelength"""
    # Interpolate onto uniform grid
    vgrid = np.interp(wavgrid, waves.to_numpy(), values.to_numpy())
    if window_size is not None:
        vgrid = pd.Series(vgrid).rolling(window=window_size, center=True).mean()
    return vgrid

def finesse(waves, values, window_size):
    """Fix up a quantity as a function of wavelength"""
    return values.rolling(window=window_size, center=True).mean()

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
    fig, axes = plt.subplots(3,  1, figsize=(8, 8), sharex=True)
    figfile = f"nearest-neighbor-stats-{id_label}.pdf"
    for _type, type_data in info["types"].items():
        if not _type in types:
            continue
        cmap = sns.dark_palette(
            tuple(type_data["husl"]),
            input="husl",
            as_cmap=True,
        )
        df = pd.read_csv(
            f"all-lines-{id_label}/nearest-neighbors-{_type}.csv"
        ).set_index("Index").sort_values(by="wave")
        # Make window-averaged arrays to plot
        wave = finesse(df.wave, df.wave, window_size)
        density = finesse(df.wave, 1 / df.dwave_mean, window_size)
        ratio_nn = finesse(df.wave, df.dwave_nn / df.dwave_mean, window_size)
        reciprocity = finesse(df.wave, df.mutual.astype(float), window_size)

        axes[0].plot(wave, density, color=cmap(1.0), label=_type)
        axes[1].plot(wave, ratio_nn, color=cmap(1.0))
        axes[2].plot(wave, reciprocity, color=cmap(1.0))



    axes[-1].set(
        xlabel=r"Wavelength, Angstrom",
    )

    axes[0].set(
        ylim=[0.0, None],
        yscale="linear",
        ylabel="Spectral density:\n" + r"$\rho_\lambda$, Lines / Å",
    )

    axes[1].axhline(1/2, color="b", linestyle="dashed")
    axes[1].set(
        ylim=[0.0, 1.05],
        ylabel="Nearest-neighbor factor:\n" + r"$\rho_\lambda \times \langle \Delta\lambda_{\mathrm{NN}} \rangle$",
    )

    axes[2].axhline(2/3, color="b", linestyle="dashed")
    axes[2].set(
        ylim=[0.0, 1.05],
        ylabel="Degree of reciprocity",
    )


    axes[0].legend()
    for ax in axes:
        ax.grid(axis="x", color="k", linestyle="dotted")
    sns.despine()
    fig.savefig(figfile)
    print(figfile, end="")





if __name__ == "__main__":
    typer.run(main)
