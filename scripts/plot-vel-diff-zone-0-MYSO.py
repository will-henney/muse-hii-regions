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
from astropy.stats import sigma_clipped_stats

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value

REPLACEMENTS = {
    "Deep": "Deep Neutral",
    "Fe": "Fe-Ni-Ca-Si",
}
TYPES = ["Unidentified", "Deep", "Neutral", "Low", "Medium", "Fe"]
META_TYPES = {
    "Unidentified": ["Unidentified"],
    "Neutral-Low": ["Deep", "Neutral", "Low", "Fe"],
    "Medium": ["Medium"], 
}
def main(
        id_label: str,
        minimum_signal_noise: float=1.0,
        signal_noise_log_range: tuple[float, float]=(0.0, 2.0),
        species_file: str="species.yaml",
        zone_file: str="zones.yaml",
):
    """Histograms of velocity difference between Zone 0 and MYSO"""

    # Load the species and line type database
    with open(species_file) as f:
        info = yaml.safe_load(f)
    # Load the zone database
    with open(zone_file) as f:
        zones = yaml.safe_load(f)


    # Load the data for the two zones
    df0 = pd.read_csv(f"all-lines-{id_label}/zone-0-velocities.csv").set_index("Index")
    dfy = pd.read_csv(f"all-lines-{id_label}/zone-MYSO-velocities.csv").set_index("Index")

    # Make a new frame with the columns that we need
    df = df0[["ID", "flux", "s_n"]].assign(
        # Note that we use dfy.wave in denominator instead of dfy.wave0
        # since we do not have the latter for the UILs. The resultant
        # inaccuracy is very small (order v/c ~ 0.001)
        d_vel=LIGHT_SPEED_KMS * (df0.wave - dfy.wave) / dfy.wave,
        # Do not use the original type field, instead make our own
        type="Unknown",
        species="Unknown",
    )

    # Use the species info to fill in the types column
    for species in reversed(info["species"]):
        prefix = species["name"]
        if not species["name"] == "UIL":
            prefix += " "
        mask = df.ID.str.startswith(prefix)
        df.type[mask] = species["type"]
        df.species[mask] = species["name"]

    # Brightness thresholds for partitioning lines
    thresholds = [0.0, 0.1, 1.0, 10.0]
    df = df.assign(fcat=0.0)
    for thresh in thresholds:
        mask = df.flux > thresh
        df.fcat[mask] = thresh

    # The thresholds correspond to rows
    n_rows = len(thresholds)
    # The line types correspond to columns
    n_cols = len(META_TYPES)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        sharex=True, sharey=True,
        figsize=(8, 7),
    )
    figfile = f"vel-diff-zone-0-MYSO-{id_label}.pdf"
    vrange = [-105, 105]
    # Iterate over the flux categories
    d_value = 0.5 / n_rows
    value = 1 -  n_rows * d_value
    for jrow, thresh in enumerate(thresholds):
        mask1 = df.fcat == thresh
        if jrow + 1 < n_rows:
            bstring = f"{thresh} to {thresholds[jrow + 1]}"
        else:
            bstring = f"> {thresh}"
        axes[jrow, -1].text(
            1.3, 0.5, "Brightness\n" + bstring,
            transform=axes[jrow, -1].transAxes,
            ha="center",
        )
        # Iterate over the line types
        for icol, (label, line_types) in enumerate(META_TYPES.items()):
            mask = (
                mask1
                & df.type.str.startswith(tuple(line_types))
                & (df.s_n > minimum_signal_noise)
            )
            data = df[mask].sort_values(by="flux")
            if len(data) == 0:
                # If we have no lines, then skip to the next
                continue
            type_data = info["types"][line_types[0]]
            color = sns.light_palette(
                tuple(type_data["husl"]),
                input="husl",
                as_cmap=True,
            )(value)
            ax = axes[jrow, icol]
            ax.hist(data.d_vel, bins=31, range=vrange, color=color)
            mean, median, std = sigma_clipped_stats(data.d_vel, sigma=2)
            stats_string = fr"${mean:.1f} \pm {std:.1f}$"
            stats_string += "\n" + fr"$N = {len(data)}$"
            ax.text(
                0.05, 0.95, stats_string,
                transform=ax.transAxes,
                ha="left", va="top",
                color=color,
            )

        value += d_value

    # Add labels for the line types at the top of each column
    for itype, meta_type in  enumerate(META_TYPES):
        ax = axes[0, itype]
        ax.text(
            0.5, 1.2, meta_type, transform=ax.transAxes,
            ha="center", va="center",
        )
    # And a title for these
    fig.text(0.5, 0.99, "Emission line type", ha="center", va="top")


    axes[-1, 1].set(
        xlabel=r"Velocity difference: Zone 0 - Zone MYSO (km/s)",
    )
    axes[1, 0].set(
        ylabel=r"Number of emission lines",
    )

    sns.despine()
    #fig.tight_layout()
    #fig.subplots_adjust(left=0.15, bottom=0.22, top=0.93)
    fig.savefig(
        figfile,
        bbox_inches="tight",
    )
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
