from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys


UNWANTED_ZONES = ["zone-S"]
UNWANTED_TYPES = ["Unidentified"]

def main(
        id_label: str,
        minimum_signal_noise: float=3.0,
        species_file: str="species.yaml",
        zone_file: str="zones.yaml",
):
    """Plot of line velocity versus flux for all zones"""

    # Load the species and line type database
    with open(species_file) as f:
        info = yaml.safe_load(f)
    # Load the zone database
    with open(zone_file) as f:
        # But drop zones we do not want
        zones = [_ for _ in yaml.safe_load(f) if _["label"] not in UNWANTED_ZONES]

    # The zones will correspond to the rows
    n_rows = len(zones)
    # And the line types to the columns (but miss out unwanted ones)
    line_types = [_ for _ in info["types"] if _ not in UNWANTED_TYPES]
    n_cols = len(line_types)
    # Reverse mapping from line type to index in the list
    line_type_index = {line_types[i]: i for i in range(n_cols)}

    fig, axes = plt.subplots(
        n_rows, n_cols,
        sharex=True, sharey=True,
        figsize=(8, 6),
    )


    figfile = f"velocity-vs-flux-by-zone-{id_label}.pdf"

    # Iterate over the zones
    for jzone, zone in enumerate(zones):
        # Read in the velocity table
        df = pd.read_csv(f"all-lines-{id_label}/{zone['label']}-velocities.csv")

        # So long as so not want to fade out the colors for the lower
        # flux lines, we can use a single call to scatter and errorbar for each species

        # Iterate over the species
        for species in reversed(info["species"]):
            if species["type"] in UNWANTED_TYPES:
                continue
            label = species["name"]
            itype = line_type_index[species["type"]]

            # Which one of the axes do we plot into?
            ax = axes[jzone, itype]

            type_data = info["types"][species["type"]]
            mask = df.ID.str.startswith(label) & (df.s_n > minimum_signal_noise)
            data = df[mask].sort_values(by="s_n")
            cmap = sns.light_palette(
                tuple(type_data["husl"]),
                input="husl",
                as_cmap=True,
            )
            size = species["size"]
            marker = species["marker"]
            # Log brightness in this line's most favorable zone
            log_f = np.log10(data.flux)
            # Remap to range [0, 1]
            norm_value = (
                (log_f - type_data["log_min"])
                / (type_data["log_max"] - type_data["log_min"])
            )
            if type_data["husl"][2] <= 50:
                edgecolors = "w"
            else:
                edgecolors =  matplotlib.colormaps["gray_r"](norm_value)

            ax.scatter(
                "flux", "vel",
                data=data,
                alpha=1.0, s=size,
                c=norm_value,
                vmin=0.0, vmax=1.0,
                cmap=cmap,
                edgecolors=edgecolors,
                linewidths=0.1, marker=marker,
            )

    ax = axes[-1, 0]
    ax.set_xscale("log")
    ax.set(
        xlim=[2e-3,  1e3],
        ylim=[80.0, 240.0],
    )
    axes[-1, 2].set(
        xlabel=r"Line Intensity (H$\beta$ = 100)",
    )
    axes[2, 0].set(
        ylabel="Heliocentric Velocity, km/s",
    )
    for ax in axes.flat:
        ax.axhspan(
            171.1 - 2.7, 171.1 + 2.7,
            alpha=0.2, zorder=0, color="k", linewidth=0,
        )
    for itype, line_type in  enumerate(line_types):
        axes[0, itype].set_title(line_type)
    for jzone, zone in enumerate(zones):
        label = zone["label"].replace("zone-", "Zone\n")
        ax = axes[jzone, -1]
        ax.text(
            1.3, 0.5, label, transform=ax.transAxes,
            ha="center", va="center",
        )
        ax = axes[jzone, 0]
        ax.text(
            -1.0, 0.5, label, transform=ax.transAxes,
            ha="center", va="center",
        )

    # Now do one more loop over the species to make the legend
    handles = []
    for species in info["species"]:
        if species["type"] in UNWANTED_TYPES:
            continue
        label = species["name"]
        size = species["size"]
        marker = species["marker"]
        type_data = info["types"][species["type"]]
        cmap = sns.light_palette(
            tuple(type_data["husl"]),
            input="husl",
            as_cmap=True,
        )
        # Get darkest color from the color map
        color = cmap(1.0)
        # Edge color is white or black, whichever gives best contrast with the key color
        if type_data["husl"][2] <= 50:
            edgecolor = "w"
        else:
            edgecolor = "k"
        handles.append(
            matplotlib.lines.Line2D(
                [], [],
                linestyle="none",
                label=label,
                marker=marker,
                markersize=np.sqrt(size),
                color=color,
                markeredgecolor=edgecolor,
                markeredgewidth=0.1,
            )
        )

    fig.legend(
        handles=handles,
        title="Line Type",
        fontsize="small",
        loc="lower center",
        ncol=8,
    )
    sns.despine()
    #fig.tight_layout()
    fig.subplots_adjust(left=0.15, bottom=0.25, top=0.95)
    fig.savefig(
        figfile,
        # bbox_inches="tight",
    )
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
