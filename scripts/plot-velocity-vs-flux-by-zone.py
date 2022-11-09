from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys


# UNWANTED_ZONES = ["zone-S"]
UNWANTED_ZONES = []
UNWANTED_TYPES = ["Unidentified"]

REPLACEMENTS = {
    "Deep": "Deep Neutral",
    "Fe": "Fe-Ni-Ca-Si",
}
BEST_TYPES = {
    "zone-0": ["Deep", "Neutral", "Low", "Medium"],
    "zone-I": ["Neutral", "Low", "Medium"],
    "zone-II": ["Low", "Medium"],
    "zone-III": ["Medium"],
    "zone-IV": ["Medium", "High"],
    "zone-MYSO": ["Deep", "Neutral", "Low", "Medium", "Fe"],
    "zone-S": ["Medium"],
}
def main(
        id_label: str,
        minimum_signal_noise: float=3.0,
        signal_noise_log_range: tuple[float, float]=(0.0, 2.0),
        species_file: str="species.yaml",
        zone_file: str="zones.yaml",
        v_sys: float=171.1,
        d_v_sys: float=2.7,
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
        figsize=(8, 7),
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
            # Remap to range [0, 1]
            # Switch to using s/n instead of flux
            norm_value = (
                (np.log10(data.s_n) - signal_noise_log_range[0])
                / (signal_noise_log_range[1] - signal_noise_log_range[0])
            )
            # log_f = np.log10(data.flux)
            # norm_value = (
            #     (log_f - type_data["log_min"])
            #     / (type_data["log_max"] - type_data["log_min"])
            # )

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
        xticks=[0.01, 1, 100],
        xticklabels=["0.01", "1", "100"],
        ylim=[80.0, 240.0],
    )
    axes[-1, 2].set(
        xlabel=" " * 24 + r"Line Intensity (H$\beta$ = 100)",
    )
    axes[3, 0].set(
        ylabel="Heliocentric Velocity, km/s",
    )
    # Add a line for the systemic velocity 
    for ax in axes.flat:
        ax.axhspan(
            v_sys - d_v_sys, v_sys + d_v_sys,
            alpha=0.2, zorder=0, color="k", linewidth=0,
        )
    # And a label for it
    axes[4, 1].text(
        0.8e-3, v_sys - 5,
        # fr"$V_\mathrm{{sys}} = {v_sys:.1f} \pm {d_v_sys:.1f}$ km/s",
        fr"$V_{{\odot}} = {v_sys:.0f}$ km/s",
        ha="right", va="center", fontsize="small",
        bbox={"facecolor": "white", "edgecolor": "none"},
    )

    # Add labels for the line types at the top of each column
    for itype, line_type in  enumerate(line_types):
        ax = axes[0, itype]
        type_data = info["types"][line_type]
        label = REPLACEMENTS.get(line_type, line_type)
        # Use a darker version of the color for the labels
        color = sns.dark_palette(
            tuple(type_data["husl"]),
            input="husl",
            as_cmap=True,
        )(0.5)
        ax.text(
            0.5, 1.2, label, transform=ax.transAxes,
            ha="center", va="center",
            color=color,
        )
    # And a title for these
    fig.text(0.5, 0.99, "Emission line type", ha="center", va="top")

    # Label the zones down the left and right sides
    for jzone, zone in enumerate(zones):
        label = zone["label"].replace("zone-", "Zone\n")
        ax = axes[jzone, -1]
        # Use a darker version of the color for the labels
        color = sns.dark_palette(
            tuple(zone["husl"]),
            input="husl",
            as_cmap=True,
        )(0.8)
        ax.text(
            1.3, 0.5, label, transform=ax.transAxes,
            ha="center", va="center",
            color=color,
        )
        ax = axes[jzone, 0]
        ax.text(
            -1.0, 0.5, label, transform=ax.transAxes,
            ha="center", va="center",
            color=color,
        )
        # Add an asterisk in the panels of line types that are
        # brightest for this zone
        for line_type in BEST_TYPES[zone["label"]]:
            itype = line_type_index[line_type]
            ax = axes[jzone, itype]
            ax.text(
                0.95, 0.05, "[*]", transform=ax.transAxes,
                ha="right", va="bottom",
                color=color, fontsize="medium", fontweight="bold",
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
        title="Species",
        fontsize="small",
        loc="lower center",
        ncol=8,
    )
    sns.despine()
    #fig.tight_layout()
    fig.subplots_adjust(left=0.15, bottom=0.22, top=0.93)
    fig.savefig(
        figfile,
        # bbox_inches="tight",
    )
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
