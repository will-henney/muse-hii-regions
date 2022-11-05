from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys

def plot_one_ratio_ratio_point(
        ax: matplotlib.axes.Axes,
        data: pd.Series,
        species_data: dict,
        type_data: dict,
        cmap: matplotlib.colors.Colormap,
        x_zone_ratio: str="II / 0",
        y_zone_ratio: str="I / 0",
):
    """Plot the zone ratio vs zone ratio for a single emission line

      Deal with the error bars and possible upper limits in any of the
      4 zones.
    """
    assert " / " in x_zone_ratio and " / " in y_zone_ratio
    assert type(data) == pd.Series, f"{type(data)}: {data}"

    # Fluxes in each zone, assuming ratio is B / A
    xB, xA  = [data[f"zone-{zone}"] for zone in x_zone_ratio.split(" / ")]
    yB, yA  = [data[f"zone-{zone}"] for zone in y_zone_ratio.split(" / ")]

    # Corresponding uncertainties
    xBe, xAe  = [data[f"zone-{zone}_sig"] for zone in x_zone_ratio.split(" / ")]
    yBe, yAe  = [data[f"zone-{zone}_sig"] for zone in y_zone_ratio.split(" / ")]

    size = species_data["size"]
    marker = species_data["marker"]
    # Log brightness in this line's most favorable zone
    log_f = np.log10(data[type_data["zone"]])
    # Remap to range [0, 1]
    norm_value = (log_f - type_data["log_min"]) / (type_data["log_max"] - type_data["log_min"])

    # Fluxes that are negative or smaller than the uncertainty are
    # replaced with 2-sigma upper limits
    x_upper = x_lower = y_lower = y_upper = False
    if yB < yBe:
        y_upper = True
        yB = max(yB, 0.0) + 2 * yBe
    if yA < yAe:
        y_lower = True
        yA = max(yA, 0.0) + 2 * yAe
    if y_upper and y_lower:
        # If both A and B fluxes are upper limits, then the ratio is indeterminate
        return

    if xB < xBe:
        x_upper = True
        xB = max(xB, 0.0) + 2 * xBe
    if xA < xAe:
        x_lower = True
        xA = max(xA, 0.0) + 2 * xAe
    if x_upper and x_lower:
        # If both A and B fluxes are upper limits, then the ratio is indeterminate
        return

    try: 
        xratio = xB / xA
        yratio = yB / yA
    except ZeroDivisionError:
        return

    # For the upper/lower limits, use a fixed fraction of the ratio to
    # give a constant length on a log scale. Also, make sure
    # upper/lower limits are penalized in brightness
    if x_upper or x_lower:
        xerr = 0.3 * xratio
        if species_data["name"] not in ["[C I]"]:
            norm_value /= 2
    else:
        xerr = np.hypot(xBe / xB, xAe / xA) * xratio
    if y_upper or y_lower:
        yerr = 0.3 * yratio
        if species_data["name"] not in ["[C I]"]:
            norm_value /= 2
    else:
        yerr = np.hypot(yBe / yB, yAe / yA) * yratio

    # Get the color for this point from the color map and normalized brightness
    color = cmap(norm_value)
    # Edge color is white or black, whichever gives best contrast with the key color
    if type_data["husl"][2] <= 50:
        edgecolor = "w"
    else:
        # In the black case, make it fade along with the face color
        edgecolor = matplotlib.colormaps["gray_r"](norm_value)

    ax.scatter(
        xratio, yratio,
        alpha=1.0, s=size, color=color, edgecolors=edgecolor,
        linewidths=0.1, marker=marker,
    )

    ax.errorbar(
        xratio,
        yratio,
        xerr=xerr,
        yerr=yerr,
        lolims=y_lower,
        uplims=y_upper,
        xlolims=x_lower,
        xuplims=x_upper,
        ecolor=color,
        fmt="none",
        elinewidth=0.3,
        alpha=1.0,
        zorder=0,
        errorevery=1,
        capsize=2.0,
        capthick=0.0,
    )


def main(
        id_label: str,
        species_file: str="species.yaml",
):
    """Plot of line fluxes from zones 0 and I"""

    with open(species_file) as f:
        info = yaml.safe_load(f)

    fig, axes = plt.subplots(2,  2, figsize=(8, 7))
    figfile = f"ratios-vs-ratios-by-zone-{id_label}.pdf"

    df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv").set_index("Index")
    df_sig = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv").set_index("Index")
    df = df.join(df_sig, rsuffix="_sig")

    xmin, xmax = 0.01, 100.0

    for ax, [x_zone_ratio, y_zone_ratio], exclude in zip(
            axes.flat,
            [
                ["II / 0", "I / 0"],
                ["II / 0", "MYSO / 0"],
                ["II / 0", "III / II"],
                ["III / II", "IV / III"],
            ],
            [
                ["IV"],
                ["IV"],
                [],
                ["0", "I"],                
            ],
    ):
        for species in reversed(info["species"]):
            label = species["name"]
            type_data = info["types"][species["type"]]
            if type_data["zone"].split("-")[-1] in exclude:
                continue
            mask = df.ID.str.startswith(label) & (df[type_data["zone"]] > 0.0)
            data = df[mask].sort_values(by=type_data["zone"])
            cmap = sns.light_palette(
                tuple(type_data["husl"]),
                input="husl",
                as_cmap=True,
            )
            for idx, row in data.iterrows():
                plot_one_ratio_ratio_point(
                    ax, row, species, type_data, cmap,
                    x_zone_ratio, y_zone_ratio
                )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(1.0,  color="k", linestyle="dashed", linewidth=0.5, zorder=-100)
        ax.axvline(1.0,  color="k", linestyle="dashed", linewidth=0.5, zorder=-100)
        ax.set(
            xlabel=x_zone_ratio,
            ylabel=y_zone_ratio,
            xlim=[xmin, xmax],
            ylim=[xmin, xmax],
        )
        ax.set_aspect("equal")
    # Now do one more loop over the species to mak the legend
    handles = []
    for species in info["species"]:
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

    fig.legend(handles=handles, title="Line Type", fontsize="small", loc="center right")
    #fig.tight_layout()
    #fig.subplots_adjust(right=0.8)
    fig.savefig(figfile, bbox_inches="tight")
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
