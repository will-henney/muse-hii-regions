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

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value

XY_BY_TYPE = {
    "Unidentified": (7000, 0.15),
    "Deep": (8730, 1.15),
    "Neutral": (8200, 1.35),
    "Low": (7300, 1.15),
    "Medium": (6563, 1.35),
    "High": (5000, 1.35),
    "Fe": (5500, 1.15),
}
REPLACEMENTS = {
    "Deep": "Deep Neutral",
    "Low": "Low Ionization",
    "Medium": "Medium Ionization",
    "High": "High Ionization",
    "Fe": "Fe-Ni-Si-Ca",
}
def main(
        id_label: str,
        species_file: str="species.yaml",
        vsys: float=171.1,
):
    """Plot of line type distribution over wavelength"""

    with open(species_file) as f:
        info = yaml.safe_load(f)

    fig, ax = plt.subplots(1,  1, figsize=(10, 3))
    figfile = f"line-type-wave-distro-{id_label}.pdf"

    df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv").set_index("Index")
    df_wave = pd.read_csv(f"all-lines-{id_label}/line-gauss-waves.csv").set_index("Index")
    df = df.join(df_wave,  rsuffix="_wave")

    n_type = {_t: 0 for _t in info["types"]}

    for species in reversed(info["species"]):
        label = species["name"]
        type_data = info["types"][species["type"]]
        mask = df.ID.str.startswith(label) & (df[type_data["zone"]] > 0.0)
        data = df[mask].sort_values(by=type_data["zone"])
        n_type[species["type"]] += len(data)
        cmap = sns.light_palette(
            tuple(type_data["husl"]),
            input="husl",
            as_cmap=True,
        )
        y0 = 0 if label == "UIL" else 1
        for idx, row in data.iterrows():
            # Log brightness in this line's most favorable zone
            log_f = np.log10(row[type_data["zone"]])
            # Corresponding wavelengtn
            wave = row[type_data["zone"] + "_wave"]
            wave /= 1.0 + vsys / LIGHT_SPEED_KMS
            # Remap to range [0, 1]
            # norm_value = (log_f - type_data["log_min"]) / (type_data["log_max"] - type_data["log_min"])
            norm_value = (log_f - (-2)) / (3 - (-2))
            # Get the color for this point from the color map and normalized brightness
            color = cmap(1.0)
            try:
                ax.plot([wave, wave], [y0, y0 + 0.9], lw=0.2 * (2.5 + log_f), color=color)
            except:
                ...
    for _type, type_data in info["types"].items():
        label = f"{REPLACEMENTS.get(_type, _type)}: $N = {n_type[_type]}$"
        x0, y0 = XY_BY_TYPE[_type]
        color = sns.light_palette(
            tuple(type_data["husl"]),
            input="husl",
            as_cmap=True,
        )(1.0)
        ax.text(
            x0, y0, label,
            color=color, ha="center", va="center",
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "boxstyle": "Round, pad=0.1"},
        )

    ax.set(
        xlabel="Wavelength, Angstrom",
        ylabel="",
        yticks=[],
    )
    sns.despine(left=True)
    fig.savefig(figfile, bbox_inches="tight")
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
