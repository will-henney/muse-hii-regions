from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import yaml
import typer

def main(
        id_label: str,
):
    """Plot of line fluxes from zones 0 and I"""
    fig, axes = plt.subplots(2,  2, figsize=(10, 10))
    figfile = f"fluxes-by-zone-{id_label}.pdf"

    df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv")
    df_sig = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv")

    # m_deep = df.Type.str.startswith("Deep")
    m_deep = df.ID.str.startswith("UIL")
    m_neut = df.Type.str.startswith("Neutral")
    m_low = df.Type.str.startswith("Low")
    m_med = df.Type.str.startswith("Med Neb")
    m_high = df.Type.str.startswith("High")

    xmin, xmax = 5e-4, 2000.0
    for ax, [zone_A,  zone_B] in zip(
            axes.flat,
            [
                ["0", "I"],
                ["0", "II"],
                ["0", "MYSO"],
                ["III", "IV"],
            ]
    ):

        ax.scatter(f"zone-{zone_A}", f"zone-{zone_B}", alpha=0.4, data=df[m_deep], label="Deep")
        ax.scatter(f"zone-{zone_A}", f"zone-{zone_B}", alpha=0.4, data=df[m_neut], label="Neutral")
        ax.scatter(f"zone-{zone_A}", f"zone-{zone_B}", alpha=0.4, data=df[m_low], label="Low Neb")
        ax.scatter(f"zone-{zone_A}", f"zone-{zone_B}", alpha=0.4, data=df[m_med], label="Medium Neb")
        ax.scatter(f"zone-{zone_A}", f"zone-{zone_B}", alpha=0.4, data=df[m_high], label="High Neb")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(0.0,  color="k")
        ax.axvline(0.0,  color="k")
        ax.plot([xmin, xmax], [xmin, xmax], color="k")
        ax.set(
            xlabel=f"Zone {zone_A} Intensity",
            ylabel=f"Zone {zone_B} Intensity",
            xlim=[xmin, xmax],
            ylim=[xmin, xmax],
        )
        ax.set_aspect("equal")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
