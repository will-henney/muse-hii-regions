from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import yaml
import typer

def main(
        id_label: str,
        use_gauss: bool=False,
):
    """Plot of line fluxes from zones 0 and I"""
    fig, axes = plt.subplots(2,  2, figsize=(7, 6))
    figfile = f"ratios-by-zone-{id_label}.pdf"

    if use_gauss:
        df = pd.read_csv(f"all-lines-{id_label}/line-gauss-fluxes.csv")
        figfile = figfile.replace(".pdf", "-gauss.pdf")
    else:
        df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv")
    df_sig = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv")

    m_deep = df.Type.str.startswith("Deep") & (df["zone-0"] > 0.0)
    m_neut = df.Type.str.startswith("Neutral")
    m_low = df.Type.str.startswith("Low")
    m_med = df.Type.str.startswith("Med Neb")
    m_high = df.Type.str.startswith("High")
    m_fe = df.Type.str.startswith("Fe")
    m_si = df.Type.str.startswith("Med Perm")

    xmin, xmax = 1e-4, 2000.0
    for ax, [zone_A,  zone_B] in zip(
            axes.flat,
            [
                ["0", "I"],
                ["0", "II"],
                ["0", "MYSO"],
                ["III", "IV"],
            ]
    ):

        for data, label, size in [
                [df[m_deep], "Deep", 10],
                [df[m_neut], "Neutral", 15],
                [df[m_low], "Low Neb", 30],
                [df[m_med], "Medium Neb", 8],
                [df[m_high], "High Neb", 25],
                [df[m_fe], "Fe, Ni", 8],
                [df[m_si], "Si II", 20],
        ]:
            A = data[f"zone-{zone_A}"]
            B = data[f"zone-{zone_B}"]
            ax.scatter(A, B / A, alpha=1.0, s=size, label=label)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(1.0,  color="k", linestyle="dashed", zorder=-100)
        ax.set(
            xlabel=f"Zone {zone_A} Intensity (Hβ = 100)",
            ylabel=f"Intensity Ratio: Zone {zone_B} / Zone {zone_A}",
            xlim=[xmin, xmax],
            ylim=[2e-3, 5e2],
        )
        ax.set_aspect("equal")
    ax.legend(ncol=2, fontsize="xx-small")
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
