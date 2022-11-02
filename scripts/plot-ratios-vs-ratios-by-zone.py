from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer

def main(
        id_label: str,
):
    """Plot of line fluxes from zones 0 and I"""
    fig, axes = plt.subplots(2,  2, figsize=(6, 5))
    figfile = f"ratios-vs-ratios-by-zone-{id_label}.pdf"

    df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv")
    df_sig = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv")
    df.loc[:, "zone-I":"zone-II"] = np.where(
        df.loc[:, "zone-I":"zone-II"] > 0.0,
        df.loc[:, "zone-I":"zone-II"],
        2 * df_sig.loc[:, "zone-I":"zone-II"]
    )

    m_deep = df.Type.str.startswith("Deep")  & (df["zone-0"] > 0.0)
    m_neut = df.Type.str.startswith("Neutral") & (df["zone-I"] > 0.0)
    m_low = df.Type.str.startswith("Low") & (df["zone-II"] > 0.0)
    m_med = df.Type.str.startswith("Med Neb") & (df["zone-III"] > 0.0)
    m_high = df.Type.str.startswith("High") & (df["zone-IV"] > 0.0)
    m_fe = df.Type.str.startswith("Fe") & (df["zone-MYSO"] > 0.0)
    m_si = df.Type.str.startswith("Med Perm") & (df["zone-MYSO"] > 0.0)
    m_ci = df.ID.str.startswith("[C I]")

    xmin, xmax = 0.01, 100.0
    for ax, [[zone_A,  zone_B], [zone_AA,  zone_BB]] in zip(
            axes.flat,
            [
                # [["0", "I"]   , ["I", "II"]], 
                [["0", "II"]   , ["0", "I"]], 
                [["I", "II"]  , ["II", "III"]], 
                [["II", "III"], ["III", "IV"]], 
                [["0", "II"], ["0", "MYSO"]], 
            ]
    ):

        for data, label, size, marker in [
                [df[m_deep], "Unidentified", 25, "*"],
                [df[m_neut], "NO I", 12, "s"],
                [df[m_low], "[O I], [NOSCl II]", 25, "p"],
                [df[m_med], "HHe I, [SArClO III]", 10, "P"],
                [df[m_high], "[ClArK IV], He II", 15, "d"],
                [df[m_fe], "[FeNi II,III]", 8, "o"],
                [df[m_si], "Si II", 20, "X"],
                [df[m_ci], "[C I]", 100, "*"],
        ]:
            A = data[f"zone-{zone_A}"]
            B = data[f"zone-{zone_B}"]
            AA = data[f"zone-{zone_AA}"]
            BB = data[f"zone-{zone_BB}"]
            ax.scatter(B / A, BB / AA, alpha=1.0,
                       s=1.5 * size, label=label,
                       edgecolors="w", linewidths=0.1, marker=marker,
                       )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(1.0,  color="k", linestyle="dashed", zorder=-100)
        ax.axvline(1.0,  color="k", linestyle="dashed", zorder=-100)
        ax.set(
            xlabel=f"{zone_B} / {zone_A}",
            ylabel=f"{zone_BB} / {zone_AA}",
            xlim=[xmin, xmax],
            ylim=[xmin, xmax],
        )
        ax.set_aspect("equal")
    axes[1, 0].legend(ncol=3, fontsize="xx-small")
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
