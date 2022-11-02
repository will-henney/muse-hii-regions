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
    fig, ax = plt.subplots(figsize=(8, 8))
    figfile = f"fluxes-0-I-{id_label}.pdf"

    df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv")
    df_sig = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv")

    m_deep = df["Type"].map(lambda x: x.lower().startswith("deep"))
    m_neut = df["Type"].map(lambda x: x.lower().startswith("neutral"))
    m_low = df["Type"].map(lambda x: x.lower().startswith("low"))

    ax.scatter("zone-0", "zone-MYSO", data=df[m_deep], label="Deep")
    ax.scatter("zone-0", "zone-MYSO", data=df[m_neut], label="Neutral")
    ax.scatter("zone-0", "zone-MYSO", data=df[m_low], label="Low Neb", color="r")
    ax.set_xscale("symlog", linthresh=0.01, linscale=1.0)
    ax.set_yscale("symlog", linthresh=0.01, linscale=1.0)

    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
