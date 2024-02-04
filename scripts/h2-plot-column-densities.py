import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns

sns.set_context("talk")

def main(
            data_dir: str=".",
            ortho_para_ratio: float=3.0,
):
      df = pd.read_csv(Path(data_dir) / "h2-excitation-table.csv").set_index("Index")
      ortho = (df.Jhi % 2) == 1
      x = df.Excit_hi_K
      y = np.log(
            (df["F(0)"] / df.g_u_h_nu_Aul)
            / (df.loc[3548, "F(0)"] / df.loc[3548, "g_u_h_nu_Aul"])
      )
      c = df.Vhi
      fig, ax = plt.subplots(figsize=(12, 4))
      ax.scatter(
            x[ortho], y[ortho], c=c[ortho],
            s=10, cmap="tab20", vmin=2, vmax=13,
      )
      ax.scatter(
            x[~ortho], y[~ortho] + np.log(ortho_para_ratio/3), c=c[~ortho],
            s=10, marker="v", cmap="tab20", vmin=2, vmax=13,
      )
      ax.set(
            xlabel="Excitation energy, $E/k$, K",
            ylabel="Column density, $\ln(N/g)$",
            xlim=[0, 55000],
            ylim=[-8, 5],
      )
      figfile = "h2-excitation.pdf"
      fig.savefig(figfile, bbox_inches="tight")
      print(figfile, end="")




if __name__ == "__main__":
    typer.run(main)
