import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
def main(
            data_path: Path=Path.cwd(),
):
      df = pd.read_csv(data_path / "h2-reddening-pairs-table.csv")
      fig, ax = plt.subplots()
      x = df.wav_ratio - 1.0
      y = np.log(df.obs_theo_ratio)
      ey = np.log(1.0 + df.rel_error)
      size = 100 * df.bright_product
      ax.scatter(x, y, s=size)
      # errorbar does not accept array of lienwidths, so do one at a time
      for _x, _y, _ey in zip(x, y, ey):
            ax.errorbar(_x, _y, yerr=_ey, fmt="none", elinewidth=0.1/_ey)
      ax.set(
            xlim=[-0.04, 0.39],
            xlabel="($\lambda_2$ / $\lambda_1$) - 1",
            ylabel="log [(Observed/Predicted)$_1$ / (Observed/Predicted)$_2$]",
      )

      sns.despine()
      figfile = "h2-reddening.pdf"
      fig.savefig(figfile, bbox_inches="tight")
      print(figfile, end="")

if __name__ == "__main__":
    typer.run(main)
