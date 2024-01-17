import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys
from matplotlib import pyplot as plt
from astropy.stats import SigmaClip
import matplotlib
import seaborn as sns

sns.set_context("talk")

def main(
            data_path: Path=Path.cwd(),
            clip_sigma: float=2.0,
):
      df = pd.read_csv(data_path / "h2-reddening-pairs-table.csv")

      # Clip the outliers
      clip = SigmaClip(sigma=clip_sigma, maxiters=5)
      clipped = clip(df.obs_theo_ratio, masked=True)
      df = df.loc[~clipped.mask]

      fig, ax = plt.subplots()
      x = df.wav_ratio
      y = np.log(df.obs_theo_ratio)
      ey = np.log(1.0 + df.rel_error)
      size = 100 * df.bright_product
      ax.scatter(x, y, s=size)
      # errorbar does not accept array of lienwidths, so do one at a time
      for _x, _y, _ey in zip(x, y, ey):
            ax.errorbar(_x, _y, yerr=_ey, fmt="none", elinewidth=0.1/_ey)

      # Plot the theoretical curve
      xmin, xmax = 0.94, 1.39
      xpts = np.linspace(xmin, xmax)
      # Mean wavelength of long line in each pair, micron
      wav_2 = 0.87
      # Central wavelength of V and K bands 
      wav_V = 0.55
      wav_K = 2.0

      for Av, beta, linestyle in [[1.0, 1.4, "solid"], [1.0, 2.0, "dashed"], [0.3, 1.4, "dashdot"]]:
            tau_2 = Av / (1.0857 * (wav_2 / wav_V) ** beta)
            Ak = Av * (wav_V / wav_K) ** beta
            label = fr"$A_V = {Av:.1f}$, $\beta = {beta:.1f}$ ($A_K = {Ak:.2f}$)"
            ax.plot(xpts, tau_2 * (1.0 - xpts ** beta), label=label, linestyle=linestyle)
      ax.axhline(0.0, linewidth=0.5, color="k", linestyle="dotted")
      ax.axvline(1.0, linewidth=0.5, color="k", linestyle="dotted")
      ax.legend(loc="lower left", fontsize="x-small")
      ax.set(
            xlim=[xmin, xmax],
            ylim=[-0.49, 0.49],
            xlabel="Wavelength ratio: $\lambda_2$ / $\lambda_1$",
            ylabel=r"$\ln \left[\dfrac{(\mathrm{Observed\ /\ Predicted})_1}{(\mathrm{Observed\ /\ Predicted})_2}\right]$",
      )

      sns.despine()
      figfile = "h2-reddening.pdf"
      fig.savefig(figfile, bbox_inches="tight")
      print(figfile, end="")

if __name__ == "__main__":
    typer.run(main)
