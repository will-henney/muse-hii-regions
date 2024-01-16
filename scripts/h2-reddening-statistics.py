import numpy as np
import yaml
from pathlib import Path
import pandas as pd
from statsmodels.stats.weightstats import DescrStatsW
from astropy.stats import SigmaClip
import typer
import sys
from matplotlib import pyplot as plt
import matplotlib

import seaborn as sns


def main(
            data_path: Path=Path.cwd(),
            min_wav_ratio: float=1.15,
            good_bright_product: float=1.5,
            clip_sigma: float=2.0,
):
      df = pd.read_csv(data_path / "h2-reddening-pairs-table.csv");
      # Weight is UNITY for pairs where bright_product is large
      # enough. These correspond to a "good" observation. Weight may
      # not exceed 1.0, but may be lower for less trustworthy pairs
      df["weights"] = np.minimum(1.0, df.bright_product / good_bright_product)
      df0 = df[df.wav_ratio < min_wav_ratio]
      df1 = df[df.wav_ratio >= min_wav_ratio]
      cols = ["wav_ratio", "obs_theo_ratio", "long_wav"]
      for df, descrip in [[df0, "CLOSE PAIRS"], [df1, "FAR PAIRS"]]:
            print("\n", descrip, "\n")
            print("Weighted statistics", "N = ", len(df))
            wstats = DescrStatsW(df[cols].to_numpy(), weights=df.weights)
            print("Column:\t", (cols))
            print("Mean:\t", np.round(wstats.mean, 2))
            print("Std:\t", np.round(wstats.std, 2))
            print("SEM:\t", np.round(wstats.std_mean, 2))
            print("N effective:\t", np.round(wstats.nobs, 2))
            print("CI 95:\n", np.round(wstats.tconfint_mean(alpha=0.05), 2))
            print("CI 68:\n", np.round(wstats.tconfint_mean(alpha=0.32), 2))
            print()

            print("Sigma-clipped statistics\n")
            clip = SigmaClip(sigma=clip_sigma, maxiters=5)
            clipped = clip(df.obs_theo_ratio, masked=True)
            dfc = df.loc[~clipped.mask]
            print(dfc.describe())
            print()

            print("Weighted AND sigma-clipped statistics", "N = ", len(dfc))
            wstats = DescrStatsW(dfc[cols].to_numpy(), weights=dfc.weights)
            print("Column:\t", (cols))
            print("Mean:\t", np.round(wstats.mean, 2))
            print("Std:\t", np.round(wstats.std, 2))
            print("SEM:\t", np.round(wstats.std_mean, 2))
            print("N effective:\t", np.round(wstats.nobs, 2))
            print("CI 95:\n", np.round(wstats.tconfint_mean(alpha=0.05), 2))
            print("CI 68:\n", np.round(wstats.tconfint_mean(alpha=0.32), 2))
            print()

if __name__ == "__main__":
    typer.run(main)
