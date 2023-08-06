import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys
COLS = ["H2_line", "wl_lab", "wl_obs", "Notes"]
def main(
            data_dir: str="n346-h2-lines",
            orig_data_dir: str="n346-lines/all-lines-c007-chop-mean",
):
      df1 = pd.read_csv(Path(orig_data_dir) / "uil-final-table.csv")
      df2 = pd.read_csv(Path(data_dir) / "h2-line-ids.csv")
      df2 = df2[COLS].sort_values("wl_obs").dropna(subset=["wl_obs"])
      df0 = pd.merge_asof(df1, df2, left_on="wave0", right_on="wl_obs", direction="nearest", tolerance=1.0)
      df0.to_csv(Path(data_dir) / "h2-final-table.csv")
      print(df0)

if __name__ == "__main__":
    typer.run(main)
