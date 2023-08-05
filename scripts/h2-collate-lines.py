import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys

def main(data_dir: str="n346-h2-lines"):
      # Next, get all the lines into a big list of dicts
      line_files = sorted(Path(data_dir).glob("*.yaml"))
      data = [
            yaml.safe_load(path.open()) for path in line_files
      ]
      df0 = pd.DataFrame(data)
      df0.to_csv(Path(data_dir) / "h2-line-ids.csv")

if __name__ == "__main__":
    typer.run(main)
