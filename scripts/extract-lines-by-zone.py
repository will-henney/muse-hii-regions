from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer

INDEX_PATTERN = "[0-9]" * 4

def main(
        zone_file: str="zones.yaml",
        data_dir: str="all-lines-c007-chop-mean",
):
    # First, get the zones
    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    # Next, get all the lines
    line_files = sorted(Path(data_dir).glob(f"{INDEX_PATTERN}-*.yaml"))
    data = [
        yaml.safe_load(path.open()) for path in line_files
    ]
    print(pd.DataFrame(
        {**{k: row[k] for k in ("Index", "ID")}, **row["Zones Strength"]}
        for row in data if row["Type"] and "Deep" in row["Type"]
    ))


if __name__ == "__main__":
    typer.run(main)
