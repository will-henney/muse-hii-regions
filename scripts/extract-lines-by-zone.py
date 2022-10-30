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
        orig_data_dir: str="all-lines-orig",
):
    # First, get the zones
    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    # Next, get all the lines
    line_files = sorted(Path(orig_data_dir).glob(f"{INDEX_PATTERN}.yaml"))
    data = [
        yaml.safe_load(path.open()) for path in line_files
    ]
    print(pd.DataFrame(data))


if __name__ == "__main__":
    typer.run(main)
