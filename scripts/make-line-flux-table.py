from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys

INDEX_PATTERN = "[0-9]" * 4

def main(
        zone_file: str="zones.yaml",
        data_dir: str="all-lines-c007-chop-mean",
):
    # First, get the zones
    with open(zone_file) as f:
        zones = yaml.safe_load(f)
    zone_labels = [z["label"] for z in zones]

    # Next, get all the lines into a big list of dicts
    line_files = sorted(Path(data_dir).glob(f"{INDEX_PATTERN}-*.yaml"))
    data = [
        yaml.safe_load(path.open()) for path in line_files
    ]

    # # Now make some dataframes, trying different techniques

    # This one gives a flat dataframe with columns for all the nested
    # dicts, but with some long column names like 'zone-III.Strength'
    df0 = pd.json_normalize(data).set_index("Index")
    unwanted_types = ["sky",  "telluric", "noise", "nan"]
    # Filter out unwanted rows. Note use of fancy .str accessor to use
    # regular string methods elementwise on a Series. Also use of @
    # inside the query string to reference python variables
    df0 = df0.query('Type.str.rstrip("?").str.lower() not in @unwanted_types')
    df0 = df0.fillna(value={"Type": "None", "ID": "None"})
    print(df0)

    # Make a separate frame for Strength, renaming the columns to just the zone labels
    df = df0[
        ["Type", "ID"] + [f"{zlabel}.Strength" for zlabel in zone_labels]
    ].rename(columns=lambda c: c.split(".")[0])
    # And do the same for Sigma
    df_sig = df0[
        ["Type", "ID"] + [f"{zlabel}.Sigma" for zlabel in zone_labels]
    ].rename(columns=lambda c: c.split(".")[0])
    # And the same for Gaussian fit amplitude
    df_g = df0[
        ["Type", "ID"] + [f"{zlabel}.Gauss Fit.Amplitude" for zlabel in zone_labels]
    ].rename(columns=lambda c: c.split(".")[0])
    # And for the mean wavelength
    df_wave = df0[
        ["Type", "ID"] + [f"{zlabel}.Gauss Fit.Mean Wave" for zlabel in zone_labels]
    ].rename(columns=lambda c: c.split(".")[0])
    # And for the line width
    df_width = df0[
        ["Type", "ID"] + [f"{zlabel}.Gauss Fit.RMS Width" for zlabel in zone_labels]
    ].rename(columns=lambda c: c.split(".")[0])



    # And put on scale of H beta = 100 (H beta has Index = 211)
    # We have to do the uncertainties first, otherwise Hb has already changed!
    df_sig.loc[:, zone_labels] *= 100 / df.loc[211, zone_labels]
    df.loc[:, zone_labels] *= 100 / df.loc[211, zone_labels]
    # For the gaussian amplitudes, we use Hb from the same. Note that
    # we do not multiply by the widths, since they have much larger
    # uncertainties than the amplitudes.  It makes more sense to
    # assume that all lines have approximately the same width
    df_g.loc[:, zone_labels] *= 100 / df_g.loc[211, zone_labels]

    print(df_g)

    print(df[df0.sky_blend])

    df.to_csv(Path(data_dir) / "line-fluxes.csv")
    df_sig.to_csv(Path(data_dir) / "line-uncertainties.csv")
    df_g.to_csv(Path(data_dir) / "line-gauss-fluxes.csv")
    df_wave.to_csv(Path(data_dir) / "line-gauss-waves.csv")
    df_width.to_csv(Path(data_dir) / "line-gauss-widths.csv")

if __name__ == "__main__":
    typer.run(main)
