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

    # Next, get all the lines into a big list of dicts
    line_files = sorted(Path(data_dir).glob(f"{INDEX_PATTERN}-*.yaml"))
    data = [
        yaml.safe_load(path.open()) for path in line_files
    ]

    # # Now make some dataframes, trying different techniques

    # # This one gives a flat dataframe with columns for all the nested
    # # dicts, but with long ugly column names like 'Zones Strength.zone-III'
    # df0 = pd.json_normalize(data).set_index("Index")
    # # Switch the long names to tuples. We cannot use a MultiIndex
    # # since that would require all columns to be tuples, which they
    # # are not
    # df0 = df0.rename(
    #     columns={
    #         s: tuple(x.replace("Zones ", "") for x in s.split(".")[::-1])
    #         for s in df0.columns if "." in s
    #     }
    # )
    # # Actually, I end up not using this. Leave here now for future reference


    # Better to select just the columns I need at the moment of creation of the DataFrame
    unwanted_types = ["sky",  "telluric", "noise", "nan"]
    # Make a seperate frame for Strength, BG, and BG sig
    df, df_bg, df_sig = [
        pd.DataFrame(
            {**{k: row[k] for k in ("Index", "Type", "ID")}, **row[f"Zones {label}"]}
            for row in data if row["Type"] and row["Type"].rstrip("?").lower() not in unwanted_types
        ).set_index("Index")
        for label in ["Strength", "BG", "BG sig"]
    ]

    # Now correct for the BG:
    print(df)
    print('Subtracting BG')
    zone_labels = [z["label"] for z in zones]
    df.loc[:, zone_labels] -=  df_bg.loc[:, zone_labels]

    # And put on scale of H beta = 100 (H beta has Index = 211)
    # We have to do the uncrtainties first, otherwise Hb has already changed!
    df_sig.loc[:, zone_labels] *= 100 / df.loc[211, zone_labels]
    df.loc[:, zone_labels] *= 100 / df.loc[211, zone_labels]

    print(df)
    df.to_csv(Path(data_dir) / "line-fluxes.csv")
    df_sig.to_csv(Path(data_dir) / "line-uncertainties.csv")

if __name__ == "__main__":
    typer.run(main)
