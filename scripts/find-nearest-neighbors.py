from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value


def main(
        id_label: str,
        species_file: str="species.yaml",
):
    """Get tables of wavelengths by type of line"""

    with open(species_file) as f:
        info = yaml.safe_load(f)

    df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv").set_index("Index")
    df_wav = pd.read_csv(f"all-lines-{id_label}/line-gauss-waves.csv").set_index("Index")
    df = df.join(df_wav, rsuffix="_wav")

    # First, make separate tables for each line type
    # And another for all the identified lines
    table_types = list(info["types"].keys()) + ["Identified"]
    # Initialize a dict to hold the tables
    tables = {k: [] for k in table_types}
    # Go accumulating rows in the tables by iterating over the species
    for species in reversed(info["species"]):
        _type = species["type"]
        type_data = info["types"][_type]
        zone = type_data["zone"]
        prefix = species["name"]
        if not species["name"] == "UIL":
            prefix += " "
        mask = df.ID.str.startswith(prefix) & (df[zone] > 0.0)
        data = df[mask]
        # Stripped down to only the data we need
        newdata = data[["Type", "ID"]].assign(
            flux=data[zone],
            wave=data[f"{zone}_wav"],
        )
        tables[_type].append(newdata)
        if not _type == "Unidentified":
            tables["Identified"].append(newdata)
    for _type, tab in tables.items():
        # Reassemble each table into a dataframe
        tab = pd.concat(tab)
        # Calculate wave numbers and sort in ascending wn
        tab = tab.assign(
            wn=1.0 / (tab.wave * (u.Angstrom.to(u.cm))),
        ).sort_values(by="wn")
        # Calculate backward and forward differences
        tab = tab.assign(
            dwn_b=np.ediff1d(tab.wn, to_begin=np.inf),
            dwn_f=np.ediff1d(tab.wn, to_end=np.inf),
        )
        # For nearest neighbor, choose minimum of the two and drop originals
        tab = tab.assign(
            dwn=np.minimum(tab.dwn_b, tab.dwn_f),
        ).drop(columns=["dwn_b", "dwn_f"])
        # Take second order differences squared to check reciprocity
        tab = tab.assign(
            dd_b=np.ediff1d(tab.dwn, to_begin=10)**2,
            dd_f=np.ediff1d(tab.dwn, to_end=10)**2,
        )
        # Reciprocal if at least one of these is zero
        tab = tab.assign(
            mutual=np.minimum(tab.dd_b, tab.dd_f) == 0.0,
        ).drop(columns=["dd_b", "dd_f"])
        print()
        print(_type, "Degree of reciprocity =", tab.mutual.sum() / len(tab))
        print(tab.describe())
        tab.to_csv(f"all-lines-{id_label}/nearest-neighbors-{_type}.csv")





if __name__ == "__main__":
    typer.run(main)
