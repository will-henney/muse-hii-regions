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
    """Make tables of nearest neighbor distances"""

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
        # Reassemble each table into a dataframe and sort on wavelength
        tab = pd.concat(tab).sort_values(by="wave")
        # Calculate backward and forward wave differences, padding the
        # ends with the mean separation to fake extra lines outside the range
        if len(tab) > 1:
            fill_value = (tab.wave.max() - tab.wave.min()) / (len(tab) - 1)
        else:
            fill_value = 0.0
        tab = tab.assign(
            dwave_b=np.ediff1d(tab.wave, to_begin=fill_value),
            dwave_f=np.ediff1d(tab.wave, to_end=fill_value),
        )
        # For nearest neighbor, choose minimum of the two, whereas
        # mean spacing is average of the two, then drop originals
        tab = tab.assign(
            dwave_nn=np.minimum(tab.dwave_b, tab.dwave_f),
            dwave_mean=(tab.dwave_b + tab.dwave_f) / 2,
        ).drop(columns=["dwave_b", "dwave_f"])
        # Calculate wave numbers and sort in ascending wn
        tab = tab.assign(
            wn=1.0 / (tab.wave * (u.Angstrom.to(u.cm))),
        ).sort_values(by="wn")
        # Calculate backward and forward wn differences
        if len(tab) > 1:
            fill_value = (tab.wn.max() - tab.wn.min()) / (len(tab) - 1)
        else:
            fill_value = 0.0
        tab = tab.assign(
            dwn_b=np.ediff1d(tab.wn, to_begin=fill_value),
            dwn_f=np.ediff1d(tab.wn, to_end=fill_value),
        )
        # For nearest neighbor, choose minimum of the two and drop originals
        tab = tab.assign(
            dwn_nn=np.minimum(tab.dwn_b, tab.dwn_f),
            dwn_mean=(tab.dwn_b + tab.dwn_f) / 2,
        ).drop(columns=["dwn_b", "dwn_f"])
        # Take second order differences squared to check reciprocity
        tab = tab.assign(
            dd_b=np.ediff1d(tab.dwn_nn, to_begin=10)**2,
            dd_f=np.ediff1d(tab.dwn_nn, to_end=10)**2,
        )
        # Reciprocal if at least one of these is zero
        tab = tab.assign(
            mutual=np.minimum(tab.dd_b, tab.dd_f) == 0.0,
        ).drop(columns=["dd_b", "dd_f"])
        print()
        # Ignore first and last points when calculating reciprocity, since their status is unknown
        print(_type, "Degree of reciprocity =", tab[1:-1].mutual.sum() / len(tab[1:-1]))
        print(tab.describe())
        tab.to_csv(f"all-lines-{id_label}/nearest-neighbors-{_type}.csv")





if __name__ == "__main__":
    typer.run(main)
