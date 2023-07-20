from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import typer
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value
UNWANTED_ZONES = ["zone-S"]
UNWANTED_TYPES = ["Unidentified"]

REPLACEMENTS = {
    "Deep": "Deep Neutral",
    "Fe": "Fe-Ni-Ca-Si",
}
BEST_TYPES = {
    "zone-0": ["Deep", "Neutral", "Low", "Medium"],
    "zone-I": ["Neutral", "Low", "Medium"],
    "zone-II": ["Low", "Medium"],
    "zone-III": ["Medium"],
    "zone-IV": ["Medium", "High"],
    "zone-MYSO": ["Deep", "Neutral", "Low", "Medium", "Fe"],
    "zone-S": ["Medium"],
}

def main(
        id_label: str,
        debug: bool=False,
        minimum_signal_noise: float=3.0,
        species_file: str="species.yaml",
        zone_file: str="zones.yaml",
        vsys: float=171.1,
        d_vsys: float=5.0,
):
    """Table of all identified lines"""

    # Load the zone database
    with open(zone_file) as f:
        # But drop zones we do not want
        zones = [_ for _ in yaml.safe_load(f) if _["label"] not in UNWANTED_ZONES]

    # Iterate over the zones
    for jzone, zone in enumerate(zones):
        # Read in the velocity table
        df = pd.read_csv(f"all-lines-{id_label}/{zone['label']}-velocities.csv").set_index("Index")
        if zone == zones[0]:
            # Initialize the output table
            df0 = df[["Type", "ID", "blend", "wave", "e_wave"]]
        eflux = df.flux / df.s_n
        flabel = f"F({zone['label'].split('-')[1]})"
        elabel = f"E({zone['label'].split('-')[1]})"
        df0 = df0.assign(**{flabel: df.flux, elabel: eflux})
        # Ignore values with s/n that is too low
        mask = df.s_n < minimum_signal_noise
        df0.loc[mask, flabel] = np.nan
        df0.loc[mask, elabel] = np.nan
    print(df0)
    df0.to_csv(f"all-lines-{id_label}/known-lines-final-table.csv")


if __name__ == "__main__":
    typer.run(main)
