from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import pyneb as pn
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
    "zone-0": ["Deep"],
    "zone-I": ["Neutral"],
    "zone-II": ["Low"],
    "zone-III": ["Med"],
    "zone-IV": ["High"],
    "zone-MYSO": ["Deep", "Fe"],
    "zone-S": ["Med"],
}
ACCEPTABLE_TYPES = {
    "zone-0": ["Deep", "Neutral", "Low", "Med", "Fe"],
    "zone-I": ["Deep", "Neutral", "Low", "Med", "Fe"],
    "zone-II": ["Deep", "Neutral", "Low", "Med", "Fe"],
    "zone-III": ["Low", "Med", "Fe", "High"],
    "zone-IV": ["Med", "High"],
    "zone-MYSO": ["Deep", "Neutral", "Low", "Med", "Fe"],
    "zone-S": ["Med"],
}
H2_COLS = ["H2_line", "wl_lab", "wl_obs", "Notes", "H2_index"]

# Find intrinsic Balmer decrement
hi = pn.RecAtom("H", 1)
tem, den = 12500, 100
R0 = hi.getEmissivity(tem, den, wave=6563) / hi.getEmissivity(tem, den, wave=4861)
# Set up reddening law for SMC
rc = pn.RedCorr()
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
rc.law = "F99"


def main(
    id_label: str,
    debug: bool = False,
    minimum_signal_noise: float = 3.0,
    zone_file: str = "zones.yaml",
    h2_data_dir: str = "n346-h2-lines",
    orig_data_dir: str = "n346-lines",
    vsys: float = 171.1,
    d_vsys: float = 5.0,
):
    """Table of all identified lines"""

    # Load the zone database
    with open(Path(orig_data_dir) / zone_file) as f:
        # But drop zones we do not want
        zones = [_ for _ in yaml.safe_load(f) if _["label"] not in UNWANTED_ZONES]

    # Iterate over the zones
    wstrings = []
    fstrings = []
    for jzone, zone in enumerate(zones):
        # Read in the velocity table
        df = pd.read_csv(
            Path(orig_data_dir) / f"all-lines-{id_label}/{zone['label']}-velocities.csv"
        ).set_index("Index")
        if zone == zones[0]:
            # Initialize the output table
            df0 = df[["Type", "ID", "blend"]]
            # And a parallel shadow table for using only the best zones for determing the waves
            df00 = df[["Type", "ID", "blend"]]
        # Columns of wavelength and flux for each zone, with their respective errors
        eflux = df.flux / df.s_n
        zstring = zone["label"].split("-")[1]
        flabel = f"F({zstring})"
        elabel = f"E({zstring})"
        wlabel = f"W({zstring})"
        dwlabel = f"dW({zstring})"
        # Save a list of the wavelength columns
        wstrings.append(wlabel)
        fstrings.append(flabel)
        # Calculate the reddening correction from H alpha, which is channel 1573
        obs_decrement = df.loc[1573].flux / 100.0
        rc.setCorr(obs_decrement / R0, wave1=6563, wave2=4861)
        correction = rc.getCorr(df.wave) / rc.getCorr(df.loc[211].wave)
        # Add the 2 flux columns to the output table
        df0 = df0.assign(
            **{
                flabel: correction * df.flux,
                elabel: correction * eflux,
            }
        )
        # And the fluxes and waves to thje shadow table
        df00 = df00.assign(
            **{
                wlabel: df.wave / (1 + vsys / LIGHT_SPEED_KMS),
                dwlabel: df.e_wave,
                flabel: correction * df.flux,
                elabel: correction * eflux,
            }
        )
        # Ignore values with s/n that is too low or that are not acceptable types for this zone
        mask = (df.s_n < minimum_signal_noise) | ~df.Type.str.startswith(
            tuple(ACCEPTABLE_TYPES[zone["label"]])
        )
        for label in [flabel, elabel]:
            df0.loc[mask, label] = np.nan
        # Repeat for the shadow table, but only allowing the BEST types
        mask00 = (df.s_n < minimum_signal_noise) | ~df.Type.str.startswith(
            tuple(BEST_TYPES[zone["label"]])
        )
        for label in [wlabel, dwlabel, flabel, elabel]:
            df00.loc[mask00, label] = np.nan
    # Second pass: consolidate to a single wavelength column
    dwstrings = ["d" + _ for _ in wstrings]
    # Weighted average over valid zones
    df0.insert(
        4,
        "wave",
        np.nanmean(df00[wstrings].to_numpy(), axis=1),
    )
    # And use the smallest for wave error
    df0.insert(5, "e_wave", np.nanmin(df00[dwstrings], axis=1))
    # And drop lines with no wavelength
    df0 = df0[np.isfinite(df0.wave)]

    #
    # Add in the H2 lines
    #

    df2 = pd.read_csv(Path(h2_data_dir) / "h2-line-ids.csv")
    df2 = df2[H2_COLS].sort_values("wl_obs").dropna(subset=["wl_obs"])
    df0 = pd.merge_asof(
        df0, df2, left_on="wave", right_on="wl_obs", direction="nearest", tolerance=2.0
    ).set_index(df0.index)
    # Select the lines formerly unidentified
    h2mask = df0.ID.str.startswith("UIL")

    # Construct the ID description for the H2 lines
    df0.loc[h2mask, "ID"] = (
        "H_2 " + df0.loc[h2mask, "H2_line"] + " " + df0.loc[h2mask, "wl_lab"]
    )
    # And remove the extraneous columns
    df0 = df0.drop(columns=["H2_line", "wl_lab"])

    # And the entries with empty ID fields
    df0 = df0.dropna(subset="ID")

    print(df0.loc[h2mask])
    df0.to_csv(Path(h2_data_dir) / "known-lines-with-h2-final-table.csv")


if __name__ == "__main__":
    typer.run(main)
