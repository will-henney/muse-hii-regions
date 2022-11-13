from pathlib import Path
import pandas as pd
import numpy as np
import typer
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value

def main(
        id_label: str,
        debug: bool=False,
        vsys: float=171.1,
        d_vsys: float=5.0,
):
    """Table of all unidentified lines"""

    # Load the data for zone 0 ... 
    df0 = pd.read_csv(f"all-lines-{id_label}/zone-0-velocities.csv").set_index("Index")
    # ... for the wavelengths and velocities ...
    dfd = pd.read_csv(f"all-lines-{id_label}/interzone-0-MYSO-vel-diffs.csv").set_index("Index")
    # ... and for the zone ratios
    dfr = pd.read_csv(f"all-lines-{id_label}/ratios-vs-ratios-by-zone.csv").set_index("Index")

    # Select only the unidentified lines with credible detection
    mask = df0.ID.str.startswith("UIL") & (df0.s_n > 1.0)
    # Make a new frame combining selected columns from each
    df = pd.concat(
        [
            df0[["Type", "ID", "blend", "wave", "e_wave", "flux", "s_n"]],
            dfd[["sigma_d_vel"]],
            dfr[[col for ratio in ["I / 0", "II / 0", "MYSO / 0"]
                 for col in [ratio, f"E({ratio})"] ]],
        ],
        axis="columns",
    )[mask]
    # Add new columns
    df = df.assign(
       wave0=df.wave / (1.0 + vsys / LIGHT_SPEED_KMS),
       sig_flux=df.flux / df.s_n, 
       sig_wave0=df.wave * np.hypot(df.sigma_d_vel, d_vsys) / LIGHT_SPEED_KMS
    ).drop(columns=["ID", "wave", "e_wave", "s_n", "sigma_d_vel"])
    priority = ["wave0", "sig_wave0", "flux", "sig_flux"]
    remainder = [c for c in df.columns if c not in priority]
    df = df[priority + remainder]
    print(df)
    df.to_csv(f"all-lines-{id_label}/uil-final-table.csv")


if __name__ == "__main__":
    typer.run(main)
