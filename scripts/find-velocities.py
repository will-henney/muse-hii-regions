from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value
DWAVE = 1.25                    # Wavelength pixel size in angstrom

def extract_wave(s: str):
    """Extract rest wavelength from the ID field if possible"""
    wave_string = s.rstrip("+").split(" ")[-1]
    try:
        wave = float(wave_string)
    except ValueError:
        wave = np.nan
    return wave


def main(
        id_label: str,
        species_file: str="species.yaml",
        zones_file: str="zones.yaml",
):
    """Find Doppler velocities for all identified lines"""

    with open(species_file) as f:
        info = yaml.safe_load(f)
    with open(zones_file) as f:
        zones = yaml.safe_load(f)

    common_cols = ["Type", "ID"]
    df_f = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv").set_index("Index")
    df_g = pd.read_csv(f"all-lines-{id_label}/line-gauss-fluxes.csv").set_index("Index")
    df_wave = pd.read_csv(f"all-lines-{id_label}/line-gauss-waves.csv").set_index("Index")
    df_width = pd.read_csv(f"all-lines-{id_label}/line-gauss-widths.csv").set_index("Index")
    df_err = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv").set_index("Index")

    df0 = df_f[common_cols]

    df1 = pd.concat(
        [df.drop(columns=common_cols) for df in [df_f, df_g, df_wave, df_width, df_err]],
        axis=1,
        keys=["flux", "gflux", "wave", "width", "err"],
    )
    df = pd.concat([df0, df1], axis=1)

    # Get rest wavelengths and blend status
    df = df.assign(
        blend=df.ID.str.endswith("+"),
        wave0=[extract_wave(_) for _ in df.ID],
    )

    print(df[common_cols + ["blend", "wave0"]])

    for zone in zones:
        zlabel = zone["label"]
        # Use greater of Gaussian fit and 3-pixel flux
        flux = np.maximum(df[("flux", zlabel)], df[("gflux", zlabel)])
        # Signal to noise
        s_n = flux / df[("err", zlabel)]
        # Is it a significant measurement?
        msig = s_n > 1.0
        # Estimate error on wavelength
        e_wave = DWAVE / s_n
        wave = np.where(msig, df[("wave", zlabel)], np.nan)
        d_wave = np.where(msig, df[("width", zlabel)], np.nan)
        vel = LIGHT_SPEED_KMS * (wave - df.wave0) / df.wave0
        d_vel = LIGHT_SPEED_KMS * d_wave / df.wave0

        dfz = df[common_cols + ["blend", "wave0"]].assign(
            wave=wave,
            e_wave=e_wave,
            d_wave=d_wave,
            vel=vel,
            d_vel=d_vel,
            flux=flux,
            s_n=s_n,
        )
        dfz.to_csv(f"all-lines-{id_label}/{zlabel}-velocities.csv")

if __name__ == "__main__":
    typer.run(main)
