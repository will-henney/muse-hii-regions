from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys
from textwrap import dedent

MATCH_LINES = {
    71: ("He  II", 4685.68),
    114: ("[Ar IV]", 4740.17),
    328: ("[O  III]", 5006.84),
    737: ("[Cl III]", 5517.71),
    1023: ("He  I", 5875.62),
    1363: ("[O  I]", 6300.30),
    1372: ("[S  III]", 6312.06),
    1573: ("H   I", 6562.79),
    1590: ("[N  II]", 6583.45),
    2032: ("[Ar III]", 7135.78),
    2524: ("[Ar III]", 7751.10),
    3306: ("[C  I]", 8727.13),
    3579: ("[S  III]", 9068.90),
}

def main(
        id_label: str,
        species_file: str="species.yaml",
        zone_file: str="zones.yaml",
        v_sys: float=171.1,
        d_v_sys: float=2.7,
):
    """Write input files for EMILI"""

    # Load the species and line type database
    with open(species_file) as f: 
        info = yaml.safe_load(f)
    # Load the zone database
    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    for zone in zones:
        prefix = f"all-lines-{id_label}/{zone['label']}"
        oprefix = f"emili/{zone['label']}"
        # Read file and drop rows with no observed wavelength
        df = pd.read_csv(
            f"{prefix}-velocities.csv"
        ).set_index("Index").dropna(subset="wave")
        # Make sure that wavelength error is not below 0.3 Angstrom
        e_wave = np.hypot(0.3, df.e_wave)
        # Select the columns we want and in the order we want
        dff = df[["wave"]].assign(
            ewm=-e_wave,
            ewp=e_wave,
            flux=df.flux / 100.0,
            fwhm=3e5 * df.d_wave / df.wave,
            sn=df.s_n,
            label=df.ID,
        )
        # Write to fixed width columns file .in
        with open(f"{oprefix}-emi.in", "w") as f:
            np.savetxt(
                f,
                dff.values,
                fmt="%7.2f %5.2f %5.2f %10.2e %6.2f %9.2f %20s",
            )
        # And write out the match lines to .match too
        with open(f"{oprefix}-emi.match", "w") as f:
            for idx, (ion, wave0) in MATCH_LINES.items():
                try:
                    f.write(
                        f"{df.loc[idx].wave:.2f} {wave0:.2f} {ion:10s} {df.loc[idx].flux / 100:.2e}\n"
                    )
                except KeyError:
                    ...

        # And finally write out the .cmd file
        with open(f"{oprefix}-emi.cmd", "w") as f:
            zlabel = zone["label"]
            f.write(
                dedent(
                    f"""\
                    A abun.dat
                    M {zlabel}-emi.match
                    O {zlabel}-emi.out
                    D {zlabel}-emi.dat
                    T 10000
                    N 100
                    I 100
                    L {zlabel}-emi.in
                    Z COMMENT vel 171 171 171 171 171
                    vel+
                    icf+
                    """
                )
            )
            if "MYSO" not in zlabel:
                f.write(
                    "\n".join([f"deplete {element} 50" for element in ["Fe", "Ni", "Si", "Ca"]])
                )

if __name__ == "__main__":
    typer.run(main)
