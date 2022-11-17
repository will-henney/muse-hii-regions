from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore


def main(
        id_label: str,
        species_file: str="species.yaml",
):
    """Make tables of cross nearest neighbor distances from UILs to Identified"""

    dfi = pd.read_csv(
        f"all-lines-{id_label}/nearest-neighbors-Identified.csv"
    ).set_index("Index").sort_index()
    dfu = pd.read_csv(
        f"all-lines-{id_label}/nearest-neighbors-Unidentified.csv"
    ).set_index("Index").sort_index()

    # Copy over the first few columns from the Unidentified table
    df = dfu.loc[:, "Type":"wave"].assign(
        # Use appply to vectorize over the wavelengths
        dwave_nn=dfu.wave.apply(
            # Find the minimum absolute separation with the Identified lines
            lambda x: (dfi.wave - x).abs().min()
        ),
        index_other=dfu.wave.apply(
            # And the index of that closest line in the Identified table
            lambda x: (dfi.wave - x).abs().idxmin()
        ),
    )
    # Now use that index to obtain other info. I had to resort to list
    # comprehensions here since pandas will not let you use repeated
    # values in the index
    df = df.assign(
        # The mean separation in the Identified lines
        dwave_mean=[dfi.dwave_mean[_i] for _i in df.index_other],
        # Whether the Identified line is closer to the UIL than to its
        # own nearest neighbor among the other Identifieds
        mutual=[this < dfi.dwave_nn[_i] for this, _i in zip(df.dwave_nn, df.index_other)],
    )
    df.to_csv(f"all-lines-{id_label}/nearest-neighbors-Cross.csv")





if __name__ == "__main__":
    typer.run(main)
