from matplotlib import pyplot as plt
import matplotlib
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import typer
import sys

def find_ratios_with_errors(
        df: pd.DataFrame,
        zone_ratio: str="II / 0",
) -> pd.DataFrame:
    """
    Find the zone ratios and errors for all emission lines

    Returns new dataframe with added columns.  This uses vector
    version of the same logic as in the plotting program, but is in
    some ways simpler since it only deals with one ratio, not two.
    """
    assert " / " in zone_ratio
    assert type(df) == pd.DataFrame

    # Fluxes in each zone, assuming ratio is B / A
    xB, xA  = [df[f"zone-{zone}"].to_numpy(copy=True) for zone in zone_ratio.split(" / ")]
    # Corresponding uncertainties
    xBe, xAe  = [df[f"zone-{zone}_sig"].to_numpy(copy=True) for zone in zone_ratio.split(" / ")]

    # Fluxes that are negative or smaller than the uncertainty are
    # replaced with 2-sigma upper limits
    #
    # For the numerator B, this gives upper limit to the ratio
    m_upper = xB < xBe
    xB = np.where(
        m_upper,
        np.maximum(xB, 0.0) + 2 * xBe,
        xB
    )
    # For the denominator A, it gives lower limit to the ratio
    m_lower = xA < xAe
    xA = np.where(
        m_lower,
        np.maximum(xA, 0.0) + 2 * xAe,
        xA
    )
    # Then calculate the ratio
    xratio = np.where(
        # Ratio is undefined if both fluxes are upper limits
        m_upper & m_lower,
        np.nan,
        xB / xA
    )
    # And the error
    xerr = np.hypot(xBe / xB, xAe / xA) * xratio
    e_zone_ratio = f"E({zone_ratio})"
    _df = pd.DataFrame(
        {
            zone_ratio: xratio,
            e_zone_ratio: xerr,
        },
        index=df.index,
    )
    # Patch in indicators of upper and lower limits in err column
    _df.loc[m_upper, e_zone_ratio] = "<"
    _df.loc[m_lower, e_zone_ratio] = ">"
    return _df

def main(
        id_label: str,
        species_file: str="species.yaml",
        use_gauss: bool=False,
):
    """Write table of ratios between the zones"""

    outfile = f"all-lines-{id_label}/ratios-vs-ratios-by-zone.csv"

    if use_gauss:
        df = pd.read_csv(f"all-lines-{id_label}/line-gauss-fluxes.csv")
        outfile = outfile.replace(".csv", "-gauss.csv")
    else:
        df = pd.read_csv(f"all-lines-{id_label}/line-fluxes.csv")
    df = df.set_index("Index")
    df_sig = pd.read_csv(f"all-lines-{id_label}/line-uncertainties.csv").set_index("Index")
    df = df.join(df_sig, rsuffix="_sig")

    dfr = df[["ID"]]
    ratio_frames = [
        find_ratios_with_errors(df, zone_ratio)
        for zone_ratio in [
                "I / 0", "II / 0", "MYSO / 0", "III / II", "IV / III",
        ]
    ]
    dfr = pd.concat([dfr] + ratio_frames, axis="columns")

    dfr.to_csv(outfile)



if __name__ == "__main__":
    typer.run(main)
