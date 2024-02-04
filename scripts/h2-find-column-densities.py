import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys
COLS = [
      "Index", "H2_line",  "Vhi",  "Jhi",  "Vlo",  "Jlo",
      "wl_mic", "wl_lab", "I_Inorm",  "Excit_hi_K", "g_u_h_nu_Aul", 
      "F(0)", "E(0)", "F(MYSO)", "E(MYSO)",
]

def main(
            data_dir: str=".",
            min_snr: float=2.0,
):
      df1 = pd.read_csv(Path(data_dir) / "known-lines-with-h2-final-table.csv")
      df1 = (df1       
             # removing low s/n lines
             [df1["F(0)"] > min_snr * df1["E(0)"]]
             # removing blends
             [~df1["Notes"].str.contains("blend", na=False)]
             [~df1["blend"]]
             # And other dubious lines
             [~df1["Notes"].str.contains("elam = 2 sig", na=False)]
             )
      df2 = pd.read_csv(Path(data_dir) / "h2-line-ids.csv")
      # Merge the two tables on the name of the H2 lines
      merged = (df2.merge(df1, how="inner", on="H2_index")[COLS]
                 # sorted by energy, then wavelength
                 .sort_values(by=["Excit_hi_K", "wl_mic"])
                 )
      merged.to_csv(Path(data_dir) / "h2-excitation-table.csv")


if __name__ == "__main__":
    typer.run(main)
