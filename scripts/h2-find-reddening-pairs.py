import numpy as np
import yaml
from pathlib import Path
import pandas as pd
import typer
import sys
COLS = ["Index", "H2_line",  "Vhi",  "Jhi",  "Vlo",  "Jlo", "wl_mic", "wl_lab", "I_Inorm",  "Excit_hi_K", "g_u_h_nu_Aul", 
        "F(0)", "E(0)", "F(MYSO)", "E(MYSO)"]

def main(
            data_dir: str="n346-h2-lines",
            min_snr: float=2.0,
):
      df1 = pd.read_csv(Path(data_dir) / "known-lines-with-h2-final-table.csv")
      df1 = (df1       
             # removing low s/n lines
             [df1["F(0)"] > min_snr * df1["E(0)"]]
             # removing blends
             [~df1["Notes"].str.contains("blend", na=False)]
             #[~df1["blend"]]
             # And other dubious lines
             [~df1["Notes"].str.contains("elam = 2 sig", na=False)]
             )
      df2 = pd.read_csv(Path(data_dir) / "h2-line-ids.csv")
      # Merge the two tables on the name of the H2 lines
      grouped = (df2.merge(df1, how="inner", on="H2_index")[COLS]
                 # sorted by energy, then wavelength
                 .sort_values(by=["Excit_hi_K", "wl_mic"])
                 .groupby("Excit_hi_K")
                 )
      # select groups of two or more lines with the same upper level and save to csv
      grouped.filter(lambda x: len(x) > 1).to_csv(Path(data_dir) / "h2-reddening-lines-table.csv")

      # now reorganise the date with one row per pair
      red_pairs = []
      for name, group in grouped:
            if len(group) > 1:
                  print(name)
                  group["cloudy_theo"] = 1e-16 * 100 * group["I_Inorm"] / group["g_u_h_nu_Aul"]
                  group["obs_theo"] = 1e-16 * group["F(0)"] / group["g_u_h_nu_Aul"]
                  group["s_n"] = group["F(0)"] / group["E(0)"]
                  group["wav_ratio"] = group.iloc[-1]["wl_mic"] / group["wl_mic"]

                  group = group.drop(columns=["Excit_hi_K", "g_u_h_nu_Aul", "I_Inorm"])
                  print(group)
                  # take ratio of the shortest wave with each longer wave in turn
                  short = group.iloc[0, :]
                  for _, long in group.iloc[1:, :].iterrows():
                        if short.H2_line == long.H2_line:
                              continue
                        red_pairs.append(
                              dict(
                                    lines = f"{short.H2_line} / {long.H2_line}",
                                    waves = f"{short.wl_lab.strip('A')} / {long.wl_lab.strip('A')}",
                                    obs_theo_ratio = short.obs_theo / long.obs_theo,
                                    wav_ratio = long.wl_mic / short.wl_mic,
                                    Texcit = name,
                                    long_wav = long.wl_mic,
                                    rel_error = np.hypot(1.0 / short.s_n, 1.0 / long.s_n),
                                    bright_product = 100 * short["F(0)"] * long["F(0)"],
                              )
                        )
      # Save list of pairs
      df0 = pd.DataFrame(red_pairs).sort_values("wav_ratio")
      df0.to_csv(Path(data_dir) / "h2-reddening-pairs-table.csv")
      print(df0)


if __name__ == "__main__":
    typer.run(main)
