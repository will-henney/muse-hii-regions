import pandas as pd
import numpy as np


df = pd.read_csv("uil-final-table.csv")

# For Zone 0, sum all the relative fluxes, however weak
sum_0 = df["flux"].sum()
median_0 = df["flux"].median()
mean_0 = df["flux"].mean()

# For MYSO, first back out the ratio to get relative fluxes
flux_myso = df["MYSO / 0"] * df["flux"]
# And then mask out lines that are blends or upper limits before summing
invalid = (df["E(MYSO / 0)"] == "<") | df["blend"]
flux_myso[invalid] = np.nan
sum_myso = np.nansum(flux_myso)
median_myso = np.nanmedian(flux_myso)
mean_myso = np.nanmean(flux_myso)

# Report statistics of sums and counts
print(f"Sum of fluxes: Zone 0 = {sum_0:.2f}, MYSO = {sum_myso:.2f}")
print(f"Number of lines: Zone 0 = {len(df)}, MYSO = {np.sum(~invalid)}")
print(f"Mean brightness: Zone 0 = {mean_0:.3f}, MYSO = {mean_myso:.3f}")
print(f"Median brightness: Zone 0 = {median_0:.3f}, MYSO = {median_myso:.3f}")
