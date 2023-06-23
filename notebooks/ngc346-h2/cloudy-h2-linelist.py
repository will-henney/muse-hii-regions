# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py,md
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Make list of optical H_2 lines
#
# For comparison with our Deep Red Lines

from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from astropy import units as u
from specutils.utils.wcs_utils import vac_to_air, air_to_vac
from astropy.table import QTable

cloudy_path = Path.home() / "Work/PNe/Peter/herschel-helix-models/models/static/"
save_path = Path.cwd().parent.parent / "data-atomic"

# Read in the molecular hydrogen line list.  This is a Helix Nebula model from an abandoned project from 2013. Eventually, I should replace it with a better model for NGC 346

df = pd.read_csv(cloudy_path / "static_n20big.h2l", sep="\t")

# Take the second half of the file and remove the Ehi, Elo columns that we do not need

df = df.loc[618:]
del df['Ehi']
del df['Elo']

# Select all lines with wavelengths less than 0.93 micron and sort by wavelength

df = df.query("`wl(mic)` < 0.93").sort_values("wl(mic)")
df

# It turns out that the wavelengths from Cloudy are already on an air scale

df = df.assign(
    wav_air=(df["wl(mic)"].to_numpy() * u.micron).to(u.Angstrom).value,
)

df

df.to_csv(save_path / "cloudy-h2-lines.csv")

# Now sort by upper level 

df = df.sort_values(["Vhi", "Jhi"])
df

df.to_csv(save_path / "cloudy-h2-lines-sort-upper.csv")

    
