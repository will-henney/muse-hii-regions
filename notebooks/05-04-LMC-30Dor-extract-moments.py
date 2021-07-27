# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
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

# # Extract line moments to image files for 30 Dor fields

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import regions
import sys
sys.path.append("../lib")
import moments
import extract

sns.set_context("talk")
sns.set_color_codes()
# -

labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-46-55-contsub.fits") 
    for label in labels
}

csub

from dataclasses import dataclass

# +
C_KMS = 3e5

@dataclass
class EmissionLine:
    """Emission line"""
    name: str
    wav0: float
    vlim: tuple[float, float] = (-50.0, 500.0)
    
    def save_moments(self, cube: Cube, prefix: str):
        wav1 = self.wav0 * (1.0 + self.vlim[0] / C_KMS)
        wav2 = self.wav0 * (1.0 + self.vlim[1] / C_KMS)
        moments.save_moments_to_fits(
            moments.find_moments(
                cube.select_lambda(wav1, wav2)
            ),
            label=self.name,
            flabel=prefix,
            restwav=self.wav0,
            irange=None,
            vrange=None,
            srange=None,
        )
            


# -

moments.SAVEPATH = Path("../data")

emlines = [
    EmissionLine("heii-4686", 4685.68),
    EmissionLine("ariv-4740", 4740.17),
    EmissionLine("ni-5198-5200", 5199.12),
    EmissionLine("feiii-5270", 5270.40),
    EmissionLine("feiii-4658", 4658.10),
    EmissionLine("hi-4861", 4861.32),
    EmissionLine("hei-4922", 4921.93),
    EmissionLine("oiii-5007", 5006.84),
    EmissionLine("si-ii-5041", 5041.03),
]

for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")

emlines = [
    EmissionLine("oii-4650", 4650.0, (-600.0, 600.0)),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")


