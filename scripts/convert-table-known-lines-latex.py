import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import re
import string

LETTERS = list(string.ascii_lowercase) + [f"{a}{b}" for a in string.ascii_lowercase for b in string.ascii_lowercase]

def format_blend(value, dp=2):
    if not value or not np.isfinite(value):
        return "0.0:"
    try:
        return str(round(float(value), dp)) + ":"
    except TypeError:
        return "0.0:"


def format_float(value, dp):
    return f"{round(float(value), dp):.{dp}f}"


def format_pair(value, uncertainty, dp=2):
    if not value or not np.isfinite(value):
        return ""
    try:
        svalue = format_float(value, dp)
    except TypeError:
        svalue = ""

    if uncertainty in  ["<", ">"]:
        return f"{uncertainty} {svalue}"

    # Special case of H beta has no uncertainty
    if svalue == "100.00":
        return f"{svalue}"

    try:
        suncertainty = format_float(uncertainty, dp)
    except TypeError:
        suncertainty = ""

    return fr"{svalue} \pm {suncertainty}"


def format_ion(ion):
    """Convert naive representation to latex version"""
    # Short circuit for sky lines
    if "OH" in ion:
        return r"Sky OH"

    # Strip off any surrounding brackets
    prefix, suffix = "", ""
    if ion.startswith("["):
        ion = ion[1:]
        prefix = "["
    if ion.endswith("]"):
        ion = ion[:-1]
        suffix = "]"
    # And just in case of extra brackets
    ion = ion.strip("[]")
    # Split into element, stage
    element, stage = ion.rsplit(maxsplit=1)
    if " " in element.strip():
        # Something has gone wrong if there are still internal spaces
        raise ValueError(f"Cannot parse ion {ion}")
    # Convert stage to arabic numerals
    stage = stage.replace("IV", "4").replace("III", "3").replace("II", "2").replace("I", "1")
    return fr"{prefix}\ion{{{element}}}{{{stage}}}{suffix}"


def extra_e_wave(row, zones):
    """Extra uncertainty in wavelength, based on Fig A2 of paper"""
    max_flux = np.nanmax([row[f"F({zone})"] for zone in zones])
    if max_flux > 10.0:
        return row["wave"] * 1.0 / 3e5
    elif max_flux > 1.0:
        return row["wave"] * 3.0 / 3e5
    elif max_flux > 0.1:
        return row["wave"] * 5.0 / 3e5
    else:
        return row["wave"] * 15.0 / 3e5

def extract_blends(text):
    """Extract blend information from notes"""
    if not text:
        return []
    if "plus" in text:
        blends = re.split("plus|and maybe", text)
    elif text.startswith("Blend with"):
        blends = text.replace("Blend with", "").split(",")
    elif text.startswith("Doublet with components"):
        blends = [text.split(",")[-1]]
    else:
        blends = [text]
    return blends

df = pd.read_csv("known-lines-final-table.csv")
table = []
zones = ["0", "I", "II", "III", "IV", "MYSO"]
# Special case for [Cl II] and He I lines that are sort of blended, but both have flux measurements
# And also He I triplet
skip_these_blends = ("8578", "8582", "8776", "8045", "8216", "8230", "8306")
for _, row in df.iterrows():
    ion, wavrest = row["ID"].rsplit(maxsplit=1)
    # Special case for [N I], where I unwisely put the mean doublet wavelength in the spreadsheet
    wavrest = wavrest.replace("5199.00", "5197.98")
    if row["blend"] and not wavrest.startswith(skip_these_blends):
        blend_label = LETTERS.pop(0)
    else:
        blend_label = ""
    table.append(
        {
            r"\lambda(\text{obs})": format_pair(row["wave"], np.hypot(row["e_wave"], extra_e_wave(row, zones))),
            "Ion": format_ion(ion.strip()),
            r"\lambda(\text{rest})": format_float(wavrest.strip("+?"), dp=2),
            "Blend": blend_label,
            **{
                fr"I(\text{{{zone}}})": format_pair(row[f"F({zone})"], row[f"E({zone})"] + 0.005, dp=2)
                for zone in zones
            },
        }
    )
    if blend_label:
        datafile = list(Path.cwd().glob(f"{row['Index']:04d}*.yaml"))[0]
        data = yaml.safe_load(datafile.read_text())
        note = data["Notes"]["ID"][0]
        blends = extract_blends(note)
        for blend in blends:
            try:
                _ion, _wavrest = blend.rsplit(maxsplit=1)
            except ValueError:
                # A bare wavelength is assumed to be the same ion
                _ion = ion
                _wavrest = blend
            if _wavrest.lower() == "sky":
                _ion, _wavrest = _ion.rsplit(maxsplit=1)
            try:
                table.append(
                    {
                        r"\lambda(\text{obs})": "",
                        "Ion": format_ion(_ion.strip()),
                        r"\lambda(\text{rest})": format_float(_wavrest.strip("+?"), dp=2),
                        "Blend": blend_label,
                        **{fr"I(\text{{{zone}}})": "" for zone in zones},
                    }
                )
            except:
                pass

    dff = pd.DataFrame(table).sort_values(by=r"\lambda(\text{rest})")

s = dff.style.hide()



with open("known-lines-final-table.tex", "w") as f:
    f.write(
        s.to_latex(
            hrules=True,
            siunitx=False,
            environment="longtable",
            column_format="RrRr " + "R" * len(zones),
        )
    )
