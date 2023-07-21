import pandas as pd
import numpy as np

def format_blend(value, dp=2):
    if not value or not np.isfinite(value):
        return "0.0:"
    try:
        return str(round(float(value), dp)) + ":"
    except TypeError:
        return "0.0:"

def format_pair(value, uncertainty, dp=2):
    if not value or not np.isfinite(value):
        return "0.0"
    try:
        svalue = str(round(float(value), dp))
    except TypeError:
        svalue = "0.0"

    if uncertainty in  ["<", ">"]:
        return f"{uncertainty} {svalue}"

    try:
        suncertainty = str(round(float(uncertainty), dp))
    except TypeError:
        suncertainty = "0.0"

    return fr"{svalue} \pm {suncertainty}"


def format_ion(ion):
    """Convert naive representation to latex version"""
    # Strip off any surrounding brackets
    prefix, suffix = "", ""
    if ion.startswith("["):
        ion = ion[1:]
        prefix = "["
    if ion.endswith("]"):
        ion = ion[:-1]
        suffix = "]"
    # Split into element, stage
    element, stage = ion.rsplit(maxsplit=1)
    # Convert stage to arabic numerals
    stage = stage.replace("IV", "4").replace("III", "3").replace("II", "2").replace("I", "1")
    return fr"{prefix}\ion{{{element}}}{{{stage}}}{suffix}"


df = pd.read_csv("known-lines-final-table.csv")
table = []
zones = ["0", "I", "II", "III", "IV", "MYSO"]
for _, row in df.iterrows():
    ion, wavrest = row["ID"].rsplit(maxsplit=1)
    table.append(
        {
            r"{\(\lambda(\text{obs})\)}": format_pair(row["wave"], row["e_wave"]),
            "Ion": format_ion(ion),
            r"{\(\lambda(\text{rest})\)}": wavrest.strip("+"),
            **{
                f"F({zone})": format_pair(row[f"F({zone})"], row[f"E({zone})"])
                for zone in zones
            },
            # r"Notes": format_notes(row["flux"], row["sig_flux"], 3),
        }
    )
    dff = pd.DataFrame(table)

s = dff.style.hide()



with open("known-lines-final-table.tex", "w") as f:
    f.write(
        s.to_latex(
            hrules=True,
            siunitx=True,
            environment="longtable",
            column_format="SSS" + "S" * len(zones),
        )
    )
