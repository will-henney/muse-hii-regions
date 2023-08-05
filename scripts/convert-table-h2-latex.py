import pandas as pd
import numpy as np


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
        return "0.0"
    try:
        svalue = str(round(float(value), dp))
    except TypeError:
        svalue = "0.0"

    if uncertainty in ["<", ">"]:
        return f"{uncertainty} {svalue}"

    try:
        suncertainty = str(round(float(uncertainty), dp))
    except TypeError:
        suncertainty = "0.0"

    return rf"{svalue} \pm {suncertainty}"


df = pd.read_csv("h2-final-table.csv")
table = []
for _, row in df.iterrows():
    if row["blend"]:
        table.append(
            {
                r"\lambda(\text{obs})": format_blend(row["wave0"]),
                r"Transition": row["H2_line"],
                r"\lambda(\text{rest})": format_float(str(row["wl_lab"]).rstrip("A"), dp=2),
                r"$I(\hb = 100)$": format_blend(row["flux"], 3),
                r"$\mathrm{I / 0}$": "",
                r"$\mathrm{II / 0}$": "",
                r"$\mathrm{MYSO / 0}$": "",
            }
        )
    else:
        table.append(
            {
                r"\lambda(\text{obs})": format_blend(row["wave0"]),
                r"Transition": row["H2_line"],
                r"\lambda(\text{rest})": format_float(str(row["wl_lab"]).rstrip("A"), dp=2),
                r"$I(\hb = 100)$": format_pair(row["flux"], row["sig_flux"], 3),
                r"$\mathrm{I / 0}$": format_pair(row["I / 0"], row["E(I / 0)"]),
                r"$\mathrm{II / 0}$": format_pair(row["II / 0"], row["E(II / 0)"]),
                r"$\mathrm{MYSO / 0}$": format_pair(
                    row["MYSO / 0"], row["E(MYSO / 0)"]
                ),
            }
        )
    dff = pd.DataFrame(table)  

s = dff.style.hide()


with open("h2-final-table.tex", "w") as f:
    f.write(
        s.to_latex(
            hrules=True,
            siunitx=False,
            environment="longtable",
            column_format="RrR RRRR",
        )
    )
