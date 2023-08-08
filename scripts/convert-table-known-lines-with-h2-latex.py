import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import re
import string

# List of labels for blends
LETTERS = list(string.ascii_lowercase) + [
    f"{a}{b}" for a in string.ascii_lowercase for b in string.ascii_lowercase
]
# Reversed so we can use it as a stack
LETTERS.reverse()

ORIG_DATA_PATH = Path.cwd().parent / "n346-lines/all-lines-c007-chop-mean"

VSYS = 171.1

MINIMUM_H2_BLEND_STRENGTH = 0.1


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

    if uncertainty in ["<", ">"]:
        return f"{uncertainty} {svalue}"

    # Special case of H beta has no uncertainty
    if svalue == "100.00":
        return f"{svalue}"

    try:
        suncertainty = format_float(uncertainty, dp)
    except TypeError:
        suncertainty = ""

    return rf"{svalue} \pm {suncertainty}"


def format_ion(ion):
    """Convert naive representation to latex version"""
    # Short circuit for sky lines
    if "OH" in ion:
        return r"Sky OH"

    if "O_2" in ion:
        return r"Sky \chem{O_2}"

    if ion.startswith("H_2"):
        molecule, transition = ion.split(maxsplit=1)
        return rf"\chem{{{molecule}}} {transition}"

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
    stage = (
        stage.replace("IV", "4")
        .replace("III", "3")
        .replace("II", "2")
        .replace("I", "1")
    )
    return rf"{prefix}\ion{{{element}}}{{{stage}}}{suffix}"


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


def extract_blends(notes):
    """Extract blend information from notes"""
    blends = []
    for text in notes:
        if not text:
            continue
        if "plus" in text:
            blends += re.split("plus|and maybe|,|and", text)
        elif text.startswith("Blend with"):
            blends += text.replace("Blend with", "").split(",")
        elif text.startswith("Doublet with components"):
            blends += [text.split(",")[-1]]
        else:
            blends += [text]
    return blends


df = pd.read_csv("known-lines-with-h2-final-table.csv")
table = []
zones = ["0", "I", "II", "III", "IV", "MYSO"]
skip_these_blends = (
    # New ones from the H2 lines
    "6270",
    "6529",
    "7803",
    "7837",
    "8694",
    "8304",
    "8680",
    # Original ones from the atomic list - Special case for [Cl II] and
    # He I lines that are sort of blended, but both have flux
    # measurements And also He I triplet
    "8578",
    "8582",
    "8776",
    "8045",
    "8216",
    "8230",
    "8306",
)
for _, row in df.iterrows():
    ion, wavrest = row["ID"].rsplit(maxsplit=1)
    # H_2 lines need the A stripping off
    wavrest = wavrest.rstrip("A")
    # Special case for [N I], where I unwisely put the mean doublet wavelength in the spreadsheet
    wavrest = wavrest.replace("5199.00", "5197.98")
    # For H_2 lines we look for extra blend info
    is_h2blend = False
    if ion.startswith("H_2"):
        h2_index = int(row["H2_index"])
        h2datafile = list(Path.cwd().glob(f"{h2_index:04d}*.yaml"))[0]
        h2data = yaml.safe_load(h2datafile.read_text())
        is_h2blend = h2data["Notes"] and h2data["Notes"].lower().startswith("blend")

    if (row["blend"] or is_h2blend) and not wavrest.startswith(skip_these_blends):
        blend_label = LETTERS.pop()
    else:
        blend_label = ""
    table.append(
        {
            r"\lambda(\text{obs})": format_pair(
                row["wave"], np.hypot(row["e_wave"], extra_e_wave(row, zones))
            ),
            "Species": format_ion(ion.strip()),
            r"\lambda(\text{rest})": format_float(wavrest.strip("+?"), dp=2),
            r"Wav sort": format_float(wavrest.strip("+?"), dp=2),
            "Blend": blend_label,
            **{
                rf"\text{{Zone {zone}}}": format_pair(
                    row[f"F({zone})"], row[f"E({zone})"] + 0.005, dp=2
                )
                for zone in zones
            },
        }
    )

    # Add extra entries for all the blends with this line
    #
    if blend_label:
        datafile = list(ORIG_DATA_PATH.glob(f"{row['Index']:04d}*.yaml"))[0]
        data = yaml.safe_load(datafile.read_text())
        try:
            notes = data["Notes"]["ID"]
            blends = extract_blends(notes)
        except KeyError:
            # Guard against missing Notes data
            blends = []

        # Also look for H_2 blends
        if is_h2blend:
            # Could be previous or next in sequence of lab wavelengths
            for j in h2_index - 1, h2_index + 1:
                bdatafiles = list(Path.cwd().glob(f"{j:04d}*.yaml"))
                if len(bdatafiles) == 1:
                    bdata = yaml.safe_load(bdatafiles[0].read_text())
                    # Check that this line has been tagged as a blend
                    # of the main line that we are looking at, and
                    # also that the predicted intensity is at least a
                    # certain fraction of the main line
                    if (
                        bdata.get("blend_index") == h2_index
                        and bdata["I_Inorm"]
                        > MINIMUM_H2_BLEND_STRENGTH * h2data["I_Inorm"]
                    ):
                        blends.append(
                            f"H_2 {bdata['H2_line']} {bdata['wl_lab'].strip('A')}"
                        )

        print(blend_label, blends)
        if len(blends) == 0:
            # No blends found, so we must hand back the label we took
            LETTERS.append(blend_label)
            print("Returning label", blend_label, "to the pool")
            # And also remove label from the main entry
            table[-1]["Blend"] = ""
            # And no point carrying on
            continue
        # What to assume when bare wavelength is given
        implicit_ion = ion
        for blend in blends:
            if not blend.strip():
                # Skip case of empty or blank string
                continue
            try:
                _ion, _wavrest = blend.strip().rsplit(maxsplit=1)
            except ValueError:
                # A bare wavelength is assumed to be the last ion that was mentioned
                _ion = implicit_ion
                _wavrest = blend
            if _wavrest.lower() == "sky":
                _ion, _wavrest = _ion.rsplit(maxsplit=1)
            _wavrest = _wavrest.strip("+?,").strip()
            if _wavrest.startswith(skip_these_blends):
                # Another chance to bail out
                print("Skipping", _ion, _wavrest)
                continue
            try:
                ion_string = format_ion(_ion.strip())
                wavrest_string = format_float(_wavrest, dp=2)
                if ion_string.startswith("Sky"):
                    # Make use of the lambda(obs) column to put rest
                    # wavelength of sky line in frame of the nebula
                    _wavobs = float(_wavrest) / (1.0 + VSYS / 3.0e5)
                    # And use this wavelength for sorting
                    wavsort = format_float(_wavobs, dp=2)
                    wavobs_string = r"\mathit{" + wavsort + "}"
                else:
                    # Other lines just use the lab rest wavelength
                    wavsort = wavrest_string
                    wavobs_string = ""
            except:
                print("Failed to format blend: ion =", _ion, "wav =", _wavrest)
                continue
            # Since this was a successfully identified blend, update
            # the ion that we use in the case that following blends have
            # bare wavelength
            implicit_ion = _ion
            table.append(
                {
                    "Wav sort": wavsort,
                    r"\lambda(\text{obs})": wavobs_string,
                    "Species": ion_string,
                    r"\lambda(\text{rest})": wavrest_string,
                    "Blend": blend_label,
                    **{rf"\text{{Zone {zone}}}": "" for zone in zones},
                }
            )
    dff = pd.DataFrame(table).sort_values(by="Wav sort").drop(columns=["Wav sort"])

s = dff.style.hide()


with open("known-lines-with-h2-final-table.tex", "w") as f:
    f.write(
        s.to_latex(
            hrules=True,
            siunitx=False,
            environment="longtable",
            column_format="RrRr " + "R" * len(zones),
        )
    )
