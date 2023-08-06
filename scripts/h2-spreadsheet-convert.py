import pandas as pd
import sys
from pathlib import Path
import typer
import openpyxl
import yaml
import slugify

def main(
        excel_file: str,
        out_folder: str="n346-lines",
):
    """Convert excel spreadsheet of emission lines to YAML files, one per row

    Preserves Notes and Comments on each cell
    """
    # Read in the spreadsheet
    workbook = openpyxl.load_workbook(excel_file, data_only=True)
    # And select the first sheet
    sheet = workbook.active

    # Make a list of row data from the sheet
    values_array = list(sheet.values)

    # Make sure the output folder exists
    out_path = Path(out_folder)
    out_path.mkdir(parents=True, exist_ok=True)

    # Column headers are in first row
    kwds = [
        # Try to make sure headers are valid identifiers
        slugify.slugify(str(x), lowercase=False, separator="_", replacements=[["λ", "lambda"]])
        for x in values_array[0]
        # And skip empty columns
        if x
    ]
    # sys.exit(str(kwds))

    # Loop over all the following rows
    for index, values in enumerate(values_array[1:], start=2):
        if not any(values):
            # Skip any blank rows
            continue
        # Make a dict of the data from this row
        data = dict(zip(kwds, values))
        # Save row number in the Google Sheet for cross-referencing of blends
        data["H2_index"] = index
        # For any lines that do not have an observed counterpart ...
        if not data["wl_obs"]:
            # ... check if they are blends with next or previous
            if data["Notes"] and data["Notes"].startswith("blend with"):
                # save pointers  to the line they are blended with
                if "blend with prev" in data["Notes"]:
                    data["blend_index"] = index - 1
                elif "blend with next" in data["Notes"]:
                    data["blend_index"] = index + 1
                else:
                    continue
            else:
                # ... othewise skip
                continue
        # rovib quantum numbers of upper and lower states
        vhi, vlo = int(data["Vhi"]), int(data["Vlo"])
        jhi, jlo = int(data["Jhi"]), int(data["Jlo"])
        # Make a string for the transition
        transition = f"{vhi:02d}_{jhi:02d}-{vlo:02d}_{jlo:02d}"
        wavstring = slugify.slugify(data["wl_lab"], separator="")
        # Save the data to a YAML file
        with open(out_path / f"{index:04d}-{transition}-{wavstring}.yaml", "w") as f:
            yaml.dump(data, f, allow_unicode=True, sort_keys=False, default_flow_style=False)



if __name__ == "__main__":
    typer.run(main)
