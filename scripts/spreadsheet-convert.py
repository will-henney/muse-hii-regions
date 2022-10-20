import pandas as pd
import sys
from pathlib import Path
import typer
import openpyxl
import yaml
import slugify

def unpack_notes_from_string(s):
    """Extract list of notes from string

    Filter out the separators and author bylines
    """
    return [
        note for note in s.split("\n")
        if not note.startswith(("----", "\t-"))
    ]


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
    #sys.exit(str(kwds))

    # Notes to each cell are called comments in the API
    notes_array = [[x.comment.content if x.comment else None for x in row] for row in sheet.rows]

    # Loop over all the following rows
    for values, notes in zip(values_array[1:], notes_array[1:]):
        if not any(values):
            # Skip any blank rows
            continue
        # Make a dict of the data from this row
        data = dict(zip(kwds, values))
        # Add the notes, but only where they exist
        if any(notes):
            data["Notes"] = {k: unpack_notes_from_string(x) for k, x in zip(kwds, notes) if x}

        # We use the Index column padded to 4 digits to construct the file stem
        index = data["Index"] = int(data["Index"])
        stem = f"{index:04d}"
        # Save the data to a JSON file
        with open(out_path / f"{stem}.yaml", "w") as f:
            yaml.dump(data, f, allow_unicode=True, sort_keys=False, default_flow_style=False)



if __name__ == "__main__":
    typer.run(main)
