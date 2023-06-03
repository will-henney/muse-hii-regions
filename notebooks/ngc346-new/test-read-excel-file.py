# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# + [markdown] pycharm={"name": "#%% md\n"}
# # Tests of reading in excel spreadsheet
#
# This has been exported from google sheets. I am particularly interested in being able to get hold of the notes

# + pycharm={"name": "#%%\n"}
import pandas as pd
from pathlib import Path

# + pycharm={"name": "#%%\n"}
datapath = Path.cwd().parent.parent / "data" / "spec1d"

# + pycharm={"name": "#%%\n"}
df = pd.read_excel(datapath / "All-Lines-MUSE-NGC-346.xlsx")

# + pycharm={"name": "#%%\n"}
df

# + [markdown] pycharm={"name": "#%% md\n"}
# So, that gave me the table, but I do not see any of the notes.

# + pycharm={"name": "#%%\n"}
set(str(x).rstrip('?') for x in df.Type)

# + pycharm={"name": "#%%\n"}
import openpyxl

# + pycharm={"name": "#%%\n"}
workbook = openpyxl.load_workbook(datapath / "All-Lines-MUSE-NGC-346.xlsx", data_only=True)
sheet = workbook.active
sheet

# + pycharm={"name": "#%%\n"}
pd.DataFrame(list(sheet.values)[1:], columns=list(sheet.values)[0])


# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}
cell = sheet['A1']

# + pycharm={"name": "#%%\n"}
cell.comment.text, cell.comment.author

# + [markdown] pycharm={"name": "#%% md\n"}
# Aha, so it turns out that the notes are called "Comments" in the excel file. OK, that is easy then.

# + pycharm={"name": "#%%\n"}
pd.DataFrame(
    [
        [
            # Take content from non-empty comments
            x.comment.content if x.comment else ""
            for x in row
            # Only use the columns we want
            if x.column_letter in "ABCDEFGHI"
        ]
        for row in sheet.rows
        # And only use rows that have at least some data
        if any(x.value for x in row)
    ]
)

# + pycharm={"name": "#%%\n"}
# pd.read_excel??

# + pycharm={"name": "#%%\n"}
any(x.value for x in list(sheet.rows)[468])

# + [markdown] pycharm={"name": "#%% md\n"}
# So, the above works for automatically extracting all the notes in the same shape as the original table

# + [markdown] pycharm={"name": "#%% md\n"}
# What happens when it is a Sheets "comment" rather than "note"?

# + pycharm={"name": "#%%\n"}
sheet['A468'].comment


# + [markdown] pycharm={"name": "#%% md\n"}
# The author is included in the content field, while the author field is still `None`

# + [markdown] pycharm={"name": "#%% md\n"}
# What happens when we have several comments on the same cell?

# + pycharm={"name": "#%%\n"}
x = sheet['F455'].comment
x

# + pycharm={"name": "#%%\n"}
x.content

# + [markdown] pycharm={"name": "#%% md\n"}
# They all get concatenated into the `.content` field

# + pycharm={"name": "#%%\n"}

