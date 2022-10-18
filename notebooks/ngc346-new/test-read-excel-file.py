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
import openpyxl

# + pycharm={"name": "#%%\n"}
workbook = openpyxl.load_workbook(datapath / "All-Lines-MUSE-NGC-346.xlsx")
sheet = workbook.active
sheet

# + pycharm={"name": "#%%\n"}
pd.DataFrame(sheet.values)


# + pycharm={"name": "#%%\n"}
cf = sheet.conditional_formatting

# + pycharm={"name": "#%%\n"}
cell = sheet['A1']

# + pycharm={"name": "#%%\n"}
cell.comment.text, cell.comment.author

# + [markdown] pycharm={"name": "#%% md\n"}
# Aha, so it turns out that the notes are called "Comments" in the excel file. OK, that is easy then.

# + pycharm={"name": "#%%\n"}
pd.DataFrame([[x.comment for x in row] for row in sheet.rows])

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

