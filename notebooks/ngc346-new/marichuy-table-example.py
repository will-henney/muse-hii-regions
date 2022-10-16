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

# + pycharm={"name": "#%%\n"}
from astropy.table import Table

# + pycharm={"name": "#%%\n"}
data_dict = {
    "Sources": ["A", "B", "C"],
    "CO (3-2)": [1.0, 2.0, 3.0],
    "CO (2-1)": [99.0, 50.0, 0.0],
}

# + pycharm={"name": "#%%\n"}
tab = Table(data_dict)
tab

# + pycharm={"name": "#%%\n"}
tab.write("sources.votable", format="votable")

# + pycharm={"name": "#%%\n"}

