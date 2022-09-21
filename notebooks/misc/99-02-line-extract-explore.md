---
jupyter:
  jupytext:
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Exporatory notebook for enhancing/replacing 05-04

We want to be able to improve on the previous line extraction in two ways:

1. Remove the sky component where present
2. Use a better velocity window that includes the entire line

We will try and develop some functions here, which will eventually be incoporated into whispy. 

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import regions
import sys
import warnings
import whispy as wp

sns.set_context("talk")
sns.set_color_codes()
```

```python

```
