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

# Issue with multiplying two cubes with different units



Reported as issue on MPDAF github site 2021-05-10

https://github.com/musevlt/mpdaf/issues/20

```python
import numpy as np
from mpdaf.obj import Cube
import astropy.units as u
a = Cube(data=np.ones((1, 1, 1)), unit=u.s)
b = Cube(data=np.ones((1, 1, 1)), unit=u.s)
c = Cube(data=np.ones((1, 1, 1)), unit=u.cm)
# Multiplying cubes with the same units is OK
a * b
# Multiplying cubes with different units yields an error
a * c
```

```python

```
