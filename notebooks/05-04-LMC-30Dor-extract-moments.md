---
jupyter:
  jupytext:
    encoding: '# -*- coding: utf-8 -*-'
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Extract line moments to image files for 30 Dor fields

```python
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import regions
import sys
import warnings
sys.path.append("../lib")
import moments
import extract

sns.set_context("talk")
sns.set_color_codes()
```

## Blue lines

```python
labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-46-55-contsub.fits") 
    for label in labels
}
```

```python
csub
```

```python
from dataclasses import dataclass
```

```python
C_KMS = 3e5

@dataclass
class EmissionLine:
    """Emission line"""
    name: str
    wav0: float
    vlim: tuple[float, float] = (100.0, 400.0)
    
    def save_moments(self, cube: Cube, prefix: str):
        wav1 = self.wav0 * (1.0 + self.vlim[0] / C_KMS)
        wav2 = self.wav0 * (1.0 + self.vlim[1] / C_KMS)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            moments.save_moments_to_fits(
                moments.find_moments(
                    cube.select_lambda(wav1, wav2)
                ),
                label=self.name,
                flabel=prefix,
                restwav=self.wav0,
                irange=None,
                vrange=[100.0, 400.0],
                srange=None,
            )
            


```

```python
moments.SAVEPATH = Path("../data")
```

The He II 4686 line is mainly a broad scattered component

```python tags=[]
emlines = [
    EmissionLine("nii-4607", 4607.12),
    EmissionLine("oii-4642", 4641.81),
    EmissionLine("oii-4650", 4649.985),
    EmissionLine("heii-4686", 4685.68, (-700, 900)),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("ariv-4740", 4740.17),
    EmissionLine("ni-5198-5200", 5199.12),
    EmissionLine("feiii-5270", 5270.40),
    EmissionLine("feiii-4658", 4658.10),
    EmissionLine("hi-4861", 4861.32),
    EmissionLine("hei-4922", 4921.93),
    EmissionLine("oiii-4959", 4958.91),
    EmissionLine("oiii-5007", 5006.84),
    EmissionLine("si-ii-5041", 5041.03),
]
```

```python
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("oii-4650", 4650.0, (-600.0, 600.0)),
    EmissionLine("feiii-4702", 4701.62),
    EmissionLine("ariv-4711", 4711.37),
    EmissionLine("feiii-4734", 4733.93),
    EmissionLine("feiii-4755", 4754.83),
    EmissionLine("feiii-4770", 4769.60),    
    EmissionLine("feiii-4881", 4881.073),
    EmissionLine("oiii-4931", 4931.32),
    EmissionLine("feiii-4987", 4987.20),    
    EmissionLine("hei-5016", 5015.68),
    EmissionLine("hei-5048", 5047.74),
    EmissionLine("si-ii-5056", 5055.98),
    EmissionLine("feii-5159", 5158.81),
    EmissionLine("ariii-5192", 5191.82),
    EmissionLine("feii-5262", 5261.61),
    EmissionLine("feiii-5412", 5412.00), 
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

## Green lines

```python tags=[]
labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-54-63-contsub.fits") 
    for label in labels
}
```

```python
emlines = [
    EmissionLine("cliii-5517", 5517.71), 
    EmissionLine("cliii-5538", 5537.88),
    EmissionLine("nii-5755", 5755.08),    
    EmissionLine("hei-5876", 5875.62),
    EmissionLine("si-ii-5958", 5957.56),
    EmissionLine("si-ii-5979", 5978.93),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python tags=[]
emlines = [
    EmissionLine("feiii-5412", 5412.00),
    EmissionLine("ni-iv-5820", 5820.10), 
    EmissionLine("wr-5840", 5840.0),
    EmissionLine("ni-iii-6000", 6000.20), 
    EmissionLine("oi-6046", 6046.23),
    EmissionLine("ca-i-6104", 6104.0),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("c-ii-6151", 6151.43),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

## Red lines

```python tags=[]
labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-62-71-contsub.fits") 
    for label in labels
}
```

```python
emlines = [
    EmissionLine("oi-6300", 6300.30), 
    EmissionLine("siii-6312", 6312.06),
    EmissionLine("si-ii-6347", 6347.11),
    EmissionLine("oi-6363", 6363.78),    
    EmissionLine("si-ii-6371", 6371.36),
    EmissionLine("cii-6462", 6461.95),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python tags=[]
emlines = [
    EmissionLine("nii-6548", 6548.05),
    EmissionLine("nii-6583", 6583.45), 
    EmissionLine("hei-6678", 6678.15),
    EmissionLine("sii-6716", 6716.44), 
    EmissionLine("sii-6731", 6730.816),
    EmissionLine("oi-7002", 7001.92),
    EmissionLine("hei-7065", 7065.28)
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-70-79-contsub.fits") 
    for label in labels
}
```

```python
emlines = [
    EmissionLine("ariii-7136", 7135.78),
    EmissionLine("feii-7155", 7155.14), 
    EmissionLine("hei-7160", 7160.13),
    EmissionLine("cii-7231", 7231.34), 
    EmissionLine("cii-7236", 7236.42),
    EmissionLine("oi-7254", 7254.15),
    EmissionLine("hei-7281", 7281.35)
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")

 emlines = [
     EmissionLine("hei-7281", 7281.35),
     EmissionLine("hei-7298", 7298.050),
     EmissionLine("oii-7318", 7318.39),
     EmissionLine("oii-7330", 7329.66), 
     EmissionLine("ni-ii-7378", 7377.83),
     EmissionLine("feii-7388", 7388.16),
     EmissionLine("ni-ii-7412", 7411.61),
     EmissionLine("ni-7442", 7442.30),
     EmissionLine("fe-ii-7453", 7452.54),
     EmissionLine("ni-7468", 7468.31),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")   
```

```python
emlines = [
    EmissionLine("ariv-7171", 7170.5),
    EmissionLine("ariv-7237", 7237.4), 
    EmissionLine("ariv-7263", 7262.7),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("cliv-7531", 7530.8),
    EmissionLine("ariii-7751", 7751.10),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("hei-7499", 7499.85),
    EmissionLine("feii-7516", 7515.83),
    EmissionLine("cai-7890", 7890.07, vlim=(150.0, 350.0)),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-78-87-contsub.fits") 
    for label in labels
}
```

```python
emlines = [
    EmissionLine("cr-ii-8000", 8000.08),
    EmissionLine("cliv-8046", 8045.62),
    EmissionLine("cai-8125", 8125.31),
    EmissionLine("xxx-8152", 8151.7, vlim=(200.0, 300.0)),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("he-i-8156", 8155.66),
    EmissionLine("crii-8188", 8188.11),
    EmissionLine("ni-8216", 8216.34),
    EmissionLine("ni-8223", 8223.14),
    EmissionLine("oi-8446", 8446.36),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
emlines = [
    EmissionLine("hi-8467", 8467.25),
    EmissionLine("hi-8545", 8545.38),
    EmissionLine("hi-8598", 8598.39),
    EmissionLine("feii-8617", 8616.95),
    EmissionLine("hi-8665", 8665.02),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python
labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-86-94-contsub.fits") 
    for label in labels
}
```

```python
emlines = [
    EmissionLine("CONT-8730", 8730.0, vlim=(-50.0, 50.0)),
    EmissionLine("ci-8727", 8727.13, vlim=(200.0, 300.0)),
    EmissionLine("CONT-8749", 8749.0, vlim=(-50.0, 50.0)),
    EmissionLine("hi-8750", 8750.47),
    EmissionLine("hi-8863", 8862.79),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

There are some misterious very neutral/molecular lines at 9114 and 9145 Å.  The 9114 line is similar to 8152 and seems to come from the highest extinction regions. 

The 9145 line is more diffuse, but also from deeper than the C I line.

The 9219 line may not be Ca I since we don't see the 9259 line that would be expected.  But then I don't know what it is. 

```python
emlines = [
    EmissionLine("hi-9015", 9014.91),
    EmissionLine("siii-9069", 9068.90),
    EmissionLine("cai-9095", 9095.09),
    EmissionLine("CONT-9112", 9112.5, vlim=(-50.0, 50.0)),
    EmissionLine("xxx-9114", 9113.8, vlim=(200.0, 300.0)),
    EmissionLine("clii-9124", 9123.6),
    EmissionLine("CONT-9145", 9145.0, vlim=(-50.0, 50.0)),
    EmissionLine("xxx-9145", 9144.9, vlim=(200.0, 300.0)),
    EmissionLine("cai-9204", 9204.09),
    EmissionLine("hei-9210", 9210.28),
    EmissionLine("cai-9219", 9218.52),
    EmissionLine("hi-9229", 9229.01),
    EmissionLine("cai-9244", 9244.31),
]
for label in "ABCD":
    for em in emlines:
        em.save_moments(csub[label], f"lmc-30dor-{label}")
```

```python

```
