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

I am having trouble with the `rc.setCorr` method, so I want to try it out in a simpler context. 

```python
import pyneb as pn
```

```python
rc = pn.RedCorr(E_BV = 1.2, R_V = 3.2, law = 'F99')
# SMC parameters from Gordon et al 2003
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
```

```python
rc.setCorr(obs_over_theo=6.5 / 2.86, wave1=6563., wave2=4861.)
```

```python
rc.E_BV
```

So that works, but if I do exactly the same with `F99-like` it doesn't.

```python
rc2 = pn.RedCorr(E_BV = 1.2, R_V = 3.2, law = 'F99-like')
# SMC parameters from Gordon et al 2003
rc2.R_V = 2.74
rc2.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
```

On the other hand, `F99` seems to work perfectly well.

```python

```
