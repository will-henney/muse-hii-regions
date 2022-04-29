# cloudytab
Ultra-lightweight python wrapper for Cloudy output files

## Example of usage

```python
from cloudytab import CloudyModel
m = CloudyModel("myfolder/mymodel")
```

`m.files` contains a list of all the files that were found: `['myfolder/mymodel.in', 'myfolder/mymodel.ovr', ...]`

`m.data` contains dict of `astropy.Table`, one for each save file: `{'ovr': <Table length=289> ..., ...}`

`m.io['in']` and `m.io['out']` contain the input and output streams

`m.skipped` contains a dict of each extension that was skipped with its reason.
