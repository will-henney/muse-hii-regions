# wcsfile
> Implements `wcsfile.read()` for reading World Coordinate System information from a file.


Read World Cooordinate System parameters from a `.wcs` file, such as those written by SAOImageDS9.

## Install

`pip install wcsfile`

(Currently aspirational only, since it has not been submitted to PyPI yet.)

## How to use

Use `wcsfile.read()` to read a WCS file into a python dict. 

```python
import wcsfile
wcsdict = wcsfile.read("testdata/mosaic-1996-HH204-align-robberto.wcs")
wcsdict
```




    {'DATE-OBS': '21/03/95',
     'EQUINOX': 2000,
     'CTYPE1': 'RA---TAN',
     'CRPIX1': 388.39887,
     'CRVAL1': 83.8433551,
     'CDELT1': 0.0,
     'CTYPE2': 'DEC--TAN',
     'CRPIX2': 567.14165,
     'CRVAL2': -5.4192964,
     'CDELT2': 0.0,
     'CD1_1': -2.76664e-05,
     'CD1_2': 3.804223e-08,
     'CD2_1': 4.913188e-09,
     'CD2_2': 2.764793e-05}



The dict can be used to initialize a FITS header object.

```python
from astropy.io import fits

hdr = fits.Header(wcsdict)
hdr
```




    DATE-OBS= '21/03/95'                                                            
    EQUINOX =                 2000                                                  
    CTYPE1  = 'RA---TAN'                                                            
    CRPIX1  =            388.39887                                                  
    CRVAL1  =           83.8433551                                                  
    CDELT1  =                  0.0                                                  
    CTYPE2  = 'DEC--TAN'                                                            
    CRPIX2  =            567.14165                                                  
    CRVAL2  =           -5.4192964                                                  
    CDELT2  =                  0.0                                                  
    CD1_1   =         -2.76664E-05                                                  
    CD1_2   =         3.804223E-08                                                  
    CD2_1   =         4.913188E-09                                                  
    CD2_2   =         2.764793E-05                                                  



Or it can be used to initialize a WCS object.

```python
from astropy.wcs import WCS

w = WCS(wcsdict)
w
```

    WARNING: FITSFixedWarning: 'datfix' made the change 'Set MJD-OBS to 49797.000000 from DATE-OBS.
    Changed DATE-OBS from '21/03/95' to '1995-03-21''. [astropy.wcs.wcs]





    WCS Keywords
    
    Number of WCS axes: 2
    CTYPE : 'RA---TAN'  'DEC--TAN'  
    CRVAL : 83.8433551  -5.4192964  
    CRPIX : 388.39887  567.14165  
    CD1_1 CD1_2  : -2.76664e-05  3.804223e-08  
    CD2_1 CD2_2  : 4.913188e-09  2.764793e-05  
    NAXIS : 0  0



```python
assert w.wcs.lattyp == "DEC"

```
