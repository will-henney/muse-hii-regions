# wcsfile
> Implements `wcsfile.read()` for reading World Coordinate System information from a file.


Read World Cooordinate System parameters from a `.wcs` file, such as those written by SAOImageDS9.

## Install

`pip install wcsfile`

(Currently aspirational only, since it has not been submitted to PyPI yet.)

## How to use

Use `wcsfile.read()` to read a WCS file into a python dict. 

```python
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



The dict can be used to initialize a WCS object.

```python
from astropy.wcs import WCS

w = WCS(wcsdict)
w.wcs
```




           flag: 137
          naxis: 2
          crpix: 0x7ffbf283a420
                   388.40       567.14    
             pc: 0x7ffbf2854210
        pc[0][]:  -2.7666e-05   3.8042e-08
        pc[1][]:   4.9132e-09   2.7648e-05
          cdelt: 0x7ffbf28682b0
                   1.0000       1.0000    
          crval: 0x7ffbf28682c0
                   83.843      -5.4193    
          cunit: 0x7ffbf2870f00
                 "deg"
                 "deg"
          ctype: 0x7ffbf285a200
                 "RA---TAN"
                 "DEC--TAN"
        lonpole: 180.000000
        latpole: -5.419296
        restfrq: 0.000000
        restwav: 0.000000
            npv: 0
         npvmax: 0
             pv: 0x0
            nps: 0
         npsmax: 0
             ps: 0x0
             cd: 0x7ffbf2874820
        cd[0][]:  -2.7666e-05   3.8042e-08
        cd[1][]:   4.9132e-09   2.7648e-05
          crota: 0x7ffbf28682d0
                   0.0000       0.0000    
         altlin: 2
         velref: 0
            alt: ' '
         colnum: 0
          colax: 0x7ffbf2874840
                     0      0
          cname: 0x7ffbf285a290
                 UNDEFINED
                 UNDEFINED
          crder: 0x7ffbf285a320
                   UNDEFINED    UNDEFINED
          csyer: 0x7ffbf2874850
                   UNDEFINED    UNDEFINED
          czphs: 0x7ffbf2880320
                   UNDEFINED    UNDEFINED
          cperi: 0x7ffbf2880330
                   UNDEFINED    UNDEFINED
        wcsname: UNDEFINED
        timesys: UNDEFINED
        trefpos: UNDEFINED
        trefdir: UNDEFINED
        plephem: UNDEFINED
       timeunit: UNDEFINED
        dateref: UNDEFINED
         mjdref:        UNDEFINED       UNDEFINED
       timeoffs: UNDEFINED
        dateobs: "1995-03-21"
        datebeg: UNDEFINED
        dateavg: UNDEFINED
        dateend: UNDEFINED
         mjdobs:  49797.000000000
         mjdbeg: UNDEFINED
         mjdavg: UNDEFINED
         mjdend: UNDEFINED
         jepoch: UNDEFINED
         bepoch: UNDEFINED
         tstart: UNDEFINED
          tstop: UNDEFINED
        xposure: UNDEFINED
        telapse: UNDEFINED
        timsyer: UNDEFINED
        timrder: UNDEFINED
        timedel: UNDEFINED
       timepixr: UNDEFINED
         obsgeo:        UNDEFINED       UNDEFINED       UNDEFINED
                        UNDEFINED       UNDEFINED       UNDEFINED
       obsorbit: UNDEFINED
        radesys: "FK5"
        equinox:   2000.000000000
        specsys: UNDEFINED
        ssysobs: UNDEFINED
        velosys: UNDEFINED
        zsource: UNDEFINED
        ssyssrc: UNDEFINED
        velangl: UNDEFINED
           ntab: 0
            tab: 0x0
           nwtb: 0
            wtb: 0x0
          types: 0x7ffbf287f5f0
                2200 2201
         lngtyp: "RA"
         lattyp: "DEC"
            lng: 0
            lat: 1
           spec: -1
       cubeface: -1
            err: 0x0
            lin: (see below)
            cel: (see below)
            spc: (see below)
         m_flag: 137
        m_naxis: 2
        m_crpix: 0x7ffbf283a420  (= crpix)
           m_pc: 0x7ffbf2854210  (= pc)
        m_cdelt: 0x7ffbf28682b0  (= cdelt)
        m_crval: 0x7ffbf28682c0  (= crval)
        m_cunit: 0x7ffbf2870f00  (= cunit)
        m_ctype: 0x7ffbf285a200  (= ctype)
           m_pv: 0x0  (= pv)
           m_ps: 0x0  (= ps)
           m_cd: 0x7ffbf2874820  (= cd)
        m_crota: 0x7ffbf28682d0  (= crota)
    
        m_colax: 0x7ffbf2874840  (= colax)
        m_cname: 0x7ffbf285a290  (= cname)
        m_crder: 0x7ffbf285a320  (= crder)
        m_csyer: 0x7ffbf2874850  (= csyer)
        m_czphs: 0x7ffbf2880320  (= czphs)
        m_cperi: 0x7ffbf2880330  (= cperi)
          m_tab: 0x0  (= tab)
          m_wtb: 0x0  (= wtb)
    
       lin.*
           flag: 137
          naxis: 2
          crpix: 0x7ffbf283a420
                   388.40       567.14    
             pc: 0x7ffbf2854210
        pc[0][]:  -2.7666e-05   3.8042e-08
        pc[1][]:   4.9132e-09   2.7648e-05
          cdelt: 0x7ffbf28682b0
                   1.0000       1.0000    
         dispre: 0x0
         disseq: 0x0
    piximg[0][]:  -2.7666e-05   3.8042e-08
    piximg[1][]:   4.9132e-09   2.7648e-05
    imgpix[0][]:  -36145.       49.734    
    imgpix[1][]:   6.4231       36169.    
        i_naxis: 2
          unity: 0
         affine: 1
         simple: 0
            err: 0x0
         tmpcrd: 0x7ffbf28315b0
         m_flag: 0
        m_naxis: 0
        m_crpix: 0x0
           m_pc: 0x0
        m_cdelt: 0x0
       m_dispre: 0x0
       m_disseq: 0x0
    
       cel.*
          flag: 137
         offset: 0
           phi0:  0.000000
         theta0: 90.000000
            ref:   83.843      -5.4193       180.00      -5.4193    
            prj: (see below)
          euler:   83.843       95.419       180.00      -0.094444     0.99553   
        latpreq: 0 (not required)
         isolat: 0
            err: 0x0
    
       prj.*
           flag: 103
           code: "TAN"
             r0: 57.295780
             pv: (not used)
           phi0:  0.000000
         theta0: 90.000000
         bounds: 7
    
           name: "gnomonic"
       category: 1 (zenithal)
        pvrange: 0
      simplezen: 1
      equiareal: 0
      conformal: 0
         global: 0
      divergent: 1
             x0: 0.000000
             y0: 0.000000
            err: 0x0
            w[]:   0.0000       0.0000       0.0000       0.0000       0.0000    
                   0.0000       0.0000       0.0000       0.0000       0.0000    
              m: 0
              n: 0
         prjx2s: 0x0118b74f20
         prjs2x: 0x0118b75380
    
       spc.*
           flag: 0
           type: "    "
           code: "   "
          crval: UNDEFINED
        restfrq: 0.000000
        restwav: 0.000000
             pv: (not used)
              w:   0.0000       0.0000       0.0000      (remainder unused)
        isGrism: 0
            err: 0x0
         spxX2P: 0x0
         spxP2S: 0x0
         spxS2P: 0x0
         spxP2X: 0x0



Or it can be used to initialize a FITS header object.

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


