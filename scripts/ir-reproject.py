filetab=[[5, "IRAC", "~Spitzer/r4384256/ch3/pbcd/SPITZER_I3_4384256_0000_7_E8758299_maic.fits~"], [8, "IRAC", "~Spitzer/r4384256/ch4/pbcd/SPITZER_I4_4384256_0000_7_E8758329_maic.fits~"], [12, "WISE", "~WISE/0145m727_ac51-w3-int-3_ra14.756957499999999_dec-72.17516_asec600.000.fits~"], [24, "MIPS", "~Spitzer/r4384512/ch1/pbcd/SPITZER_M1_4384512_0000_10_E6046561_maic.fits~"], [60, "MIPS", "~Spitzer/r10743808/ch2/pbcd/SPITZER_M2_10743808_0000_10_E6429330_maic.fits~"], [100, "PACS", "~Herschel/science/0001_14.75696000_-72.17516000_SMC.HERITAGE.PACS100.img.fits~"]]
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import reproject

indir = Path("~/Work/Muse-Hii-Data/SMC-NGC-346").expanduser()
outdir = Path("../data")
c0 = SkyCoord.from_name("Cl* NGC 346 W 3")


NY, NX = 200, 200
w0 = WCS(naxis=2)
w0.wcs.crpix = [NX / 2, NY / 2]
w0.wcs.crval = [c0.ra.deg, c0.dec.deg]
w0.wcs.cdelt = np.ones((2,)) / 3600.0
w0.wcs.ctype = ["RA---TAN", "DEC--TAN"]
print(w0.to_header())

for wav, cam, pathstring in filetab:
    pathstring = str(indir / pathstring.strip("~"))
    hdulist = fits.open(pathstring)
    hdu = hdulist[0]
    newdata, footprint = reproject.reproject_interp(hdu, w0, (NY, NX))
    newfile = outdir / f"ngc346-ir-{10*wav:04d}-{cam}.fits"
    fits.PrimaryHDU(data=newdata, header=w0.to_header()).writeto(newfile, overwrite=True)
