filetab=[[3.6, "IRAC1", 1.0, "~Spitzer/r4384256/ch1/pbcd/SPITZER_I1_4384256_0000_7_E8758509_maic.fits~"], [4.5, "IRAC2", 1.0, "~Spitzer/r4384256/ch2/pbcd/SPITZER_I2_4384256_0000_7_E8758310_maic.fits~"], [5.731, "IRAC3", 1.0, "~Spitzer/r4384256/ch3/pbcd/SPITZER_I3_4384256_0000_7_E8758299_maic.fits~"], [8.0, "IRAC4", 1.0, "~Spitzer/r4384256/ch4/pbcd/SPITZER_I4_4384256_0000_7_E8758329_maic.fits~"], [8.276, "MSX-A", 7133000.0, "~MSX/SMCA.FIT~"], [12.082, "WISE3", 0.04123, "~WISE/0145m727_ac51-w3-int-3_ra14.756957499999999_dec-72.17516_asec600.000.fits~"], [22.194, "WISE4", 1.176, "~WISE/0145m727_ac51-w4-int-3_ra14.756957499999999_dec-72.17516_asec600.000.fits~"], [12.126, "MSX-C", 28630000.0, "~MSX/SMCC.FIT~"], [14.649, "MSX-D", 32160000.0, "~MSX/SMCD.FIT~"], [21.411, "MSX-E", 24760000.0, "~MSX/SMCE.FIT~"], [23.68, "MIPS1", 1.0, "~Spitzer/r4384512/ch1/pbcd/SPITZER_M1_4384512_0000_10_E6046561_maic.fits~"], [71.42, "MIPS2", 1.0, "~Spitzer/r10743808/ch2/pbcd/SPITZER_M2_10743808_0000_10_E6429330_maic.fits~"], [100, "PACS-B", 1.0, "~Herschel/science/0001_14.75696000_-72.17516000_SMC.HERITAGE.PACS100.img.fits~"], [160, "PACS-R", 1.0, "~Herschel/science/0001_14.75696000_-72.17516000_SMC.HERITAGE.PACS160.img.fits~"], [250, "SPIRE250", 1.0, "~Herschel/science/0001_14.75696000_-72.17516000_SMC.HERITAGE.SPIRE250.img.fits~"], [350, "SPIRE350", 1.0, "~Herschel/science/0001_14.75696000_-72.17516000_SMC.HERITAGE.SPIRE350.img.fits~"], [500, "SPIRE500", 1.0, "~Herschel/science/0001_14.75696000_-72.17516000_SMC.HERITAGE.SPIRE500.img.fits~"]]
import sys
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import reproject
sys.path.append("../lib")
from wcsfile import wcsfile


indir = Path("~/Work/Muse-Hii-Data/SMC-NGC-346").expanduser()
outdir = Path("../data")
c0 = SkyCoord.from_name("Cl* NGC 346 W 3")

wmsx = WCS(wcsfile.read(outdir / "ngc346-msx-correct2.wcs"))

NY, NX = 5 * 300, 5 * 300
w0 = WCS(naxis=2)
w0.wcs.crpix = [NX / 2, NY / 2]
w0.wcs.crval = [c0.ra.deg, c0.dec.deg]
w0.wcs.cdelt = np.array([-0.2, 0.2]) / 3600.0
w0.wcs.ctype = ["RA---TAN", "DEC--TAN"]

for wav, cam, norm, pathstring in filetab:
    print(wav, cam)
    pathstring = str(indir / pathstring.strip("~"))
    hdulist = fits.open(pathstring)
    hdu = hdulist[0]
    if "MSX" in cam.upper():
        # Remove all trace of 3rd dimension
        del hdu.header["*3"]
        hdu.header["NAXIS"] = 2
        hdu.data = hdu.data[0, :, :]
        # Small shift to alignment
        hdu.header.update(wmsx.to_header())
    print(WCS(hdu.header))
    newdata, footprint = reproject.reproject_interp(
        hdu,
        w0,
        (NY, NX),
        order="nearest-neighbor",
    )
    newfile = outdir / f"ngc346-ir-{int(10*wav):04d}-{cam}.fits"
    fits.PrimaryHDU(
        data=newdata * norm,
        header=w0.to_header()
    ).writeto(newfile, overwrite=True)
