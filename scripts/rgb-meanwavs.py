TABLE=[["ci 8727.13", "ci", 8732.24, -3, 35], ["", "ci", 8731.19, -3, 30], ["", "ci", 8729.94, -2, 15], ["cliv 8045.62", "cliv", 8051.19, -2, 10], ["", "cliv", 8049.94, -2, 20], ["", "cliv", 8048.69, -2, 10], ["oiii 5006.84", "oiii", 5011.22, -8500, 25000], ["", "oiii", 5009.97, -17500, 65000], ["", "oiii", 5008.72, -14000, 60000], ["xxx 8037.0", "xxx", 8043.69, -1, 10], ["", "xxx", 8042.44, -1, 12], ["", "xxx", 8041.19, -1, 10], ["xxx 8151.3424", "xxx", 8156.19, -1, 25], ["", "xxx", 8154.94, -2, 35], ["", "xxx", 8153.69, -1, 25], ["oi 6363.78", "oi", 6368.69, -13, 60], ["", "oi", 6367.44, -25, 150], ["", "oi", 6366.19, -18, 150], ["ariv 4740.17", "ariv", 4743.72, -10, 70], ["", "ariv", 4742.47, -15, 90], ["", "ariv", 4741.22, -10, 40], ["oi 8446.36", "oi", 8452.44, 0, 30], ["", "oi", 8451.19, 0, 40], ["", "oi", 8449.94, 0, 30], ["oiii 4958.91", "oiii", 4962.47, -5200, 15000], ["", "oiii", 4961.22, -5500, 20000], ["", "oiii", 4959.97, -2600, 8000], ["siii 9068.90", "siii", 9074.94, -320, 2000], ["", "siii", 9073.69, -475, 3500], ["", "siii", 9072.44, -250, 2000]]
from pathlib import Path
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.constants import c as light_speed
import json

DATADIR = Path.cwd().parent / "data"
data = {}
for maybe_lineid, ion, wav, fmin, fmax in TABLE:
     # Reorganize input data
     wavstring = f"{wav:.2f}".replace(".", "_")
     filename = f"ngc346-slice-{ion}-{wavstring}"
     if maybe_lineid:
          lineid = maybe_lineid
          wav0 = float(lineid.split()[-1])
          prefix = f"{ion}-{round(wav0)}"
          data[lineid] = {
               "ion": ion,
               "prefix": prefix,
               "wavs": [wav],
               "wav0": wav0,
               "fmins": [fmin],
               "fmaxs": [fmax],
               "filenames": [filename],
          }
     else:
          data[lineid]["wavs"].append(wav)
          data[lineid]["fmins"].append(fmin)
          data[lineid]["fmaxs"].append(fmax)
          data[lineid]["filenames"].append(filename)

with open("rgb-channels-data.json", "w") as f:
     json.dump(data, f, indent=3)

for lineid, d in data.items():
     # Get list of HDUs, one for each channel
     hdus = [
          fits.open(DATADIR / (fname + ".fits"))[0]
          for fname in d["filenames"]
     ]
     # Stack the three data arrays in a cube
     imstack = np.stack([hdu.data for hdu in hdus], axis=0)
     # Make matching stacks of the min fluxes and wavelengths
     fmins = np.array(d["fmins"]).reshape((3, 1, 1))
     wavs = np.array(d["wavs"]).reshape((3, 1, 1))
     wav0 = d["wav0"]
     # Correct zeropoint in each channel
     imstack -= fmins
     # Calculate maps of total flux and mean wavelength
     imsum = np.sum(imstack, axis=0)
     imwav = np.sum(imstack * wavs, axis=0) / imsum
     imvel = (((imwav - wav0) / wav0) * light_speed).to(u.km / u.s).value
     bfrac = imstack[2] / imsum
     rfrac = imstack[0] / imsum
     imm1 = rfrac - bfrac
     imm2 = bfrac + rfrac
     # Write out the results
     header = hdus[0].header
     prefix = "ngc346-rgbchan-" + d["prefix"] 
     fits.PrimaryHDU(header=header, data=imsum).writeto(
          DATADIR / f"{prefix}-sum.fits", overwrite=True,
     )
     fits.PrimaryHDU(header=header, data=imwav).writeto(
          DATADIR / f"{prefix}-mean-wav.fits", overwrite=True,
     )
     fits.PrimaryHDU(header=header, data=imvel).writeto(
          DATADIR / f"{prefix}-mean-vel.fits", overwrite=True,
     )
     fits.PrimaryHDU(header=header, data=bfrac).writeto(
          DATADIR / f"{prefix}-3wav-b.fits", overwrite=True,
     )
     fits.PrimaryHDU(header=header, data=rfrac).writeto(
          DATADIR / f"{prefix}-3wav-r.fits", overwrite=True,
     )
     fits.PrimaryHDU(header=header, data=imm1).writeto(
          DATADIR / f"{prefix}-3wav-m1.fits", overwrite=True,
     )
     fits.PrimaryHDU(header=header, data=imm2).writeto(
          DATADIR / f"{prefix}-3wav-m2.fits", overwrite=True,
     )
