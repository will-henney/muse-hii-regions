from __future__ import print_function
import sys
import os
from distutils.dep_util import newer, newer_group
import numpy as np
from rebin_utils import downsample, oversample
from astropy.io import fits

nlist = [1, 2, 4, 8, 16, 32, 64, 128, 256]
mingoods = [2, 2, 2, 2, 2, 2, 2, 2, 2]

def pad_array(a, n):
    """Pad 2d array `a` to nearest multiple of `n` in each dimension"""
    newshape = (n*np.ceil(np.array(a.shape).astype(float)/n)).astype(int)
    b = np.zeros(newshape, dtype=a.dtype)
    b[:a.shape[0], :a.shape[1]] = a
    return b


try: 
    infile = sys.argv[1]
except:
    sys.exit('Usage: {} FITSFILE'.format(sys.argv[0]))


hdu = fits.open(infile)[0]
if hdu.data is None:
    hdu = fits.open(infile)[1]
hdr = hdu.header
# Maximum binning
nmax = nlist[-1]

# Pad arrays to nearest multiple of nmax
im = pad_array(hdu.data, nmax)

basename = os.path.basename(infile)
map_type = basename.split('-')[0]

if map_type in ['mean', 'sigma']:
    # For the mean velocity and sigma width maps, weight by brightness
    wfile = infile.replace(map_type, 'linesum')
    if infile.endswith('-patfixx.fits'):
        # Strip off the pattern fix prefix if present
        wfile = wfile.replace('-patfixx', '')
    whdu = fits.open(wfile)[0]
    if whdu.data is None:
        # try second HDU if first has no data
        whdu = fits.open(wfile)[1]
    w = pad_array(whdu.data, nmax)
else:
    # Otherwise, just natural weighting
    w = np.ones_like(im)

continuum = fits.open('muse-hr-image-wfc3-f547m.fits')['DATA'].data
starmask = continuum > 30
m =  np.isfinite(hdu.data) & (~starmask)
m = pad_array(m, nmax)

for n, mingood in zip(nlist, mingoods):
    im[~m] = 0.0
    outfile = infile.replace('.fits', '-bin{:03d}.fits'.format(n))
    if n == nlist[0]:
        # Do dependency checking on the first iteration
        if not newer(infile, outfile):
            # Bail out if dependency not newer than target
            sys.exit(outfile + ' is already up to date.')
    print('Saving', outfile)
    # Save both the scaled image and the weights, but at the full resolution
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=oversample(im, n), header=hdr, name='scaled'),
        fits.ImageHDU(data=oversample(w, n), header=hdr, name='weight'),
    ]).writeto(outfile, clobber=True)
    # Now do the rebinning by a factor of two
    [im,], m, w = downsample([im,], m, weights=w, mingood=mingood)
