"""Functions to extract emission lines from MUSE cubes

Author: Will Henney, IRyA-UNAM, 2021
"""
import numpy as np
from mpdaf.obj import Cube
from numpy.polynomial import Chebyshev as T
import itertools


def _am_i_special(i, j):
    return i == 0 and j % 10 == 0


def fit_continuum(cube, wav_ranges, deg=1, median=True, verbose=True):
    """Fit the continuum pixel-by-pixel to the MUSE datacube `cube`

    Returns a cube of the same shape as `cube` that continas the
    fitted continuum.

    Fitting is done by a polynomial of degree `deg` (default: 1)

    Wavelength ranges to be treated as (largely) continuum are listed
    in `wav_ranges`, which should be a sequence of one or more pairs
    of (min, max) wavelength.  E.g., [(6445.0, 6690.0), (6710, 6745)]

    If `median=True` (default), then only the median value of each
    (min, max) range is used in the fitting.  This should be (1)
    faster and (2) allows the range to contain some weak lines without
    affecting the fit.  If `median=False`, then all wavelength steps
    within each range are used in the fit.

    I found the MPDAF `Cube.iter_spe()` method to be too slow (this
    iterates over a set of`Spectrum` objects, one for each pixel).
    Instead, I iterate over slices of the `Cube.data` array, which is
    much faster.
    """
    # Extract basic data from the cube
    nv, ny, nx = cube.shape
    wavs = cube.wave.coord()
    # Make a new cube to receive the fitted continuum
    cont_cube = cube.clone(data_init=np.empty, var_init=np.zeros)
    # Loop over all spatial pixels in the cube
    for j, i in itertools.product(range(ny), range(nx)):
        if verbose and _am_i_special(i, j):
            print("extract.fit_continuum: row", j)
        spec = cube.data.data[:, j, i]
        mask = cube.mask[:, j, i]
        x, y = [], []
        # For each wave range, collect the data to fit
        for w1, w2 in wav_ranges:
            # Logical array to select all wave points that are (1)
            # unmasked in the data cube and (2) are within the wave
            # range
            m = ~mask & (wavs >= w1) & (wavs <= w2)
            if median:
                # Either a single point representative of this range
                x.append(np.mean(wavs[m]))
                y.append(np.median(spec[m]))
            else:
                # Or all the unmasked points in the range
                x += wavs[m]
                y += spec[m]
        try:
            p = T.fit(x, y, deg=deg)
            cont_cube.data[:, j, i] = p(wavs)
        except:
            cont_cube.data[:, j, i] = 0.0
            cont_cube.mask[:, j, i] = True
    return cont_cube
