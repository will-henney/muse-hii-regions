"""Functions to extract emission lines from MUSE cubes and PV images

Author: Will Henney, IRyA-UNAM, 2021
"""
from collections.abc import Sequence
from typing import Optional, Union
import numpy as np
from mpdaf.obj import Cube  # type: ignore
from numpy.polynomial import Chebyshev as T
import itertools
from astropy.wcs import WCS  # type: ignore
from astropy.io import fits  # type: ignore


def _am_i_special(i: int, j: int) -> bool:
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
    (min, max) range is used in the fitting.  This allows the range to
    contain some weak lines without affecting the fit.  On the other
    hand, it is generally slower.  If `median=False`, then all
    wavelength steps within each range are used in the fit.

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
                x += list(wavs[m])
                y += list(spec[m])
        try:
            p = T.fit(x, y, deg=deg)
            cont_cube.data[:, j, i] = p(wavs)
        except:
            cont_cube.data[:, j, i] = 0.0
            cont_cube.mask[:, j, i] = True
    return cont_cube


def pv_fit_continuum(
    pvim: np.ndarray,
    wcs: WCS,
    wav_ranges: Sequence[tuple[float, float]],
    deg: int = 1,
    median: bool = True,
    verbose: bool = False,
) -> np.ndarray:
    """Fit the continuum pixel-by-pixel to the PV image `pvim`.

    Second argument `wcs` must contain the pixel-to-wavelength
    conversion.  The wavelength axis is assumed to be the
    WCS-convention x pixel axis, which is the SECOND python axis.  The
    WCS-convention y pixel axis of `wcs` is not used.

    Returns an image of the same shape as `pvim` that contains the
    fitted continuum.

    Fitting is done by a polynomial of degree `deg` (default: 1)

    Wavelength ranges to be treated as (largely) continuum are listed
    in `wav_ranges`, which should be a sequence of one or more pairs
    of (min, max) wavelength.  E.g., [(6445.0, 6690.0), (6710, 6745)]

    If `median=True` (default), then only the median value of each
    (min, max) range is used in the fitting.  This allows the range to
    contain some weak lines without affecting the fit.  On the other
    hand, it is generally slower.  If `median=False`, then all
    wavelength steps within each range are used in the fit.
    """
    # Extract basic data from the cube
    ny, nx = pvim.shape
    # Assume that wavelength is along the x axis
    wavs, _ = wcs.pixel_to_world_values(np.arange(nx), [0] * nx)
    # Make a new image to receive the fitted continuum
    cont_pvim = np.empty_like(pvim)
    # Loop over the spatial (y) axis of the image
    for j in range(ny):
        if verbose and _am_i_special(0, j):
            print("extract.pv_fit_continuum: row", j)
        spec = pvim[j, :]
        x, y = [], []
        # For each wave range, collect the data to fit
        #
        # FIXME: This is just a simplified version of the algorithm I
        # used for the MPDAF cube, even though that is not really
        # necessary here since we do not have a mask to account for.
        # There may be a more efficient way of doing it without a
        # loop.
        for w1, w2 in wav_ranges:
            # Logical array to select all wave points that are within
            # the wave range
            m = (wavs >= w1) & (wavs <= w2)
            if median:
                # Either a single point representative of this range
                x.append(np.mean(wavs[m]))
                y.append(np.median(spec[m]))
            else:
                # Or all the unmasked points in the range
                x += list(wavs[m])
                y += list(spec[m])
        try:
            # Fit the polynomial to the selected wavelength points
            p = T.fit(x, y, deg=deg)
            # Fill in corresponding row of the output image
            cont_pvim[j, :] = p(wavs)
        except:
            # If we fail for any reason, then just set zeros
            cont_pvim[j, :] = 0.0
    return cont_pvim


def pvslice(im, w, wavrange, posrange):
    """
    Return the (image, wcs) tuple of a sub-image of the PV image `im`
    with WCS `w` for the wavelength range `wavrange` and the position
    range `posrange`
    """
    ny, nx = im.shape
    if posrange is None:
        ylim = np.array([0, ny])
    else:
        _, ylim = w.world_to_pixel_values([0, 0], posrange)
    if wavrange is None:
        xlim = np.array([0, nx])
    else:
        xlim, _ = w.world_to_pixel_values(wavrange, [0, 0])

    # Force pixel limits to be in bounds
    ylim[0] = max(0, ylim[0])
    ylim[1] = min(ny - 1, ylim[1])
    xlim[0] = max(0, xlim[0])
    xlim[1] = min(nx - 1, xlim[1])

    xslice, yslice = slice(*xlim.astype(int)), slice(*ylim.astype(int))
    return im[yslice, xslice], w.slice((yslice, xslice))


from dataclasses import dataclass
from tetrabloks import rebin_utils


@dataclass
class EmissionLine:
    """
    Class for keeping track of an emission line

    Minimum functionality: name and rest wavelength

    May also include other stuff like spatial profiles
    """

    name: str
    wav0: float
    vlim: tuple[float, float] = (-50.0, 500.0)


@dataclass
class SlitProfile:
    """
    A profile as a function of position along the slit
    """

    position: np.ndarray
    data: np.ndarray

    def multibin(
        self,
        kmax: int = 5,
        mask: Optional[np.ndarray] = None,
        weights: Optional[np.ndarray] = None,
    ) -> dict[int, "SlitProfile"]:
        """
        Calculate a sequence of rebinnings: x 2, x 4, x 8, ..., x 2^`kmax`

        Each rebinning is oversampled back to the original resolution (with padding).
        Returns a dict: {1: profile@1x, 2: profile@2x}
        """
        nlist = 2 ** np.arange(kmax)
        nmax = nlist[-1]
        x = self.position[:]
        y = self.data[:]
        # Check that array length is multiple of maximum n
        remainder = len(x) % nmax
        npad = (nmax - remainder) % nmax
        if npad > 0:
            # If not, then pad the right edge to the correct size
            x = np.pad(x, (0, npad), mode="empty")
            y = np.pad(y, (0, npad), mode="empty")
        assert len(x) % nmax == 0, f"Does not divide: {len(x)}, {nmax}"
        # Save the padded positions at original resolution
        x0 = x[:]

        # Make an all-true mask if none was provided
        if mask is None:
            m = np.ones_like(y).astype(bool)
        else:
            # Make sure mask conforms to the padded array
            m = np.pad(mask, (0, npad), mode="constant", constant_values=False)
        # Save the padded mask at original resolution
        m0 = m[:]
        # Policy is to set values where mask is false to NaNs
        x0[~m0] = np.nan

        # Use constant weights if none were provided
        if weights is None:
            w = np.ones_like(y)
        else:
            w = weights

        rslt = {}
        for n in nlist:
            assert x.shape == y.shape == m.shape == w.shape
            # Upsample the data array back to the original resolution
            y0 = rebin_utils.oversample1d(y, n)
            # And set all values invalid that are outside the original mask
            y0[~m0] = np.nan
            assert y0.shape == x0.shape, (
                f"Different shapes: y0 = {y0.shape}, " f"x0 = {x0.shape} (n = {n})"
            )
            # Save the profile at current binning level
            rslt[n] = SlitProfile(x0, y0)
            # Downsample to get the next coarser binning level
            [x, y], m, w = rebin_utils.downsample1d([x, y], m, weights=w)
        return rslt


@dataclass
class PositionVelocityImage:
    """
    Class to contain data and accompanying WCS of a position-velocity image

    The wavelength axis must be the second python axis
    """

    data: np.ndarray
    wcs: WCS

    @classmethod
    def from_hdu(
        cls, hdu: Union[fits.PrimaryHDU, fits.ImageHDU]
    ) -> "PositionVelocityImage":
        return PositionVelocityImage(hdu.data, WCS(hdu))

    def slit_profile(self, emline: EmissionLine) -> SlitProfile:
        """
        Calculate integrated line flux as function of position
        """
        wavrange = emline.wav0 * (1.0 + np.array(emline.vlim) / 3e5)
        pvwin = self.slice(wavrange, None)
        ny, nx = pvwin.data.shape
        _, positions = pvwin.wcs.pixel_to_world_values(
            [0] * ny,
            np.arange(ny),
        )
        return SlitProfile(positions, np.atleast_1d(pvwin.data.sum(axis=1)))

    def slit_ew_profile(
        self, emline: EmissionLine, continuum: "PositionVelocityImage"
    ) -> SlitProfile:
        """
        Calculate equivalent width (angstrom units) of line as
        function of position
        """
        wavrange = emline.wav0 * (1.0 + np.array(emline.vlim) / 3e5)
        pvwin = self.slice(wavrange, None)
        contwin = continuum.slice(wavrange, None)
        ny, nx = pvwin.data.shape
        _, positions = pvwin.wcs.pixel_to_world_values(
            [0] * ny,
            np.arange(ny),
        )
        dwav = pvwin.wcs.wcs.cdelt[0]
        ew = dwav * (pvwin.data / contwin.data).sum(axis=1)
        return SlitProfile(positions, ew)

    def slice(
        self,
        wavrange: Optional[Sequence[float]],
        posrange: Optional[Sequence[float]],
    ) -> "PositionVelocityImage":
        """
        Return a slice (sub-image) of the position-velocity image for
        the wavelength range `wavrange` and the position range
        `posrange`.  If either range is None, then the full extent
        along that axis is used.
        """
        ny, nx = self.data.shape
        if posrange is None:
            ylim = np.array([0, ny])
        else:
            _, ylim = self.wcs.world_to_pixel_values([0, 0], posrange)
        if wavrange is None:
            xlim = np.array([0, nx])
        else:
            xlim, _ = self.wcs.world_to_pixel_values(wavrange, [0, 0])

        # Force pixel limits to be in bounds. This allows one to use,
        # for example, a very high upper limit on wavrange or posrange
        # to guarantee getting all the way up to the end of that axis
        ylim[0] = max(0, ylim[0])
        ylim[1] = min(ny - 1, ylim[1])
        xlim[0] = max(0, xlim[0])
        xlim[1] = min(nx - 1, xlim[1])

        xslice, yslice = slice(*xlim.astype(int)), slice(*ylim.astype(int))
        return PositionVelocityImage(
            self.data[yslice, xslice],
            self.wcs.slice((yslice, xslice)),
        )
