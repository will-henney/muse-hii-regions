"""Utilities for fitting and removing sky lines

Will Henney 2022-04-27 - Initial version is for MUSE spectra
"""

from .linetools import EmissionLine, SpectralRange, VelocityScale
import numpy as np
from astropy.modeling import models, fitting
from mpdaf.obj import Cube  # type: ignore

FITTER = fitting.LevMarLSQFitter()


class Sky:
    """
    An emission line with a sky component in a MUSE cube.  The sky
    component is assumed spatially constant over the field and is
    fitted by a Gaussian
    """

    def __init__(
        self,
        cube: Cube,
        em: EmissionLine,
        full_vlim: tuple[float, float] = (-500.0, 500.0),
        sky_vlim: tuple[float, float] = (-50.0, 50.0),
    ):
        self.em = em
        # The full spectral range is only for visualization purposes
        self.full_range = SpectralRange(self.em.wav0, full_vlim)
        # The sky spectral range is where we expect to find the sky line
        self.sky_range = SpectralRange(self.em.wav0, sky_vlim)
        # Take a sub-cube of the full spectral range
        self.cube = cube.select_lambda(*self.full_range.wavlim)
        # Take median over all spaxels to get the spectrum
        self.spec = self.cube.median(axis=(1, 2))
        # Wavelength corresponding to the spectrum pixels
        self.wavs = self.spec.wave.coord()
        # And the corresponding velocities
        self.vels = VelocityScale(self.em.wav0).wav2vel(self.wavs)
        # Mask to select pixels in the spectrum that are in the sky
        # component spectral range
        self.sky_mask = (self.wavs >= self.sky_range.wavlim[0]) & (
            self.wavs <= self.sky_range.wavlim[1]
        )
        # And we might as well fit sky now since it is cheap
        # Fit a single Gaussian to the sky profile
        self.gauss = FITTER(
            models.Gaussian1D(
                amplitude=self.spec.data[self.sky_mask].max(),
                mean=self.em.wav0,
                stddev=1.0,
            ),
            self.wavs[self.sky_mask],
            self.spec.data[self.sky_mask],
        )


def remove_sky(sky: Sky) -> Cube:
    """
    Subtract the Gaussian model that was fitted to the sky component
    and return the result as a new cube
    """
    skyspec = sky.gauss(sky.wavs)
    newcube = sky.cube.copy()
    newcube.data -= skyspec[:, None, None]
    return newcube
