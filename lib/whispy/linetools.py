"""Classes and functions to represent emission lines
"""
from dataclasses import dataclass
import astropy.units as u  # type: ignore
import astropy.constants as const  # type: ignore

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value


class VelocityScale:
    """
    A spectral velocity scale in km/s around a reference wavelength
    `wav0` in Angstrom that is defined with the optical VOPT
    convention.
    """

    def __init__(self, wav0: float):
        self.wav0 = wav0

    def wav2vel(self, wav):
        """Convert wavelength in Angstrom to velocity in km/s"""
        return LIGHT_SPEED_KMS * (wav - self.wav0) / self.wav0

    def vel2wav(self, vel):
        """Convert velocity in km/s to wavelength in Angstrom"""
        return self.wav0 * (1.0 + vel / LIGHT_SPEED_KMS)


@dataclass
class SpectralRange:
    """
    A spectral range around a reference wavelength `wav0` in Angstrom
    that is given as a velocity interval `vlim` in km/s (tuple of min, max)
    """

    wav0: float
    vlim: tuple[float, float] = (-50.0, 500.0)

    @property
    def wavlim(self):
        """Wavelength limits, synchronised with `vlim`"""
        return tuple(map(VelocityScale(self.wav0).vel2wav, self.vlim))

    @wavlim.setter
    def wavlim(self, value: tuple[float, float]):
        self.vlim = tuple(map(VelocityScale(self.wav0).wav2vel, value))


@dataclass
class EmissionLine:
    """
    Class for keeping track of an emission line

    Minimum functionality: name and rest wavelength

    May also include other stuff like spatial profiles
    """

    name: str
    wav0: float
