import pytest
import numpy as np
from .linetools import SpectralRange, EmissionLine, VelocityScale, LIGHT_SPEED_KMS


def test_light_speed():
    """
    I thought this test was sure to fail with an assertion of
    equality, but amazingly it passes!  This means that Emacs's
    constants.el and astropy.constants have the exact same value of c
    to machine precision.  I am leaving the test in, just to see if
    speed of light ever changes ...
    """
    light_speed_emacs = 2.99792458e10 / 1.0e5
    assert LIGHT_SPEED_KMS == light_speed_emacs
    # assert np.isclose(LIGHT_SPEED_KMS, light_speed_emacs, rtol=0.0, atol=0.0)


@pytest.fixture()
def vscale():
    """Example velocity scale for tests"""
    yield VelocityScale(6300.30)
    # No clean-up required


@pytest.fixture()
def specrange():
    """Example spectral range for tests"""
    yield SpectralRange(6300.30, (-50.0, 500.0))
    # No clean-up required


def test_wav2vel2wav(vscale):
    wav = 6303.0
    wav_roundtrip = vscale.vel2wav(vscale.wav2vel(wav))
    assert np.isclose(wav, wav_roundtrip)


def test_vel2wav2vel(vscale):
    """Test round-tripping to a ridiculously high tolerance"""
    vel = -1000 / 3
    vel_roundtrip = vscale.wav2vel(vscale.vel2wav(vel))
    assert np.isclose(vel, vel_roundtrip, rtol=1.0e-12, atol=0.0)


def test_wavlim_read(specrange):
    wavlim_expected = 6299.24922, 6310.80776
    assert np.allclose(specrange.wavlim, wavlim_expected)


def test_wavlim_write(specrange):
    specrange.wavlim = 6299.24922, 6310.80776
    vlim_expected = -50.0, 500.0
    assert np.allclose(specrange.vlim, vlim_expected)


def test_wavlim_linear(specrange):
    """Check that adding a large offset to wavlim preserves velocity delta"""
    dv_orig = specrange.vlim[1] - specrange.vlim[0]
    wavoffset = 200.0
    # Add the same offset to lower and upper wav limits
    specrange.wavlim = specrange.wavlim[0] + wavoffset, specrange.wavlim[1] + wavoffset
    dv_new = specrange.vlim[1] - specrange.vlim[0]
    assert np.isclose(dv_new, dv_orig)


def test_vlim_scaling(specrange):
    """Check that vlim changes by the expected amount when adding an offset to wavlim"""
    vmax_orig = specrange.vlim[1]
    wavoffset = 200.0
    specrange.wavlim = specrange.wavlim[0], specrange.wavlim[1] + wavoffset
    vmax_new = specrange.vlim[1]
    expected_voffset = LIGHT_SPEED_KMS * wavoffset / specrange.wav0
    assert np.isclose(vmax_new, vmax_orig + expected_voffset)


def test_wavlim_scaling(specrange):
    """Check that wavlim changes by the expected amount when adding an offset to vlim"""
    wavmin_orig = specrange.wavlim[0]
    voffset = -1.0e4
    specrange.vlim = specrange.vlim[0] + voffset, specrange.vlim[1]
    wavmin_new = specrange.wavlim[0]
    expected_wavoffset = specrange.wav0 * voffset / LIGHT_SPEED_KMS
    assert np.isclose(wavmin_new, wavmin_orig + expected_wavoffset)
