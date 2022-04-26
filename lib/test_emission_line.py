import pytest
import numpy as np
from extract import EmissionLine, LIGHT_SPEED_KMS


@pytest.fixture()
def emline():
    """Example emision line for tests"""
    yield EmissionLine("[O I] 6300", 6300.30)
    # No clean-up required


def test_wav2vel2wav(emline):
    wav = 6303.0
    wav_roundtrip = emline.vel2wav(emline.wav2vel(wav))
    assert np.isclose(wav, wav_roundtrip)


def test_vel2wav2vel(emline):
    """Test round-tripping to a ridiculously high tolerance"""
    vel = -1000 / 3
    vel_roundtrip = emline.wav2vel(emline.vel2wav(vel))
    assert np.isclose(vel, vel_roundtrip, rtol=1.0e-12, atol=0.0)


def test_wavlim_read(emline):
    wavlim_expected = 6299.24922, 6310.80776
    assert np.allclose(emline.wavlim, wavlim_expected)


def test_wavlim_write(emline):
    emline.wavlim = 6299.24922, 6310.80776
    vlim_expected = -50.0, 500.0
    assert np.allclose(emline.vlim, vlim_expected)


def test_wavlim_linear(emline):
    """Check that adding a large offset to wavlim preserves velocity delta"""
    dv_orig = emline.vlim[1] - emline.vlim[0]
    wavoffset = 200.0
    # Add the same offset to lower and upper wav limits
    emline.wavlim = emline.wavlim[0] + wavoffset, emline.wavlim[1] + wavoffset
    dv_new = emline.vlim[1] - emline.vlim[0]
    assert np.isclose(dv_new, dv_orig)


def test_vlim_scaling(emline):
    """Check that vlim changes by the expected amount when adding an offset to wavlim"""
    vmax_orig = emline.vlim[1]
    wavoffset = 200.0
    emline.wavlim = emline.wavlim[0], emline.wavlim[1] + wavoffset
    vmax_new = emline.vlim[1]
    expected_voffset = LIGHT_SPEED_KMS * wavoffset / emline.wav0
    assert np.isclose(vmax_new, vmax_orig + expected_voffset)


def test_wavlim_scaling(emline):
    """Check that wavlim changes by the expected amount when adding an offset to vlim"""
    wavmin_orig = emline.wavlim[0]
    voffset = -1.0e4
    emline.vlim = emline.vlim[0] + voffset, emline.vlim[1]
    wavmin_new = emline.wavlim[0]
    expected_wavoffset = emline.wav0 * voffset / LIGHT_SPEED_KMS
    assert np.isclose(wavmin_new, wavmin_orig + expected_wavoffset)
