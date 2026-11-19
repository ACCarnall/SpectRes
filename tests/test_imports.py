import numpy as np
import pytest

import spectres
from spectres.spectral_resampling import spectres as python_spectres


def test_default_resampler_stays_uncompiled():
    assert spectres.spectres is python_spectres


def test_optional_numba_resampler_is_exported():
    pytest.importorskip("numba")
    wavs = np.arange(1., 6.)
    fluxes = wavs + 0.25
    actual = spectres.spectres_numba(wavs[1:-1], wavs, fluxes)
    np.testing.assert_allclose(actual, fluxes[1:-1])
