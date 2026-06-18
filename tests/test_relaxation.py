"""Tests for epyr.relaxation: T1/T2 decay and recovery models and fitting."""

import numpy as np
import pytest

from epyr.relaxation import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)


@pytest.mark.smoke
class TestModelFunctions:
    def test_mono_exponential_at_zero(self):
        v = mono_exponential(0.0, amplitude=5.0, T=10.0, offset=1.0)
        assert v == pytest.approx(6.0)

    def test_mono_exponential_decays_to_offset(self):
        t = np.array([0.0, 1e6])
        y = mono_exponential(t, amplitude=5.0, T=10.0, offset=1.0)
        assert y[-1] == pytest.approx(1.0, abs=1e-6)

    def test_stretched_exponential_reduces_to_mono_at_beta_one(self):
        t = np.linspace(0, 50, 25)
        y_stretched = stretched_exponential(
            t, amplitude=3.0, T=8.0, beta=1.0, offset=0.5
        )
        y_mono = mono_exponential(t, amplitude=3.0, T=8.0, offset=0.5)
        np.testing.assert_allclose(y_stretched, y_mono)

    def test_biexponential_at_zero(self):
        v = biexponential(
            0.0, amplitude1=2.0, tau1=5.0, amplitude2=3.0, tau2=50.0, offset=1.0
        )
        assert v == pytest.approx(6.0)

    def test_inversion_recovery_at_zero_is_negative_amplitude(self):
        v = inversion_recovery(0.0, amplitude=4.0, T1=20.0, offset=0.5)
        assert v == pytest.approx(0.5 - 4.0)

    def test_inversion_recovery_long_time_recovers_amplitude(self):
        v = inversion_recovery(np.array([1e6]), amplitude=4.0, T1=20.0, offset=0.5)
        assert v[0] == pytest.approx(4.5, abs=1e-6)

    def test_saturation_recovery_at_zero_is_offset(self):
        v = saturation_recovery(0.0, amplitude=4.0, T1=20.0, offset=0.5)
        assert v == pytest.approx(0.5)

    def test_gamma_gaussian_decay_at_zero(self):
        v = gamma_gaussian_decay(0.0, amplitude=2.0, Gamma0=0.1, GammaG=0.05, offset=0.2)
        assert v == pytest.approx(2.2)

    def test_gamma_gaussian_decay_array_shape_and_dtype(self):
        t = np.linspace(0, 100, 50)
        y = gamma_gaussian_decay(t, amplitude=1.0, Gamma0=0.02, GammaG=0.01)
        assert y.shape == t.shape
        assert y.dtype == np.float64
