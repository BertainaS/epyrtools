"""Tests for epyr.relaxation: T1/T2 decay and recovery models and fitting."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from epyr.relaxation import (
    RelaxationFitResult,
    biexponential,
    fit_relaxation,
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
        v = gamma_gaussian_decay(
            0.0, amplitude=2.0, Gamma0=0.1, GammaG=0.05, offset=0.2
        )
        assert v == pytest.approx(2.2)

    def test_gamma_gaussian_decay_array_shape_and_dtype(self):
        t = np.linspace(0, 100, 50)
        y = gamma_gaussian_decay(t, amplitude=1.0, Gamma0=0.02, GammaG=0.01)
        assert y.shape == t.shape
        assert y.dtype == np.float64


@pytest.mark.smoke
class TestFitRelaxationMonoExponential:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=15.0, offset=1.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.parameters["amplitude"] == pytest.approx(5.0, rel=1e-3)
        assert result.parameters["T"] == pytest.approx(15.0, rel=1e-3)
        assert result.parameters["offset"] == pytest.approx(1.0, abs=1e-3)
        assert result.r_squared > 0.999

    def test_returns_relaxation_fit_result(self):
        t = np.linspace(0, 50, 100)
        y = mono_exponential(t, amplitude=2.0, T=10.0, offset=0.0)
        result = fit_relaxation(t, y, plot=False)
        assert isinstance(result, RelaxationFitResult)
        assert result.model == "mono_exponential"

    def test_plot_true_does_not_raise(self):
        t = np.linspace(0, 50, 80)
        y = mono_exponential(t, amplitude=3.0, T=12.0, offset=0.5)
        result = fit_relaxation(t, y, plot=True)
        assert result.success
        plt.close("all")


@pytest.mark.standard
class TestFitRelaxationNoisy:
    def test_noisy_mono_exponential(self):
        rng = np.random.default_rng(0)
        t = np.linspace(0, 100, 300)
        y_true = mono_exponential(t, amplitude=5.0, T=20.0, offset=1.0)
        y = y_true + 0.02 * 5.0 * rng.standard_normal(t.size)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.r_squared > 0.99
        assert result.parameters["T"] == pytest.approx(20.0, rel=0.05)


@pytest.mark.smoke
class TestFitRelaxationValidation:
    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError):
            fit_relaxation(np.array([1, 2, 3, 4]), np.array([1, 2, 3]), plot=False)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError):
            fit_relaxation(np.array([1, 2, 3]), np.array([1, 2, 3]), plot=False)

    def test_unsupported_model_raises(self):
        t = np.linspace(0, 10, 20)
        y = np.ones_like(t)
        with pytest.raises(ValueError):
            fit_relaxation(t, y, model="not_a_model", plot=False)

    def test_mask_excludes_points(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=15.0, offset=1.0)
        mask = np.ones_like(t, dtype=bool)
        mask[50:60] = False
        result = fit_relaxation(t, y, mask=mask, plot=False)
        assert result.success
        assert result.t_fit.size == 190
        assert result.parameters["T"] == pytest.approx(15.0, rel=1e-3)


@pytest.mark.smoke
class TestFitRelaxationStretchedExponential:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0.1, 100, 200)
        y = stretched_exponential(t, amplitude=4.0, T=25.0, beta=0.7, offset=0.5)
        result = fit_relaxation(t, y, model="stretched_exponential", plot=False)
        assert result.success
        assert result.parameters["T"] == pytest.approx(25.0, rel=1e-3)
        assert result.parameters["beta"] == pytest.approx(0.7, rel=1e-3)
        assert result.r_squared > 0.999


@pytest.mark.smoke
class TestFitRelaxationGammaGaussianDecay:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = gamma_gaussian_decay(t, amplitude=3.0, Gamma0=0.05, GammaG=0.02, offset=0.2)
        result = fit_relaxation(t, y, model="gamma_gaussian_decay", plot=False)
        assert result.success
        assert result.parameters["Gamma0"] == pytest.approx(0.05, rel=1e-3)
        assert result.parameters["GammaG"] == pytest.approx(0.02, rel=1e-3)
        assert result.r_squared > 0.999
