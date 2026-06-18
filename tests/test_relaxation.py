"""Tests for epyr.relaxation: T1/T2 decay and recovery models and fitting."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from epyr.relaxation import (
    RelaxationFitResult,
    biexponential,
    fit_multiple_decays,
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


@pytest.mark.smoke
class TestFitRelaxationBiexponential:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 200, 300)
        y = biexponential(
            t, amplitude1=3.0, tau1=8.0, amplitude2=2.0, tau2=80.0, offset=0.3
        )
        result = fit_relaxation(t, y, model="biexponential", plot=False)
        assert result.success
        assert result.parameters["tau1"] == pytest.approx(8.0, rel=1e-2)
        assert result.parameters["tau2"] == pytest.approx(80.0, rel=1e-2)
        assert result.r_squared > 0.999


@pytest.mark.smoke
class TestFitRelaxationInversionRecovery:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = inversion_recovery(t, amplitude=4.0, T1=20.0, offset=0.5)
        result = fit_relaxation(t, y, model="inversion_recovery", plot=False)
        assert result.success
        assert result.parameters["T1"] == pytest.approx(20.0, rel=1e-3)
        assert result.parameters["amplitude"] == pytest.approx(4.0, rel=1e-3)
        assert result.r_squared > 0.999


@pytest.mark.smoke
class TestFitRelaxationSaturationRecovery:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = saturation_recovery(t, amplitude=4.0, T1=20.0, offset=0.5)
        result = fit_relaxation(t, y, model="saturation_recovery", plot=False)
        assert result.success
        assert result.parameters["T1"] == pytest.approx(20.0, rel=1e-3)
        assert result.parameters["amplitude"] == pytest.approx(4.0, rel=1e-3)
        assert result.r_squared > 0.999


@pytest.mark.standard
class TestFitMultipleDecays:
    def test_picks_mono_exponential_when_data_is_mono_exponential(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0, 100, 300)
        y_true = mono_exponential(t, amplitude=5.0, T=20.0, offset=1.0)
        y = y_true + 0.01 * (y_true.max() - y_true.min()) * rng.standard_normal(t.size)

        results = fit_multiple_decays(t, y, plot=False)
        best = min(results, key=lambda k: results[k].chi_squared)
        assert best == "mono_exponential"

    def test_picks_stretched_exponential_when_data_is_stretched(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0.1, 100, 300)
        y_true = stretched_exponential(t, amplitude=4.0, T=25.0, beta=0.6, offset=0.5)
        y = y_true + 0.01 * (y_true.max() - y_true.min()) * rng.standard_normal(t.size)

        results = fit_multiple_decays(t, y, plot=False)
        best = min(results, key=lambda k: results[k].chi_squared)
        assert best == "stretched_exponential"

    def test_picks_biexponential_when_data_is_biexponential(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0, 200, 300)
        y_true = biexponential(
            t, amplitude1=3.0, tau1=8.0, amplitude2=2.0, tau2=80.0, offset=0.3
        )
        y = y_true + 0.01 * (y_true.max() - y_true.min()) * rng.standard_normal(t.size)

        results = fit_multiple_decays(t, y, plot=False)
        best = min(results, key=lambda k: results[k].chi_squared)
        assert best == "biexponential"

    def test_plot_true_does_not_raise(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=20.0, offset=1.0)
        results = fit_multiple_decays(t, y, plot=True)
        assert all(r.success for r in results.values())
        plt.close("all")


@pytest.mark.deep
class TestFitRelaxationEdgeCases:
    def test_decay_time_much_shorter_than_time_range(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=0.5, offset=1.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.parameters["T"] == pytest.approx(0.5, rel=1e-2)

    def test_decay_time_much_longer_than_time_range(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=10000.0, offset=1.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        # The decay is nearly linear over this range, so individual
        # parameters are not well constrained; only require convergence
        # and a good fit, not exact parameter recovery.
        assert result.success
        assert result.r_squared > 0.999

    def test_beta_near_lower_bound(self):
        t = np.linspace(0.1, 100, 200)
        y = stretched_exponential(t, amplitude=4.0, T=25.0, beta=0.1, offset=0.5)
        result = fit_relaxation(t, y, model="stretched_exponential", plot=False)
        assert result.success
        assert result.parameters["beta"] == pytest.approx(0.1, rel=1e-2)

    def test_beta_above_one_compressed_exponential(self):
        t = np.linspace(0.1, 100, 200)
        y = stretched_exponential(t, amplitude=4.0, T=25.0, beta=2.0, offset=0.5)
        result = fit_relaxation(t, y, model="stretched_exponential", plot=False)
        assert result.success
        assert result.parameters["beta"] == pytest.approx(2.0, rel=1e-2)

    def test_negative_offset(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=15.0, offset=-3.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.parameters["offset"] == pytest.approx(-3.0, rel=1e-2)


@pytest.mark.scientific
class TestFitRelaxationAnalyticCrossCheck:
    def test_mono_exponential_matches_log_linear_regression(self):
        rng = np.random.default_rng(2)
        t = np.linspace(0, 100, 200)
        T_true = 18.0
        amplitude_true = 6.0
        offset_true = 0.8
        y_true = mono_exponential(
            t, amplitude=amplitude_true, T=T_true, offset=offset_true
        )
        y = y_true + 0.01 * amplitude_true * rng.standard_normal(t.size)

        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success

        # Closed-form solution: with the true offset known, log(|y - offset|)
        # is linear in t with slope -1 / T.
        slope, _ = np.polyfit(t, np.log(np.abs(y - offset_true)), 1)
        T_analytic = -1.0 / slope

        assert result.parameters["T"] == pytest.approx(T_analytic, rel=0.08)
