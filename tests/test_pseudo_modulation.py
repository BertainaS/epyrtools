"""Tests for epyr.lineshapes.pseudo_modulation."""

import numpy as np
import pytest

from epyr.lineshapes.gaussian import gaussian
from epyr.lineshapes.pseudo_modulation import pseudo_modulation


@pytest.mark.smoke
class TestPseudoModulationBasic:
    def test_basic_shape_and_finiteness(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)

        assert y_pm.shape == y.shape
        assert np.all(np.isfinite(y_pm))
        assert np.isrealobj(y_pm)

    def test_harmonic_1_real_input_has_nonzero_signal(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)

        assert np.max(np.abs(y_pm)) > 1e-6


@pytest.mark.standard
class TestPseudoModulationValidation:
    def test_harmonic_2_runs_and_differs_from_harmonic_1(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm_1 = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)
        y_pm_2 = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=2)

        assert y_pm_2.shape == y.shape
        assert np.all(np.isfinite(y_pm_2))
        assert not np.allclose(y_pm_1, y_pm_2)

    def test_pad_false_same_shape(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1, pad=False)

        assert y_pm.shape == y.shape
        assert np.all(np.isfinite(y_pm))

    def test_complex_input_stays_complex(self):
        x = np.linspace(-50, 50, 500)
        y = gaussian(x, 0.0, 8.0).astype(complex) * (1 + 0.3j)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)

        assert np.iscomplexobj(y_pm)
        assert y_pm.shape == y.shape

    def test_nonuniform_x_raises(self):
        x = np.concatenate([np.linspace(-50, -10, 500), np.linspace(-9, 50, 300)])
        y = gaussian(x, 0.0, 8.0)

        with pytest.raises(ValueError, match="uniformly spaced"):
            pseudo_modulation(x, y, mod_amplitude=2.0)

    @pytest.mark.parametrize("bad_amplitude", [0.0, -1.0])
    def test_invalid_mod_amplitude_raises(self, bad_amplitude):
        x = np.linspace(-50, 50, 500)
        y = gaussian(x, 0.0, 8.0)

        with pytest.raises(ValueError, match="mod_amplitude"):
            pseudo_modulation(x, y, mod_amplitude=bad_amplitude)

    @pytest.mark.parametrize("bad_harmonic", [0, 3, -1])
    def test_invalid_harmonic_raises(self, bad_harmonic):
        x = np.linspace(-50, 50, 500)
        y = gaussian(x, 0.0, 8.0)

        with pytest.raises(ValueError, match="harmonic"):
            pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=bad_harmonic)

    def test_shape_mismatch_raises(self):
        x = np.linspace(-50, 50, 500)
        y = np.zeros(400)

        with pytest.raises(ValueError, match="Shape mismatch"):
            pseudo_modulation(x, y, mod_amplitude=2.0)


@pytest.mark.scientific
class TestPseudoModulationAnalyticLimit:
    def test_first_harmonic_small_amplitude_matches_derivative(self):
        x = np.linspace(-50, 50, 4000)
        center, width = 0.0, 8.0
        y = gaussian(x, center, width, derivative=0)
        dy_analytic = gaussian(x, center, width, derivative=1)

        eps = 1e-3
        y_pm = pseudo_modulation(x, y, mod_amplitude=eps, harmonic=1)
        scaled = y_pm / (eps / 4)

        # Stay away from the outermost few points where the FFT treats
        # the (non-periodic) signal as if it wrapped around.
        mask = np.abs(x) < 30
        max_dy = np.max(np.abs(dy_analytic[mask]))

        assert np.allclose(scaled[mask], dy_analytic[mask], atol=0.02 * max_dy)


@pytest.mark.smoke
def test_pseudo_modulation_importable_from_package():
    from epyr.lineshapes import pseudo_modulation as exported

    assert exported is pseudo_modulation
