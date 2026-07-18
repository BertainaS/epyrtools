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
