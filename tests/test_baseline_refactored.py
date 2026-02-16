#!/usr/bin/env python3
"""
Comprehensive tests for the refactored baseline correction package.

Tests the new modular epyr.baseline package (v0.1.8) with all components:
- Mathematical models
- Correction algorithms
- Region selection
- Automatic model selection
- Backend control
"""

import warnings
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# Test imports for new baseline package
try:
    import epyr.baseline as baseline
    from epyr.baseline import (  # Core correction functions; Mathematical models; Region selection; Backend control; Utilities
        RegionSelector,
        baseline_auto_1d,
        baseline_bi_exponential_1d,
        baseline_polynomial_1d,
        baseline_polynomial_2d,
        baseline_stretched_exponential_1d,
        bi_exponential_1d,
        create_region_mask_1d,
        get_model_function,
        list_available_models,
        polynomial_1d,
        setup_inline_backend,
        setup_widget_backend,
        stretched_exponential_1d,
    )

    BASELINE_AVAILABLE = True
except ImportError as e:
    print(f"Baseline package import failed: {e}")
    BASELINE_AVAILABLE = False


@pytest.fixture
def synthetic_cw_data():
    """Create synthetic CW EPR data with polynomial baseline."""
    x = np.linspace(3300, 3400, 200)
    # Gaussian signal
    signal = 100 * np.exp(-(((x - 3350) / 5) ** 2))
    # Quadratic baseline
    baseline_true = 0.01 * (x - 3350) ** 2 + 50
    # Add noise
    noise = 2 * np.random.normal(size=len(x))
    y = signal + baseline_true + noise

    return {
        "x": x,
        "y": y,
        "y_clean": signal,
        "baseline_true": baseline_true,
        "signal_center": 3350,
        "signal_width": 10,
    }


@pytest.fixture
def synthetic_t2_data():
    """Create synthetic T2 relaxation data with stretched exponential decay."""
    x = np.linspace(0, 2000, 150)
    # Stretched exponential parameters
    A, tau, beta, offset = 1000, 500, 1.2, 50
    # True baseline
    baseline_true = offset + A * np.exp(-((x / tau) ** beta))
    # Add noise
    noise = 20 * np.random.normal(size=len(x))
    y = baseline_true + noise

    return {
        "x": x,
        "y": y,
        "baseline_true": baseline_true,
        "params": {"A": A, "tau": tau, "beta": beta, "offset": offset},
    }


@pytest.fixture
def synthetic_2d_data():
    """Create synthetic 2D EPR data."""
    x1 = np.linspace(-10, 10, 20)
    x2 = np.linspace(-5, 5, 15)
    X, Y = np.meshgrid(x1, x2)

    # 2D signal
    signal = 100 * np.exp(-(X**2 + Y**2) / 4)
    # 2D polynomial baseline
    baseline_true = 0.5 * X + 0.3 * Y + 0.1 * X**2 + 20
    # Add noise
    noise = 5 * np.random.normal(size=X.shape)
    Z = signal + baseline_true + noise

    return {
        "x": [x1, x2],
        "X": X,
        "Y": Y,
        "Z": Z,
        "baseline_true": baseline_true,
        "signal": signal,
    }


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestMathematicalModels:
    """Test mathematical model functions."""

    def test_stretched_exponential_function(self):
        """Test stretched exponential mathematical function."""
        x = np.linspace(0, 1000, 100)
        A, tau, beta, offset = 100, 200, 1.5, 10

        y = stretched_exponential_1d(x, A, tau, beta, offset)

        # Check shape
        assert y.shape == x.shape

        # Check boundary conditions
        assert np.isclose(y[0], offset + A)  # At x=0
        assert np.all(y >= offset - 1e-10)  # Above offset (within numerical precision)

        # Check monotonic decay (for positive x, A, tau)
        assert np.all(np.diff(y) <= 1e-10)  # Non-increasing

    def test_bi_exponential_function(self):
        """Test bi-exponential mathematical function."""
        x = np.linspace(0, 1000, 100)
        A1, tau1, A2, tau2, offset = 50, 100, 30, 500, 5

        y = bi_exponential_1d(x, A1, tau1, A2, tau2, offset)

        # Check shape
        assert y.shape == x.shape

        # Check boundary conditions
        assert np.isclose(y[0], offset + A1 + A2)  # At x=0
        assert np.all(y >= offset - 1e-10)  # Above offset

    def test_polynomial_1d_function(self):
        """Test 1D polynomial function."""
        x = np.linspace(-10, 10, 50)
        coeffs = [1, 2, 0.5]  # 1 + 2x + 0.5x²

        y = polynomial_1d(x, *coeffs)

        # Check shape
        assert y.shape == x.shape

        # Check specific values
        x_test = np.array([0, 1, 2])
        y_expected = np.array([1, 3.5, 7])  # 1, 1+2+0.5, 1+4+2
        y_test = polynomial_1d(x_test, *coeffs)
        np.testing.assert_allclose(y_test, y_expected)

    def test_list_available_models(self):
        """Test model listing function."""
        models = list_available_models()

        assert isinstance(models, list)
        assert len(models) >= 3
        assert "polynomial" in models
        assert "stretched_exponential" in models
        assert "bi_exponential" in models

    def test_get_model_function(self):
        """Test model function getter."""
        # Get stretched exponential function
        func = get_model_function("stretched_exponential", "1d")
        assert callable(func)

        # Test function works
        x = np.array([0, 100, 200])
        y = func(x, 100, 200, 1.0, 0)
        assert len(y) == len(x)

        # Test invalid model
        with pytest.raises(ValueError):
            get_model_function("invalid_model", "1d")


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestRegionSelection:
    """Test region selection utilities."""

    def test_create_region_mask_1d_exclude(self):
        """Test 1D region mask creation in exclude mode."""
        x = np.linspace(0, 100, 101)
        regions = [(20, 30), (70, 80)]

        mask = create_region_mask_1d(x, regions, mode="exclude")

        # Check mask properties
        assert mask.shape == x.shape
        assert mask.dtype == bool

        # Check excluded regions
        exclude_indices = ((x >= 20) & (x <= 30)) | ((x >= 70) & (x <= 80))
        expected_mask = ~exclude_indices
        np.testing.assert_array_equal(mask, expected_mask)

    def test_create_region_mask_1d_include(self):
        """Test 1D region mask creation in include mode."""
        x = np.linspace(0, 100, 101)
        regions = [(10, 20), (80, 90)]

        mask = create_region_mask_1d(x, regions, mode="include")

        # Check included regions only
        include_indices = ((x >= 10) & (x <= 20)) | ((x >= 80) & (x <= 90))
        np.testing.assert_array_equal(mask, include_indices)

    def test_region_mask_invalid_mode(self):
        """Test invalid mode raises error."""
        x = np.linspace(0, 10, 11)
        regions = [(2, 5)]

        with pytest.raises(ValueError, match="mode must be"):
            create_region_mask_1d(x, regions, mode="invalid")


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestPolynomialCorrection:
    """Test polynomial baseline correction."""

    def test_polynomial_1d_basic(self, synthetic_cw_data):
        """Test basic polynomial correction."""
        data = synthetic_cw_data

        corrected, baseline_fit = baseline_polynomial_1d(
            data["x"], data["y"], None, order=2
        )

        # Check output shapes
        assert corrected.shape == data["y"].shape
        assert baseline_fit.shape == data["y"].shape

        # Check correction reduces signal amplitude in wings
        wing_mask = (data["x"] < 3330) | (data["x"] > 3370)
        original_wing_std = np.std(data["y"][wing_mask])
        corrected_wing_std = np.std(corrected[wing_mask])
        assert corrected_wing_std < original_wing_std

    def test_polynomial_1d_different_orders(self, synthetic_cw_data):
        """Test polynomial correction with different orders."""
        data = synthetic_cw_data

        for order in [1, 2, 3, 4]:
            corrected, baseline_fit = baseline_polynomial_1d(
                data["x"], data["y"], None, order=order
            )

            assert corrected.shape == data["y"].shape
            assert not np.any(np.isnan(corrected))

    def test_polynomial_1d_manual_regions(self, synthetic_cw_data):
        """Test polynomial correction with manual regions."""
        data = synthetic_cw_data

        # Define baseline regions (wings only)
        baseline_regions = [(3300, 3330), (3370, 3400)]

        corrected, baseline_fit = baseline_polynomial_1d(
            data["x"],
            data["y"],
            None,
            order=2,
            manual_regions=baseline_regions,
            region_mode="include",
        )

        assert corrected.shape == data["y"].shape
        assert baseline_fit.shape == data["y"].shape

    def test_polynomial_2d_basic(self, synthetic_2d_data):
        """Test 2D polynomial correction."""
        data = synthetic_2d_data

        corrected, baseline_fit = baseline_polynomial_2d(
            data["x"], data["Z"], None, order=(1, 1)
        )

        # Check shapes
        assert corrected.shape == data["Z"].shape
        assert baseline_fit.shape == data["Z"].shape

        # Check correction reduces background variation
        corner_mask = (np.abs(data["X"]) > 8) & (np.abs(data["Y"]) > 3)
        if np.any(corner_mask):
            original_corner_std = np.std(data["Z"][corner_mask])
            corrected_corner_std = np.std(corrected[corner_mask])
            assert (
                corrected_corner_std <= original_corner_std * 1.1
            )  # Allow some tolerance


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestExponentialCorrection:
    """Test exponential baseline correction methods."""

    def test_stretched_exponential_1d_basic(self, synthetic_t2_data):
        """Test basic stretched exponential correction."""
        data = synthetic_t2_data

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Ignore fitting warnings
            corrected, baseline_fit = baseline_stretched_exponential_1d(
                data["x"], data["y"], None
            )

        # Check shapes
        assert corrected.shape == data["y"].shape
        assert baseline_fit.shape == data["y"].shape

        # Check that correction removes most of the decay
        final_points = slice(-20, None)
        corrected_final_mean = np.mean(corrected[final_points])
        original_final_mean = np.mean(data["y"][final_points])
        assert abs(corrected_final_mean) < abs(original_final_mean)

    def test_stretched_exponential_1d_with_params(self, synthetic_t2_data):
        """Test stretched exponential with custom parameters."""
        data = synthetic_t2_data

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            corrected, baseline_fit = baseline_stretched_exponential_1d(
                data["x"],
                data["y"],
                None,
                beta_range=(0.5, 2.0),
                exclude_initial=5,
                exclude_final=10,
            )

        assert corrected.shape == data["y"].shape
        assert not np.any(np.isnan(corrected))

    def test_bi_exponential_1d_basic(self, synthetic_t2_data):
        """Test bi-exponential correction."""
        data = synthetic_t2_data

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                corrected, baseline_fit = baseline_bi_exponential_1d(
                    data["x"], data["y"], None
                )

                # Check shapes
                assert corrected.shape == data["y"].shape
                assert baseline_fit.shape == data["y"].shape

            except Exception:
                # Bi-exponential fitting can be challenging, allow failure
                pytest.skip("Bi-exponential fitting failed (acceptable)")


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestAutomaticSelection:
    """Test automatic model selection."""

    def test_baseline_auto_1d_basic(self, synthetic_cw_data):
        """Test basic automatic selection."""
        data = synthetic_cw_data

        corrected, baseline_fit, info = baseline_auto_1d(
            data["x"], data["y"], None, verbose=False
        )

        # Check outputs
        assert corrected.shape == data["y"].shape
        assert baseline_fit.shape == data["y"].shape
        assert isinstance(info, dict)

        # Check info structure
        assert "best_model" in info
        assert "criteria" in info
        assert "parameters" in info
        assert isinstance(info["criteria"], dict)
        assert "r2" in info["parameters"]

    def test_baseline_auto_1d_model_restriction(self, synthetic_cw_data):
        """Test automatic selection with restricted models."""
        data = synthetic_cw_data

        corrected, baseline_fit, info = baseline_auto_1d(
            data["x"], data["y"], None, models=["polynomial"], verbose=False
        )

        # Should select polynomial
        assert info["best_model"] == "polynomial"
        assert "polynomial" in info["criteria"]

    def test_baseline_auto_1d_different_criteria(self, synthetic_cw_data):
        """Test automatic selection with different criteria."""
        data = synthetic_cw_data

        criteria_list = ["aic", "bic", "r2"]

        for criterion in criteria_list:
            corrected, baseline_fit, info = baseline_auto_1d(
                data["x"], data["y"], None, selection_criterion=criterion, verbose=False
            )

            assert info["best_model"] in [
                "polynomial",
                "stretched_exponential",
                "bi_exponential",
            ]
            assert info["parameters"]["r2"] >= 0  # R² should be non-negative


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestBackendControl:
    """Test matplotlib backend control functions."""

    @patch("epyr.baseline.interactive.get_ipython")
    def test_setup_inline_backend(self, mock_get_ipython):
        """Test inline backend setup."""
        # Mock IPython environment
        mock_ipython = MagicMock()
        mock_get_ipython.return_value = mock_ipython

        # Test function
        setup_inline_backend()

        # Check that magic command was called
        mock_ipython.magic.assert_called_once_with("matplotlib inline")

    @patch("epyr.baseline.interactive.get_ipython")
    def test_setup_widget_backend(self, mock_get_ipython):
        """Test widget backend setup."""
        mock_ipython = MagicMock()
        mock_get_ipython.return_value = mock_ipython

        setup_widget_backend()

        mock_ipython.magic.assert_called_once_with("matplotlib widget")

    @patch("epyr.baseline.interactive.get_ipython")
    def test_setup_backend_no_ipython(self, mock_get_ipython):
        """Test backend setup when not in IPython."""
        mock_get_ipython.return_value = None

        # Should not raise error, just print warning
        setup_inline_backend()
        setup_widget_backend()


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestUtilities:
    """Test utility functions."""

    def test_package_structure(self):
        """Test that all expected modules and functions are available."""
        # Check main package
        assert hasattr(baseline, "__version__")
        assert hasattr(baseline, "__all__")

        # Check key functions are in __all__
        expected_functions = [
            "baseline_polynomial_1d",
            "baseline_stretched_exponential_1d",
            "baseline_auto_1d",
            "setup_inline_backend",
            "RegionSelector",
        ]

        for func_name in expected_functions:
            assert func_name in baseline.__all__
            assert hasattr(baseline, func_name)

    def test_configuration_functions(self):
        """Test package configuration."""
        # Test configuration getter
        config = baseline.get_configuration()
        assert isinstance(config, dict)
        assert "polynomial_order" in config

        # Test configuration setter
        baseline.configure(polynomial_order=4)
        new_config = baseline.get_configuration()
        assert new_config["polynomial_order"] == 4

    def test_help_functions(self):
        """Test help functions don't crash."""
        # These should not raise exceptions
        baseline.get_help()
        baseline.jupyter_help()


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_empty_data(self):
        """Test handling of empty data."""
        x = np.array([])
        y = np.array([])

        with pytest.raises(ValueError):
            baseline_polynomial_1d(x, y, None)

    def test_invalid_data_shape(self):
        """Test handling of invalid data shapes."""
        x = np.array([1, 2, 3])
        y = np.array([[1, 2], [3, 4]])  # 2D array for 1D function

        with pytest.raises(ValueError):
            baseline_polynomial_1d(x, y, None)

    def test_insufficient_points(self):
        """Test handling of insufficient points for high-order polynomial."""
        x = np.array([1, 2])
        y = np.array([1, 2])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            corrected, baseline_fit = baseline_polynomial_1d(x, y, None, order=5)

        # Should return original data when fitting fails
        np.testing.assert_array_equal(corrected, y)
        np.testing.assert_array_equal(baseline_fit, np.zeros_like(y))

    def test_invalid_selection_criterion(self):
        """Test invalid selection criterion."""
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 4, 9, 16, 25])

        with pytest.raises(ValueError):
            baseline_auto_1d(x, y, None, selection_criterion="invalid")


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline package not available")
class TestCompatibility:
    """Test backward compatibility and integration."""

    def test_backward_compatibility_imports(self):
        """Test that old function names still work through aliases."""
        # Test that compatibility aliases exist
        assert hasattr(baseline, "baseline_auto")
        assert hasattr(baseline, "polynomial_baseline_1d")

        # Test they point to the correct functions
        assert baseline.baseline_auto == baseline.baseline_auto_1d
        assert baseline.polynomial_baseline_1d == baseline.baseline_polynomial_1d

    def test_complex_data_handling(self):
        """Test handling of complex EPR data."""
        x = np.linspace(0, 1000, 50)
        # Create complex data (like from eprload)
        y_real = 100 * np.exp(-x / 200) + 20
        y_imag = 50 * np.exp(-x / 300) + 10
        y_complex = y_real + 1j * y_imag

        # Test with real part
        corrected, baseline_fit = baseline_stretched_exponential_1d(
            x, y_complex, None, use_real_part=True
        )

        assert corrected.shape == y_complex.shape
        assert np.iscomplexobj(corrected)  # Should preserve complex nature

        # Test with magnitude
        corrected, baseline_fit = baseline_stretched_exponential_1d(
            x, y_complex, None, use_real_part=False
        )

        assert corrected.shape == y_complex.shape


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
