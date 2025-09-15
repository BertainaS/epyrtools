import numpy as np
import pytest

try:
    from epyr.baseline import (
        baseline_constant_offset,
        baseline_mono_exponential,
        baseline_polynomial,
        baseline_polynomial_2d,
        baseline_stretched_exponential,
    )

    BASELINE_AVAILABLE = True
except ImportError:
    BASELINE_AVAILABLE = False


@pytest.mark.skipif(not BASELINE_AVAILABLE, reason="Baseline module not available")
class TestBaseline:
    """Test suite for baseline correction functions."""

    def test_baseline_polynomial_basic(self, baseline_test_data):
        """Test basic polynomial baseline correction."""
        data = baseline_test_data
        x = data["x"]
        y = data["y_with_baseline"]
        signal_regions = data["signal_regions"]

        # Perform baseline correction excluding signal regions
        y_corrected, baseline_fit = baseline_polynomial(
            y, x_data=x, poly_order=1, exclude_regions=signal_regions
        )

        # Check that baseline was removed
        assert len(y_corrected) == len(y)
        assert len(baseline_fit) == len(y)

        # The corrected data should have much smaller baseline slope
        corrected_slope = np.polyfit(x, y_corrected, 1)[0]
        original_slope = np.polyfit(x, y, 1)[0]
        assert abs(corrected_slope) < abs(original_slope)

        # Baseline should be approximately the known baseline shape
        correlation = np.corrcoef(baseline_fit, data["true_baseline"])[0, 1]
        assert correlation > 0.8  # Should be highly correlated

    def test_baseline_polynomial_zero_order(self):
        """Test constant baseline correction (order 0)."""
        # Create test data with constant offset
        np.random.seed(42)  # For reproducible tests
        x = np.linspace(0, 100, 50)
        baseline_offset = 10.0
        y = np.ones_like(x) * baseline_offset + np.random.normal(0, 0.1, len(x))

        y_corrected, baseline_fit = baseline_polynomial(y, x_data=x, poly_order=0)

        # Mean of corrected data should be close to zero
        assert abs(np.mean(y_corrected)) < 1.0
        # Baseline fit should be approximately constant
        assert np.std(baseline_fit) < 1.0
        assert abs(np.mean(baseline_fit) - baseline_offset) < 1.0

    def test_baseline_polynomial_input_validation(self):
        """Test input validation for baseline correction."""
        x = np.linspace(0, 10, 10)
        y = np.ones(10)

        # Test with invalid polynomial order
        with pytest.raises((ValueError, TypeError)):
            baseline_polynomial(y, x_data=x, poly_order=-1)

        # Test with mismatched x and y lengths
        with pytest.raises(ValueError):
            baseline_polynomial(y, x_data=x[:-1], poly_order=1)

        # Test with invalid exclude regions format
        with pytest.raises((ValueError, TypeError)):
            baseline_polynomial(y, x_data=x, exclude_regions=["invalid"])

    def test_baseline_polynomial_higher_orders(self, baseline_test_data):
        """Test polynomial baseline correction with different orders."""
        data = baseline_test_data
        x = data["x"]
        y = data["y_with_baseline"]

        # Test different polynomial orders
        for order in [0, 1, 2, 3]:
            y_corrected, baseline_fit = baseline_polynomial(
                y, x_data=x, poly_order=order, exclude_regions=data["signal_regions"]
            )

            assert len(y_corrected) == len(y)
            assert len(baseline_fit) == len(y)
            # Higher order should generally fit better (lower residuals)
            residuals = y - baseline_fit
            assert np.std(residuals) > 0  # Should have some variance
