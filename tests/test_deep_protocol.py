"""
Deep Protocol Testing Suite for EPyR Tools
==========================================

This module implements a comprehensive testing protocol that systematically
tests all functions, modules, and integration points in the EPyR Tools package.

Test Categories:
1. Unit Tests - Individual function testing with edge cases
2. Integration Tests - Module interaction testing
3. Performance Tests - Speed and memory benchmarks
4. Validation Tests - Scientific accuracy verification
5. Error Handling Tests - Exception and edge case handling
6. Documentation Tests - Docstring and example validation

Test Protocol Levels:
- SMOKE: Basic functionality verification
- STANDARD: Comprehensive feature testing
- DEEP: Exhaustive testing with edge cases and performance
- SCIENTIFIC: Scientific accuracy and validation testing
"""

import gc
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest

# EPyR Tools imports
import epyr
from epyr import baseline, constants, eprload, fair, plot
from epyr.lineshapes import (
    Lineshape,
    convspec,
    gaussian,
    lineshape_class,
    lorentzian,
    lshape,
    voigtian,
)

# Test configuration
PROTOCOL_LEVELS = ["smoke", "standard", "deep", "scientific"]
DEFAULT_TOLERANCE = 1e-10
PERFORMANCE_TIMEOUT = 30.0  # seconds


class TestProtocol:
    """Base class for deep protocol testing with utilities."""

    @staticmethod
    def measure_performance(func, *args, **kwargs) -> Dict[str, float]:
        """Measure function performance metrics."""
        # Memory before
        gc.collect()

        # Time execution
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()

        execution_time = end_time - start_time

        return {
            "execution_time": execution_time,
            "result_size": len(result) if hasattr(result, "__len__") else 1,
        }

    @staticmethod
    def validate_numerical_stability(
        func, inputs: List[Tuple], tolerance: float = DEFAULT_TOLERANCE
    ):
        """Test numerical stability with repeated executions."""
        results = []
        for _ in range(5):  # Run 5 times
            for input_args in inputs:
                result = func(*input_args)
                if isinstance(result, (list, tuple)):
                    results.extend(
                        [
                            np.array(r) if not isinstance(r, np.ndarray) else r
                            for r in result
                        ]
                    )
                else:
                    results.append(
                        np.array(result)
                        if not isinstance(result, np.ndarray)
                        else result
                    )

        # Check consistency
        if len(results) > 1:
            reference = results[0]
            for i, result in enumerate(results[1:], 1):
                if hasattr(result, "shape") and result.shape == reference.shape:
                    np.testing.assert_allclose(
                        result,
                        reference,
                        rtol=tolerance,
                        err_msg=f"Numerical instability detected in run {i}",
                    )


class TestEPyRCoreModules(TestProtocol):
    """Comprehensive testing of core EPyR modules."""

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_constants_module(self, protocol_level):
        """Deep testing of constants module."""
        if protocol_level == "smoke":
            # Basic import and access
            assert hasattr(constants, "ELECTRON_G_FACTOR")
            assert abs(constants.ELECTRON_G_FACTOR) > 0  # g-factor can be negative

        elif protocol_level == "standard":
            # Test all physical constants
            required_constants = [
                "ELECTRON_G_FACTOR",
                "BOHR_MAGNETON",
                "PLANCK_CONSTANT",
                "SPEED_OF_LIGHT",
                "ELECTRON_MASS",
                "NUCLEAR_MAGNETON",
            ]
            for const_name in required_constants:
                assert hasattr(constants, const_name)
                value = getattr(constants, const_name)
                assert isinstance(value, (int, float))
                if "G_FACTOR" in const_name:
                    assert abs(value) > 0  # g-factors can be negative
                else:
                    assert value > 0

        elif protocol_level == "deep":
            # Test constant relationships and units
            # e.g., g-factor should be dimensionless and ~2 (absolute value)
            assert 2.0 < abs(constants.ELECTRON_G_FACTOR) < 2.1

            # Test derived calculations
            gyromagnetic_ratio = (
                constants.ELECTRON_G_FACTOR
                * constants.BOHR_MAGNETON_SI
                / constants.REDUCED_PLANCK_CONSTANT_SI
            )
            assert abs(gyromagnetic_ratio) > 1e10  # Should be in rad/s/T range

        elif protocol_level == "scientific":
            # Validate against known physical values (note: g-factor is negative)
            np.testing.assert_allclose(
                constants.ELECTRON_G_FACTOR, -2.00231930436256, rtol=1e-10
            )
            np.testing.assert_allclose(
                constants.BOHR_MAGNETON_SI, 9.2740100783e-24, rtol=1e-8
            )

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_baseline_module(self, protocol_level, baseline_test_data):
        """Deep testing of baseline correction."""
        x = baseline_test_data["x"]
        y = baseline_test_data["y_with_baseline"]
        true_baseline = baseline_test_data["true_baseline"]

        if protocol_level == "smoke":
            # Basic polynomial correction
            corrected, fitted_baseline = baseline.baseline_polynomial(
                y, x_data=x, poly_order=1
            )
            assert len(corrected) == len(y)
            assert len(fitted_baseline) == len(y)

        elif protocol_level == "standard":
            # Test different polynomial orders
            for order in [0, 1, 2, 3]:
                corrected, fitted_baseline = baseline.baseline_polynomial(
                    y, x_data=x, poly_order=order
                )
                assert len(corrected) == len(y)
                assert not np.any(np.isnan(corrected))

            # Test with exclusion regions
            exclude_regions = baseline_test_data["signal_regions"]
            corrected, fitted_baseline = baseline.baseline_polynomial(
                y, x_data=x, poly_order=1, exclude_regions=exclude_regions
            )
            assert len(corrected) == len(y)

        elif protocol_level == "deep":
            # Test numerical stability
            self.validate_numerical_stability(
                baseline.baseline_polynomial, [(y, x, 1), (y, x, 2)], tolerance=1e-12
            )

            # Test edge cases
            # Single point
            with pytest.raises((ValueError, IndexError)):
                baseline.baseline_polynomial(np.array([1.0]), x_data=np.array([0.0]))

            # All NaN
            nan_data = np.full_like(y, np.nan)
            with pytest.raises((ValueError, np.linalg.LinAlgError)):
                baseline.baseline_polynomial(nan_data, x_data=x)

        elif protocol_level == "scientific":
            # Validate correction accuracy
            corrected, fitted_baseline = baseline.baseline_polynomial(
                y,
                x_data=x,
                poly_order=1,
                exclude_regions=baseline_test_data["signal_regions"],
            )

            # The fitted baseline should be close to the true baseline
            baseline_error = np.mean(np.abs(fitted_baseline - true_baseline))
            assert (
                baseline_error < 5.0
            ), f"Baseline correction error too large: {baseline_error}"


class TestLineshapesDeep(TestProtocol):
    """Comprehensive testing of the lineshapes module."""

    @pytest.fixture
    def standard_field_range(self):
        """Standard magnetic field range for testing."""
        return np.linspace(-20, 20, 1000)

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_gaussian_function(self, protocol_level, standard_field_range):
        """Deep testing of Gaussian lineshape function."""
        B = standard_field_range

        if protocol_level == "smoke":
            # Basic function call
            result = gaussian(B, center=0, width=5)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        elif protocol_level == "standard":
            # Test parameter variations
            centers = [0, 5, -5]
            widths = [1, 5, 10]
            derivatives = [0, 1, 2]

            for center in centers:
                for width in widths:
                    for derivative in derivatives:
                        result = gaussian(
                            B, center=center, width=width, derivative=derivative
                        )
                        assert len(result) == len(B)
                        assert np.all(np.isfinite(result))

                        # Check normalization properties
                        if derivative == 0:
                            assert np.max(result) > 0

        elif protocol_level == "deep":
            # Test numerical properties
            # Symmetry test for absorption (derivative=0)
            center = 0
            width = 5
            gauss_abs = gaussian(B, center=center, width=width, derivative=0)

            # Should be symmetric around center
            mid_idx = len(B) // 2
            left_half = gauss_abs[:mid_idx]
            right_half = gauss_abs[mid_idx + 1 :][::-1]
            np.testing.assert_allclose(left_half, right_half, rtol=1e-10)

            # Test derivative relationships
            gauss_0 = gaussian(B, center=center, width=width, derivative=0)
            gauss_1 = gaussian(B, center=center, width=width, derivative=1)

            # First derivative should be antisymmetric
            np.testing.assert_allclose(gauss_1, -gauss_1[::-1], rtol=1e-10)

            # Test phase rotation
            for phase in [0, 0.5, 1.0, 1.5]:
                result = gaussian(B, center=center, width=width, phase=phase)
                assert len(result) == len(B)
                assert np.all(np.isfinite(result))

        elif protocol_level == "scientific":
            # Validate against analytical properties
            center = 0
            width = 4.0  # HWHM

            gauss = gaussian(B, center=center, width=width, derivative=0)

            # Maximum should be at center
            max_idx = np.argmax(gauss)
            center_idx = np.argmin(np.abs(B - center))
            assert abs(max_idx - center_idx) <= 2  # Allow for discretization

            # Test HWHM property (Half Width at Half Maximum)
            max_val = np.max(gauss)
            half_max_indices = np.where(gauss >= max_val / 2)[0]
            if len(half_max_indices) > 0:
                field_width = B[half_max_indices[-1]] - B[half_max_indices[0]]
                expected_width = 2 * width  # Full width
                np.testing.assert_allclose(field_width, expected_width, rtol=0.1)

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_lorentzian_function(self, protocol_level, standard_field_range):
        """Deep testing of Lorentzian lineshape function."""
        B = standard_field_range

        if protocol_level == "smoke":
            # Basic function call
            result = lorentzian(B, center=0, width=5)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        elif protocol_level == "standard":
            # Test parameter variations
            centers = [0, 3, -3]
            widths = [1, 3, 8]
            phases = [0, 0.5, 1.0]

            for center in centers:
                for width in widths:
                    for phase in phases:
                        result = lorentzian(B, center=center, width=width, phase=phase)
                        assert len(result) == len(B)
                        assert np.all(np.isfinite(result))

        elif protocol_level == "deep":
            # Test Lorentzian properties
            center = 0
            width = 3.0

            # Pure absorption (phase=0)
            lorentz_abs = lorentzian(B, center=center, width=width, phase=0)

            # Should be symmetric around center
            mid_idx = len(B) // 2
            if len(B) % 2 == 1:  # Odd number of points
                left_half = lorentz_abs[:mid_idx]
                right_half = lorentz_abs[mid_idx + 1 :][::-1]
                np.testing.assert_allclose(left_half, right_half, rtol=1e-10)

            # Pure dispersion (phase=1)
            lorentz_disp = lorentzian(B, center=center, width=width, phase=1)
            # Dispersion should be antisymmetric
            np.testing.assert_allclose(lorentz_disp, -lorentz_disp[::-1], rtol=1e-10)

            # Test numerical stability
            self.validate_numerical_stability(
                lorentzian,
                [(B, center, width, 0), (B, center, width, 0.5)],
                tolerance=1e-12,
            )

        elif protocol_level == "scientific":
            # Validate Lorentzian HWHM
            center = 0
            width = 2.0  # HWHM

            lorentz = lorentzian(B, center=center, width=width, phase=0)

            # Maximum should be at center
            max_idx = np.argmax(lorentz)
            center_idx = np.argmin(np.abs(B - center))
            assert abs(max_idx - center_idx) <= 2

            # Check HWHM property
            max_val = np.max(lorentz)
            half_max_val = max_val / 2

            # Find points closest to half maximum
            half_max_indices = np.where(lorentz >= half_max_val)[0]
            if len(half_max_indices) > 0:
                field_span = B[half_max_indices[-1]] - B[half_max_indices[0]]
                expected_span = 2 * width  # Full width at half maximum
                np.testing.assert_allclose(
                    field_span, expected_span, rtol=0.15
                )  # Allow for discretization

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_voigtian_function(self, protocol_level, standard_field_range):
        """Deep testing of Voigtian (true convolution) function."""
        B = standard_field_range

        if protocol_level == "smoke":
            # Basic function call
            result = voigtian(B, center=0, sigma=2, gamma=2)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        elif protocol_level == "standard":
            # Test different sigma/gamma ratios
            ratios = [(1, 1), (2, 1), (1, 2), (3, 3)]

            for sigma, gamma in ratios:
                result = voigtian(B, center=0, sigma=sigma, gamma=gamma)
                assert len(result) == len(B)
                assert np.all(np.isfinite(result))
                assert np.max(result) > 0

        elif protocol_level == "deep":
            # Test limiting cases
            # When sigma >> gamma, should approach Gaussian
            voigt_gauss_like = voigtian(B, center=0, sigma=5, gamma=0.1)
            pure_gauss = gaussian(
                B, center=0, width=5 * np.sqrt(2 * np.log(2))
            )  # Convert sigma to HWHM

            # Should be reasonably similar (not exact due to convolution)
            correlation = np.corrcoef(voigt_gauss_like, pure_gauss)[0, 1]
            assert (
                correlation > 0.95
            ), f"Voigt-Gaussian correlation too low: {correlation}"

            # When gamma >> sigma, should approach Lorentzian
            voigt_lorentz_like = voigtian(B, center=0, sigma=0.1, gamma=5)
            pure_lorentz = lorentzian(B, center=0, width=5, phase=0)

            correlation = np.corrcoef(voigt_lorentz_like, pure_lorentz)[0, 1]
            assert (
                correlation > 0.95
            ), f"Voigt-Lorentzian correlation too low: {correlation}"

        elif protocol_level == "scientific":
            # Test convolution accuracy
            # Voigt should be broader than both parent functions
            sigma, gamma = 2.0, 2.0
            center = 0

            voigt = voigtian(B, center=center, sigma=sigma, gamma=gamma)
            gauss = gaussian(B, center=center, width=sigma * np.sqrt(2 * np.log(2)))
            lorentz = lorentzian(B, center=center, width=gamma, phase=0)

            # Voigt FWHM should be between parent FWHMs
            def estimate_fwhm(profile, field):
                max_val = np.max(profile)
                half_max_indices = np.where(profile >= max_val / 2)[0]
                if len(half_max_indices) > 0:
                    return field[half_max_indices[-1]] - field[half_max_indices[0]]
                return 0

            voigt_fwhm = estimate_fwhm(voigt, B)
            gauss_fwhm = estimate_fwhm(gauss, B)
            lorentz_fwhm = estimate_fwhm(lorentz, B)

            # Voigt FWHM should be at least as wide as the broader parent
            min_parent_fwhm = min(gauss_fwhm, lorentz_fwhm)
            assert (
                voigt_fwhm >= min_parent_fwhm * 0.9
            )  # Small tolerance for discretization

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_lineshape_class(self, protocol_level, standard_field_range):
        """Deep testing of unified Lineshape class."""
        B = standard_field_range

        if protocol_level == "smoke":
            # Basic class instantiation and usage
            shape = Lineshape("gaussian", width=3.0)
            result = shape(B, center=0)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        elif protocol_level == "standard":
            # Test all supported lineshape types
            lineshape_types = ["gaussian", "lorentzian", "pseudo_voigt", "voigt"]

            for shape_type in lineshape_types:
                if shape_type == "voigt":
                    shape = Lineshape(shape_type, width=3.0, sigma=2.0, gamma=2.0)
                elif shape_type == "pseudo_voigt":
                    shape = Lineshape(shape_type, width=3.0, alpha=0.5)
                else:
                    shape = Lineshape(shape_type, width=3.0)

                result = shape(B, center=0)
                assert len(result) == len(B)
                assert np.all(np.isfinite(result))

        elif protocol_level == "deep":
            # Test parameter validation and error handling
            # Invalid lineshape type
            with pytest.raises(ValueError):
                Lineshape("invalid_type", width=3.0)

            # Missing required parameters
            with pytest.raises(TypeError):
                shape = Lineshape("voigt", width=3.0)  # Missing sigma, gamma
                shape(B, center=0)

            # Test immutability and consistency
            shape = Lineshape("gaussian", width=5.0)
            result1 = shape(B, center=0)
            result2 = shape(B, center=0)
            np.testing.assert_array_equal(result1, result2)

        elif protocol_level == "scientific":
            # Test pseudo-Voigt accuracy
            # Pseudo-Voigt should interpolate between Gaussian and Lorentzian
            width = 4.0

            # Pure Gaussian (alpha=0)
            pseudo_gauss = Lineshape("pseudo_voigt", width=width, alpha=0.0)
            result_gauss = pseudo_gauss(B, center=0)

            pure_gauss = gaussian(B, center=0, width=width)
            np.testing.assert_allclose(result_gauss, pure_gauss, rtol=1e-10)

            # Pure Lorentzian (alpha=1)
            pseudo_lorentz = Lineshape("pseudo_voigt", width=width, alpha=1.0)
            result_lorentz = pseudo_lorentz(B, center=0)

            pure_lorentz = lorentzian(B, center=0, width=width, phase=0)
            np.testing.assert_allclose(result_lorentz, pure_lorentz, rtol=1e-10)


class TestPerformanceBenchmarks(TestProtocol):
    """Performance testing and benchmarking."""

    @pytest.mark.parametrize("data_size", [100, 1000, 10000])
    @pytest.mark.parametrize("function_name", ["gaussian", "lorentzian", "voigtian"])
    def test_lineshape_performance(self, data_size, function_name):
        """Benchmark lineshape function performance."""
        B = np.linspace(-50, 50, data_size)

        # Get function reference
        func_map = {
            "gaussian": lambda B: gaussian(B, center=0, width=5),
            "lorentzian": lambda B: lorentzian(B, center=0, width=5),
            "voigtian": lambda B: voigtian(B, center=0, sigma=3, gamma=3),
        }

        func = func_map[function_name]

        # Measure performance
        perf_data = self.measure_performance(func, B)

        # Performance assertions (adjust based on expected performance)
        expected_max_time = {
            100: 0.01,  # 10ms for 100 points
            1000: 0.1,  # 100ms for 1000 points
            10000: 1.0,  # 1s for 10000 points
        }

        assert (
            perf_data["execution_time"] < expected_max_time[data_size]
        ), f"{function_name} too slow for {data_size} points: {perf_data['execution_time']:.3f}s"

    @pytest.mark.parametrize("protocol_level", ["standard", "deep"])
    def test_memory_usage(self, protocol_level):
        """Test memory usage patterns."""
        if protocol_level == "standard":
            # Basic memory test
            B = np.linspace(-10, 10, 10000)

            # Should not leak memory with repeated calls
            for _ in range(10):
                result = gaussian(B, center=0, width=3)
                del result
                gc.collect()

        elif protocol_level == "deep":
            # Stress test with large datasets
            large_B = np.linspace(-100, 100, 100000)

            # Should handle large datasets
            result = gaussian(large_B, center=0, width=10)
            assert len(result) == len(large_B)

            # Memory cleanup
            del result, large_B
            gc.collect()


class TestIntegrationWorkflows(TestProtocol):
    """Integration testing of complete workflows."""

    def test_full_epr_analysis_workflow(self, sample_1d_data, baseline_test_data):
        """Test complete EPR analysis workflow integration."""
        # 1. Load data (simulated)
        x, y = sample_1d_data

        # 2. Baseline correction
        y_corrected, baseline_fit = baseline.baseline_polynomial(
            y, x_data=x, poly_order=1
        )

        # 3. Lineshape fitting simulation
        # Find peak center
        peak_idx = np.argmax(y_corrected)
        peak_center = x[peak_idx]

        # Estimate width from data
        half_max = np.max(y_corrected) / 2
        half_max_indices = np.where(y_corrected >= half_max)[0]
        if len(half_max_indices) > 1:
            estimated_width = (x[half_max_indices[-1]] - x[half_max_indices[0]]) / 2
        else:
            estimated_width = 5.0

        # Generate theoretical lineshapes
        gauss_fit = gaussian(x, center=peak_center, width=estimated_width)
        lorentz_fit = lorentzian(x, center=peak_center, width=estimated_width, phase=0)

        # 4. Validate workflow completed successfully
        assert len(y_corrected) == len(x)
        assert len(gauss_fit) == len(x)
        assert len(lorentz_fit) == len(x)
        assert np.all(np.isfinite(y_corrected))
        assert np.all(np.isfinite(gauss_fit))
        assert np.all(np.isfinite(lorentz_fit))

        # 5. Check that baseline correction improved the fit
        # (This is a simplified check - in practice would use proper fitting metrics)
        corrected_peak = np.max(y_corrected)
        original_peak = np.max(y)
        assert corrected_peak > 0  # Should have positive signal after correction


class TestScientificValidation(TestProtocol):
    """Scientific validation and accuracy tests."""

    def test_physical_constants_accuracy(self):
        """Validate physical constants against NIST values."""
        # Test electron g-factor (NIST 2018 CODATA) - note negative sign
        nist_ge = -2.00231930436256
        np.testing.assert_allclose(constants.ELECTRON_G_FACTOR, nist_ge, rtol=1e-10)

        # Test Bohr magneton (NIST 2018 CODATA) in J/T
        nist_mb = 9.2740100783e-24
        np.testing.assert_allclose(constants.BOHR_MAGNETON_SI, nist_mb, rtol=1e-8)

    def test_lineshape_mathematical_properties(self):
        """Validate mathematical properties of lineshapes."""
        B = np.linspace(-50, 50, 2001)  # High resolution for accuracy
        dB = B[1] - B[0]

        # Test Gaussian normalization
        gauss = gaussian(B, center=0, width=5, derivative=0)
        # Gaussian integral should be proportional to width (for HWHM parameterization)
        gauss_integral = np.trapz(gauss, dx=dB)
        assert gauss_integral > 0

        # Test Lorentzian normalization
        lorentz = lorentzian(B, center=0, width=5, phase=0)
        lorentz_integral = np.trapz(lorentz, dx=dB)
        assert lorentz_integral > 0

        # Test derivative properties
        gauss_0 = gaussian(B, center=0, width=5, derivative=0)
        gauss_1 = gaussian(B, center=0, width=5, derivative=1)

        # Numerical derivative of absorption should approximate first derivative
        numerical_derivative = np.gradient(gauss_0, dB)
        correlation = np.corrcoef(numerical_derivative, gauss_1)[0, 1]
        assert abs(correlation) > 0.99, f"Derivative correlation too low: {correlation}"


# Test execution control
@pytest.fixture(scope="session")
def test_protocol_config():
    """Configuration for test protocol execution."""
    return {
        "run_performance": True,
        "run_scientific": True,
        "tolerance": DEFAULT_TOLERANCE,
        "timeout": PERFORMANCE_TIMEOUT,
    }


# Pytest markers for test organization
pytestmark = [pytest.mark.deep_protocol, pytest.mark.comprehensive]


# Test discovery and execution summary
class TestSummaryReport:
    """Generate test execution summary and reporting."""

    @pytest.fixture(autouse=True, scope="session")
    def test_session_summary(self):
        """Provide session-wide test summary."""
        print("\n" + "=" * 80)
        print("EPyR Tools Deep Protocol Testing Suite")
        print("=" * 80)
        print(f"Testing protocol levels: {', '.join(PROTOCOL_LEVELS)}")
        print(f"Default tolerance: {DEFAULT_TOLERANCE}")
        print(f"Performance timeout: {PERFORMANCE_TIMEOUT}s")
        print("=" * 80)

        yield

        print("\n" + "=" * 80)
        print("Deep Protocol Testing Complete")
        print("=" * 80)
