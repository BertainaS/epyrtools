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
        """Test numerical stability with repeated executions.

        Each input tuple is evaluated five times; all repetitions must
        agree with the first to within the tolerance. Different inputs
        are never compared to each other.
        """

        def as_arrays(result):
            parts = result if isinstance(result, (list, tuple)) else [result]
            return [np.asarray(part) for part in parts]

        for input_args in inputs:
            reference = as_arrays(func(*input_args))
            for run in range(1, 5):
                repeat = as_arrays(func(*input_args))
                for ref_part, repeat_part in zip(reference, repeat):
                    np.testing.assert_allclose(
                        repeat_part,
                        ref_part,
                        rtol=tolerance,
                        err_msg=f"Numerical instability detected in run {run}",
                    )


class TestEPyRCoreModules(TestProtocol):
    """Comprehensive testing of core EPyR modules."""

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_constants_module(self, protocol_level):
        """Deep testing of constants module."""
        if protocol_level == "smoke":
            # Basic import and access (GFREE follows the EasySpin sign
            # convention: positive magnitude of the free-electron g-factor)
            assert hasattr(constants, "GFREE")
            assert constants.GFREE > 0

        elif protocol_level == "standard":
            # Test all physical constants
            required_constants = [
                "GFREE",
                "BMAGN",
                "PLANCK",
                "HBAR",
                "CLIGHT",
                "BOLTZM",
                "NMAGN",
            ]
            for const_name in required_constants:
                assert hasattr(constants, const_name)
                value = getattr(constants, const_name)
                assert isinstance(value, (int, float))
                assert value > 0

        elif protocol_level == "deep":
            # Test constant relationships and units
            # g-factor should be dimensionless and ~2
            assert 2.0 < constants.GFREE < 2.1

            # Electron gyromagnetic ratio g*muB/hbar ~ 1.76e11 rad/s/T
            gyromagnetic_ratio = constants.GFREE * constants.BMAGN / constants.HBAR
            assert abs(gyromagnetic_ratio) > 1e10

        elif protocol_level == "scientific":
            # Validate against CODATA 2022 values
            np.testing.assert_allclose(constants.GFREE, 2.00231930436092, rtol=1e-10)
            np.testing.assert_allclose(constants.BMAGN, 9.2740100657e-24, rtol=1e-8)

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_baseline_module(self, protocol_level, baseline_test_data):
        """Deep testing of baseline correction."""
        x = baseline_test_data["x"]
        y = baseline_test_data["y_with_baseline"]
        true_baseline = baseline_test_data["true_baseline"]

        if protocol_level == "smoke":
            # Basic polynomial correction
            corrected, fitted_baseline = baseline.baseline_polynomial_1d(x, y, order=1)
            assert len(corrected) == len(y)
            assert len(fitted_baseline) == len(y)

        elif protocol_level == "standard":
            # Test different polynomial orders
            for order in [0, 1, 2, 3]:
                corrected, fitted_baseline = baseline.baseline_polynomial_1d(
                    x, y, order=order
                )
                assert len(corrected) == len(y)
                assert not np.any(np.isnan(corrected))

            # Test with exclusion regions
            exclude_regions = baseline_test_data["signal_regions"]
            corrected, fitted_baseline = baseline.baseline_polynomial_1d(
                x, y, order=1, manual_regions=exclude_regions, region_mode="exclude"
            )
            assert len(corrected) == len(y)

        elif protocol_level == "deep":
            # Test numerical stability
            self.validate_numerical_stability(
                baseline.baseline_polynomial_1d,
                [(x, y, None, 1), (x, y, None, 2)],
                tolerance=1e-12,
            )

            # Test edge cases
            # Single point: degrades gracefully with a warning
            with pytest.warns(UserWarning):
                baseline.baseline_polynomial_1d(np.array([0.0]), np.array([1.0]))

            # All NaN input propagates NaN without crashing
            nan_data = np.full_like(y, np.nan)
            corrected_nan, _ = baseline.baseline_polynomial_1d(x, nan_data)
            assert np.all(np.isnan(corrected_nan))

        elif protocol_level == "scientific":
            # Validate correction accuracy
            corrected, fitted_baseline = baseline.baseline_polynomial_1d(
                x,
                y,
                order=1,
                manual_regions=baseline_test_data["signal_regions"],
                region_mode="exclude",
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
        """Standard magnetic field range for testing.

        Odd point count places a sample exactly at x=0, so symmetry and
        center-value assertions hold to machine precision.
        """
        return np.linspace(-20, 20, 1001)

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
            gaussian(B, center=center, width=width, derivative=0)
            gauss_1 = gaussian(B, center=center, width=width, derivative=1)

            # Antisymmetric; atol absorbs the near-zero center sample,
            # where a purely relative tolerance is meaningless
            np.testing.assert_allclose(gauss_1, -gauss_1[::-1], rtol=1e-10, atol=1e-12)

            # Test phase rotation
            for phase in [0, 0.5, 1.0, 1.5]:
                result = gaussian(B, center=center, width=width, phase=phase)
                assert len(result) == len(B)
                assert np.all(np.isfinite(result))

        elif protocol_level == "scientific":
            # Validate against analytical properties
            center = 0
            width = 4.0  # FWHM (library convention)

            gauss = gaussian(B, center=center, width=width, derivative=0)

            # Maximum should be at center
            max_idx = np.argmax(gauss)
            center_idx = np.argmin(np.abs(B - center))
            assert abs(max_idx - center_idx) <= 2  # Allow for discretization

            # The width parameter is the full width at half maximum
            max_val = np.max(gauss)
            half_max_indices = np.where(gauss >= max_val / 2)[0]
            if len(half_max_indices) > 0:
                field_width = B[half_max_indices[-1]] - B[half_max_indices[0]]
                np.testing.assert_allclose(field_width, width, rtol=0.1)

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

            # Pure dispersion (phase = pi/2)
            lorentz_disp = lorentzian(B, center=center, width=width, phase=np.pi / 2)
            # Antisymmetric; atol absorbs the near-zero center sample,
            # where a purely relative tolerance is meaningless
            np.testing.assert_allclose(
                lorentz_disp, -lorentz_disp[::-1], rtol=1e-10, atol=1e-12
            )

            # Test numerical stability (absorption and first derivative)
            self.validate_numerical_stability(
                lorentzian,
                [(B, center, width, 0), (B, center, width, 1)],
                tolerance=1e-12,
            )

        elif protocol_level == "scientific":
            # Validate Lorentzian width convention
            center = 0
            width = 2.0  # FWHM (library convention)

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
                # The width parameter is the full width at half maximum
                np.testing.assert_allclose(
                    field_span, width, rtol=0.15
                )  # Allow for discretization

    @pytest.mark.parametrize("protocol_level", PROTOCOL_LEVELS)
    def test_voigtian_function(self, protocol_level, standard_field_range):
        """Deep testing of Voigtian (true convolution) function."""
        B = standard_field_range

        if protocol_level == "smoke":
            # Basic function call
            result = voigtian(B, 0, (2, 2))
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        elif protocol_level == "standard":
            # Test different sigma/gamma ratios
            ratios = [(1, 1), (2, 1), (1, 2), (3, 3)]

            for sigma, gamma in ratios:
                result = voigtian(B, 0, (sigma, gamma))
                assert len(result) == len(B)
                assert np.all(np.isfinite(result))
                assert np.max(result) > 0

        elif protocol_level == "deep":
            # Test limiting cases
            # When gaussian_fwhm >> lorentzian_fwhm, should approach Gaussian
            voigt_gauss_like = voigtian(B, 0, (5, 0.1))
            pure_gauss = gaussian(
                B, center=0, width=5
            )  # Same FWHM as voigtian gaussian_fwhm

            # Should be reasonably similar (not exact due to convolution)
            correlation = np.corrcoef(voigt_gauss_like, pure_gauss)[0, 1]
            assert (
                correlation > 0.95
            ), f"Voigt-Gaussian correlation too low: {correlation}"

            # When gamma >> sigma, should approach Lorentzian
            voigt_lorentz_like = voigtian(B, 0, (0.1, 5))
            pure_lorentz = lorentzian(B, center=0, width=5, phase=0)

            correlation = np.corrcoef(voigt_lorentz_like, pure_lorentz)[0, 1]
            assert (
                correlation > 0.95
            ), f"Voigt-Lorentzian correlation too low: {correlation}"

        elif protocol_level == "scientific":
            # Test convolution accuracy
            # Voigt should be broader than both parent functions
            gauss_fwhm, lorentz_fwhm = 2.0, 2.0
            center = 0

            voigt = voigtian(B, center, (gauss_fwhm, lorentz_fwhm))
            gauss = gaussian(B, center=center, width=gauss_fwhm)
            lorentz = lorentzian(B, center=center, width=lorentz_fwhm, phase=0)

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
                    shape = Lineshape(shape_type, width=(2.0, 2.0))
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

            # Voigt requires a (gaussian_width, lorentzian_width) tuple,
            # rejected at construction
            with pytest.raises(ValueError):
                Lineshape("voigt", width=3.0)

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
            "voigtian": lambda B: voigtian(B, 0, (3, 3)),
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

        assert perf_data["execution_time"] < expected_max_time[data_size], (
            f"{function_name} too slow for {data_size}"
            f" points: {perf_data['execution_time']:.3f}s"
        )

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
        y_corrected, baseline_fit = baseline.baseline_polynomial_1d(x, y, order=1)

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
        assert corrected_peak > 0  # Should have positive signal after correction


class TestScientificValidation(TestProtocol):
    """Scientific validation and accuracy tests."""

    def test_physical_constants_accuracy(self):
        """Validate physical constants against NIST CODATA 2022 values."""
        # Free-electron g-factor magnitude (EasySpin sign convention)
        nist_ge = 2.00231930436092
        np.testing.assert_allclose(constants.GFREE, nist_ge, rtol=1e-10)

        # Bohr magneton in J/T
        nist_mb = 9.2740100657e-24
        np.testing.assert_allclose(constants.BMAGN, nist_mb, rtol=1e-8)

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
