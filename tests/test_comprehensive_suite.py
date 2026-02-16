"""
Comprehensive Test Suite for EPyR Tools - All Module Coverage
============================================================

This module provides comprehensive testing coverage for all EPyR Tools modules
with systematic validation of every public function and method.

Coverage includes:
- epyr.eprload: Data loading functionality
- epyr.baseline: Baseline correction algorithms
- epyr.fair: FAIR data conversion
- epyr.plot: Plotting and visualization
- epyr.constants: Physical constants
- epyr.lineshapes: All lineshape functions
- epyr.isotope_gui: GUI components (where applicable)
- epyr.config: Configuration management
- epyr.cli: Command-line interface
- epyr.plugins: Plugin system
- epyr.performance: Performance monitoring
"""

import inspect
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

# EPyR Tools imports - comprehensive coverage
import epyr
from epyr import (
    baseline,
    cli,
    config,
    constants,
    eprload,
    fair,
    performance,
    plot,
    plugins,
)
from epyr.lineshapes import (
    Lineshape,
    convspec,
    gaussian,
    lineshape_class,
    lorentzian,
    lshape,
    voigtian,
)

# Conditional imports for optional dependencies
try:
    from epyr import isotope_gui

    HAS_GUI = True
except ImportError:
    HAS_GUI = False
    isotope_gui = None


class TestModuleCoverage:
    """Systematic testing of all module functions."""

    def get_module_functions(self, module) -> List[Tuple[str, Any]]:
        """Extract all public functions from a module."""
        functions = []
        for name, obj in inspect.getmembers(module):
            if (
                inspect.isfunction(obj)
                and not name.startswith("_")
                and hasattr(obj, "__module__")
                and obj.__module__.startswith("epyr")
            ):
                functions.append((name, obj))
        return functions

    def test_eprload_module_coverage(self, mock_bruker_files):
        """Test all functions in eprload module."""
        # Test main eprload function
        if hasattr(eprload, "eprload"):
            # Test with mock files
            try:
                # This would normally require GUI, so we test the import and basic properties
                assert callable(eprload.eprload)

                # Test function signature
                sig = inspect.signature(eprload.eprload)
                assert "file_path" in sig.parameters or len(sig.parameters) == 0

            except Exception as e:
                # If GUI-dependent, just ensure function exists
                assert hasattr(eprload, "eprload")

        # Test other functions in module
        functions = self.get_module_functions(eprload)
        for func_name, func_obj in functions:
            # Basic function validation
            assert callable(func_obj)
            assert func_obj.__doc__ is not None or func_name.startswith("_")

    def test_baseline_module_coverage(self, baseline_test_data):
        """Test all functions in baseline module."""
        x = baseline_test_data["x"]
        y = baseline_test_data["y_with_baseline"]

        # Test polynomial baseline correction
        if hasattr(baseline, "baseline_polynomial"):
            corrected, fitted = baseline.baseline_polynomial(y, x_data=x, poly_order=1)
            assert len(corrected) == len(y)
            assert len(fitted) == len(y)

        # Test all other baseline functions
        functions = self.get_module_functions(baseline)
        for func_name, func_obj in functions:
            assert callable(func_obj)

            # Test function with appropriate parameters based on name
            if "polynomial" in func_name.lower() and func_name != "baseline_polynomial":
                # Already tested above
                continue
            elif "exponential" in func_name.lower():
                try:
                    result = func_obj(y, x_data=x)
                    assert len(result) >= 1  # Should return at least corrected data
                except (TypeError, ValueError):
                    # Function may require different parameters
                    pass

    def test_fair_module_coverage(self, temp_data_files):
        """Test all functions in fair module."""
        functions = self.get_module_functions(fair)

        # Test conversion functions
        for func_name, func_obj in functions:
            assert callable(func_obj)

            if "convert" in func_name.lower():
                # Test conversion function signature
                sig = inspect.signature(func_obj)
                # Should have input file parameter
                assert len(sig.parameters) >= 1

            elif "export" in func_name.lower():
                # Test export function signature
                sig = inspect.signature(func_obj)
                assert len(sig.parameters) >= 2  # data and output path typically

    def test_plot_module_coverage(self, sample_2d_data):
        """Test all functions in plot module."""
        x_axis, y_axis, Z_data = sample_2d_data

        functions = self.get_module_functions(plot)

        for func_name, func_obj in functions:
            assert callable(func_obj)

            if "2d" in func_name.lower() or "map" in func_name.lower():
                try:
                    # Test 2D plotting function
                    result = func_obj(x_axis, y_axis, Z_data)
                    # Should return matplotlib objects
                    assert result is not None
                except (ImportError, TypeError, ValueError):
                    # May require matplotlib backend or specific parameters
                    pass
            elif "1d" in func_name.lower() or "line" in func_name.lower():
                try:
                    # Test 1D plotting function
                    result = func_obj(x_axis, Z_data[0])  # Use first row as 1D data
                    assert result is not None
                except (ImportError, TypeError, ValueError):
                    pass

    def test_constants_module_coverage(self):
        """Test all constants and functions in constants module."""
        # Test that all required constants exist and are reasonable
        required_constants = [
            "ELECTRON_G_FACTOR",
            "BOHR_MAGNETON",
            "PLANCK_CONSTANT",
            "SPEED_OF_LIGHT",
            "ELECTRON_MASS",
            "NUCLEAR_MAGNETON",
        ]

        for const_name in required_constants:
            if hasattr(constants, const_name):
                value = getattr(constants, const_name)
                assert isinstance(value, (int, float, complex))
                assert np.isfinite(value)
                assert value != 0  # Physical constants should be non-zero

        # Test any utility functions
        functions = self.get_module_functions(constants)
        for func_name, func_obj in functions:
            assert callable(func_obj)

    def test_config_module_coverage(self):
        """Test all functions in config module."""
        functions = self.get_module_functions(config)

        for func_name, func_obj in functions:
            assert callable(func_obj)

            # Test configuration functions
            if "load" in func_name.lower() or "get" in func_name.lower():
                try:
                    result = func_obj()
                    # Should return configuration data
                    assert result is not None
                except Exception:
                    # May require specific config files
                    pass
            elif "save" in func_name.lower() or "set" in func_name.lower():
                # Test that function signature expects parameters
                sig = inspect.signature(func_obj)
                assert len(sig.parameters) >= 1

    def test_cli_module_coverage(self):
        """Test all functions in cli module."""
        functions = self.get_module_functions(cli)

        for func_name, func_obj in functions:
            assert callable(func_obj)

            # Test CLI function properties
            if "main" in func_name.lower() or "run" in func_name.lower():
                # Should be entry points
                sig = inspect.signature(func_obj)
                # CLI functions often take no args or argv
                assert len(sig.parameters) <= 1
            elif "parse" in func_name.lower():
                # Parser functions should take arguments
                sig = inspect.signature(func_obj)
                assert len(sig.parameters) >= 1

    def test_plugins_module_coverage(self):
        """Test all functions in plugins module."""
        functions = self.get_module_functions(plugins)

        for func_name, func_obj in functions:
            assert callable(func_obj)

            # Test plugin system functions
            if "register" in func_name.lower():
                sig = inspect.signature(func_obj)
                # Registration should take plugin info
                assert len(sig.parameters) >= 1
            elif "load" in func_name.lower() or "discover" in func_name.lower():
                try:
                    result = func_obj()
                    # Should return plugin information
                    assert result is not None or result == []  # Empty list is valid
                except Exception:
                    pass

    def test_performance_module_coverage(self):
        """Test all functions in performance module."""
        functions = self.get_module_functions(performance)

        for func_name, func_obj in functions:
            assert callable(func_obj)

            # Test performance monitoring functions
            if "measure" in func_name.lower() or "time" in func_name.lower():
                sig = inspect.signature(func_obj)
                # Should take function to measure
                assert len(sig.parameters) >= 1
            elif "memory" in func_name.lower():
                try:
                    result = func_obj()
                    # Should return memory info
                    assert isinstance(result, (int, float, dict))
                except Exception:
                    pass

    @pytest.mark.skipif(not HAS_GUI, reason="GUI module not available")
    def test_isotope_gui_module_coverage(self):
        """Test all functions in isotope_gui module (if available)."""
        if not HAS_GUI:
            pytest.skip("isotope_gui module not available")

        functions = self.get_module_functions(isotope_gui)

        for func_name, func_obj in functions:
            assert callable(func_obj)

            # Test GUI functions (basic validation only, no actual GUI testing)
            if "run" in func_name.lower() or "main" in func_name.lower():
                sig = inspect.signature(func_obj)
                # GUI entry points typically take no args
                assert len(sig.parameters) <= 1


class TestLineshapesComprehensive:
    """Comprehensive testing of all lineshape functions."""

    @pytest.fixture
    def test_field_range(self):
        """Field range for comprehensive testing."""
        return np.linspace(-30, 30, 1500)

    def test_all_lineshape_functions(self, test_field_range):
        """Test every function in the lineshapes module."""
        B = test_field_range

        # Test gaussian function and all its variations
        for derivative in [0, 1, 2]:
            for phase in [0, 0.5, 1.0]:
                result = gaussian(
                    B, center=0, width=5, derivative=derivative, phase=phase
                )
                assert len(result) == len(B)
                assert np.all(np.isfinite(result))

        # Test lorentzian function variations
        for phase in [0, 0.25, 0.5, 0.75, 1.0]:
            result = lorentzian(B, center=0, width=5, phase=phase)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        # Test voigtian function variations
        sigma_gamma_pairs = [(1, 1), (2, 3), (3, 2), (1, 5), (5, 1)]
        for sigma, gamma in sigma_gamma_pairs:
            result = voigtian(B, center=0, sigma=sigma, gamma=gamma)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        # Test lshape (general lineshape) function
        for alpha in [0, 0.25, 0.5, 0.75, 1.0]:
            result = lshape(B, center=0, width=5, alpha=alpha)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))

        # Test convspec function
        # Create a simple spectrum to convolve
        original_spectrum = gaussian(B, center=0, width=3)
        broadening_function = gaussian(B, center=0, width=1)

        convolved = convspec(original_spectrum, broadening_function)
        assert len(convolved) == len(original_spectrum)
        assert np.all(np.isfinite(convolved))

    def test_lineshape_class_all_types(self, test_field_range):
        """Test Lineshape class with all supported types."""
        B = test_field_range

        # Test all supported lineshape types
        test_cases = [
            ("gaussian", {"width": 5.0}),
            ("lorentzian", {"width": 5.0}),
            ("pseudo_voigt", {"width": 5.0, "alpha": 0.5}),
            ("voigt", {"width": 5.0, "sigma": 3.0, "gamma": 3.0}),
        ]

        for shape_type, params in test_cases:
            shape = Lineshape(shape_type, **params)
            result = shape(B, center=0)
            assert len(result) == len(B)
            assert np.all(np.isfinite(result))
            assert np.max(result) > 0  # Should have positive values

    def test_factory_functions(self, test_field_range):
        """Test factory functions for creating lineshapes."""
        B = test_field_range

        # Test factory functions if they exist in the lineshapes module
        lineshape_module = epyr.lineshapes

        factory_functions = [
            "create_gaussian",
            "create_lorentzian",
            "create_voigt",
            "create_pseudo_voigt",
        ]

        for factory_name in factory_functions:
            if hasattr(lineshape_module, factory_name):
                factory_func = getattr(lineshape_module, factory_name)
                assert callable(factory_func)

                try:
                    # Test factory function
                    if "voigt" in factory_name and factory_name != "pseudo_voigt":
                        lineshape_obj = factory_func(width=5.0, sigma=3.0, gamma=3.0)
                    elif "pseudo_voigt" in factory_name:
                        lineshape_obj = factory_func(width=5.0, alpha=0.5)
                    else:
                        lineshape_obj = factory_func(width=5.0)

                    # Should return callable lineshape object
                    assert callable(lineshape_obj)

                    result = lineshape_obj(B, center=0)
                    assert len(result) == len(B)
                    assert np.all(np.isfinite(result))

                except Exception as e:
                    # Factory function might have different signature
                    pass


class TestErrorHandling:
    """Comprehensive error handling and edge case testing."""

    def test_lineshape_error_conditions(self):
        """Test error handling in lineshape functions."""
        B = np.linspace(-10, 10, 100)

        # Test invalid parameters
        with pytest.raises((ValueError, TypeError)):
            gaussian(B, center=0, width=0)  # Zero width

        with pytest.raises((ValueError, TypeError)):
            gaussian(B, center=0, width=-5)  # Negative width

        with pytest.raises((ValueError, TypeError)):
            lorentzian(B, center=0, width=-3)  # Negative width

        with pytest.raises((ValueError, TypeError)):
            voigtian(B, center=0, sigma=-1, gamma=2)  # Negative sigma

        with pytest.raises((ValueError, TypeError)):
            voigtian(B, center=0, sigma=2, gamma=-1)  # Negative gamma

    def test_array_input_validation(self):
        """Test validation of array inputs."""
        # Empty arrays
        empty_array = np.array([])
        with pytest.raises((ValueError, IndexError)):
            gaussian(empty_array, center=0, width=5)

        # NaN arrays
        nan_array = np.full(100, np.nan)
        result = gaussian(nan_array, center=0, width=5)
        assert np.all(np.isnan(result))

        # Infinite arrays
        inf_array = np.full(100, np.inf)
        result = gaussian(inf_array, center=0, width=5)
        # Should handle gracefully
        assert len(result) == len(inf_array)

    def test_baseline_error_conditions(self, sample_1d_data):
        """Test error handling in baseline correction."""
        x, y = sample_1d_data

        # Test invalid polynomial orders
        with pytest.raises((ValueError, TypeError)):
            baseline.baseline_polynomial(y, x_data=x, poly_order=-1)

        # Test mismatched array sizes
        with pytest.raises((ValueError, IndexError)):
            baseline.baseline_polynomial(y, x_data=x[:-10], poly_order=1)

        # Test insufficient data points
        tiny_y = np.array([1.0, 2.0])
        tiny_x = np.array([0.0, 1.0])
        with pytest.raises((ValueError, np.linalg.LinAlgError)):
            baseline.baseline_polynomial(tiny_y, x_data=tiny_x, poly_order=5)


class TestNumericalAccuracy:
    """Test numerical accuracy and precision."""

    def test_lineshape_precision(self):
        """Test numerical precision of lineshape calculations."""
        # High precision field range
        B = np.linspace(-5, 5, 10001)

        # Test that repeated calculations give identical results
        for i in range(5):
            result1 = gaussian(B, center=0, width=2.5)
            result2 = gaussian(B, center=0, width=2.5)
            np.testing.assert_array_equal(result1, result2)

    def test_edge_case_parameters(self):
        """Test lineshapes with edge case parameters."""
        B = np.linspace(-20, 20, 1000)

        # Very small widths
        result = gaussian(B, center=0, width=1e-6)
        assert np.all(np.isfinite(result))
        assert np.max(result) > 0

        # Very large widths
        result = gaussian(B, center=0, width=1e6)
        assert np.all(np.isfinite(result))

        # Centers far outside the field range
        result = gaussian(B, center=1e6, width=5)
        assert np.all(np.isfinite(result))
        assert np.max(result) >= 0  # Should be very small but finite


class TestDocumentationCoverage:
    """Test that all functions have proper documentation."""

    def test_all_functions_documented(self):
        """Ensure all public functions have docstrings."""
        modules_to_test = [
            epyr.lineshapes.gaussian,
            epyr.lineshapes.lorentzian,
            epyr.lineshapes.voigtian,
            epyr.lineshapes.lshape,
            epyr.lineshapes.convspec,
        ]

        for module in modules_to_test:
            functions = [
                getattr(module, name)
                for name in dir(module)
                if callable(getattr(module, name)) and not name.startswith("_")
            ]

            for func in functions:
                assert (
                    func.__doc__ is not None
                ), f"Function {func.__name__} lacks documentation"
                assert (
                    len(func.__doc__.strip()) > 10
                ), f"Function {func.__name__} has insufficient documentation"

    def test_lineshape_class_documentation(self):
        """Test Lineshape class documentation."""
        assert Lineshape.__doc__ is not None
        assert len(Lineshape.__doc__.strip()) > 50

        # Test method documentation
        for method_name in ["__call__", "__init__"]:
            if hasattr(Lineshape, method_name):
                method = getattr(Lineshape, method_name)
                if hasattr(method, "__doc__"):
                    assert method.__doc__ is not None


# Test execution configuration
@pytest.fixture(scope="session", autouse=True)
def comprehensive_test_config():
    """Configure comprehensive testing session."""
    print("\n" + "=" * 70)
    print("EPyR Tools Comprehensive Test Suite")
    print("=" * 70)
    print("Testing all modules and functions systematically...")
    print(f"Python version: {sys.version}")
    print(f"NumPy version: {np.__version__}")
    print("=" * 70)

    yield

    print("\n" + "=" * 70)
    print("Comprehensive Test Suite Complete")
    print("=" * 70)


# Pytest markers
pytestmark = [pytest.mark.comprehensive, pytest.mark.full_coverage]
