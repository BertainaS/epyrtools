"""
Core functionality test suite for EPyR Tools release validation.

This module contains the essential tests that must pass for a stable release.
These tests validate the core user-facing functionality without requiring
the full test suite to pass.

Test Categories:
- Data loading (eprload)
- Basic plotting (eprplot)
- FAIR conversion
- Plugin system
"""

from pathlib import Path

import numpy as np
import pytest

# Core functionality imports that should always work
from epyr import eprload, plot_1d
from epyr.fair import convert_bruker_to_fair
from epyr.plugins import plugin_manager


class TestCoreDataLoading:
    """Test core data loading functionality."""

    def test_eprload_import(self):
        """Test that eprload can be imported."""
        from epyr import eprload

        assert callable(eprload)

    def test_eprload_file_detection(self):
        """Test file extension detection by trying actual extensions."""
        # This tests the internal logic by seeing what error we get
        from epyr.eprload import eprload

        # Test that different extensions are recognized as different formats
        # We expect different error messages for supported vs unsupported extensions
        try:
            eprload("test.DSC")
        except FileNotFoundError:
            pass  # Expected - file doesn't exist but extension is recognized
        except Exception as e:
            assert (
                "extension" not in str(e).lower()
            )  # Should not complain about extension

    def test_eprload_invalid_file_handling(self):
        """Test graceful handling of invalid files."""
        from epyr.eprload import eprload

        # Test non-existent file
        with pytest.raises(FileNotFoundError):
            eprload("nonexistent_file.DSC")

    def test_eprload_unsupported_extension(self):
        """Test handling of unsupported file extensions."""
        from epyr.eprload import eprload

        # Create a temporary file with unsupported extension
        temp_file = Path("test_file.xyz")
        temp_file.touch()
        try:
            # Test that unsupported extension is handled gracefully
            result = eprload(str(temp_file))
            # If no exception, result should be None or empty
            assert result is None or (isinstance(result, tuple) and len(result) == 4)
        except (ValueError, RuntimeError) as e:
            # Should get an error about unsupported extension
            assert "extension" in str(e).lower() or "unsupported" in str(e).lower()
        finally:
            temp_file.unlink()


class TestCorePlotting:
    """Test core plotting functionality."""

    def test_plot_1d_import(self):
        """Test that plot_1d can be imported."""
        from epyr import plot_1d

        assert callable(plot_1d)

    def test_plot_1d_basic(self):
        """Test basic 1D plotting functionality."""
        import matplotlib

        from epyr import plot_1d

        matplotlib.use("Agg")  # Non-interactive backend

        # Create test data
        x = np.linspace(0, 10, 100)
        y = np.sin(x)

        # Should not raise an exception
        fig, ax = plot_1d(x, y, title="Test Plot")
        assert fig is not None
        assert ax is not None


class TestCoreFAIR:
    """Test core FAIR conversion functionality."""

    def test_fair_import(self):
        """Test that FAIR module can be imported."""
        from epyr.fair import convert_bruker_to_fair

        assert callable(convert_bruker_to_fair)

    def test_fair_parameter_mapping(self):
        """Test that parameter mapping is available."""
        from epyr.fair.parameter_mapping import BRUKER_PARAM_MAP

        assert isinstance(BRUKER_PARAM_MAP, dict)
        assert len(BRUKER_PARAM_MAP) > 100  # Should have many parameters

    def test_fair_invalid_input_handling(self):
        """Test FAIR conversion with invalid input."""
        from epyr.fair import convert_bruker_to_fair

        # Test with non-existent file
        result = convert_bruker_to_fair(
            input_file="nonexistent.DSC", output_dir="./test_output", formats=["json"]
        )
        # Should return False for failed conversion
        assert result is False


class TestCorePlugins:
    """Test core plugin system functionality."""

    def test_plugin_manager_import(self):
        """Test that plugin manager can be imported."""
        from epyr.plugins import plugin_manager

        assert hasattr(plugin_manager, "register_plugin")
        assert hasattr(plugin_manager, "get_export_plugin")

    def test_plugin_discovery(self):
        """Test plugin discovery functionality."""
        from epyr.plugins import plugin_manager

        # Test that we can access plugin collections
        assert hasattr(plugin_manager, "export_plugins")
        assert hasattr(plugin_manager, "file_format_plugins")
        assert isinstance(plugin_manager.export_plugins, dict)


class TestCoreConfiguration:
    """Test core configuration functionality."""

    def test_config_import(self):
        """Test that config can be imported."""
        from epyr.config import config

        assert hasattr(config, "get")
        assert hasattr(config, "set")

    def test_config_basic_operations(self):
        """Test basic config operations."""
        from epyr.config import config

        # Test getting a default value
        value = config.get("plotting.dpi", default=300)
        assert isinstance(value, (int, float))

        # Test setting and getting a value
        config.set("test.value", 42)
        assert config.get("test.value") == 42


class TestCoreLineshapes:
    """Test core lineshape functionality."""

    def test_lineshape_imports(self):
        """Test that basic lineshapes can be imported."""
        from epyr.lineshapes import gaussian, lorentzian

        assert callable(gaussian)
        assert callable(lorentzian)

    def test_gaussian_basic(self):
        """Test basic Gaussian function."""
        from epyr.lineshapes import gaussian

        x = np.linspace(-5, 5, 100)
        y = gaussian(x, center=0, width=1)

        assert len(y) == len(x)
        assert np.isfinite(y).all()
        assert y[50] > y[0]  # Peak should be higher than edges


def run_core_tests():
    """Run the core test suite and return results."""
    import subprocess
    import sys

    # Run this specific test file
    result = subprocess.run(
        [sys.executable, "-m", "pytest", __file__, "-v", "--tb=short"],
        capture_output=True,
        text=True,
    )

    return {
        "return_code": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "success": result.returncode == 0,
    }


if __name__ == "__main__":
    # Allow running this test file directly
    pytest.main([__file__, "-v"])
