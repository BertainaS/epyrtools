"""Tests for plot module (specialized EPR plotting functionality).

plot_2d_map follows the eprload data convention: the first argument is
the axis data (a list of two 1D arrays for 2D datasets), the second is
the 2D signal array of shape (ny, nx), and the third is the params dict
used for axis labels and units.
"""

from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pytest

from epyr import eprplot
from epyr.eprplot import plot_2d_map


class TestEPRPlotting:
    """Test suite for EPR plotting functions."""

    def test_eprplot_module_exists(self):
        """Test that eprplot module has expected functions."""
        # Check that module has core plotting functions
        assert hasattr(eprplot, "plot_1d")
        assert hasattr(eprplot, "plot_2d_map")
        assert hasattr(eprplot, "plot_2d_waterfall")
        assert callable(eprplot.plot_1d)
        assert callable(eprplot.plot_2d_map)

    @patch("matplotlib.pyplot.show")
    def test_plot_2d_map_basic(self, mock_show, sample_2d_data):
        """Test basic 2D mapping functionality."""
        x_axis, y_axis, z_data = sample_2d_data

        # Test basic 2D plot
        fig, ax = plot_2d_map([x_axis, y_axis], z_data)

        # Check that figure and axis were created
        assert fig is not None
        assert ax is not None

        # Check that the plot has the expected properties
        assert ax.get_xlabel() != ""  # Should have labels
        assert ax.get_ylabel() != ""

        # Clean up
        plt.close(fig)

    @patch("matplotlib.pyplot.show")
    def test_plot_2d_map_with_units(self, mock_show, sample_2d_data):
        """Test 2D mapping with axis names and units from params."""
        x_axis, y_axis, z_data = sample_2d_data

        params = {
            "XAXIS_NAME": "Magnetic Field",
            "XAXIS_UNIT": "mT",
            "YAXIS_NAME": "Frequency",
            "YAXIS_UNIT": "GHz",
        }
        fig, ax = plot_2d_map([x_axis, y_axis], z_data, params)

        # Check labels include units
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()

        assert "mT" in xlabel or "Magnetic Field" in xlabel
        assert "GHz" in ylabel or "Frequency" in ylabel

        plt.close(fig)

    @patch("matplotlib.pyplot.show")
    def test_plot_2d_map_colormap_options(self, mock_show, sample_2d_data):
        """Test 2D mapping with different colormaps."""
        x_axis, y_axis, z_data = sample_2d_data

        colormaps = ["viridis", "plasma", "inferno", "coolwarm"]

        for cmap in colormaps:
            try:
                fig, ax = plot_2d_map([x_axis, y_axis], z_data, cmap=cmap)

                # Check that colorbar was created
                assert len(fig.axes) >= 2  # Main plot + colorbar

                plt.close(fig)

            except Exception as e:
                pytest.fail(f"Failed with colormap {cmap}: {e}")

    @patch("matplotlib.pyplot.show")
    def test_plot_2d_map_complex_data(self, mock_show, sample_2d_data):
        """Test 2D mapping with complex data (real part is plotted)."""
        x_axis, y_axis, z_real = sample_2d_data

        # Create complex data
        z_complex = z_real + 1j * z_real * 0.5

        fig, ax = plot_2d_map([x_axis, y_axis], z_complex)

        assert fig is not None
        assert ax is not None

        plt.close(fig)

    @patch("matplotlib.pyplot.show")
    def test_plot_2d_map_input_validation(self, mock_show):
        """Test input validation for 2D plotting."""
        x_axis = np.linspace(0, 10, 10)
        y_axis = np.linspace(0, 5, 5)
        z_data = np.random.random((5, 10))  # Correct shape (ny, nx)

        # Test correct input
        fig, ax = plot_2d_map([x_axis, y_axis], z_data)
        plt.close(fig)

        # 1D data is rejected
        with pytest.raises(ValueError, match="Expected 2D data"):
            plot_2d_map([x_axis, y_axis], np.random.random(10))

        # Mismatched axis lengths fall back to index coordinates
        z_other = np.random.random((3, 8))
        fig, ax = plot_2d_map([x_axis, y_axis], z_other)
        assert "Index" in ax.get_xlabel()
        plt.close(fig)

    def test_plot_module_integration(self, sample_2d_data):
        """Test integration of eprplot module with plotting functions."""
        x_axis, y_axis, z_data = sample_2d_data

        with patch("matplotlib.pyplot.show"):
            # Test that plotting works
            fig, ax = plot_2d_map([x_axis, y_axis], z_data)

            # Check figure properties
            figsize = fig.get_size_inches()
            assert len(figsize) == 2  # Should have width and height

            plt.close(fig)

    @patch("matplotlib.pyplot.show")
    def test_plot_2d_map_edge_cases(self, mock_show):
        """Test 2D plotting with edge cases."""
        # Very small data
        x_small = np.array([1, 2])
        y_small = np.array([1, 2])
        z_small = np.array([[1, 2], [3, 4]])

        fig, ax = plot_2d_map([x_small, y_small], z_small)
        assert fig is not None
        plt.close(fig)

        # Single point (degenerate case)
        x_single = np.array([1])
        y_single = np.array([1])
        z_single = np.array([[1]])

        try:
            fig, ax = plot_2d_map([x_single, y_single], z_single)
            plt.close(fig)
        except (ValueError, IndexError, TypeError):
            # Expected for degenerate cases
            pass

    @patch("matplotlib.pyplot.show")
    def test_plot_styling_options(self, mock_show, sample_2d_data):
        """Test title and color scale options."""
        x_axis, y_axis, z_data = sample_2d_data

        fig, ax = plot_2d_map(
            [x_axis, y_axis], z_data, title="Test map", vmin=0, vmax=10
        )

        assert ax.get_title() == "Test map"

        plt.close(fig)

    def test_plot_function_returns(self, sample_2d_data):
        """Test that plotting functions return expected objects."""
        x_axis, y_axis, z_data = sample_2d_data

        with patch("matplotlib.pyplot.show"):
            result = plot_2d_map([x_axis, y_axis], z_data)

            # Should return figure and axis objects
            assert len(result) == 2
            fig, ax = result

            # Check object types
            assert hasattr(fig, "savefig")  # Figure-like object
            assert hasattr(ax, "plot")  # Axis-like object

            plt.close(fig)

    @patch("matplotlib.pyplot.show")
    def test_plot_save_functionality(self, mock_show, sample_2d_data):
        """Test plot saving functionality."""
        import tempfile
        from pathlib import Path

        x_axis, y_axis, z_data = sample_2d_data

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp_path = Path(tmp.name)

        try:
            fig, ax = plot_2d_map([x_axis, y_axis], z_data)

            # Test saving
            fig.savefig(tmp_path, dpi=150, bbox_inches="tight")

            # Check that file was created
            assert tmp_path.exists()
            assert tmp_path.stat().st_size > 0

            plt.close(fig)

        finally:
            if tmp_path.exists():
                tmp_path.unlink()


class TestPlotUtilities:
    """Test utility functions for plotting."""

    def test_colormap_validation(self):
        """Test colormap validation if implemented."""
        valid_cmaps = ["viridis", "plasma", "inferno", "coolwarm", "RdYlBu"]

        for cmap in valid_cmaps:
            # Should not raise error for valid colormaps
            try:
                import matplotlib.cm as cm

                cm.get_cmap(cmap)
            except (ValueError, AttributeError):
                pytest.skip(f"Colormap {cmap} not available")

    def test_plot_data_preprocessing(self, sample_2d_data):
        """Test data preprocessing for plotting."""
        x_axis, y_axis, z_data = sample_2d_data

        # Test with NaN values
        z_with_nan = z_data.copy()
        z_with_nan[0, 0] = np.nan

        with patch("matplotlib.pyplot.show"):
            # pcolormesh masks NaN cells, so this must not raise
            fig, ax = plot_2d_map([x_axis, y_axis], z_with_nan)
            plt.close(fig)

    def test_axis_range_calculation(self):
        """Test automatic axis range calculation."""
        # Test with various data ranges
        ranges = [
            (np.linspace(-100, 100, 50), np.linspace(0, 1, 50)),
            (np.linspace(1e-6, 1e-3, 50), np.linspace(1e9, 1e10, 50)),
            (np.array([1, 2, 3]), np.array([10, 20, 30])),
        ]

        for x_range, y_range in ranges:
            z_test = np.random.random((len(y_range), len(x_range)))

            with patch("matplotlib.pyplot.show"):
                fig, ax = plot_2d_map([x_range, y_range], z_test)

                # Axis limits must cover the data range (pcolormesh
                # extends half a cell beyond the outermost points)
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()

                assert xlim[0] <= x_range.min() <= xlim[1]
                assert xlim[0] <= x_range.max() <= xlim[1]
                assert ylim[0] <= y_range.min() <= ylim[1]
                assert ylim[0] <= y_range.max() <= ylim[1]

                plt.close(fig)
