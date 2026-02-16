from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from epyr.eprload import eprload


class TestEPRLoad:
    """Test suite for eprload module."""

    def test_eprload_with_invalid_file(self):
        """Test that eprload handles invalid file paths gracefully."""
        with pytest.raises(FileNotFoundError):
            eprload("nonexistent_file.dta", plot_if_possible=False)

    def test_eprload_plotting_disabled(self, temp_data_files):
        """Test that plotting can be disabled."""
        test_file = temp_data_files["test.dta"]

        # Mock the loading functions to return dummy data
        with patch("epyr.sub.loadBES3T.load") as mock_load:
            mock_load.return_value = (
                np.array([1, 2, 3]),  # y data
                np.array([0, 1, 2]),  # x data
                {"MWFQ": 9.4e9},  # params
            )

            with patch("matplotlib.pyplot.show") as mock_show:
                result = eprload(str(test_file), plot_if_possible=False)
                mock_show.assert_not_called()
                assert result[0] is not None  # x data
                assert result[1] is not None  # y data

    def test_scaling_parameter_validation(self, temp_data_files):
        """Test that scaling parameter is handled correctly."""
        test_file = temp_data_files["test.dta"]

        # Test with invalid scaling characters
        with pytest.raises(
            ValueError, match="Scaling string contains invalid characters"
        ):
            eprload(str(test_file), scaling="XYZ", plot_if_possible=False)

        # Test with valid scaling
        valid_scaling = "nPGTc"
        with patch("epyr.sub.loadBES3T.load") as mock_load:
            mock_load.return_value = (
                np.array([1, 2, 3]),
                np.array([0, 1, 2]),
                {"MWFQ": 9.4e9},
            )
            result = eprload(
                str(test_file), scaling=valid_scaling, plot_if_possible=False
            )
            assert result[0] is not None

    def test_file_extension_detection(self):
        """Test file extension detection logic."""
        test_files = ["data.dta", "data.dsc", "data.spc", "data.par"]

        for filename in test_files:
            # Each should raise FileNotFoundError (file doesn't exist)
            # but NOT ValueError (extension is valid)
            with pytest.raises(FileNotFoundError):
                eprload(filename, plot_if_possible=False)

    def test_unsupported_extension(self, temp_data_files):
        """Test that unsupported file extensions raise ValueError."""
        # Create a file with unsupported extension
        test_file = Path(str(temp_data_files["test.dta"])).parent / "data.txt"
        test_file.write_text("dummy content")

        with pytest.raises(ValueError, match="Unsupported file extension"):
            eprload(str(test_file), plot_if_possible=False)

    def test_eprload_with_mock_bruker_files(self, mock_bruker_files):
        """Test loading with properly mocked Bruker files."""
        base_path = mock_bruker_files["base_path"]

        # Mock the actual loading functions with correct import path
        with patch("epyr.sub.loadBES3T.load") as mock_load:
            mock_load.return_value = (
                np.array([1.0, 2.0, 3.0, 4.0]),  # y data
                np.array([320.0, 320.01, 320.02, 320.03]),  # x data
                {"MWFQ": 9.4e9, "MWPW": 20.0},  # parameters
            )

            result = eprload(str(base_path) + ".dsc", plot_if_possible=False)

            # Check that we got valid results
            x, y, params, file_path = result
            assert x is not None
            assert y is not None
            assert params is not None
            assert file_path is not None
            assert len(y) == 4
            assert params["MWFQ"] == 9.4e9

    def test_return_type_parameter(self, temp_data_files):
        """Test return_type parameter with real, imag, and default options."""
        test_file = temp_data_files["test.dta"]

        # Create complex data for testing
        complex_y = np.array([1 + 2j, 3 + 4j, 5 + 6j])
        real_x = np.array([0, 1, 2])

        with patch("epyr.sub.loadBES3T.load") as mock_load:
            mock_load.return_value = (complex_y, real_x, {"MWFQ": 9.4e9})

            # Test default (should return complex data as-is)
            x, y, params, file_path = eprload(
                str(test_file), plot_if_possible=False, return_type="default"
            )
            assert np.iscomplexobj(y)
            assert np.array_equal(y, complex_y)

            # Test real (should return only real part)
            mock_load.return_value = (complex_y.copy(), real_x, {"MWFQ": 9.4e9})
            x, y, params, file_path = eprload(
                str(test_file), plot_if_possible=False, return_type="real"
            )
            assert np.isrealobj(y)
            assert np.array_equal(y, np.real(complex_y))

            # Test imag (should return only imaginary part)
            mock_load.return_value = (complex_y.copy(), real_x, {"MWFQ": 9.4e9})
            x, y, params, file_path = eprload(
                str(test_file), plot_if_possible=False, return_type="imag"
            )
            assert np.isrealobj(y)
            assert np.array_equal(y, np.imag(complex_y))

            # Test invalid return_type
            mock_load.return_value = (complex_y.copy(), real_x, {"MWFQ": 9.4e9})
            with pytest.raises(ValueError, match="Invalid return_type"):
                eprload(str(test_file), plot_if_possible=False, return_type="invalid")

    def test_return_type_with_real_data(self, temp_data_files):
        """Test that return_type works correctly with already-real data."""
        test_file = temp_data_files["test.dta"]

        # Real data should remain real regardless of return_type
        real_y = np.array([1.0, 2.0, 3.0])
        real_x = np.array([0, 1, 2])

        with patch("epyr.sub.loadBES3T.load") as mock_load:
            # Test with real data and return_type="real"
            mock_load.return_value = (real_y.copy(), real_x, {"MWFQ": 9.4e9})
            x, y, params, file_path = eprload(
                str(test_file), plot_if_possible=False, return_type="real"
            )
            assert np.isrealobj(y)
            assert np.array_equal(y, real_y)

            # Test with real data and return_type="imag" (should give zeros)
            mock_load.return_value = (real_y.copy(), real_x, {"MWFQ": 9.4e9})
            x, y, params, file_path = eprload(
                str(test_file), plot_if_possible=False, return_type="imag"
            )
            assert np.isrealobj(y)
            assert np.allclose(y, np.zeros_like(real_y))
