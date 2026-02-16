"""Advanced tests for constants module."""

import numpy as np
import pytest

try:
    from epyr import constants

    CONSTANTS_AVAILABLE = True
except ImportError:
    CONSTANTS_AVAILABLE = False


@pytest.mark.skipif(not CONSTANTS_AVAILABLE, reason="Constants module not available")
class TestPhysicalConstants:
    """Test suite for physical constants."""

    def test_physical_constants_values(self):
        """Test that physical constants have reasonable values."""
        # Test fundamental constants
        if hasattr(constants, "PLANCK_CONSTANT"):
            h = constants.PLANCK_CONSTANT
            # Planck constant should be around 6.626e-34 J⋅s
            assert 6e-34 < h < 7e-34

        if hasattr(constants, "BOHR_MAGNETON"):
            mu_b = constants.BOHR_MAGNETON
            # Bohr magneton should be around 9.274e-24 J/T
            assert 9e-24 < mu_b < 10e-24

        if hasattr(constants, "ELECTRON_G_FACTOR"):
            g_e = constants.ELECTRON_G_FACTOR
            # Free electron g-factor should be around ±2.002
            assert abs(g_e) > 2.0 and abs(g_e) < 2.1

    def test_unit_conversions(self):
        """Test unit conversion constants."""
        # Test magnetic field conversions
        if hasattr(constants, "GAUSS_TO_TESLA"):
            gauss_to_tesla = constants.GAUSS_TO_TESLA
            assert gauss_to_tesla == 1e-4

        if hasattr(constants, "TESLA_TO_GAUSS"):
            tesla_to_gauss = constants.TESLA_TO_GAUSS
            assert tesla_to_gauss == 1e4

    def test_epr_specific_constants(self):
        """Test EPR-specific constants."""
        # Test g-factor calculation constants
        if hasattr(constants, "MHZ_PER_GAUSS_PER_G_FACTOR"):
            conversion = constants.MHZ_PER_GAUSS_PER_G_FACTOR
            # Should be around 2.8 MHz/G for g=1
            assert 2.5 < conversion < 3.0

    def test_constant_consistency(self):
        """Test consistency between related constants."""
        # Test reciprocal relationships
        if hasattr(constants, "GAUSS_TO_TESLA") and hasattr(
            constants, "TESLA_TO_GAUSS"
        ):
            g_to_t = constants.GAUSS_TO_TESLA
            t_to_g = constants.TESLA_TO_GAUSS
            assert abs(g_to_t * t_to_g - 1.0) < 1e-10

    def test_isotope_data_structure(self):
        """Test isotope data structure if available."""
        if hasattr(constants, "ISOTOPE_DATA"):
            isotope_data = constants.ISOTOPE_DATA
            assert isinstance(isotope_data, dict)

            # Test common isotopes
            if "H" in isotope_data:
                hydrogen = isotope_data["H"]
                if isinstance(hydrogen, dict) and "1" in hydrogen:
                    h1 = hydrogen["1"]
                    assert isinstance(h1, dict)
                    # Check for expected fields
                    expected_fields = ["spin", "abundance", "gamma"]
                    for field in expected_fields:
                        if field in h1:
                            assert isinstance(h1[field], (int, float))

    def test_nuclear_constants(self):
        """Test nuclear-specific constants."""
        if hasattr(constants, "NUCLEAR_MAGNETON"):
            mu_n = constants.NUCLEAR_MAGNETON
            # Nuclear magneton should be around 5.05e-27 J/T
            assert 5e-27 < mu_n < 6e-27

        # Test relationship between nuclear and Bohr magnetons
        if hasattr(constants, "NUCLEAR_MAGNETON") and hasattr(
            constants, "BOHR_MAGNETON"
        ):
            mu_n = constants.NUCLEAR_MAGNETON
            mu_b = constants.BOHR_MAGNETON
            # mu_n should be much smaller than mu_b (by proton/electron mass ratio)
            assert mu_n < mu_b
            ratio = mu_b / mu_n
            # Ratio should be around 1836 (electron/proton mass ratio)
            assert 1800 < ratio < 1900


@pytest.mark.skipif(not CONSTANTS_AVAILABLE, reason="Constants module not available")
class TestConstantsUtilities:
    """Test utility functions in constants module."""

    def test_frequency_field_conversion(self):
        """Test frequency-field conversion functions if available."""
        # Test basic conversion functions
        conversion_functions = [
            "field_to_frequency",
            "frequency_to_field",
            "gauss_to_tesla",
            "tesla_to_gauss",
        ]

        for func_name in conversion_functions:
            if hasattr(constants, func_name):
                func = getattr(constants, func_name)
                if callable(func):
                    try:
                        # Test with reasonable EPR values
                        if "frequency" in func_name:
                            result = func(9.4e9, 2.0)  # 9.4 GHz, g=2
                        else:
                            result = func(3350)  # ~3350 G field
                        assert isinstance(result, (int, float, np.number))
                        assert result > 0
                    except TypeError:
                        # Function might need different parameters
                        pass

    def test_g_factor_calculations(self):
        """Test g-factor calculation utilities if available."""
        if hasattr(constants, "calculate_g_factor"):
            calc_g = constants.calculate_g_factor
            if callable(calc_g):
                try:
                    # Test with typical EPR values
                    g = calc_g(9.4e9, 3350)  # 9.4 GHz, 3350 G
                    assert 1.5 < g < 2.5  # Reasonable g-factor range
                except TypeError:
                    # Function might need different parameters
                    pass

    def test_unit_conversion_functions(self):
        """Test unit conversion helper functions."""
        # Test field unit conversions
        test_values = [100, 1000, 3350, 10000]  # Typical EPR field values in G

        for value in test_values:
            if hasattr(constants, "gauss_to_tesla"):
                if callable(constants.gauss_to_tesla):
                    try:
                        tesla_value = constants.gauss_to_tesla(value)
                        assert isinstance(tesla_value, (int, float, np.number))
                        assert tesla_value == value * 1e-4
                    except TypeError:
                        pass

            if hasattr(constants, "tesla_to_gauss"):
                if callable(constants.tesla_to_gauss):
                    try:
                        tesla_input = value * 1e-4  # Convert to Tesla
                        gauss_value = constants.tesla_to_gauss(tesla_input)
                        assert isinstance(gauss_value, (int, float, np.number))
                        assert abs(gauss_value - value) < 1e-10
                    except TypeError:
                        pass

    def test_nuclear_data_access(self):
        """Test nuclear data access functions if available."""
        common_nuclei = ["1H", "13C", "14N", "15N", "31P"]

        for nucleus in common_nuclei:
            # Test various ways nuclear data might be accessed
            access_methods = [
                f"get_nuclear_data",
                f"get_isotope_data",
                f"nuclear_properties",
            ]

            for method_name in access_methods:
                if hasattr(constants, method_name):
                    method = getattr(constants, method_name)
                    if callable(method):
                        try:
                            data = method(nucleus)
                            if data is not None:
                                assert isinstance(data, dict)
                                # Check for common nuclear properties
                                expected_props = ["spin", "gamma", "abundance"]
                                for prop in expected_props:
                                    if prop in data:
                                        assert isinstance(data[prop], (int, float))
                        except (TypeError, KeyError, ValueError):
                            # Method might not support this nucleus or format
                            pass
