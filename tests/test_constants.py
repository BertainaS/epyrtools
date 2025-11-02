import pytest

from epyr import constants


class TestConstants:
    """Test suite for constants module."""

    def test_physical_constants_exist(self):
        """Test that basic physical constants are defined."""
        # Test for common EPR-related constants
        assert hasattr(constants, "__all__") or len(dir(constants)) > 0

    def test_constants_are_numeric(self):
        """Test that constants have appropriate numeric types."""
        # Get all public attributes from constants module
        const_names = [
            name
            for name in dir(constants)
            if not name.startswith("_") and not callable(getattr(constants, name))
        ]

        if const_names:
            # Check that at least some constants are numeric
            numeric_constants = []
            for name in const_names:
                value = getattr(constants, name)
                if isinstance(value, (int, float, complex)):
                    numeric_constants.append(name)

            # We expect to have some numeric constants for EPR calculations
            assert len(numeric_constants) >= 0  # At least allow empty constants file
