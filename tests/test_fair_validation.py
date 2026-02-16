"""
Tests for FAIR validation functionality.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from epyr.fair.validation import (
    ValidationResult,
    create_validation_report,
    validate_data_integrity,
    validate_epr_parameters,
    validate_fair_dataset,
    validate_fair_metadata,
    validate_file_format,
)


class TestValidationResult:
    """Test ValidationResult class."""

    def test_validation_result_initialization(self):
        """Test ValidationResult initialization."""
        result = ValidationResult()

        assert result.is_valid is True
        assert result.errors == []
        assert result.warnings == []
        assert result.info == []

    def test_add_error(self):
        """Test adding error messages."""
        result = ValidationResult()

        result.add_error("Test error")

        assert result.is_valid is False
        assert "Test error" in result.errors
        assert len(result.errors) == 1

    def test_add_warning(self):
        """Test adding warning messages."""
        result = ValidationResult()

        result.add_warning("Test warning")

        assert result.is_valid is True  # Warnings don't invalidate
        assert "Test warning" in result.warnings
        assert len(result.warnings) == 1

    def test_add_info(self):
        """Test adding info messages."""
        result = ValidationResult()

        result.add_info("Test info")

        assert result.is_valid is True
        assert "Test info" in result.info
        assert len(result.info) == 1

    def test_get_summary(self):
        """Test validation summary generation."""
        result = ValidationResult()

        result.add_error("Error 1")
        result.add_warning("Warning 1")
        result.add_info("Info 1")

        summary = result.get_summary()

        assert summary["valid"] is False
        assert summary["error_count"] == 1
        assert summary["warning_count"] == 1
        assert "Error 1" in summary["errors"]
        assert "Warning 1" in summary["warnings"]
        assert "Info 1" in summary["info"]


class TestValidateFairMetadata:
    """Test FAIR metadata validation."""

    def test_valid_complete_metadata(self):
        """Test validation with complete valid metadata."""
        metadata = {
            "title": "EPR Spectrum of Sample X",
            "description": "Continuous wave EPR measurement at X-band",
            "creator": "Dr. Jane Smith",
            "date_created": "2023-09-06T10:30:00Z",
            "identifier": "EPR-2023-001",
            "format": "Bruker BES3T",
            "license": "CC-BY-4.0",
            "instrument": "Bruker EMXnano",
            "measurement_parameters": {"frequency": 9.4e9},
            "sample_information": "DPPH radical standard",
            "processing_history": "Baseline corrected",
            "units": {"x_axis": "mT", "y_axis": "a.u."},
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is True
        assert len(result.errors) == 0
        # May have some info messages but no errors/warnings for complete metadata

    def test_missing_required_fields(self):
        """Test validation with missing required fields."""
        metadata = {
            "title": "EPR Spectrum",
            # Missing: description, creator, date_created, identifier, format, license
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is False
        assert len(result.errors) >= 6  # At least 6 missing required fields

        error_text = " ".join(result.errors)
        assert "description" in error_text
        assert "creator" in error_text
        assert "date_created" in error_text

    def test_missing_recommended_fields(self):
        """Test validation with missing recommended fields."""
        metadata = {
            # All required fields present
            "title": "EPR Spectrum",
            "description": "Test measurement",
            "creator": "Test User",
            "date_created": "2023-09-06T10:30:00Z",
            "identifier": "TEST-001",
            "format": "Bruker",
            "license": "CC-BY-4.0",
            # Missing recommended fields: instrument, measurement_parameters, etc.
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is True  # Valid but incomplete
        assert len(result.warnings) > 0

        warning_text = " ".join(result.warnings)
        assert "instrument" in warning_text
        assert "measurement_parameters" in warning_text

    def test_invalid_date_format(self):
        """Test validation with invalid date format."""
        metadata = {
            "title": "EPR Spectrum",
            "description": "Test",
            "creator": "Test User",
            "date_created": "2023-09-06",  # Missing time part
            "identifier": "TEST-001",
            "format": "Bruker",
            "license": "CC-BY-4.0",
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is True  # Just a warning for date format
        warning_text = " ".join(result.warnings)
        assert "ISO 8601" in warning_text

    def test_short_identifier(self):
        """Test validation with very short identifier."""
        metadata = {
            "title": "EPR Spectrum",
            "description": "Test",
            "creator": "Test User",
            "date_created": "2023-09-06T10:30:00Z",
            "identifier": "123",  # Very short
            "format": "Bruker",
            "license": "CC-BY-4.0",
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is True
        warning_text = " ".join(result.warnings)
        assert "more descriptive" in warning_text

    def test_invalid_units_format(self):
        """Test validation with invalid units format."""
        metadata = {
            "title": "EPR Spectrum",
            "description": "Test",
            "creator": "Test User",
            "date_created": "2023-09-06T10:30:00Z",
            "identifier": "TEST-001",
            "format": "Bruker",
            "license": "CC-BY-4.0",
            "units": "mT",  # Should be dict, not string
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is True
        warning_text = " ".join(result.warnings)
        assert "dictionary" in warning_text

    def test_missing_units_axes(self):
        """Test validation with incomplete units specification."""
        metadata = {
            "title": "EPR Spectrum",
            "description": "Test",
            "creator": "Test User",
            "date_created": "2023-09-06T10:30:00Z",
            "identifier": "TEST-001",
            "format": "Bruker",
            "license": "CC-BY-4.0",
            "units": {"x_axis": "mT"},  # Missing y_axis
        }

        result = validate_fair_metadata(metadata)

        assert result.is_valid is True
        warning_text = " ".join(result.warnings)
        assert "y-axis units" in warning_text


class TestValidateDataIntegrity:
    """Test data integrity validation."""

    def test_valid_data(self):
        """Test validation with valid data."""
        x_data = np.linspace(3400, 3500, 100)
        y_data = np.sin(x_data)
        metadata = {"data_points": 100}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is True
        assert len(result.errors) == 0

    def test_empty_y_data(self):
        """Test validation with empty y-data."""
        x_data = np.array([])
        y_data = np.array([])
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is False
        assert "Y-data is empty" in " ".join(result.errors)

    def test_none_y_data(self):
        """Test validation with None y-data."""
        result = validate_data_integrity(None, None, {})

        assert result.is_valid is False
        assert "Y-data is empty or None" in " ".join(result.errors)

    def test_nan_values(self):
        """Test validation with NaN values."""
        x_data = np.array([1, 2, np.nan, 4, 5])
        y_data = np.array([1, np.nan, 3, 4, 5])
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is False
        error_text = " ".join(result.errors)
        assert "Y-data contains NaN" in error_text

        warning_text = " ".join(result.warnings)
        assert "X-data contains NaN" in warning_text

    def test_infinite_values(self):
        """Test validation with infinite values."""
        x_data = np.array([1, 2, np.inf, 4, 5])
        y_data = np.array([1, 2, 3, np.inf, 5])
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is False
        error_text = " ".join(result.errors)
        assert "Y-data contains infinite" in error_text

        warning_text = " ".join(result.warnings)
        assert "X-data contains infinite" in warning_text

    def test_mismatched_lengths(self):
        """Test validation with mismatched x/y data lengths."""
        x_data = np.array([1, 2, 3])
        y_data = np.array([1, 2, 3, 4, 5])  # Different length
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is False
        error_text = " ".join(result.errors)
        assert "X-data length" in error_text and "Y-data length" in error_text

    def test_non_monotonic_x_data(self):
        """Test validation with non-monotonic x-data."""
        x_data = np.array([1, 3, 2, 4, 5])  # Not monotonic
        y_data = np.array([1, 2, 3, 4, 5])
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is True  # Warning, not error
        warning_text = " ".join(result.warnings)
        assert "not monotonic" in warning_text

    def test_zero_range_data(self):
        """Test validation with constant y-data (zero range)."""
        x_data = np.array([1, 2, 3, 4, 5])
        y_data = np.array([1, 1, 1, 1, 1])  # Constant values
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is True
        warning_text = " ".join(result.warnings)
        assert "zero range" in warning_text

    def test_few_data_points(self):
        """Test validation with very few data points."""
        x_data = np.array([1, 2])
        y_data = np.array([3, 4])
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is True
        warning_text = " ".join(result.warnings)
        assert "Very few data points" in warning_text

    def test_large_dataset(self):
        """Test validation with large dataset."""
        x_data = np.linspace(0, 1000, 200000)  # > 100k points
        y_data = np.random.randn(200000)
        metadata = {}

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is True
        info_text = " ".join(result.info)
        assert "Large dataset" in info_text

    def test_metadata_consistency(self):
        """Test validation with metadata consistency checks."""
        x_data = np.array([1, 2, 3, 4, 5])
        y_data = np.array([1, 2, 3, 4, 5])
        metadata = {"data_points": 10}  # Wrong count

        result = validate_data_integrity(x_data, y_data, metadata)

        assert result.is_valid is True  # Warning, not error
        warning_text = " ".join(result.warnings)
        assert "Data points mismatch" in warning_text


class TestValidateEprParameters:
    """Test EPR parameter validation."""

    def test_complete_epr_parameters(self):
        """Test validation with complete EPR parameters."""
        metadata = {
            "measurement_parameters": {
                "microwave_frequency": 9.4e9,
                "magnetic_field_range": [3300, 3600],
                "modulation_frequency": 100000,
                "modulation_amplitude": 1.0,
                "temperature": 298,
            }
        }

        result = validate_epr_parameters(metadata)

        assert result.is_valid is True
        # Should find all required parameters

    def test_missing_epr_parameters(self):
        """Test validation with missing EPR parameters."""
        metadata = {
            "measurement_parameters": {
                "temperature": 298
                # Missing: microwave_frequency, field_range, etc.
            }
        }

        result = validate_epr_parameters(metadata)

        assert result.is_valid is True  # Warnings, not errors
        warning_text = " ".join(result.warnings)
        assert "Missing EPR parameter" in warning_text

    def test_parameter_variations(self):
        """Test validation with parameter name variations."""
        metadata = {
            "measurement_parameters": {
                "mw_freq": 9.4e9,  # Variation of microwave_frequency
                "b_range": [3300, 3600],  # Variation of magnetic_field_range
                "ModFreq": 100000,  # Variation of modulation_frequency
                "RMA": 1.0,  # Variation of modulation_amplitude
            }
        }

        result = validate_epr_parameters(metadata)

        assert result.is_valid is True
        info_text = " ".join(result.info)
        # Should recognize parameter variations
        assert "Found parameter variation" in info_text

    def test_unusual_frequency_values(self):
        """Test validation with unusual frequency values."""
        metadata = {
            "measurement_parameters": {
                "microwave_frequency": 2000,  # Very high frequency
                "temperature": -50,  # Below absolute zero
            }
        }

        result = validate_epr_parameters(metadata)

        assert result.is_valid is True
        warning_text = " ".join(result.warnings)
        assert "Unusual microwave frequency" in warning_text
        assert "Unusual temperature" in warning_text

    def test_nested_parameter_sources(self):
        """Test validation with parameters in different metadata locations."""
        metadata = {
            "microwave_frequency": 9.4e9,  # Top level
            "instrument_parameters": {"magnetic_field_range": [3300, 3600]},  # Nested
            "measurement_parameters": {
                "modulation_frequency": 100000  # Different nested
            },
        }

        result = validate_epr_parameters(metadata)

        assert result.is_valid is True
        # Should find parameters in different locations


class TestValidateFileFormat:
    """Test file format validation."""

    def test_nonexistent_file(self):
        """Test validation of non-existent file."""
        result = validate_file_format(Path("/nonexistent/file.csv"))

        assert result.is_valid is False
        assert "File does not exist" in " ".join(result.errors)

    def test_empty_file(self):
        """Test validation of empty file."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            empty_file = Path(f.name)

        try:
            result = validate_file_format(empty_file)

            assert result.is_valid is False
            assert "File is empty" in " ".join(result.errors)
        finally:
            empty_file.unlink()

    def test_valid_json_file(self):
        """Test validation of valid JSON file."""
        test_data = {"key": "value", "number": 42}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(test_data, f)
            json_file = Path(f.name)

        try:
            result = validate_file_format(json_file, "json")

            assert result.is_valid is True
            info_text = " ".join(result.info)
            assert "JSON format validation passed" in info_text
        finally:
            json_file.unlink()

    def test_invalid_json_file(self):
        """Test validation of invalid JSON file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write("invalid json content {")
            invalid_json_file = Path(f.name)

        try:
            result = validate_file_format(invalid_json_file, "json")

            assert result.is_valid is False
            error_text = " ".join(result.errors)
            assert "Format validation failed" in error_text
        finally:
            invalid_json_file.unlink()

    def test_valid_csv_file(self):
        """Test validation of valid CSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("field,intensity\n3400,1.5\n3401,1.8\n")
            csv_file = Path(f.name)

        try:
            result = validate_file_format(csv_file, "csv")

            assert result.is_valid is True
            info_text = " ".join(result.info)
            assert "CSV format validation passed" in info_text
        finally:
            csv_file.unlink()

    def test_poorly_delimited_csv(self):
        """Test validation of poorly delimited CSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(
                "field intensity\n3400 1.5\n3401 1.8\n"
            )  # Space-separated, not comma
            csv_file = Path(f.name)

        try:
            result = validate_file_format(csv_file, "csv")

            assert result.is_valid is True  # Warning, not error
            warning_text = " ".join(result.warnings)
            assert "not be properly delimited" in warning_text
        finally:
            csv_file.unlink()

    def test_format_extension_mismatch(self):
        """Test validation with format/extension mismatch."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("some text content")
            text_file = Path(f.name)

        try:
            result = validate_file_format(text_file, "csv")

            assert result.is_valid is True  # Warning, not error
            warning_text = " ".join(result.warnings)
            assert "Expected CSV format" in warning_text
        finally:
            text_file.unlink()

    def test_large_file_info(self):
        """Test validation with large file."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            # Write > 100 MB of data
            large_data = b"x" * (101 * 1024 * 1024)
            f.write(large_data)
            large_file = Path(f.name)

        try:
            result = validate_file_format(large_file)

            assert result.is_valid is True
            info_text = " ".join(result.info)
            assert "Large file size" in info_text
        finally:
            large_file.unlink()


class TestValidateFairDataset:
    """Test comprehensive FAIR dataset validation."""

    def test_complete_valid_dataset(self):
        """Test validation of complete valid dataset."""
        x_data = np.linspace(3400, 3500, 100)
        y_data = np.sin(x_data)
        metadata = {
            "title": "Test EPR Spectrum",
            "description": "Test measurement",
            "creator": "Test User",
            "date_created": "2023-09-06T10:30:00Z",
            "identifier": "TEST-001",
            "format": "Bruker BES3T",
            "license": "CC-BY-4.0",
            "instrument": "Test Spectrometer",
            "measurement_parameters": {
                "microwave_frequency": 9.4e9,
                "magnetic_field_range": [3400, 3500],
            },
            "units": {"x_axis": "mT", "y_axis": "a.u."},
        }

        data_dict = {"x_data": x_data, "y_data": y_data, "metadata": metadata}

        result = validate_fair_dataset(data_dict)

        assert result.is_valid is True

    def test_dataset_with_multiple_issues(self):
        """Test validation of dataset with multiple issues."""
        x_data = np.array([1, 2, np.nan])  # Contains NaN
        y_data = np.array([1, 2])  # Wrong length
        metadata = {
            "title": "Test Spectrum",
            # Missing many required fields
        }

        data_dict = {"x_data": x_data, "y_data": y_data, "metadata": metadata}

        result = validate_fair_dataset(data_dict)

        assert result.is_valid is False
        assert len(result.errors) > 1  # Multiple errors

        # Should have errors from different validation types
        error_text = " ".join(result.errors)
        assert "Metadata" in error_text
        assert "Data Integrity" in error_text

    def test_dataset_with_file_validation(self):
        """Test validation including file format validation."""
        x_data = np.array([1, 2, 3])
        y_data = np.array([4, 5, 6])
        metadata = {"title": "Test"}

        data_dict = {"x_data": x_data, "y_data": y_data, "metadata": metadata}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"test": "data"}, f)
            test_file = Path(f.name)

        try:
            result = validate_fair_dataset(data_dict, test_file)

            # Should include file format validation
            combined_text = " ".join(result.errors + result.warnings + result.info)
            assert "File Format" in combined_text

        finally:
            test_file.unlink()


class TestCreateValidationReport:
    """Test validation report generation."""

    def test_passed_report(self):
        """Test report generation for passed validation."""
        result = ValidationResult()
        result.add_info("All checks passed")

        report = create_validation_report(result)

        assert "✓ PASSED" in report
        assert "Errors: 0" in report
        assert "All checks passed" in report

    def test_failed_report(self):
        """Test report generation for failed validation."""
        result = ValidationResult()
        result.add_error("Critical error 1")
        result.add_error("Critical error 2")
        result.add_warning("Minor warning")
        result.add_info("Some info")

        report = create_validation_report(result)

        assert "✗ FAILED" in report
        assert "Errors: 2" in report
        assert "Warnings: 1" in report
        assert "Critical error 1" in report
        assert "Critical error 2" in report
        assert "Minor warning" in report
        assert "Some info" in report

    def test_report_save_to_file(self):
        """Test saving report to file."""
        result = ValidationResult()
        result.add_error("Test error")
        result.add_warning("Test warning")

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            report_file = Path(f.name)

        try:
            report = create_validation_report(result, report_file)

            # File should be created and contain report
            assert report_file.exists()
            file_content = report_file.read_text()
            assert "Test error" in file_content
            assert "Test warning" in file_content
            assert file_content == report

        finally:
            if report_file.exists():
                report_file.unlink()
