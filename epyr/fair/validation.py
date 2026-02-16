"""
Data Validation Module for FAIR Compliance
==========================================

This module provides validation functions to ensure data meets FAIR
(Findable, Accessible, Interoperable, Reusable) principles and
scientific data standards.

Features:
- Metadata completeness validation
- Data integrity checks
- Format compliance verification
- Scientific metadata standards validation
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from ..config import config
from ..logging_config import get_logger

logger = get_logger(__name__)


class ValidationResult:
    """Container for validation results."""

    def __init__(self):
        self.is_valid = True
        self.warnings: List[str] = []
        self.errors: List[str] = []
        self.info: List[str] = []

    def add_error(self, message: str):
        """Add an error message."""
        self.errors.append(message)
        self.is_valid = False

    def add_warning(self, message: str):
        """Add a warning message."""
        self.warnings.append(message)

    def add_info(self, message: str):
        """Add an info message."""
        self.info.append(message)

    def get_summary(self) -> Dict[str, Any]:
        """Get validation summary."""
        return {
            "valid": self.is_valid,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "errors": self.errors,
            "warnings": self.warnings,
            "info": self.info,
        }


def validate_fair_metadata(metadata: Dict[str, Any]) -> ValidationResult:
    """Validate metadata for FAIR compliance.

    Args:
        metadata: Metadata dictionary to validate

    Returns:
        ValidationResult with detailed findings
    """
    result = ValidationResult()

    # Required FAIR metadata fields
    required_fields = {
        "title": "Dataset title",
        "description": "Dataset description",
        "creator": "Data creator/author",
        "date_created": "Creation date",
        "identifier": "Unique identifier",
        "format": "Data format",
        "license": "Usage license",
    }

    # Check for required fields
    for field, description in required_fields.items():
        if field not in metadata or not metadata[field]:
            result.add_error(f"Missing required FAIR metadata: {description} ({field})")

    # Recommended fields for EPR data
    recommended_fields = {
        "instrument": "EPR spectrometer information",
        "measurement_parameters": "Measurement conditions",
        "sample_information": "Sample description",
        "processing_history": "Data processing steps",
        "units": "Data units information",
    }

    for field, description in recommended_fields.items():
        if field not in metadata or not metadata[field]:
            result.add_warning(
                f"Missing recommended EPR metadata: {description} ({field})"
            )

    # Validate specific field formats
    if "date_created" in metadata:
        try:
            # Basic ISO 8601 date format check
            date_str = str(metadata["date_created"])
            if len(date_str) < 10 or "T" not in date_str:
                result.add_warning(
                    "Date format should be ISO 8601 (YYYY-MM-DDTHH:MM:SS)"
                )
        except Exception:
            result.add_error("Invalid date_created format")

    if "identifier" in metadata:
        identifier = str(metadata["identifier"])
        if len(identifier) < 5:
            result.add_warning("Identifier should be more descriptive")

    # Check for units information
    if "units" in metadata:
        units = metadata["units"]
        if isinstance(units, dict):
            if "x_axis" not in units:
                result.add_warning("Missing x-axis units specification")
            if "y_axis" not in units:
                result.add_warning("Missing y-axis units specification")
        else:
            result.add_warning(
                "Units should be a dictionary with x_axis and y_axis fields"
            )

    result.add_info(f"Checked {len(metadata)} metadata fields")
    return result


def validate_data_integrity(
    x_data: Optional[np.ndarray], y_data: np.ndarray, metadata: Dict[str, Any]
) -> ValidationResult:
    """Validate data integrity and consistency.

    Args:
        x_data: X-axis data array
        y_data: Y-axis data array
        metadata: Associated metadata

    Returns:
        ValidationResult with findings
    """
    result = ValidationResult()

    # Check data arrays
    if y_data is None or len(y_data) == 0:
        result.add_error("Y-data is empty or None")
        return result

    # Check for NaN or infinite values
    if np.any(np.isnan(y_data)):
        result.add_error("Y-data contains NaN values")

    if np.any(np.isinf(y_data)):
        result.add_error("Y-data contains infinite values")

    # Check x-data consistency if present
    if x_data is not None:
        if len(x_data) != len(y_data):
            result.add_error(
                f"X-data length ({len(x_data)}) != Y-data length ({len(y_data)})"
            )

        if np.any(np.isnan(x_data)):
            result.add_warning("X-data contains NaN values")

        if np.any(np.isinf(x_data)):
            result.add_warning("X-data contains infinite values")

        # Check for monotonic x-data (expected for most EPR measurements)
        if not np.all(np.diff(x_data) >= 0) and not np.all(np.diff(x_data) <= 0):
            result.add_warning("X-data is not monotonic")

    # Validate data ranges
    y_range = np.max(y_data) - np.min(y_data)
    if y_range == 0:
        result.add_warning("Y-data has zero range (constant values)")

    # Check for reasonable data density
    data_points = len(y_data)
    if data_points < 10:
        result.add_warning("Very few data points (< 10)")
    elif data_points > 100000:
        result.add_info("Large dataset (> 100k points)")

    # Validate consistency with metadata
    if "data_points" in metadata:
        expected_points = metadata["data_points"]
        if isinstance(expected_points, (int, float)) and expected_points != data_points:
            result.add_warning(
                f"Data points mismatch: got {data_points}, expected {expected_points}"
            )

    result.add_info(f"Validated data arrays: {data_points} points")
    return result


def validate_epr_parameters(metadata: Dict[str, Any]) -> ValidationResult:
    """Validate EPR-specific measurement parameters.

    Args:
        metadata: Metadata containing EPR parameters

    Returns:
        ValidationResult with EPR-specific validation
    """
    result = ValidationResult()

    # Essential EPR parameters
    epr_required = {
        "microwave_frequency": "Microwave frequency (GHz)",
        "magnetic_field_range": "Magnetic field range",
        "modulation_frequency": "Field modulation frequency",
        "modulation_amplitude": "Field modulation amplitude",
    }

    # Look in various metadata locations
    params_sources = [metadata]
    if "measurement_parameters" in metadata:
        params_sources.append(metadata["measurement_parameters"])
    if "instrument_parameters" in metadata:
        params_sources.append(metadata["instrument_parameters"])

    found_params = set()
    for params in params_sources:
        if isinstance(params, dict):
            found_params.update(params.keys())

    for param, description in epr_required.items():
        if param not in found_params:
            # Check common variations
            variations = {
                "microwave_frequency": ["mw_freq", "frequency", "FrequencyMon"],
                "magnetic_field_range": ["field_range", "b_range", "HCF", "HSW"],
                "modulation_frequency": ["mod_freq", "ModFreq"],
                "modulation_amplitude": ["mod_amp", "ModAmp", "RMA"],
            }

            found_variation = False
            if param in variations:
                for var in variations[param]:
                    if var in found_params:
                        found_variation = True
                        result.add_info(f"Found parameter variation: {var} for {param}")
                        break

            if not found_variation:
                result.add_warning(f"Missing EPR parameter: {description}")

    # Validate parameter ranges if present
    for params in params_sources:
        if isinstance(params, dict):
            if "microwave_frequency" in params:
                freq = params["microwave_frequency"]
                if isinstance(freq, (int, float)):
                    if freq < 1 or freq > 1000:  # GHz
                        result.add_warning(f"Unusual microwave frequency: {freq} GHz")

            if "temperature" in params:
                temp = params["temperature"]
                if isinstance(temp, (int, float)):
                    if temp < 0 or temp > 1000:  # Kelvin
                        result.add_warning(f"Unusual temperature: {temp} K")

    result.add_info("Validated EPR-specific parameters")
    return result


def validate_file_format(
    file_path: Path, expected_format: Optional[str] = None
) -> ValidationResult:
    """Validate file format and structure.

    Args:
        file_path: Path to file to validate
        expected_format: Expected format (csv, json, hdf5)

    Returns:
        ValidationResult with format validation
    """
    result = ValidationResult()

    if not file_path.exists():
        result.add_error(f"File does not exist: {file_path}")
        return result

    file_size = file_path.stat().st_size
    if file_size == 0:
        result.add_error("File is empty")
        return result

    # Check file extension
    file_ext = file_path.suffix.lower()

    # Validate specific formats
    if expected_format:
        expected_format = expected_format.lower()

        if expected_format == "csv" and file_ext != ".csv":
            result.add_warning(f"Expected CSV format but file has {file_ext} extension")
        elif expected_format == "json" and file_ext != ".json":
            result.add_warning(
                f"Expected JSON format but file has {file_ext} extension"
            )
        elif expected_format == "hdf5" and file_ext not in [".h5", ".hdf5"]:
            result.add_warning(
                f"Expected HDF5 format but file has {file_ext} extension"
            )

    # Basic format validation
    try:
        if file_ext == ".json":
            with open(file_path, "r", encoding="utf-8") as f:
                json.load(f)
            result.add_info("JSON format validation passed")

        elif file_ext == ".csv":
            # Basic CSV check
            with open(file_path, "r", encoding="utf-8") as f:
                first_line = f.readline().strip()
                if "," not in first_line and "\t" not in first_line:
                    result.add_warning("CSV file may not be properly delimited")
            result.add_info("CSV format validation passed")

        elif file_ext in [".h5", ".hdf5"]:
            try:
                import h5py

                with h5py.File(file_path, "r") as f:
                    if len(f.keys()) == 0:
                        result.add_warning("HDF5 file contains no datasets")
                result.add_info("HDF5 format validation passed")
            except ImportError:
                result.add_warning("Cannot validate HDF5 format (h5py not available)")

    except Exception as e:
        result.add_error(f"Format validation failed: {e}")

    # Check file size warnings
    if file_size > 100 * 1024 * 1024:  # 100 MB
        result.add_info(f"Large file size: {file_size / (1024*1024):.1f} MB")

    return result


def validate_fair_dataset(
    data_dict: Dict[str, Any], file_path: Optional[Path] = None
) -> ValidationResult:
    """Comprehensive FAIR dataset validation.

    Args:
        data_dict: Dictionary containing data and metadata
        file_path: Optional file path for format validation

    Returns:
        Combined ValidationResult
    """
    logger.info("Starting comprehensive FAIR dataset validation")

    combined_result = ValidationResult()

    # Extract components
    x_data = data_dict.get("x_data")
    y_data = data_dict.get("y_data")
    metadata = data_dict.get("metadata", {})

    # Run all validations
    validations = []

    # Metadata validation
    meta_result = validate_fair_metadata(metadata)
    validations.append(("Metadata", meta_result))

    # Data integrity validation
    if y_data is not None:
        data_result = validate_data_integrity(x_data, y_data, metadata)
        validations.append(("Data Integrity", data_result))

    # EPR parameters validation
    epr_result = validate_epr_parameters(metadata)
    validations.append(("EPR Parameters", epr_result))

    # File format validation
    if file_path:
        format_result = validate_file_format(file_path)
        validations.append(("File Format", format_result))

    # Combine results
    for name, result in validations:
        combined_result.errors.extend([f"[{name}] {err}" for err in result.errors])
        combined_result.warnings.extend(
            [f"[{name}] {warn}" for warn in result.warnings]
        )
        combined_result.info.extend([f"[{name}] {info}" for info in result.info])

        if not result.is_valid:
            combined_result.is_valid = False

    # Summary
    logger.info(
        f"Validation complete: {len(combined_result.errors)} errors, "
        f"{len(combined_result.warnings)} warnings"
    )

    return combined_result


def create_validation_report(
    result: ValidationResult, output_path: Optional[Path] = None
) -> str:
    """Create a formatted validation report.

    Args:
        result: ValidationResult to format
        output_path: Optional path to save report

    Returns:
        Formatted report string
    """
    report_lines = []

    # Header
    report_lines.append("FAIR Dataset Validation Report")
    report_lines.append("=" * 35)
    report_lines.append("")

    # Summary
    status = "✓ PASSED" if result.is_valid else "✗ FAILED"
    report_lines.append(f"Overall Status: {status}")
    report_lines.append(f"Errors: {len(result.errors)}")
    report_lines.append(f"Warnings: {len(result.warnings)}")
    report_lines.append("")

    # Errors
    if result.errors:
        report_lines.append("ERRORS:")
        report_lines.append("-" * 8)
        for error in result.errors:
            report_lines.append(f"  • {error}")
        report_lines.append("")

    # Warnings
    if result.warnings:
        report_lines.append("WARNINGS:")
        report_lines.append("-" * 9)
        for warning in result.warnings:
            report_lines.append(f"  • {warning}")
        report_lines.append("")

    # Info
    if result.info:
        report_lines.append("INFO:")
        report_lines.append("-" * 5)
        for info in result.info:
            report_lines.append(f"  • {info}")
        report_lines.append("")

    report = "\n".join(report_lines)

    # Save to file if requested
    if output_path:
        output_path.write_text(report, encoding="utf-8")
        logger.info(f"Validation report saved to {output_path}")

    return report
