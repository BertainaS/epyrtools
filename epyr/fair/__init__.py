"""
FAIR data conversion package for EPyR Tools.

This package provides functionality to convert proprietary Bruker EPR data
files into FAIR (Findable, Accessible, Interoperable, and Reusable) formats
such as CSV, JSON, and HDF5.

The package is organized into focused modules:
- parameter_mapping: Bruker parameter to FAIR format mappings
- data_processing: Core data processing and metadata handling
- exporters: Format-specific export functions (CSV, JSON, HDF5)
- conversion: Main conversion functions and workflows
"""

# Import main public API functions
from .conversion import (
    batch_convert_directory,
    convert_bruker_to_fair,
    save_fair,
)
from .data_processing import append_fair_metadata, extract_axis_info, process_parameters
from .exporters import save_to_csv_json, save_to_hdf5
from .parameter_mapping import BRUKER_PARAM_MAP
from .validation import validate_fair_dataset, ValidationResult, create_validation_report

# Re-export for backward compatibility
__all__ = [
    "convert_bruker_to_fair",
    "save_fair",
    "append_fair_metadata",
    "batch_convert_directory",
    "process_parameters",
    "extract_axis_info",
    "save_to_csv_json",
    "save_to_hdf5",
    "BRUKER_PARAM_MAP",
    "validate_fair_dataset",
    "ValidationResult", 
    "create_validation_report",
]
