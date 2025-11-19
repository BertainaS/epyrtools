"""
Core data processing and metadata handling for FAIR format conversion.

This module contains functions for processing Bruker EPR parameters and
preparing metadata for FAIR format export.
"""

import warnings
from datetime import datetime
from typing import Any, Dict, List, Tuple, Union

import numpy as np

from ..logging_config import get_logger

logger = get_logger(__name__)

from .parameter_mapping import BRUKER_PARAM_MAP


def process_parameters(pars: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Process raw parameters using BRUKER_PARAM_MAP.

    Args:
        pars: Raw parameters dictionary from Bruker files

    Returns:
        Tuple of (fair_metadata, unmapped_parameters)
    """
    fair_metadata = {}
    unmapped_parameters = {}

    # Add conversion metadata
    fair_metadata["conversion_info"] = {
        "value": {
            "converter_script": "epyr_fair_converter",
            "conversion_timestamp": datetime.now().isoformat(),
            "epyr_version": "0.3.1",  # Current EPyR Tools version
        },
        "unit": "",
        "description": "Information about the conversion process to FAIR format.",
    }

    for key, value in pars.items():
        # Skip internal keys added by eprload
        if key.startswith("_"):
            continue

        if key in BRUKER_PARAM_MAP:
            map_info = BRUKER_PARAM_MAP[key]
            fair_key = map_info["fair_name"]

            # Determine unit: Use map unit, check for unit references
            unit = map_info["unit"]
            if unit == "refer to XUNI" and "XUNI" in pars:
                unit = pars.get("XUNI", "unknown")
            elif unit == "refer to YUNI" and "YUNI" in pars:
                yuni_val = pars.get("YUNI", "unknown")
                unit = (
                    yuni_val.split(",")[0].strip()
                    if isinstance(yuni_val, str)
                    else str(yuni_val)
                )
            elif unit == "refer to ZUNI" and "ZUNI" in pars:
                zuni_val = pars.get("ZUNI", "unknown")
                unit = (
                    zuni_val.split(",")[0].strip()
                    if isinstance(zuni_val, str)
                    else str(zuni_val)
                )

            fair_metadata[fair_key] = {
                "value": value,
                "unit": unit,
                "description": map_info["description"],
            }
        else:
            # Store unmapped parameters
            unmapped_parameters[key] = value
            logger.debug(
                f"Parameter '{key}' not found in BRUKER_PARAM_MAP. "
                f"Storing in 'unmapped_parameters'."
            )

    return fair_metadata, unmapped_parameters


def append_fair_metadata(
    data_dict: Dict[str, Any], pars: Dict[str, Any], original_file: str = ""
) -> Dict[str, Any]:
    """Append FAIR metadata to existing data dictionary.

    Args:
        data_dict: Existing data dictionary to append to
        pars: Raw parameters from Bruker file
        original_file: Original file path

    Returns:
        Updated data dictionary with FAIR metadata
    """
    fair_meta, unmapped_meta = process_parameters(pars)

    # Add file information
    if original_file:
        data_dict["original_file"] = original_file

    # Add processed metadata
    data_dict["fair_metadata"] = fair_meta

    # Add unmapped parameters if any exist
    if unmapped_meta:
        data_dict["unmapped_parameters"] = unmapped_meta

    return data_dict


def extract_axis_info(pars: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Extract axis information from parameters.

    Args:
        pars: Raw parameters dictionary

    Returns:
        Dictionary with x, y, z axis information
    """
    axis_info = {}

    # Process FAIR metadata to extract axis info
    fair_meta, _ = process_parameters(pars)

    # X-axis info
    x_info = {}
    if "x_axis_unit" in fair_meta:
        x_info["unit"] = fair_meta["x_axis_unit"]["value"]
    if "x_axis_name" in fair_meta:
        x_info["name"] = fair_meta["x_axis_name"]["value"]
    if "number_of_points_x_axis" in fair_meta:
        x_info["points"] = fair_meta["number_of_points_x_axis"]["value"]
    if "x_axis_minimum" in fair_meta:
        x_info["minimum"] = fair_meta["x_axis_minimum"]["value"]
    if "x_axis_width" in fair_meta:
        x_info["width"] = fair_meta["x_axis_width"]["value"]

    if x_info:
        axis_info["x"] = x_info

    # Y-axis info
    y_info = {}
    if "y_axis_unit" in fair_meta:
        y_info["unit"] = fair_meta["y_axis_unit"]["value"]
    if "y_axis_name" in fair_meta:
        y_info["name"] = fair_meta["y_axis_name"]["value"]
    if "number_of_points_y_axis" in fair_meta:
        y_info["points"] = fair_meta["number_of_points_y_axis"]["value"]
    if "y_axis_minimum" in fair_meta:
        y_info["minimum"] = fair_meta["y_axis_minimum"]["value"]
    if "y_axis_width" in fair_meta:
        y_info["width"] = fair_meta["y_axis_width"]["value"]

    if y_info:
        axis_info["y"] = y_info

    return axis_info


def validate_fair_metadata(metadata: Dict[str, Any]) -> List[str]:
    """Validate FAIR metadata for completeness and correctness.

    Args:
        metadata: FAIR metadata dictionary

    Returns:
        List of validation warnings/errors
    """
    warnings_list = []

    # Check for essential parameters
    essential_params = ["microwave_frequency", "x_axis_unit", "number_of_points_x_axis"]

    for param in essential_params:
        if param not in metadata:
            warnings_list.append(f"Missing essential parameter: {param}")

    # Check for reasonable values
    if "microwave_frequency" in metadata:
        freq = metadata["microwave_frequency"].get("value", 0)
        if isinstance(freq, (int, float)) and (freq <= 0 or freq > 1e12):
            warnings_list.append(f"Microwave frequency seems unreasonable: {freq} Hz")

    # Check unit consistency
    for key, info in metadata.items():
        if not isinstance(info, dict):
            continue
        if "unit" not in info:
            warnings_list.append(f"Parameter {key} missing unit information")
        if "description" not in info:
            warnings_list.append(f"Parameter {key} missing description")

    return warnings_list


def get_experiment_summary(metadata: Dict[str, Any]) -> Dict[str, Any]:
    """Extract key experimental parameters for quick overview.

    Args:
        metadata: FAIR metadata dictionary

    Returns:
        Dictionary with experiment summary
    """
    summary = {}

    # Key parameters for EPR experiments
    key_params = {
        "microwave_frequency": "MW Frequency",
        "microwave_power": "MW Power",
        "modulation_amplitude": "Modulation Amplitude",
        "receiver_time_constant": "Time Constant",
        "number_of_points_x_axis": "Data Points",
        "sample_identifier": "Sample",
        "experiment_type": "Experiment Type",
        "acquisition_date": "Date",
        "acquisition_time": "Time",
    }

    for param_key, display_name in key_params.items():
        if param_key in metadata:
            info = metadata[param_key]
            value = info.get("value", "N/A")
            unit = info.get("unit", "")

            # Format value with unit
            if unit and unit != "None":
                summary[display_name] = f"{value} {unit}"
            else:
                summary[display_name] = str(value)

    return summary
