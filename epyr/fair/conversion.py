"""
Main conversion functions and workflows for Bruker EPR to FAIR format conversion.

This module provides the high-level interface for converting Bruker EPR data
to FAIR-compliant formats using the EPyR Tools package.
"""

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from epyr.eprload import eprload

from .data_processing import append_fair_metadata
from .exporters import save_fair as _save_fair_formats


def convert_bruker_to_fair(
    input_file: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    formats: List[str] = ["csv", "json"],
    include_metadata: bool = True,
    scaling: str = "",
) -> bool:
    """Load Bruker EPR data using eprload and convert to FAIR formats.

    This is the main conversion function that handles the complete workflow
    from loading Bruker data to saving in FAIR-compliant formats with
    structured metadata.

    Args:
        input_file: Path to Bruker data file (.dta, .dsc, .spc, .par).
        output_dir: Directory to save converted files. If None, saves in
            same directory as input file.
        formats: List of formats to generate. Options: 'csv', 'json', 'hdf5'.
        include_metadata: Whether to include metadata in output files.
        scaling: Scaling options passed to eprload (e.g., 'nPGT').

    Returns:
        True if conversion successful, False otherwise.
    """
    from ..logging_config import get_logger
    logger = get_logger(__name__)
    
    try:
        logger.info("Starting FAIR conversion process...")
        input_file = Path(input_file)
        
        if not input_file.exists():
            logger.error(f"Input file not found: {input_file}")
            return False

        # Load data using eprload (disable internal plotting)
        x, y, pars, original_file_path_str = eprload(
            file_name=str(input_file),
            scaling=scaling,
            plot_if_possible=False,
        )

        # Check if loading was successful
        if y is None or pars is None or original_file_path_str is None:
            logger.error("Data loading failed. Aborting conversion.")
            return False

        logger.info(f"Successfully loaded: {original_file_path_str}")
        original_file_path = Path(original_file_path_str)

        # Determine output location and basename
        if output_dir is None:
            output_path = original_file_path.parent
        else:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

        output_basename = output_path / original_file_path.stem
        logger.info(f"Output base name: {output_basename}")

        # Perform conversions based on requested formats
        if not formats:
            logger.warning("No output formats specified. Nothing to do.")
            return False

        logger.info("Processing parameters and generating outputs...")

        # Convert formats list to legacy format for compatibility
        output_formats = []
        if 'csv' in formats and 'json' in formats:
            output_formats.append('csv_json')
        elif 'csv' in formats:
            output_formats.append('csv')
        elif 'json' in formats:
            output_formats.append('json')
        
        if 'hdf5' in formats:
            output_formats.append('hdf5')

        # Use the consolidated save function from exporters
        _save_fair_formats(
            output_basename, x, y, pars, original_file_path_str, output_formats
        )

        logger.info("FAIR conversion process finished.")
        return True
        
    except Exception as e:
        logger.error(f"Conversion failed: {e}")
        return False


def save_fair(
    output_basename: Union[str, Path],
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Dict[str, Any],
    original_file_path: str,
    output_formats: List[str] = ["csv_json", "hdf5"],
) -> None:
    """Save already-loaded EPR data to one or more FAIR formats.

    This function is useful when you have already loaded EPR data and want
    to save it in FAIR formats without going through the loading process again.

    Args:
        output_basename: Base path for output files (without extension).
        x: Abscissa data array(s) or None.
        y: Ordinate data array.
        params: Dictionary of parameters from the original file.
        original_file_path: Full path of the original loaded file.
        output_formats: List of formats to save. Options: 'csv_json', 'hdf5'.

    Returns:
        None. Files are saved to disk.
    """
    output_basename = Path(output_basename)
    output_path = output_basename.parent
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Saving data from '{original_file_path}' to FAIR formats...")

    # Use the consolidated save function from exporters
    _save_fair_formats(
        output_basename, x, y, params, original_file_path, output_formats
    )

    print("\nFAIR saving process finished.")


def batch_convert_directory(
    input_directory: Union[str, Path],
    output_directory: Optional[Union[str, Path]] = None,
    file_extensions: List[str] = [".dsc", ".spc", ".par"],
    scaling: str = "",
    output_formats: List[str] = ["csv_json", "hdf5"],
    recursive: bool = False,
) -> None:
    """Convert all Bruker EPR files in a directory to FAIR formats.

    Args:
        input_directory: Directory containing Bruker EPR files.
        output_directory: Directory to save converted files. If None, saves
            alongside original files.
        file_extensions: List of file extensions to process.
        scaling: Scaling options passed to eprload.
        output_formats: List of formats to generate.
        recursive: If True, search subdirectories recursively.

    Returns:
        None. Prints progress and summary.
    """
    input_path = Path(input_directory)
    if not input_path.is_dir():
        raise ValueError(f"Input path is not a directory: {input_directory}")

    if output_directory is not None:
        output_path = Path(output_directory)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = None

    print(f"Starting batch conversion of directory: {input_path}")
    print(f"Looking for files with extensions: {file_extensions}")

    # Find all matching files
    files_to_process = []
    for ext in file_extensions:
        if recursive:
            pattern = f"**/*{ext}"
        else:
            pattern = f"*{ext}"
        files_to_process.extend(input_path.glob(pattern))

    if not files_to_process:
        print("No matching files found.")
        return

    print(f"Found {len(files_to_process)} files to process.")

    successful_conversions = 0
    failed_conversions = 0

    for i, file_path in enumerate(files_to_process, 1):
        print(f"\n--- Processing {i}/{len(files_to_process)}: {file_path.name} ---")

        try:
            # Determine output directory for this file
            if output_path is not None:
                # Maintain relative directory structure in output
                rel_path = file_path.parent.relative_to(input_path)
                file_output_dir = output_path / rel_path
            else:
                file_output_dir = None

            convert_bruker_to_fair(
                input_file_or_dir=file_path,
                output_dir=file_output_dir,
                scaling=scaling,
                output_formats=output_formats,
            )
            successful_conversions += 1

        except Exception as e:
            print(f"Error processing {file_path}: {type(e).__name__} - {e}")
            failed_conversions += 1

    print(f"\n--- Batch conversion complete ---")
    print(f"Successfully processed: {successful_conversions} files")
    print(f"Failed to process: {failed_conversions} files")


def validate_conversion(
    fair_json_file: Union[str, Path],
    original_data_file: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Validate a FAIR conversion by checking the JSON metadata file.

    Args:
        fair_json_file: Path to the JSON metadata file from FAIR conversion.
        original_data_file: Path to original data file for comparison (optional).

    Returns:
        Dictionary with validation results including warnings and statistics.
    """
    import json

    json_path = Path(fair_json_file)

    if not json_path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    validation_results = {
        "file": str(json_path),
        "valid": True,
        "warnings": [],
        "statistics": {},
        "metadata_completeness": {},
    }

    try:
        with open(json_path, "r", encoding="utf-8") as f:
            metadata = json.load(f)

        # Check required top-level keys
        required_keys = ["original_file", "fair_metadata"]
        missing_keys = [key for key in required_keys if key not in metadata]

        if missing_keys:
            validation_results["warnings"].append(
                f"Missing required keys: {missing_keys}"
            )
            validation_results["valid"] = False

        # Analyze FAIR metadata completeness
        if "fair_metadata" in metadata:
            fair_meta = metadata["fair_metadata"]
            validation_results["statistics"]["total_fair_parameters"] = len(fair_meta)

            # Count parameters by category
            categories = {
                "measurement": [
                    "microwave_frequency",
                    "microwave_power",
                    "modulation_amplitude",
                ],
                "sample": ["sample_identifier", "sample_form"],
                "acquisition": [
                    "acquisition_date",
                    "acquisition_time",
                    "number_of_scans",
                ],
                "axes": ["x_axis_unit", "y_axis_unit", "number_of_points_x_axis"],
            }

            for category, params in categories.items():
                found = sum(1 for param in params if param in fair_meta)
                validation_results["metadata_completeness"][
                    category
                ] = f"{found}/{len(params)}"

        # Check for unmapped parameters
        if "unmapped_parameters" in metadata:
            unmapped_count = len(metadata["unmapped_parameters"])
            validation_results["statistics"]["unmapped_parameters"] = unmapped_count
            if unmapped_count > 0:
                validation_results["warnings"].append(
                    f"{unmapped_count} parameters could not be mapped to FAIR format"
                )

        # Validate conversion timestamp
        if (
            "fair_metadata" in metadata
            and "conversion_info" in metadata["fair_metadata"]
        ):
            conv_info = metadata["fair_metadata"]["conversion_info"]
            if "value" in conv_info and "conversion_timestamp" in conv_info["value"]:
                validation_results["statistics"]["conversion_timestamp"] = conv_info[
                    "value"
                ]["conversion_timestamp"]

    except Exception as e:
        validation_results["valid"] = False
        validation_results["warnings"].append(f"Error reading JSON file: {str(e)}")

    return validation_results


if __name__ == "__main__":
    print("--- Running Bruker to FAIR Converter ---")
    print("This will use eprload's dialog to select a Bruker file,")
    print("then convert it to CSV/JSON and HDF5 formats.")
    print("-" * 50)

    # Use file dialog to select input, save to same directory
    convert_bruker_to_fair()
