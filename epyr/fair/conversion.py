"""
Main conversion functions and workflows for Bruker EPR to FAIR format conversion.

This module provides the high-level interface for converting Bruker EPR data
to FAIR-compliant formats using the EPyR Tools package.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from ..eprload import eprload
from ..logging_config import get_logger

logger = get_logger(__name__)

from .exporters import save_fair as _save_fair_formats


def convert_bruker_to_fair(
    input_file: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    formats: Optional[List[str]] = None,
    include_metadata: bool = True,
    scaling: str = "",
) -> bool:
    """Load a Bruker file and convert it to FAIR-compliant formats.

    End-to-end pipeline: :func:`epyr.eprload` -> parameter normalization
    -> writing one or more output files next to the input.

    Parameters
    ----------
    input_file : str or pathlib.Path
        Path to a Bruker data file (``.dta``, ``.dsc``, ``.spc``, ``.par``).
    output_dir : str or pathlib.Path, optional
        Where to write the converted files. Defaults to the input file's
        directory. Created if missing.
    formats : list of str, optional
        Subset of ``['csv', 'json', 'hdf5', 'jpg']``. Default
        ``['csv', 'json']``.

        - ``csv``  : data array
        - ``json`` : metadata only
        - ``hdf5`` : data + metadata in one file
        - ``jpg``  : preview figure (1D plot or 2D map + waterfall)
    include_metadata : bool, optional
        Whether CSV files should carry the parameter dictionary as
        commented header lines. Default True.
    scaling : str, optional
        Scaling code passed through to :func:`epyr.eprload`. Default ``""``.

    Returns
    -------
    bool
        True on success, False if loading or writing failed (errors are
        logged, not raised).

    Examples
    --------
    >>> from epyr.fair import convert_bruker_to_fair
    >>> ok = convert_bruker_to_fair(
    ...     "examples/data/130406SB_CaWO4_Er_CW_5K_20.DSC",
    ...     output_dir="/tmp/epyr_out",
    ...     formats=["json", "hdf5"],
    ... )  # doctest: +SKIP
    >>> ok  # doctest: +SKIP
    True
    """
    if formats is None:
        formats = ["csv", "json"]
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

        # Use the consolidated save function from exporters
        # Formats are now passed directly (csv, json, hdf5)
        _save_fair_formats(output_basename, x, y, pars, original_file_path_str, formats)

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
    output_formats: Optional[List[str]] = None,
) -> None:
    """Write already-loaded EPR data to one or more FAIR formats.

    Use this when data is already in memory (e.g., after processing) and
    you want to write outputs without re-reading the Bruker file.

    Parameters
    ----------
    output_basename : str or pathlib.Path
        Base path (no extension); each format appends its own.
    x : np.ndarray, list of np.ndarray, or None
        Abscissa as returned by :func:`epyr.eprload`.
    y : np.ndarray
        Signal data.
    params : dict
        Parameter dictionary from the original file.
    original_file_path : str
        Full path of the source file, kept as provenance in the outputs.
    output_formats : list of str, optional
        Subset of ``['csv', 'json', 'hdf5', 'jpg', 'csv_json']``.
        Default ``['csv', 'json']``.

    Returns
    -------
    None
        Outputs are written to disk under ``output_basename``.

    Examples
    --------
    >>> import numpy as np
    >>> from epyr.fair import save_fair
    >>> x = np.linspace(3300, 3400, 100)
    >>> y = np.random.randn(100)
    >>> save_fair("/tmp/demo", x, y, {"MWFQ": 9.4e9},
    ...           "demo.dsc", ["json"])  # doctest: +SKIP
    """
    if output_formats is None:
        output_formats = ["csv", "json"]
    output_basename = Path(output_basename)
    output_path = output_basename.parent
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving data from '{original_file_path}' to FAIR formats...")

    # Use the consolidated save function from exporters
    _save_fair_formats(
        output_basename, x, y, params, original_file_path, output_formats
    )

    logger.info("\nFAIR saving process finished.")


def batch_convert_directory(
    input_directory: Union[str, Path],
    output_directory: Optional[Union[str, Path]] = None,
    file_extensions: Optional[List[str]] = None,
    scaling: str = "",
    output_formats: Optional[List[str]] = None,
    recursive: bool = False,
) -> None:
    """Convert every Bruker file in a directory to FAIR formats.

    Parameters
    ----------
    input_directory : str or pathlib.Path
        Directory to scan.
    output_directory : str or pathlib.Path, optional
        Where to write the converted files. Defaults to alongside each
        input file. Created if missing.
    file_extensions : list of str, optional
        Which extensions count as Bruker files. Default
        ``['.dsc', '.spc', '.par']``.
    scaling : str, optional
        Scaling code passed through to :func:`epyr.eprload`.
    output_formats : list of str, optional
        Subset of ``['csv', 'json', 'hdf5', 'jpg']``. Default
        ``['csv', 'json']``.
    recursive : bool, optional
        Also descend into subdirectories. Default False.

    Returns
    -------
    None
        Progress and a final per-file summary are logged.

    Examples
    --------
    >>> from epyr.fair import batch_convert_directory
    >>> batch_convert_directory(
    ...     "examples/data",
    ...     output_directory="/tmp/epyr_batch",
    ...     output_formats=["json"],
    ... )  # doctest: +SKIP
    """
    if file_extensions is None:
        file_extensions = [".dsc", ".spc", ".par"]
    if output_formats is None:
        output_formats = ["csv", "json"]
    input_path = Path(input_directory)
    if not input_path.is_dir():
        raise ValueError(f"Input path is not a directory: {input_directory}")

    if output_directory is not None:
        output_path = Path(output_directory)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = None

    logger.info(f"Starting batch conversion of directory: {input_path}")
    logger.info(f"Looking for files with extensions: {file_extensions}")

    # Find all matching files
    files_to_process = []
    for ext in file_extensions:
        if recursive:
            pattern = f"**/*{ext}"
        else:
            pattern = f"*{ext}"
        files_to_process.extend(input_path.glob(pattern))

    if not files_to_process:
        logger.info("No matching files found.")
        return

    logger.info(f"Found {len(files_to_process)} files to process.")

    successful_conversions = 0
    failed_conversions = 0

    for i, file_path in enumerate(files_to_process, 1):
        logger.info(
            f"\n--- Processing {i}/{len(files_to_process)}: {file_path.name} ---"
        )

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
            logger.error(f"Error processing {file_path}: {type(e).__name__} - {e}")
            failed_conversions += 1

    logger.info("\n--- Batch conversion complete ---")
    logger.info(f"Successfully processed: {successful_conversions} files")
    logger.info(f"Failed to process: {failed_conversions} files")


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
    logger.info("--- Running Bruker to FAIR Converter ---")
    logger.info("This will use eprload's dialog to select a Bruker file,")
    logger.info("then convert it to CSV/JSON and HDF5 formats.")
    logger.info("-" * 50)

    # Use file dialog to select input, save to same directory
    convert_bruker_to_fair()
