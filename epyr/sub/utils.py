# sub/utils.py
import re
import sys
import warnings
from pathlib import Path
from typing import List, Union

import numpy as np

# Regular expression to check if a string can be converted to a number
_NUMBER_RE = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")


def read_par_file(par_file_path: Path) -> dict:
    """Reads a Bruker ESP/WinEPR .par file (key-value pairs)."""
    parameters = {}
    if not par_file_path.is_file():
        raise FileNotFoundError(f"Cannot find the parameter file {par_file_path}")

    try:
        with open(par_file_path, "r", encoding="latin-1") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split(maxsplit=1)
                if len(parts) < 1:
                    continue

                key = parts[0]
                # Check if key starts with a letter (basic validation)
                if not key[0].isalpha():
                    continue

                value = parts[1].strip() if len(parts) > 1 else ""

                # Remove surrounding single quotes if present
                if len(value) >= 2 and value.startswith("'") and value.endswith("'"):
                    value = value[1:-1]

                # Basic cleaning of key for dictionary access if needed
                # Note: Bruker keys are typically valid identifiers
                # if not key.isidentifier():
                #    key = re.sub(r'\W|^(?=\d)', '_', key) # Simple sanitization

                parameters[key] = value
    except Exception as e:
        raise IOError(f"Error reading PAR file {par_file_path}: {e}") from e

    if parameters.get("JEX"):
        parameters["XAXIS_NAME"] = parameters["JEX"]
    if parameters.get("JUN"):
        parameters["XAXIS_UNIT"] = parameters["JUN"]
    if parameters.get("XXUN"):
        parameters["XAXIS_UNIT"] = parameters["XXUN"]
    if parameters.get("JEY"):
        parameters["YAXIS_NAME"] = parameters["JEY"]
    if parameters.get("XYUN"):
        parameters["YAXIS_UNIT"] = parameters["XYUN"]

    return parameters


def read_dsc_file(dsc_file_path: Path) -> dict:
    """Reads a Bruker BES3T .DSC file (key-value pairs, handles line continuation)."""
    parameters = {}
    if not dsc_file_path.is_file():
        raise FileNotFoundError(f"Cannot find the descriptor file {dsc_file_path}")

    lines = []
    try:
        with open(dsc_file_path, "r", encoding="latin-1") as f:
            lines = f.readlines()
    except Exception as e:
        raise IOError(f"Error reading DSC file {dsc_file_path}: {e}") from e

    processed_lines = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # Handle line continuation characters '\'
        while line.endswith("\\"):
            i += 1
            if i < len(lines):
                line = line[:-1] + lines[i].strip()
            else:
                line = line[:-1]  # Remove trailing '\' even if it's the last line
        processed_lines.append(line.replace("\\n", "\n"))  # Replace escaped newlines
        i += 1

    for line in processed_lines:
        if not line:
            continue

        parts = line.split(maxsplit=1)
        if len(parts) < 1:
            continue

        key = parts[0]
        # Stop if Manipulation History Layer is reached
        if key.upper() == "#MHL":
            break
        # Skip lines not starting with a letter (comments, etc.)
        if not key[0].isalpha():
            continue

        value = parts[1].strip() if len(parts) > 1 else ""

        # Remove surrounding single quotes if present
        if len(value) >= 2 and value.startswith("'") and value.endswith("'"):
            value = value[1:-1]

        # Basic cleaning of key if needed (Bruker keys usually okay)
        # if not key.isidentifier():
        #     key = re.sub(r'\W|^(?=\d)', '_', key)

        parameters[key] = value
    if parameters.get("XNAM"):
        parameters["XAXIS_NAME"] = parameters["XNAM"]
        parameters["XAXIS_UNIT"] = parameters["XUNI"]
    if parameters.get("YNAM"):
        parameters["YAXIS_NAME"] = parameters["YNAM"]
        parameters["YAXIS_UNIT"] = parameters["YUNI"]

    return parameters


def parse_field_params(parameters: dict) -> dict:
    """
    Attempts to convert string values in a dictionary to numbers (int or float).
    """
    parsed_params = {}
    for key, value in parameters.items():
        if isinstance(value, str):
            # Try converting to int
            try:
                parsed_params[key] = int(value)
                continue
            except ValueError:
                pass
            # Try converting to float
            try:
                # Use regex for more robust float check if needed,
                # but direct conversion attempt is usually fine
                if _NUMBER_RE.match(value):
                    parsed_params[key] = float(value)
                else:
                    parsed_params[key] = value  # Keep as string if not number-like
            except ValueError:
                parsed_params[key] = value  # Keep original string if conversion fails
        else:
            parsed_params[key] = value  # Keep non-string values as they are
    return parsed_params


def get_matrix(
    data_file_path: Path,
    dimensions: List[int],
    number_format_code: str,
    byte_order: str,
    is_complex: Union[bool, np.ndarray],
) -> np.ndarray:
    """
    Reads binary data from a file into a NumPy array.

    Args:
        data_file_path: Path to the data file (.DTA, .spc).
        dimensions: List of dimensions [nx, ny, nz].
        number_format: String representing numpy dtype ('int8', 'int16', etc.).
        byte_order: 'ieee-be' (big) or 'ieee-le' (little).
        is_complex: Boolean indicating if the data is complex. Can be array for multi-channel.

    Returns:
        NumPy array with the data.
    """
    if not data_file_path.is_file():
        raise FileNotFoundError(f"Data file not found: {data_file_path}")

    # Determine numpy dtype and endianness
    dt_char = ">" if byte_order == "ieee-be" else "<"
    try:
        # Construct dtype using standard codes (e.g., '>f8', '<i4')
        dtype = np.dtype(f"{dt_char}{number_format_code}")
    except TypeError as e:  # Catch potential error during dtype creation
        # Add original exception context using 'from e'
        raise ValueError(f"Unsupported number format code: {number_format_code}") from e

    # Calculate expected number of elements
    n_points_total = int(np.prod(dimensions))
    if n_points_total == 0:
        return np.array([])

    # Handle potentially complex data reading
    # For now, assume is_complex is a single boolean
    # A more complex implementation could handle mixed real/complex channels
    is_complex_flag = (
        np.any(is_complex) if isinstance(is_complex, (list, np.ndarray)) else is_complex
    )

    if is_complex_flag:
        n_values_to_read = n_points_total * 2
        actual_dtype = dtype.base  # Read underlying real type
    else:
        n_values_to_read = n_points_total
        actual_dtype = dtype

    # Read raw data from file
    try:
        raw_data = np.fromfile(
            data_file_path, dtype=actual_dtype, count=n_values_to_read
        )
    except Exception as e:
        raise IOError(f"Error reading data file {data_file_path}: {e}") from e

    # Verify number of elements read
    if raw_data.size < n_values_to_read:
        raise IOError(
            f"Could not read expected number of data points from {data_file_path}. "
            f"Expected {n_values_to_read}, got {raw_data.size}."
        )
    elif raw_data.size > n_values_to_read:
        warnings.warn(
            f"Read more data points ({raw_data.size}) than expected ({n_values_to_read}) "
            f"from {data_file_path}. Truncating."
        )
        raw_data = raw_data[:n_values_to_read]

    # Combine real and imaginary parts if complex
    if is_complex_flag:
        if raw_data.size % 2 != 0:
            raise ValueError("Read odd number of values for complex data.")
        data = raw_data[::2] + 1j * raw_data[1::2]
    else:
        data = raw_data

    # Reshape the data - NumPy uses C order (last index fastest)
    # MATLAB uses Fortran order (first index fastest)
    # BES3T/ESP files are typically C-ordered (X varies fastest)
    # Reshape to (nz, ny, nx) if 3D, (ny, nx) if 2D, (nx,) if 1D
    shape_numpy_order = [
        d for d in dimensions[::-1] if d > 1
    ]  # Reverse and remove dims of size 1
    if not shape_numpy_order:  # If all dims are 1 or empty
        shape_numpy_order = (n_points_total,)

    try:
        # Use squeeze to remove dimensions of size 1, similar to MATLAB's behavior
        data = data.reshape(shape_numpy_order).squeeze()
        # If the result is 0-dim after squeeze (single point), make it 1-dim
        if data.ndim == 0:
            data = data.reshape(1)

    except ValueError as e:
        raise ValueError(
            f"Could not reshape data with {data.size} points into desired shape {shape_numpy_order}. Original dims: {dimensions}. Error: {e}"
        ) from e

    return data


def BrukerListFiles(path, recursive=False):
    """
    List all Bruker EPR data files (.DTA, .dta, .SPC, .spc) in the given directory.

    Args:
        path (str or Path): Path to the folder containing Bruker files.
        recursive (bool, optional): If True, search subfolders recursively. Defaults to False.

    Returns:
        list[Path]: Sorted list of Path objects for found files.
    """
    exts = {".dta", ".DTA", ".spc", ".SPC"}
    path = Path(path)
    if not path.is_dir():
        raise NotADirectoryError(f"{path} is not a valid directory.")

    if recursive:
        files = [p for p in path.rglob("*") if p.suffix in exts and p.is_file()]
    else:
        files = [p for p in path.iterdir() if p.suffix in exts and p.is_file()]

    return sorted(files)
