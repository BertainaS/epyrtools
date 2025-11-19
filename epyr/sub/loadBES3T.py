# sub/loadBES3T.py
import warnings
from pathlib import Path

import numpy as np

from .utils import (  # Import from sibling utils
    get_matrix,
    parse_field_params,
    read_dsc_file,
)


def load(full_base_name: Path, file_extension: str, scaling: str) -> tuple:
    """
    Loads Bruker BES3T data (.DTA, .DSC).

    Args:
        full_base_name: Path object without extension.
        file_extension: The original file extension (e.g., '.dta', '.dsc').
        scaling: Scaling string (e.g., 'nP G').

    Returns:
        tuple: (data, abscissa, parameters)
    """
    # Determine DSC and DTA file extensions, respecting case
    dsc_extension = ".dsc"
    dta_extension = ".dta"
    if file_extension.isupper():
        dsc_extension = dsc_extension.upper()
        dta_extension = dta_extension.upper()

    # Use string concatenation instead of with_suffix() to handle filenames with multiple dots
    dsc_file = Path(str(full_base_name) + dsc_extension)
    dta_file = Path(str(full_base_name) + dta_extension)

    # Read descriptor file
    parameters = read_dsc_file(dsc_file)

    # --- Determine complexity, dimensions, byte order, format ---
    is_complex = np.array([False])  # Default to real
    n_data_values = 1
    if "IKKF" in parameters:
        parts = parameters["IKKF"].split(",")
        n_data_values = len(parts)
        is_complex = np.array([p.strip().upper() == "CPLX" for p in parts])
    else:
        warnings.warn("IKKF not found in .DSC file. Assuming IKKF=REAL.")

    # Dimensions
    try:
        nx = int(parameters.get("XPTS", 0))
        ny = int(parameters.get("YPTS", 1))
        nz = int(parameters.get("ZPTS", 1))
    except (ValueError, TypeError) as e:
        raise ValueError(
            f"Could not parse dimensions (XPTS/YPTS/ZPTS) from DSC file: {e}"
        ) from e
    if nx == 0:
        raise ValueError("XPTS is missing or zero in DSC file.")
    dimensions = [nx, ny, nz]

    # Byte Order
    byte_order = "ieee-be"  # Default big-endian
    if "BSEQ" in parameters:
        bseq_val = parameters["BSEQ"].upper()
        if bseq_val == "BIG":
            byte_order = "ieee-be"
        elif bseq_val == "LIT":
            byte_order = "ieee-le"
        else:
            raise ValueError(f"Unknown BSEQ value '{parameters['BSEQ']}' in .DSC file.")
    else:
        warnings.warn("BSEQ not found in .DSC file. Assuming BSEQ=BIG (big-endian).")

    # Number Format (assuming same for real and imag if complex)
    number_format_code = None
    if "IRFMT" in parameters:
        # For simplicity, take the first format if multiple are listed
        irfmt_val = parameters["IRFMT"].split(",")[0].strip().upper()
        fmt_map = {
            "C": "i1",  # int8   -> i1
            "S": "i2",  # int16  -> i2
            "I": "i4",  # int32  -> i4
            "F": "f4",  # float32-> f4
            "D": "f8",
        }  # float64-> f8
        if irfmt_val in fmt_map:
            number_format = fmt_map[irfmt_val]
        elif irfmt_val in ("A", "0", "N"):
            raise ValueError(
                f"Unsupported or no data format IRFMT='{irfmt_val}' in DSC."
            )
        else:
            raise ValueError(f"Unknown IRFMT value '{irfmt_val}' in .DSC file.")
    else:
        raise ValueError("IRFMT keyword not found in .DSC file.")

    if "IIFMT" in parameters and np.any(is_complex):
        iifmt_val = parameters["IIFMT"].split(",")[0].strip().upper()
        if iifmt_val != parameters["IRFMT"].split(",")[0].strip().upper():
            warnings.warn(
                "IRFMT and IIFMT differ in DSC file. Using IRFMT for reading."
            )
            # Raise error? MATLAB code enforces identity. Let's warn for now.
            # raise ValueError("IRFMT and IIFMT in DSC file must be identical.")

    # --- Construct Abscissa ---
    abscissa_list = [None] * 3  # X, Y, Z
    axis_names = ["X", "Y", "Z"]
    axis_defined = [False] * 3

    for i, axis in enumerate(axis_names):
        dim_size = dimensions[i]
        if dim_size <= 1:
            continue

        axis_type = parameters.get(f"{axis}TYP", "IDX")  # Default to linear index

        if axis_type == "IGD":  # Indirect, non-linear axis
            companion_suffix = f".{axis}GF"
            if file_extension.isupper():
                companion_suffix = companion_suffix.upper()
            # Use string concatenation instead of with_suffix() to handle filenames with multiple dots
            companion_file = Path(str(full_base_name) + companion_suffix)

            fmt_key = f"{axis}FMT"
            data_format_char = parameters.get(fmt_key, "D").upper()  # Default double
            fmt_map = {"D": "float64", "F": "float32", "I": "int32", "S": "int16"}
            if data_format_char not in fmt_map:
                warnings.warn(
                    f"Cannot read companion file format '{data_format_char}' for axis {axis}. Assuming linear."
                )
                axis_type = "IDX"  # Fallback to linear if format unknown
            else:
                companion_dtype_str = fmt_map[data_format_char]
                # Determine endianness for companion file (assume same as data)
                dt_char = ">" if byte_order == "ieee-be" else "<"
                companion_dtype = np.dtype(f"{dt_char}{companion_dtype_str}")

                if companion_file.is_file():
                    try:
                        axis_data = np.fromfile(
                            companion_file, dtype=companion_dtype, count=dim_size
                        )
                        if axis_data.size == dim_size:
                            abscissa_list[i] = axis_data
                            axis_defined[i] = True
                        else:
                            warnings.warn(
                                f"Could not read expected {dim_size} values from companion file {companion_file}. Assuming linear axis."
                            )
                            axis_type = "IDX"  # Fallback to linear
                    except Exception as e:
                        warnings.warn(
                            f"Error reading companion file {companion_file}: {e}. Assuming linear axis."
                        )
                        axis_type = "IDX"  # Fallback to linear
                else:
                    warnings.warn(
                        f"Companion file {companion_file} not found for non-linear axis {axis}. Assuming linear axis."
                    )
                    axis_type = "IDX"  # Fallback to linear

        if axis_type == "IDX":  # Linear axis
            min_key = f"{axis}MIN"
            wid_key = f"{axis}WID"
            try:
                minimum = float(parameters[min_key])
                width = float(parameters[wid_key])
                if dim_size > 1:
                    if width == 0:
                        warnings.warn(
                            f"{axis} range has zero width (WID=0). Using index range 0 to N-1."
                        )
                        # Use 0 to N-1 for indices if width is zero
                        abscissa_list[i] = np.arange(dim_size)
                    else:
                        # linspace is inclusive of endpoint
                        abscissa_list[i] = np.linspace(
                            minimum, minimum + width, dim_size
                        )
                elif dim_size == 1:
                    abscissa_list[i] = np.array([minimum])  # Single point axis
                else:  # dim_size == 0 ? Should not happen if XPTS>0
                    abscissa_list[i] = np.array([])

                axis_defined[i] = True
            except (KeyError, ValueError, TypeError):
                warnings.warn(
                    f"Could not read MIN/WID parameters for axis {axis}. Using default index."
                )
                abscissa_list[i] = np.arange(dim_size)  # Default to index
                axis_defined[i] = True  # Mark as defined (with index)

        elif axis_type == "NTUP":
            raise NotImplementedError("Cannot read data with NTUP axes.")

    # Consolidate abscissa
    defined_abscissae = [
        a for a, defined in zip(abscissa_list, axis_defined) if defined
    ]
    if len(defined_abscissae) == 1:
        abscissa = defined_abscissae[0]
    elif len(defined_abscissae) > 1:
        abscissa = defined_abscissae  # Return list for multiple dimensions
    else:
        abscissa = None  # No axes defined

    # --- Read Data Matrix ---
    # Assuming single data value type for now (n_data_values=1)
    # NOTE: Multiple data values per point (n_data_values > 1) not yet supported\n    # This would require handling interleaved data formats in some BES3T files
    if n_data_values > 1:
        warnings.warn(
            f"DSC file indicates {n_data_values} data values per point (IKKF). Only reading the first value."
        )
        # Adjust logic here if multiple channels need reading/combining

    data = get_matrix(dta_file, dimensions, number_format, byte_order, is_complex[0])

    # --- Scale Data ---
    if scaling and data is not None and data.size > 0:
        # Get experiment type and pre-scaling flag
        expt_type = parameters.get("EXPT", "CW").upper()
        is_cw = expt_type == "CW"
        # SctNorm indicates if Bruker software already applied some scaling
        data_prescaled = parameters.get("SctNorm", "false").lower() == "true"

        # Get parameters needed for scaling
        n_averages = None
        receiver_gain_db = None
        receiver_gain = None
        sampling_time_s = None
        sampling_time_ms = None
        mw_power_w = None
        mw_power_mw = None
        temperature_k = None

        try:
            n_averages = int(parameters.get("AVGS"))
        except (ValueError, TypeError, KeyError):
            pass
        try:
            receiver_gain_db = float(parameters.get("RCAG"))
        except (ValueError, TypeError, KeyError):
            pass
        try:
            sampling_time_s = float(parameters.get("SPTP"))  # Time in seconds
        except (ValueError, TypeError, KeyError):
            pass
        try:
            mw_power_w = float(parameters.get("MWPW"))  # Power in Watt
        except (ValueError, TypeError, KeyError):
            pass
        try:
            temperature_k = float(parameters.get("STMP"))  # Temperature in K
        except (ValueError, TypeError, KeyError):
            pass

        # Calculate derived values
        if receiver_gain_db is not None:
            receiver_gain = 10 ** (receiver_gain_db / 20.0)
        if sampling_time_s is not None:
            sampling_time_ms = sampling_time_s * 1000.0
        if mw_power_w is not None:
            mw_power_mw = mw_power_w * 1000.0

        # Apply scaling factors
        if "n" in scaling:
            if n_averages is not None and n_averages > 0:
                if data_prescaled:
                    # MATLAB errors here, let's warn
                    warnings.warn(
                        f"Cannot scale by number of scans ('n'): Data is already averaged (SctNorm=true, AVGS={n_averages})."
                    )
                else:
                    data = data / n_averages
            else:
                warnings.warn(
                    "Cannot scale by number of scans ('n'): AVGS missing, zero, or invalid."
                )

        if is_cw and "G" in scaling:
            if receiver_gain is not None and receiver_gain != 0:
                # Assume data not already scaled by gain, even if SctNorm=true
                # (Bruker scaling details can be complex)
                data = data / receiver_gain
            else:
                warnings.warn(
                    "Cannot scale by receiver gain ('G'): RCAG missing or invalid."
                )

        if is_cw and "c" in scaling:
            if sampling_time_ms is not None and sampling_time_ms > 0:
                # MATLAB notes Xepr scales even if SctNorm=false. Assume we should always scale if requested.
                # if data_prescaled:
                #    warnings.warn("Scaling by conversion time ('c') requested, but data may already be scaled (SctNorm=true). Applying anyway.")
                data = data / sampling_time_ms
            else:
                warnings.warn(
                    "Cannot scale by conversion time ('c'): SPTP missing, zero, or invalid."
                )

        if is_cw and "P" in scaling:
            if mw_power_mw is not None and mw_power_mw > 0:
                data = data / np.sqrt(mw_power_mw)
            else:
                warnings.warn(
                    "Cannot scale by microwave power ('P'): MWPW missing, zero, or invalid."
                )
        elif not is_cw and "P" in scaling:
            warnings.warn(
                "Microwave power scaling ('P') requested, but experiment is not CW."
            )

        if "T" in scaling:
            if (
                temperature_k is not None
            ):  # Allow T=0K ? MATLAB doesn't error but scaling makes no sense. Let's scale anyway.
                if temperature_k == 0:
                    warnings.warn(
                        "Temperature (STMP) is zero. Scaling by T will result in zero."
                    )
                data = data * temperature_k
            else:
                warnings.warn(
                    "Cannot scale by temperature ('T'): STMP missing or invalid."
                )

    # Parse string parameters to numbers where possible
    parameters = parse_field_params(parameters)

    return data, abscissa, parameters
