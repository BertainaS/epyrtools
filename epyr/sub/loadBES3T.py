"""Loader for Bruker BES3T data (.DTA, .DSC).

The public entry point is :func:`load`. Private helpers handle the four
stages: file resolution, header parsing, abscissa construction, and
optional scaling.
"""

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .utils import (
    get_matrix,
    parse_field_params,
    read_dsc_file,
)

# Format mappings (Bruker code → numpy dtype string)
_DATA_FORMAT_MAP = {
    "C": "i1",  # int8
    "S": "i2",  # int16
    "I": "i4",  # int32
    "F": "f4",  # float32
    "D": "f8",  # float64
}
_COMPANION_FORMAT_MAP = {"D": "f8", "F": "f4", "I": "i4", "S": "i2"}


def _resolve_extensions(file_extension: str) -> Tuple[str, str]:
    """Return (dsc_ext, dta_ext) with case matching ``file_extension``."""
    dsc, dta = ".dsc", ".dta"
    if file_extension.isupper():
        dsc, dta = dsc.upper(), dta.upper()
    return dsc, dta


def _parse_complexity(parameters: Dict[str, Any]) -> Tuple[np.ndarray, int]:
    """Extract complexity flags and number of data values from IKKF."""
    if "IKKF" in parameters:
        parts = parameters["IKKF"].split(",")
        is_complex = np.array([p.strip().upper() == "CPLX" for p in parts])
        return is_complex, len(parts)

    warnings.warn("IKKF not found in .DSC file. Assuming IKKF=REAL.", stacklevel=3)
    return np.array([False]), 1


def _parse_dimensions(parameters: Dict[str, Any]) -> List[int]:
    """Return ``[nx, ny, nz]`` from XPTS, YPTS, ZPTS."""
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
    return [nx, ny, nz]


def _parse_byte_order(parameters: Dict[str, Any]) -> str:
    """Return ``'ieee-be'`` or ``'ieee-le'`` from BSEQ."""
    if "BSEQ" not in parameters:
        warnings.warn(
            "BSEQ not found in .DSC file. Assuming BSEQ=BIG (big-endian).",
            stacklevel=3,
        )
        return "ieee-be"

    bseq_val = parameters["BSEQ"].upper()
    if bseq_val == "BIG":
        return "ieee-be"
    if bseq_val == "LIT":
        return "ieee-le"
    raise ValueError(f"Unknown BSEQ value '{parameters['BSEQ']}' in .DSC file.")


def _parse_number_format(parameters: Dict[str, Any], is_complex: np.ndarray) -> str:
    """Resolve IRFMT (and check IIFMT consistency) to a numpy dtype string."""
    if "IRFMT" not in parameters:
        raise ValueError("IRFMT keyword not found in .DSC file.")

    irfmt_val = parameters["IRFMT"].split(",")[0].strip().upper()
    if irfmt_val in _DATA_FORMAT_MAP:
        number_format = _DATA_FORMAT_MAP[irfmt_val]
    elif irfmt_val in ("A", "0", "N"):
        raise ValueError(f"Unsupported or no data format IRFMT='{irfmt_val}' in DSC.")
    else:
        raise ValueError(f"Unknown IRFMT value '{irfmt_val}' in .DSC file.")

    if "IIFMT" in parameters and np.any(is_complex):
        iifmt_val = parameters["IIFMT"].split(",")[0].strip().upper()
        if iifmt_val != irfmt_val:
            warnings.warn(
                "IRFMT and IIFMT differ in DSC file. Using IRFMT for reading.",
                stacklevel=3,
            )
    return number_format


def _read_companion_axis(
    companion_file: Path, dim_size: int, fmt_char: str, byte_order: str
) -> Optional[np.ndarray]:
    """Read a non-linear axis from a companion ``.XGF/.YGF/.ZGF`` file.

    Returns ``None`` if the file is missing, the format is unknown, or the
    payload size does not match ``dim_size``.
    """
    if fmt_char not in _COMPANION_FORMAT_MAP:
        warnings.warn(
            f"Cannot read companion file format '{fmt_char}'. Assuming linear.",
            stacklevel=3,
        )
        return None

    if not companion_file.is_file():
        warnings.warn(
            f"Companion file {companion_file} not found. Assuming linear axis.",
            stacklevel=3,
        )
        return None

    dt_char = ">" if byte_order == "ieee-be" else "<"
    dtype = np.dtype(f"{dt_char}{_COMPANION_FORMAT_MAP[fmt_char]}")
    try:
        axis_data = np.fromfile(companion_file, dtype=dtype, count=dim_size)
    except Exception as e:
        warnings.warn(
            f"Error reading companion file {companion_file}: {e}. Assuming linear.",
            stacklevel=3,
        )
        return None

    if axis_data.size != dim_size:
        warnings.warn(
            f"Could not read expected {dim_size} values from {companion_file}.",
            stacklevel=3,
        )
        return None
    return axis_data


def _build_linear_axis(
    parameters: Dict[str, Any], axis: str, dim_size: int
) -> np.ndarray:
    """Build a linear axis from XMIN/XWID (and equivalents for Y, Z)."""
    min_key = f"{axis}MIN"
    wid_key = f"{axis}WID"
    try:
        minimum = float(parameters[min_key])
        width = float(parameters[wid_key])
    except (KeyError, ValueError, TypeError):
        warnings.warn(
            f"Could not read MIN/WID for axis {axis}. Using default index.",
            stacklevel=3,
        )
        return np.arange(dim_size)

    if dim_size == 1:
        return np.array([minimum])
    if dim_size <= 0:
        return np.array([])
    if width == 0:
        warnings.warn(
            f"{axis} range has zero width (WID=0). Using index 0 to N-1.",
            stacklevel=3,
        )
        return np.arange(dim_size)
    return np.linspace(minimum, minimum + width, dim_size)


def _build_abscissa(
    parameters: Dict[str, Any],
    dimensions: List[int],
    byte_order: str,
    full_base_name: Path,
    file_extension: str,
) -> Optional[Union[np.ndarray, List[np.ndarray]]]:
    """Construct the abscissa for each defined axis.

    Returns a single array (1D), a list of arrays (multi-D), or ``None``
    if no axis could be constructed.
    """
    abscissa_list: List[Optional[np.ndarray]] = [None, None, None]
    upper_ext = file_extension.isupper()

    for i, axis in enumerate(["X", "Y", "Z"]):
        dim_size = dimensions[i]
        if dim_size <= 1:
            continue

        axis_type = parameters.get(f"{axis}TYP", "IDX")

        if axis_type == "IGD":
            suffix = f".{axis}GF"
            if upper_ext:
                suffix = suffix.upper()
            companion = Path(str(full_base_name) + suffix)
            fmt_char = parameters.get(f"{axis}FMT", "D").upper()
            axis_data = _read_companion_axis(companion, dim_size, fmt_char, byte_order)
            if axis_data is not None:
                abscissa_list[i] = axis_data
                continue
            axis_type = "IDX"  # Fallback

        if axis_type == "IDX":
            abscissa_list[i] = _build_linear_axis(parameters, axis, dim_size)
        elif axis_type == "NTUP":
            raise NotImplementedError("Cannot read data with NTUP axes.")

    defined = [a for a in abscissa_list if a is not None]
    if len(defined) == 0:
        return None
    if len(defined) == 1:
        return defined[0]
    return defined


def _extract_scaling_params(parameters: Dict[str, Any]) -> Dict[str, Optional[float]]:
    """Parse DSC parameters needed for scaling (returns None if missing)."""
    out: Dict[str, Optional[float]] = {}
    spec = {
        "n_averages": ("AVGS", int),
        "receiver_gain_db": ("RCAG", float),
        "sampling_time_s": ("SPTP", float),
        "mw_power_w": ("MWPW", float),
        "temperature_k": ("STMP", float),
    }
    for key, (param_key, caster) in spec.items():
        try:
            out[key] = caster(parameters.get(param_key))
        except (ValueError, TypeError, KeyError):
            out[key] = None
    return out


def _apply_scaling(
    data: np.ndarray, scaling: str, parameters: Dict[str, Any]
) -> np.ndarray:
    """Apply the requested scaling factors to ``data``.

    Each character in ``scaling`` triggers one operation: ``'n'`` (averages),
    ``'P'`` (microwave power, CW only), ``'G'`` (receiver gain, CW only),
    ``'T'`` (temperature), ``'c'`` (sampling time, CW only).
    """
    p = _extract_scaling_params(parameters)
    expt_type = parameters.get("EXPT", "CW").upper()
    is_cw = expt_type == "CW"
    data_prescaled = parameters.get("SctNorm", "false").lower() == "true"

    if "n" in scaling:
        n_avg = p["n_averages"]
        if n_avg is not None and n_avg > 0:
            if data_prescaled:
                warnings.warn(
                    f"Cannot scale by number of scans ('n'): Data already averaged"
                    f" (SctNorm=true, AVGS={n_avg}).",
                    stacklevel=3,
                )
            else:
                data = data / n_avg
        else:
            warnings.warn(
                "Cannot scale by number of scans ('n'): AVGS missing or invalid.",
                stacklevel=3,
            )

    if is_cw and "G" in scaling:
        gain_db = p["receiver_gain_db"]
        if gain_db is not None:
            gain = 10 ** (gain_db / 20.0)
            if gain != 0:
                data = data / gain
            else:
                warnings.warn(
                    "Cannot scale by receiver gain ('G'): gain is zero.", stacklevel=3
                )
        else:
            warnings.warn(
                "Cannot scale by receiver gain ('G'): RCAG missing or invalid.",
                stacklevel=3,
            )

    if is_cw and "c" in scaling:
        t_s = p["sampling_time_s"]
        if t_s is not None and t_s > 0:
            data = data / (t_s * 1000.0)
        else:
            warnings.warn(
                "Cannot scale by conversion time ('c'): SPTP missing or invalid.",
                stacklevel=3,
            )

    if "P" in scaling:
        if is_cw:
            p_w = p["mw_power_w"]
            if p_w is not None and p_w > 0:
                data = data / np.sqrt(p_w * 1000.0)
            else:
                warnings.warn(
                    "Cannot scale by microwave power ('P'): MWPW missing or invalid.",
                    stacklevel=3,
                )
        else:
            warnings.warn(
                "Microwave power scaling ('P') requested, but experiment is not CW.",
                stacklevel=3,
            )

    if "T" in scaling:
        t_k = p["temperature_k"]
        if t_k is not None:
            if t_k == 0:
                warnings.warn(
                    "Temperature (STMP) is zero. Scaling by T will result in zero.",
                    stacklevel=3,
                )
            data = data * t_k
        else:
            warnings.warn(
                "Cannot scale by temperature ('T'): STMP missing or invalid.",
                stacklevel=3,
            )
    return data


def load(full_base_name: Path, file_extension: str, scaling: str) -> tuple:
    """Load Bruker BES3T data (.DTA / .DSC).

    Parameters
    ----------
    full_base_name : pathlib.Path
        Path without extension.
    file_extension : str
        Original file extension (e.g., ``'.dta'`` or ``'.DSC'``); only its
        case is used to pick the case of the companion file extensions.
    scaling : str
        Scaling specification (e.g., ``'nP G'``). See :func:`_apply_scaling`.

    Returns
    -------
    tuple
        ``(data, abscissa, parameters)``.
    """
    dsc_ext, dta_ext = _resolve_extensions(file_extension)
    dsc_file = Path(str(full_base_name) + dsc_ext)
    dta_file = Path(str(full_base_name) + dta_ext)

    parameters = read_dsc_file(dsc_file)

    is_complex, n_data_values = _parse_complexity(parameters)
    dimensions = _parse_dimensions(parameters)
    byte_order = _parse_byte_order(parameters)
    number_format = _parse_number_format(parameters, is_complex)

    abscissa = _build_abscissa(
        parameters, dimensions, byte_order, full_base_name, file_extension
    )

    if n_data_values > 1:
        warnings.warn(
            f"DSC file indicates {n_data_values} data values per point (IKKF)."
            " Only reading the first value.",
            stacklevel=2,
        )

    data = get_matrix(dta_file, dimensions, number_format, byte_order, is_complex[0])

    if scaling and data is not None and data.size > 0:
        data = _apply_scaling(data, scaling, parameters)

    parameters = parse_field_params(parameters)
    return data, abscissa, parameters
