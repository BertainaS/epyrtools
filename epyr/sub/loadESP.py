"""Loader for Bruker ESP/WinEPR data (.spc, .par).

The public entry point is :func:`load`. Private helpers handle header
interpretation, dimension resolution, abscissa construction, and the
optional scaling pipeline.
"""

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .utils import (
    get_matrix,
    parse_field_params,
    read_par_file,
)

# Bits in the JSS flag, ESP/WinEPR convention
_JSS_BIT_COMPLEX = 1 << 4
_JSS_BIT_2D = 1 << 12

# Range parameters parsed as floats for abscissa construction
_RANGE_KEYS = ("HCF", "HSW", "GST", "GSI", "XXLB", "XXWI", "XYLB", "XYWI", "RCT")


def _resolve_extensions(file_extension: str) -> Tuple[str, str]:
    """Return (par_ext, spc_ext) with case matching ``file_extension``."""
    par, spc = ".par", ".spc"
    if file_extension.isupper():
        par, spc = par.upper(), spc.upper()
    return par, spc


def _parse_flags(parameters: Dict[str, Any]) -> Tuple[str, str, bool, bool]:
    """Resolve file type, endianness, complexity and 2D status from PAR header.

    Returns
    -------
    file_type : str
        ``'c'`` (ESP CW), ``'p'`` (pulse), or ``'w'`` (WinEPR).
    endian : str
        ``'ieee-be'`` or ``'ieee-le'``.
    is_complex : bool
    two_d : bool
    """
    file_type = "c"
    endian = "ieee-be"
    is_complex = False
    two_d = False

    if "DOS" in parameters:
        endian = "ieee-le"
        file_type = "w"

    if "JSS" in parameters:
        try:
            flags = int(parameters["JSS"])
            is_complex = bool(flags & _JSS_BIT_COMPLEX)
            two_d = bool(flags & _JSS_BIT_2D)
        except (ValueError, TypeError):
            warnings.warn(
                "Could not parse JSS flag in .par file. Assuming defaults.",
                stacklevel=3,
            )
    return file_type, endian, is_complex, two_d


def _try_int(parameters: Dict[str, Any], key: str) -> Optional[int]:
    """Parse parameters[key] as int; warn and return None on failure."""
    if key not in parameters:
        return None
    try:
        return int(parameters[key])
    except (ValueError, TypeError):
        warnings.warn(f"Could not parse {key} in .par file.", stacklevel=3)
        return None


def _apply_primary_dims(
    parameters: Dict[str, Any],
    file_type: str,
    is_complex: bool,
    two_d: bool,
) -> Tuple[int, int, str, Optional[int]]:
    """Apply ANZ, SSX, SSY (primary dimension keys). Returns (nx, ny, ftype, n_anz)."""
    nx, ny = 1024, 1

    n_anz = _try_int(parameters, "ANZ")
    if n_anz is not None and not two_d:
        if file_type == "c":
            file_type = "p"
        nx = n_anz // 2 if is_complex else n_anz

    ssx = _try_int(parameters, "SSX")
    if ssx is not None and (two_d or file_type == "p"):
        if file_type == "c":
            file_type = "p"
        nx = ssx // 2 if is_complex else ssx

    ssy = _try_int(parameters, "SSY")
    if ssy is not None and (two_d or file_type == "p"):
        if file_type == "c":
            file_type = "p"
        ny = ssy

    return nx, ny, file_type, n_anz


def _apply_legacy_dims(
    parameters: Dict[str, Any], nx: int, ny: int, two_d: bool, n_anz: Optional[int]
) -> Tuple[int, int]:
    """Apply legacy dimension keys (RES, REY, XPLS) only when primary keys absent."""
    if "RES" in parameters and not two_d and n_anz is None:
        res = _try_int(parameters, "RES")
        if res is not None:
            nx = res
    if "REY" in parameters and two_d and ny == 1:
        rey = _try_int(parameters, "REY")
        if rey is not None:
            ny = rey
    if "XPLS" in parameters and not two_d and n_anz is None and "RES" not in parameters:
        xpls = _try_int(parameters, "XPLS")
        if xpls is not None:
            nx = xpls
    return nx, ny


def _resolve_dimensions(
    parameters: Dict[str, Any], file_type: str, is_complex: bool, two_d: bool
) -> Tuple[int, int, str]:
    """Decide nx, ny, and refined file_type from dimension parameters.

    Bruker PAR files store the size under several conflicting keys: ANZ, SSX,
    SSY, RES, REY, XPLS. This function picks the right combination and warns
    on inconsistencies.
    """
    nx, ny, file_type, n_anz = _apply_primary_dims(
        parameters, file_type, is_complex, two_d
    )

    if two_d and n_anz is not None and nx * ny != n_anz:
        raise ValueError("Inconsistent 2D dimensions from ANZ, SSX, SSY in .par file.")
    if not two_d and n_anz is not None:
        expected_nx = n_anz // 2 if is_complex else n_anz
        if nx != expected_nx:
            warnings.warn(
                "ANZ conflicts with other dimension keys for 1D data. Using ANZ.",
                stacklevel=3,
            )
            nx = expected_nx

    nx, ny = _apply_legacy_dims(parameters, nx, ny, two_d, n_anz)
    return nx, ny, file_type


def _resolve_number_format(file_type: str) -> str:
    """Map file_type to the numpy dtype string for the .spc payload."""
    if file_type == "w":
        return "f4"
    if file_type in ("c", "p"):
        return "i4"
    warnings.warn(
        f"Unclear file type '{file_type}', assuming int32 (f4) data format.",
        stacklevel=3,
    )
    return "f4"


def _parse_range_params(parameters: Dict[str, Any]) -> Dict[str, Optional[float]]:
    """Convert the abscissa-related parameters to float (None on failure)."""
    out: Dict[str, Optional[float]] = {}
    for key in _RANGE_KEYS:
        if key in parameters:
            try:
                out[key] = float(parameters[key])
            except (ValueError, TypeError):
                out[key] = None
        else:
            out[key] = None
    return out


def _select_axis_source(
    params_num: Dict[str, Optional[float]], is_endor: bool, two_d: bool
) -> int:
    """Pick which parameter group to use for the X axis.

    Returns
    -------
    int
        ``1`` for GST/GSI, ``2`` for HCF/HSW, ``3`` for XXLB/XXWI, ``0`` if none.
    """

    def both(a: str, b: str) -> bool:
        return params_num.get(a) is not None and params_num.get(b) is not None

    has_gst_gsi = both("GST", "GSI")
    has_hcf_hsw = both("HCF", "HSW")
    has_xx = both("XXLB", "XXWI")

    if is_endor and has_gst_gsi:
        return 1
    if has_xx:
        return 3
    if has_hcf_hsw and has_gst_gsi:
        # Both present: MATLAB prefers GST/GSI
        return 1
    if has_hcf_hsw:
        return 2
    if has_gst_gsi:
        return 1

    # Last-chance fallback: HCF alone with implicit HSW=50G
    if params_num.get("HCF") is not None and params_num.get("HSW") is None:
        params_num["HSW"] = 50.0
        warnings.warn("HSW missing, assuming 50 G.", stacklevel=3)
        return 2
    return 0


def _build_abscissa(
    parameters: Dict[str, Any], nx: int, ny: int, two_d: bool
) -> Optional[Union[np.ndarray, List[np.ndarray]]]:
    """Construct the abscissa from PAR-file range parameters."""
    if nx <= 0:
        return None

    jex = parameters.get("JEX", "field-sweep").lower()
    is_endor = "endor" in jex
    is_time_sweep = "time-sweep" in jex

    params_num = _parse_range_params(parameters)

    if is_time_sweep and params_num.get("RCT") is not None:
        return np.arange(nx) * params_num["RCT"] / 1000.0  # seconds

    source = _select_axis_source(params_num, is_endor, two_d)
    x_axis: Optional[np.ndarray] = None
    y_axis: Optional[np.ndarray] = None

    if source == 1:
        gst, gsi = params_num["GST"], params_num["GSI"]
        x_axis = np.linspace(gst, gst + gsi, nx)
    elif source == 2:
        hcf, hsw = params_num["HCF"], params_num["HSW"]
        x_axis = np.linspace(hcf - hsw / 2, hcf + hsw / 2, nx)
    elif source == 3:
        xxlb, xxwi = params_num["XXLB"], params_num["XXWI"]
        x_axis = np.linspace(xxlb, xxlb + xxwi, nx)
        if (
            two_d
            and params_num.get("XYLB") is not None
            and params_num.get("XYWI") is not None
        ):
            xylb, xywi = params_num["XYLB"], params_num["XYWI"]
            y_axis = np.linspace(xylb, xylb + xywi, ny)

    if x_axis is None:
        warnings.warn(
            "Could not determine abscissa range from parameter file. Using indices.",
            stacklevel=3,
        )
        return np.arange(nx)
    if y_axis is not None:
        return [x_axis, y_axis]
    return x_axis


def _extract_scaling_params(parameters: Dict[str, Any]) -> Dict[str, Optional[float]]:
    """Parse PAR parameters needed for scaling (returns None on failure)."""
    spec = {
        "n_scans_done": ("JSD", int),
        "receiver_gain": ("RRG", float),
        "mw_power_mw": ("MP", float),
        "temperature_k": ("TE", float),
        "conversion_time_ms": ("RCT", float),
    }
    out: Dict[str, Optional[float]] = {}
    for key, (param_key, caster) in spec.items():
        if param_key not in parameters:
            out[key] = None
            continue
        try:
            out[key] = caster(parameters[param_key])
        except (ValueError, TypeError):
            out[key] = None
    return out


def _apply_power_sweep_scaling(
    data: np.ndarray,
    mw_power_mw: float,
    abscissa: Optional[Union[np.ndarray, List[np.ndarray]]],
) -> np.ndarray:
    """Scale a 2D power-sweep dataset by sqrt of microwave power per row.

    The Y-axis stores dB attenuation relative to ``mw_power_mw``. Rows where
    the resulting power is <= 1e-12 mW are left unscaled.
    """
    if not (isinstance(abscissa, list) and len(abscissa) == 2):
        warnings.warn(
            "Cannot apply power sweep scaling ('P'): abscissa missing or not 2D.",
            stacklevel=3,
        )
        return data

    y_axis = abscissa[1]
    if y_axis.size != data.shape[0]:
        warnings.warn(
            "Cannot apply power sweep scaling ('P'): Y-axis size mismatch.",
            stacklevel=3,
        )
        return data

    power_values_mw = mw_power_mw * (10.0 ** (-y_axis / 10.0))
    valid = power_values_mw > 1e-12
    if np.any(~valid):
        warnings.warn(
            "Some power values in power sweep are <= 0. Scaling skipped there.",
            stacklevel=3,
        )

    sqrt_power = np.sqrt(power_values_mw[valid])
    data[valid, :] = data[valid, :] / sqrt_power[:, np.newaxis]
    return data


def _apply_scaling(
    data: np.ndarray,
    scaling: str,
    parameters: Dict[str, Any],
    abscissa: Optional[Union[np.ndarray, List[np.ndarray]]],
    two_d: bool,
) -> np.ndarray:
    """Apply the requested scaling factors. See module docstring for codes."""
    p = _extract_scaling_params(parameters)
    is_power_sweep_y = "mw-power-sweep" in parameters.get("JEY", "").lower()

    if "n" in scaling:
        n = p["n_scans_done"]
        if n is not None and n > 0:
            data = data / n
        else:
            warnings.warn(
                "Cannot scale by number of scans ('n'): JSD missing or invalid.",
                stacklevel=3,
            )

    if "G" in scaling:
        gain = p["receiver_gain"]
        if gain is not None and gain != 0:
            data = data / gain
        else:
            warnings.warn(
                "Cannot scale by receiver gain ('G'): RRG missing or invalid.",
                stacklevel=3,
            )

    if "P" in scaling:
        mw = p["mw_power_mw"]
        if mw is not None and mw > 0:
            if is_power_sweep_y and two_d and data.ndim == 2:
                data = _apply_power_sweep_scaling(data, mw, abscissa)
            elif data.ndim <= 2:
                data = data / np.sqrt(mw)
            else:
                warnings.warn(
                    "Cannot apply power scaling ('P') to data with >2 dimensions.",
                    stacklevel=3,
                )
        else:
            warnings.warn(
                "Cannot scale by microwave power ('P'): MP missing or invalid.",
                stacklevel=3,
            )

    if "T" in scaling:
        t_k = p["temperature_k"]
        if t_k is not None and t_k > 0:
            data = data * t_k
        else:
            warnings.warn(
                "Cannot scale by temperature ('T'): TE missing or invalid.",
                stacklevel=3,
            )

    if "c" in scaling:
        t_ms = p["conversion_time_ms"]
        if t_ms is not None and t_ms > 0:
            data = data / t_ms
        else:
            warnings.warn(
                "Cannot scale by conversion time ('c'): RCT missing or invalid.",
                stacklevel=3,
            )
    return data


def load(full_base_name: Path, file_extension: str, scaling: str) -> tuple:
    """Load Bruker ESP/WinEPR data (.spc / .par).

    Parameters
    ----------
    full_base_name : pathlib.Path
        Path without extension.
    file_extension : str
        Original file extension (``.par`` / ``.spc``); only its case is used
        to pick the case of the companion file extensions.
    scaling : str
        Scaling specification (e.g., ``'nP G'``). See :func:`_apply_scaling`.

    Returns
    -------
    tuple
        ``(data, abscissa, parameters)``.
    """
    par_ext, spc_ext = _resolve_extensions(file_extension)
    par_file = Path(str(full_base_name) + par_ext)
    spc_file = Path(str(full_base_name) + spc_ext)

    parameters = read_par_file(par_file)

    file_type, endian, is_complex, two_d = _parse_flags(parameters)
    nx, ny, file_type = _resolve_dimensions(parameters, file_type, is_complex, two_d)
    number_format = _resolve_number_format(file_type)

    abscissa = _build_abscissa(parameters, nx, ny, two_d)

    # WinEPR sometimes saves a 2D slice as 1D: if JSS says no 2D, force ny=1.
    if not two_d and ny > 1:
        warnings.warn(
            f"Parameter file indicates ny={ny} but JSS does not. Assuming 1D data.",
            stacklevel=2,
        )
        ny = 1

    data = get_matrix(spc_file, [nx, ny, 1], number_format, endian, is_complex)

    if scaling and data is not None and data.size > 0:
        data = _apply_scaling(data, scaling, parameters, abscissa, two_d)

    parameters = parse_field_params(parameters)
    return data, abscissa, parameters
