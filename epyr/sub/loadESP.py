# sub/loadESP.py
import warnings
from pathlib import Path

import numpy as np

from .utils import (  # Import from sibling utils
    get_matrix,
    parse_field_params,
    read_par_file,
)


def load(full_base_name: Path, file_extension: str, scaling: str) -> tuple:
    """
    Loads Bruker ESP/WinEPR data (.spc, .par).

    Args:
        full_base_name: Path object without extension.
        file_extension: The original file extension (e.g., '.spc', '.par').
        scaling: Scaling string (e.g., 'nP G').

    Returns:
        tuple: (data, abscissa, parameters)
    """
    # Determine PAR and SPC file extensions, respecting case
    par_extension = ".par"
    spc_extension = ".spc"
    if file_extension.isupper():
        par_extension = par_extension.upper()
        spc_extension = spc_extension.upper()

    # Use string concatenation instead of with_suffix() to handle filenames with multiple dots
    par_file = Path(str(full_base_name) + par_extension)
    spc_file = Path(str(full_base_name) + spc_extension)

    # Read parameter file
    parameters = read_par_file(par_file)

    # --- Determine file type, endianness, dimensions, complexity ---
    file_type = "c"  # Default: ESP cw EPR
    two_d = False
    is_complex = False
    nx = 1024
    ny = 1
    endian = "ieee-be"  # Default: big-endian for older ESP

    if "DOS" in parameters:
        endian = "ieee-le"
        file_type = "w"  # Windows WinEPR

    # JSS: Job Status and Spectrometer status flags
    if "JSS" in parameters:
        try:
            flags = int(parameters["JSS"])
            is_complex = bool(flags & (1 << 4))  # Bit 5 (0-indexed)
            two_d = bool(flags & (1 << 12))  # Bit 13
        except (ValueError, TypeError):
            warnings.warn("Could not parse JSS flag in .par file. Assuming defaults.")

    # Get dimensions from various possible keys
    n_anz = None
    if "ANZ" in parameters:  # Total number of points
        try:
            n_anz = int(parameters["ANZ"])
            if not two_d:
                if file_type == "c":
                    file_type = "p"  # Assume pulse if ANZ present and not 2D
                nx = n_anz // 2 if is_complex else n_anz
            # Consistency check if 2D will happen later
        except (ValueError, TypeError):
            warnings.warn("Could not parse ANZ in .par file.")

    # SSX/SSY preferred for 2D
    if "SSX" in parameters:  # X points (often total points for 1D pulse)
        try:
            ssx_val = int(parameters["SSX"])
            if two_d or file_type == "p":
                if file_type == "c":
                    file_type = "p"  # Treat as pulse if SSX/SSY found
                nx = ssx_val // 2 if is_complex else ssx_val
        except (ValueError, TypeError):
            warnings.warn("Could not parse SSX in .par file.")

    if "SSY" in parameters:  # Y points
        try:
            if two_d or file_type == "p":
                if file_type == "c":
                    file_type = "p"
                ny = int(parameters["SSY"])
        except (ValueError, TypeError):
            warnings.warn("Could not parse SSY in .par file.")

    if two_d and n_anz is not None and nx * ny != n_anz:
        raise ValueError("Inconsistent 2D dimensions from ANZ, SSX, SSY in .par file.")
    elif not two_d and n_anz is not None:
        # ANZ overrides RES/XPLS if not 2D
        if nx != (n_anz // 2 if is_complex else n_anz):
            warnings.warn(
                "ANZ value conflicts with other dimension keys (RES/XPLS) for 1D data. Using ANZ."
            )
            nx = n_anz // 2 if is_complex else n_anz

    # Older keys (RES, REY, XPLS) - might override if SSX/SSY/ANZ absent or inconsistent
    if "RES" in parameters and not two_d and n_anz is None:  # X points (CW)
        try:
            nx = int(parameters["RES"])
        except (ValueError, TypeError):
            warnings.warn("Could not parse RES in .par file.")
    if "REY" in parameters and two_d and ny == 1:  # Y points (check if SSY was missing)
        try:
            ny = int(parameters["REY"])
        except (ValueError, TypeError):
            warnings.warn("Could not parse REY in .par file.")
    if (
        "XPLS" in parameters and not two_d and n_anz is None and "RES" not in parameters
    ):  # X points (Pulse?)
        try:
            nx = int(parameters["XPLS"])
        except (ValueError, TypeError):
            warnings.warn("Could not parse XPLS in .par file.")

    # Determine number format
    if file_type == "w":  # WinEPR/Simfonia
        number_format_code = "f4"
    elif file_type in ["c", "p"]:  # ESP CW or Pulse
        number_format_code = "i4"
    else:
        number_format_code = "f4"  # Default fallback
        warnings.warn(
            f"Unclear file type '{file_type}', assuming int32 (f4) data format."
        )

    # --- Construct Abscissa ---
    abscissa = None
    if nx > 0:  # Changed from >1 to >0 to handle single-point data
        # Get experiment type hints
        jex = parameters.get("JEX", "field-sweep").lower()
        jey = parameters.get("JEY", "").lower()
        is_endor = "endor" in jex
        is_time_sweep = "time-sweep" in jex
        is_power_sweep_y = "mw-power-sweep" in jey

        # Convert range parameters safely
        params_num = {}
        range_keys = ["HCF", "HSW", "GST", "GSI", "XXLB", "XXWI", "XYLB", "XYWI", "RCT"]
        for key in range_keys:
            if key in parameters:
                try:
                    params_num[key] = float(parameters[key])
                except (ValueError, TypeError):
                    params_num[key] = None  # Keep as None if conversion fails
            else:
                params_num[key] = None

        # Decide which parameters to use (TakeGH logic)
        take_gh = 0  # 0: Undetermined, 1: GST/GSI, 2: HCF/HSW, 3: XX/XYLB/WI
        if (
            is_endor
            and params_num.get("GST") is not None
            and params_num.get("GSI") is not None
        ):
            take_gh = 1
        elif params_num.get("XXLB") is not None and params_num.get("XXWI") is not None:
            # If XY are also present, assume 2D and use XX/XY
            if (
                params_num.get("XYLB") is not None
                and params_num.get("XYWI") is not None
                and two_d
            ):
                take_gh = 3
            # If only XX present, or if XY present but not 2D, still prefer XX for X axis
            elif not two_d or (
                params_num.get("XYLB") is None or params_num.get("XYWI") is None
            ):
                take_gh = 3  # Use XX for X-axis (might just be 1D)
        elif (
            params_num.get("HCF") is not None
            and params_num.get("HSW") is not None
            and params_num.get("GST") is not None
            and params_num.get("GSI") is not None
        ):
            take_gh = 1  # Prefer GST/GSI if all are present (MATLAB comment)
        elif params_num.get("HCF") is not None and params_num.get("HSW") is not None:
            take_gh = 2
        elif params_num.get("GST") is not None and params_num.get("GSI") is not None:
            take_gh = 1
        # Fallbacks / assumptions if primary methods fail
        elif (
            params_num.get("GSI") is None
            and params_num.get("HSW") is None
            and params_num.get("HCF") is not None
        ):
            params_num["HSW"] = 50.0  # Assume HSW=50 G if missing (like MATLAB code)
            warnings.warn("HSW missing, assuming 50 G.")
            take_gh = 2
        elif (
            params_num.get("HCF") is None
            and params_num.get("XXLB") is not None
            and params_num.get("XXWI") is not None
        ):
            take_gh = 3  # If HCF missing but XX present, use XX
        elif params_num.get("XXLB") is not None and params_num.get("XXWI") is not None:
            take_gh = 3  # Fallback to XX if other methods failed

        # Build abscissa based on decision
        if is_time_sweep and params_num.get("RCT") is not None:
            conversion_time_ms = params_num["RCT"]
            abscissa = np.arange(nx) * conversion_time_ms / 1000.0  # Time in seconds
        else:
            x_axis = None
            y_axis = None
            if (
                take_gh == 1
                and params_num.get("GST") is not None
                and params_num.get("GSI") is not None
            ):
                gst, gsi = params_num["GST"], params_num["GSI"]
                # linspace includes endpoint, MATLAB GSI*linspace(0,1,nx) implies GSI is width
                # Correct interpretation: gst is start, gsi is increment
                # x_axis = gst + gsi * np.arange(nx) # If GSI is increment
                # If GSI is total width (like HSW):
                x_axis = np.linspace(gst, gst + gsi, nx)
            elif (
                take_gh == 2
                and params_num.get("HCF") is not None
                and params_num.get("HSW") is not None
            ):
                hcf, hsw = params_num["HCF"], params_num["HSW"]
                x_axis = np.linspace(hcf - hsw / 2, hcf + hsw / 2, nx)
            elif (
                take_gh == 3
                and params_num.get("XXLB") is not None
                and params_num.get("XXWI") is not None
            ):
                xxlb, xxwi = params_num["XXLB"], params_num["XXWI"]
                x_axis = np.linspace(xxlb, xxlb + xxwi, nx)
                if (
                    two_d
                    and params_num.get("XYLB") is not None
                    and params_num.get("XYWI") is not None
                ):
                    xylb, xywi = params_num["XYLB"], params_num["XYWI"]
                    y_axis = np.linspace(xylb, xylb + xywi, ny)

            if x_axis is not None:
                if y_axis is not None:
                    abscissa = [x_axis, y_axis]
                else:
                    abscissa = x_axis
            else:
                warnings.warn(
                    "Could not determine abscissa range from parameter file. Abscissa set to None."
                )
                abscissa = np.arange(nx)  # Default to indices if range fails

    # Slice of 2D data might be saved as 1D (WinEPR behavior mentioned in MATLAB)
    # Ensure ny=1 if not explicitly 2D according to JSS flag
    if not two_d and ny > 1:
        # This might happen if REY was present but JSS didn't have 2D flag
        warnings.warn(
            f"Parameter file indicates ny={ny} (REY?) but JSS flag doesn't indicate 2D. Assuming 1D data."
        )
        ny = 1

    dimensions = [nx, ny, 1]  # Use [nx, ny, nz=1] for get_matrix

    # --- Read Data File ---

    data = get_matrix(spc_file, dimensions, number_format_code, endian, is_complex)

    # --- Scale Data ---
    if scaling and data is not None and data.size > 0:
        # Convert parameters needed for scaling, handle missing ones
        n_scans_done = None
        receiver_gain = None
        mw_power_mw = None
        temperature_k = None
        conversion_time_ms = None

        if "JSD" in parameters:
            try:
                n_scans_done = int(parameters["JSD"])
            except (ValueError, TypeError):
                pass
        if "RRG" in parameters:
            try:
                receiver_gain = float(parameters["RRG"])
            except (ValueError, TypeError):
                pass
        if "MP" in parameters:  # Power in mW
            try:
                mw_power_mw = float(parameters["MP"])
            except (ValueError, TypeError):
                pass
        if "TE" in parameters:  # Temperature in K
            try:
                temperature_k = float(parameters["TE"])
            except (ValueError, TypeError):
                pass
        if "RCT" in parameters:  # Conversion time in ms
            try:
                conversion_time_ms = float(parameters["RCT"])
            except (ValueError, TypeError):
                pass

        # Apply scaling factors
        if "n" in scaling:
            if n_scans_done is not None and n_scans_done > 0:
                data = data / n_scans_done
            else:
                warnings.warn(
                    "Cannot scale by number of scans ('n'): JSD missing, zero, or invalid."
                )
        if "G" in scaling:
            if receiver_gain is not None and receiver_gain != 0:
                data = data / receiver_gain
            else:
                warnings.warn(
                    "Cannot scale by receiver gain ('G'): RRG missing, zero, or invalid."
                )
        if "P" in scaling:
            if mw_power_mw is not None and mw_power_mw > 0:
                if is_power_sweep_y and two_d and data.ndim == 2:
                    # Assume power is along the second dimension (y-axis, index 0 in numpy)
                    if (
                        abscissa is not None
                        and isinstance(abscissa, list)
                        and len(abscissa) == 2
                    ):
                        y_axis_params = abscissa[
                            1
                        ]  # Should be dB values if XYWI/XYLB were used
                        if (
                            y_axis_params.size == data.shape[0]
                        ):  # ny should be first dim
                            power_factors_db = y_axis_params
                            # Assuming y_axis_params stores dB ATTENUATION relative to MP power
                            # Power(i) = MP * 10^(-dB(i)/10)
                            power_values_mw = mw_power_mw * (
                                10.0 ** (-power_factors_db / 10.0)
                            )
                            # Avoid division by zero or sqrt of negative
                            valid_power = (
                                power_values_mw > 1e-12
                            )  # Threshold for valid power
                            if np.any(~valid_power):
                                warnings.warn(
                                    "Some power values in power sweep are <= 0. Scaling skipped for those points."
                                )
                            # Scale each row (spectrum at specific power)
                            # Need to broadcast sqrt(power) correctly across columns (nx)
                            sqrt_power = np.sqrt(power_values_mw[valid_power])
                            # Apply scaling row-wise where power is valid
                            data[valid_power, :] = (
                                data[valid_power, :] / sqrt_power[:, np.newaxis]
                            )

                        else:
                            warnings.warn(
                                "Cannot apply power sweep scaling ('P'): Y-axis data mismatch or missing."
                            )
                    else:
                        warnings.warn(
                            "Cannot apply power sweep scaling ('P'): Abscissa data missing or not 2D."
                        )

                elif (
                    data.ndim <= 2
                ):  # Apply simple scaling for 1D or non-power-sweep 2D
                    data = data / np.sqrt(mw_power_mw)
                else:
                    warnings.warn(
                        "Cannot apply power scaling ('P') to data with >2 dimensions."
                    )
            else:
                warnings.warn(
                    "Cannot scale by microwave power ('P'): MP missing, zero, or invalid."
                )
        if "T" in scaling:
            if temperature_k is not None and temperature_k > 0:
                data = data * temperature_k
            else:
                warnings.warn(
                    "Cannot scale by temperature ('T'): TE missing, zero, or invalid."
                )
        if "c" in scaling:
            if conversion_time_ms is not None and conversion_time_ms > 0:
                data = data / conversion_time_ms
            else:
                warnings.warn(
                    "Cannot scale by conversion time ('c'): RCT missing, zero, or invalid."
                )

    # Parse string parameters to numbers where possible
    parameters = parse_field_params(parameters)

    return data, abscissa, parameters
