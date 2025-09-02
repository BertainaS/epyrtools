"""
Format-specific export functions for FAIR data conversion.

This module contains functions to export EPR data and metadata to various
FAIR-compliant formats including CSV/JSON and HDF5.
"""

import csv
import json
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Union

import numpy as np

try:
    import h5py

    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

from .data_processing import process_parameters


def save_to_csv_json(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Save data to CSV and structured metadata to JSON.

    Args:
        output_basename: Base path for output files (without extension)
        x: Abscissa data array(s) or None
        y: Intensity data array
        pars: Raw parameters dictionary
        original_file_path: Path to original data file
    """
    json_file = output_basename.with_suffix(".json")
    csv_file = output_basename.with_suffix(".csv")

    fair_meta, unmapped_meta = process_parameters(pars)

    print(f"  Saving structured metadata to: {json_file}")

    # Save metadata to JSON
    metadata_to_save = {
        "original_file": original_file_path,
        "fair_metadata": fair_meta,
    }
    if unmapped_meta:
        metadata_to_save["unmapped_parameters"] = unmapped_meta

    try:
        with open(json_file, "w", encoding="utf-8") as f:
            json.dump(metadata_to_save, f, indent=4, default=str)
    except IOError as e:
        warnings.warn(f"Could not write JSON file {json_file}: {e}")
    except TypeError as e:
        warnings.warn(
            f"Error serializing metadata to JSON for {json_file}: {e}. "
            f"Some parameters might not be saved correctly."
        )

    print(f"  Saving data to: {csv_file}")

    # Save data to CSV
    try:
        with open(csv_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)

            # Write header
            writer.writerow(["# EPR Data Export"])
            writer.writerow(["# Original File:", original_file_path])

            # Add key parameters from FAIR metadata
            mwfq_info = fair_meta.get("microwave_frequency", {})
            field_info = fair_meta.get(
                "field_center", fair_meta.get("field_sweep_start", {})
            )
            sweep_info = fair_meta.get(
                "field_sweep_width", fair_meta.get("field_sweep_increment", {})
            )

            writer.writerow(
                [
                    "# Microwave_Frequency:",
                    f"{mwfq_info.get('value', 'N/A')} {mwfq_info.get('unit', '')}".strip(),
                ]
            )
            writer.writerow(
                [
                    "# Field_Center/Start:",
                    f"{field_info.get('value', 'N/A')} {field_info.get('unit', '')}".strip(),
                ]
            )
            writer.writerow(
                [
                    "# Field_Sweep/Increment:",
                    f"{sweep_info.get('value', 'N/A')} {sweep_info.get('unit', '')}".strip(),
                ]
            )
            writer.writerow(["# Data_Shape:", str(y.shape)])
            writer.writerow(["# Data_Type:", str(y.dtype)])
            writer.writerow(["# ---"])

            # Prepare data columns
            header_row = []
            data_columns = []
            is_complex = np.iscomplexobj(y)
            is_2d = y.ndim == 2

            # Get axis units from FAIR metadata
            x_unit_val = fair_meta.get("x_axis_unit", {}).get("value", "a.u.")
            y_unit_val = fair_meta.get("y_axis_unit", {}).get("value", "a.u.")
            if isinstance(y_unit_val, str) and "," in y_unit_val:
                y_unit_val = y_unit_val.split(",")[0].strip()

            if not is_2d:  # 1D Data
                n_pts = y.shape[0]

                # Abscissa column
                if x is not None and isinstance(x, np.ndarray) and x.shape == y.shape:
                    header_row.append(f"Abscissa ({x_unit_val})")
                    data_columns.append(x)
                else:
                    header_row.append("Index")
                    data_columns.append(np.arange(n_pts))
                    if x is not None:
                        warnings.warn(
                            "Provided x-axis ignored for CSV (shape mismatch or not ndarray). Using index."
                        )

                # Intensity columns
                if is_complex:
                    header_row.extend(
                        ["Intensity_Real (a.u.)", "Intensity_Imag (a.u.)"]
                    )
                    data_columns.append(np.real(y))
                    data_columns.append(np.imag(y))
                else:
                    header_row.append("Intensity (a.u.)")
                    data_columns.append(y)

            else:  # 2D Data - use "long" format (X, Y, Value(s))
                ny, nx = y.shape
                x_coords_flat = np.arange(nx)
                y_coords_flat = np.arange(ny)
                header_row.extend([f"X_Index ({nx} points)", f"Y_Index ({ny} points)"])

                # Determine X and Y axes from input 'x'
                if isinstance(x, list) and len(x) >= 2:
                    x_axis, y_axis = x[0], x[1]
                    if isinstance(x_axis, np.ndarray) and x_axis.size == nx:
                        x_coords_flat = x_axis
                        header_row[0] = f"X_Axis ({x_unit_val})"
                    if isinstance(y_axis, np.ndarray) and y_axis.size == ny:
                        y_coords_flat = y_axis
                        header_row[1] = f"Y_Axis ({y_unit_val})"
                elif isinstance(x, np.ndarray) and x.ndim == 1 and x.size == nx:
                    x_coords_flat = x
                    header_row[0] = f"X_Axis ({x_unit_val})"

                # Create grid and flatten
                xx, yy = np.meshgrid(x_coords_flat, y_coords_flat)
                data_columns.append(xx.ravel())
                data_columns.append(yy.ravel())

                # Intensity columns
                if is_complex:
                    header_row.extend(
                        ["Intensity_Real (a.u.)", "Intensity_Imag (a.u.)"]
                    )
                    data_columns.append(np.real(y).ravel())
                    data_columns.append(np.imag(y).ravel())
                else:
                    header_row.append("Intensity (a.u.)")
                    data_columns.append(y.ravel())

            # Write data
            writer.writerow(header_row)
            rows_to_write = np.stack(data_columns, axis=-1)
            writer.writerows(rows_to_write)

    except IOError as e:
        warnings.warn(f"Could not write CSV file {csv_file}: {e}")
    except Exception as e:
        warnings.warn(f"An unexpected error occurred while writing CSV {csv_file}: {e}")


def _try_set_h5_attr(h5_object, key: str, value: Any):
    """Helper to safely set HDF5 attributes, converting to string on type error."""
    try:
        if value is None:
            h5_object.attrs[key] = "None"
        elif isinstance(value, (list, tuple)) and all(
            isinstance(i, (int, float, str, np.number, bytes)) for i in value
        ):
            try:
                h5_object.attrs[key] = value
            except TypeError:
                try:
                    h5_object.attrs[key] = np.array(value)
                except TypeError:
                    h5_object.attrs[key] = str(value)
        elif isinstance(value, Path):
            h5_object.attrs[key] = str(value)
        else:
            h5_object.attrs[key] = value

    except TypeError:
        warnings.warn(
            f"Could not store attribute '{key}' (type: {type(value)}) "
            f"directly in HDF5 attributes. Converting to string."
        )
        h5_object.attrs[key] = str(value)
    except Exception as e:
        warnings.warn(
            f"Unexpected error storing attribute '{key}': {type(e).__name__} - {e}. Skipping."
        )


def save_to_hdf5(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Save data and structured metadata to an HDF5 file.

    Args:
        output_basename: Base path for output files (without extension)
        x: Abscissa data array(s) or None
        y: Intensity data array
        pars: Raw parameters dictionary
        original_file_path: Path to original data file
    """
    if not HAS_H5PY:
        warnings.warn(
            "h5py library not found. Skipping HDF5 output. Install with 'pip install h5py'"
        )
        return

    h5_file = output_basename.with_suffix(".h5")
    fair_meta, unmapped_meta = process_parameters(pars)

    print(f"  Saving structured data and metadata to: {h5_file}")

    try:
        with h5py.File(h5_file, "w") as f:
            # Store global metadata
            f.attrs["original_file"] = original_file_path
            f.attrs[
                "description"
            ] = "FAIR representation of EPR data converted from Bruker format."
            f.attrs["conversion_timestamp"] = datetime.now().isoformat()
            f.attrs["converter_script_version"] = "epyr_fair_converter_v1.0"

            # Store structured FAIR metadata
            param_grp = f.create_group("metadata/parameters_fair")
            param_grp.attrs[
                "description"
            ] = "Mapped parameters with units and descriptions."

            for fair_key, info in fair_meta.items():
                item_grp = param_grp.create_group(fair_key)
                _try_set_h5_attr(item_grp, "value", info["value"])
                _try_set_h5_attr(item_grp, "unit", info["unit"])
                _try_set_h5_attr(item_grp, "description", info["description"])

            # Store unmapped parameters
            if unmapped_meta:
                unmap_grp = f.create_group("metadata/parameters_original")
                unmap_grp.attrs[
                    "description"
                ] = "Parameters from the original file not found in the FAIR mapping."
                for key, value in unmapped_meta.items():
                    _try_set_h5_attr(unmap_grp, key, value)

            # Store data
            data_grp = f.create_group("data")
            ds_y = data_grp.create_dataset("intensity", data=y)
            ds_y.attrs["description"] = "Experimental intensity data."
            ds_y.attrs["units"] = "a.u."
            if np.iscomplexobj(y):
                ds_y.attrs["signal_type"] = "complex"
            else:
                ds_y.attrs["signal_type"] = "real"

            # Get axis units from FAIR metadata
            x_unit_val = fair_meta.get("x_axis_unit", {}).get("value", "a.u.")
            y_unit_val = fair_meta.get("y_axis_unit", {}).get("value", "a.u.")
            if isinstance(y_unit_val, str) and "," in y_unit_val:
                y_unit_val = y_unit_val.split(",")[0].strip()
            z_unit_val = fair_meta.get("z_axis_unit", {}).get("value", "a.u.")
            if isinstance(z_unit_val, str) and "," in z_unit_val:
                z_unit_val = z_unit_val.split(",")[0].strip()

            # Store abscissa data
            axis_datasets = {}
            if x is None:
                if y.ndim >= 1:
                    nx = y.shape[-1]
                    ds_x = data_grp.create_dataset("abscissa_x", data=np.arange(nx))
                    _try_set_h5_attr(ds_x, "units", "points")
                    _try_set_h5_attr(ds_x, "description", "X axis (index)")
                    _try_set_h5_attr(ds_x, "axis_type", "index")
                    axis_datasets["x"] = ds_x
                if y.ndim >= 2:
                    ny = y.shape[-2]
                    ds_y_ax = data_grp.create_dataset("abscissa_y", data=np.arange(ny))
                    _try_set_h5_attr(ds_y_ax, "units", "points")
                    _try_set_h5_attr(ds_y_ax, "description", "Y axis (index)")
                    _try_set_h5_attr(ds_y_ax, "axis_type", "index")
                    axis_datasets["y"] = ds_y_ax

            elif isinstance(x, np.ndarray):  # 1D data
                ds_x = data_grp.create_dataset("abscissa_x", data=x)
                _try_set_h5_attr(ds_x, "units", x_unit_val)
                _try_set_h5_attr(ds_x, "description", f"X axis")
                _try_set_h5_attr(ds_x, "axis_type", "independent_variable")
                axis_datasets["x"] = ds_x

            elif isinstance(x, list):  # Multi-D data
                if len(x) >= 1 and x[0] is not None and isinstance(x[0], np.ndarray):
                    ds_x = data_grp.create_dataset("abscissa_x", data=x[0])
                    _try_set_h5_attr(ds_x, "units", x_unit_val)
                    _try_set_h5_attr(ds_x, "description", "X axis")
                    _try_set_h5_attr(ds_x, "axis_type", "independent_variable_x")
                    axis_datasets["x"] = ds_x
                if len(x) >= 2 and x[1] is not None and isinstance(x[1], np.ndarray):
                    ds_y_ax = data_grp.create_dataset("abscissa_y", data=x[1])
                    _try_set_h5_attr(ds_y_ax, "units", y_unit_val)
                    _try_set_h5_attr(ds_y_ax, "description", "Y axis")
                    _try_set_h5_attr(ds_y_ax, "axis_type", "independent_variable_y")
                    axis_datasets["y"] = ds_y_ax
                if len(x) >= 3 and x[2] is not None and isinstance(x[2], np.ndarray):
                    ds_z_ax = data_grp.create_dataset("abscissa_z", data=x[2])
                    _try_set_h5_attr(ds_z_ax, "units", z_unit_val)
                    _try_set_h5_attr(ds_z_ax, "description", "Z axis")
                    _try_set_h5_attr(ds_z_ax, "axis_type", "independent_variable_z")
                    axis_datasets["z"] = ds_z_ax

            # Link axes to data dimensions using HDF5 Dimension Scales API
            if "intensity" in data_grp:
                dims = ds_y.dims
                current_ndim = ds_y.ndim

                # Link X dimension (last dimension)
                if current_ndim >= 1 and "x" in axis_datasets:
                    x_dim_index = current_ndim - 1
                    try:
                        dims[x_dim_index].label = "x"
                        dims[x_dim_index].attach_scale(axis_datasets["x"])
                    except Exception as e:
                        warnings.warn(
                            f"Error linking X dimension scale: {type(e).__name__} - {e}"
                        )

                # Link Y dimension (second to last dimension)
                if current_ndim >= 2 and "y" in axis_datasets:
                    y_dim_index = current_ndim - 2
                    try:
                        dims[y_dim_index].label = "y"
                        dims[y_dim_index].attach_scale(axis_datasets["y"])
                    except Exception as e:
                        warnings.warn(
                            f"Error linking Y dimension scale: {type(e).__name__} - {e}"
                        )

                # Link Z dimension (third to last dimension)
                if current_ndim >= 3 and "z" in axis_datasets:
                    z_dim_index = current_ndim - 3
                    try:
                        dims[z_dim_index].label = "z"
                        dims[z_dim_index].attach_scale(axis_datasets["z"])
                    except Exception as e:
                        warnings.warn(
                            f"Error linking Z dimension scale: {type(e).__name__} - {e}"
                        )

    except IOError as e:
        warnings.warn(f"Could not write HDF5 file {h5_file}: {e}")
    except Exception as e:
        warnings.warn(
            f"An unexpected error occurred while writing HDF5 file {h5_file}: "
            f"{type(e).__name__} - {e}"
        )


def save_fair(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
    formats: List[str] = ["csv_json", "hdf5"],
) -> None:
    """Save EPR data in specified FAIR formats.

    Args:
        output_basename: Base path for output files (without extension)
        x: Abscissa data array(s) or None
        y: Intensity data array
        pars: Raw parameters dictionary
        original_file_path: Path to original data file
        formats: List of output formats ('csv_json', 'hdf5')
    """
    if "csv_json" in formats:
        save_to_csv_json(output_basename, x, y, pars, original_file_path)

    if "hdf5" in formats:
        save_to_hdf5(output_basename, x, y, pars, original_file_path)
