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
from typing import Any, Dict, List, Optional, Union

import numpy as np

try:
    import h5py

    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

from ..logging_config import get_logger

logger = get_logger(__name__)

from .data_processing import process_parameters


def save_to_json(
    output_basename: Path,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Write parameter metadata to ``<basename>.json``.

    The output contains the source-file path, the FAIR-normalized
    metadata, and any unmapped Bruker keys preserved verbatim.

    Parameters
    ----------
    output_basename : pathlib.Path
        Base path; ``.json`` is appended.
    pars : dict
        Raw Bruker parameters as returned by :func:`epyr.eprload`.
    original_file_path : str
        Source-file path kept as provenance.

    Returns
    -------
    None

    Examples
    --------
    >>> from pathlib import Path
    >>> from epyr import eprload
    >>> from epyr.fair import save_to_json
    >>> x, y, params, fp = eprload("examples/data/130406SB_CaWO4_Er_CW_5K_20.DSC")
    >>> save_to_json(Path("/tmp/demo"), params, fp)  # doctest: +SKIP
    """
    json_file = output_basename.with_suffix(".json")
    fair_meta, unmapped_meta = process_parameters(pars)

    logger.info(f"  Saving structured metadata to: {json_file}")

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
        warnings.warn(f"Could not write JSON file {json_file}: {e}", stacklevel=2)
    except TypeError as e:
        warnings.warn(
            f"Error serializing metadata to JSON for {json_file}: {e}. "
            f"Some parameters might not be saved correctly.",
            stacklevel=2,
        )


def save_to_csv(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Write EPR data and a metadata header to ``<basename>.csv``.

    The header carries microwave frequency, modulation amplitude, sample
    name, and the path to the source file, in commented (``#``) lines
    before the column data.

    Parameters
    ----------
    output_basename : pathlib.Path
        Base path; ``.csv`` is appended.
    x : np.ndarray, list of np.ndarray, or None
        Abscissa from :func:`epyr.eprload`. Lists (2D) are written as
        long-format rows.
    y : np.ndarray
        Signal array.
    pars : dict
        Raw Bruker parameters.
    original_file_path : str
        Source-file path kept as provenance in the CSV header.

    Returns
    -------
    None

    Examples
    --------
    >>> from pathlib import Path
    >>> from epyr import eprload
    >>> from epyr.fair import save_to_csv
    >>> x, y, params, fp = eprload("examples/data/130406SB_CaWO4_Er_CW_5K_20.DSC")
    >>> save_to_csv(Path("/tmp/demo"), x, y, params, fp)  # doctest: +SKIP
    """
    csv_file = output_basename.with_suffix(".csv")
    fair_meta, unmapped_meta = process_parameters(pars)

    logger.info(f"  Saving data to: {csv_file}")

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
                    (
                        f"{mwfq_info.get('value', 'N/A')} "
                        f"{mwfq_info.get('unit', '')}"
                    ).strip(),
                ]
            )
            writer.writerow(
                [
                    "# Field_Center/Start:",
                    (
                        f"{field_info.get('value', 'N/A')} "
                        f"{field_info.get('unit', '')}"
                    ).strip(),
                ]
            )
            writer.writerow(
                [
                    "# Field_Sweep/Increment:",
                    (
                        f"{sweep_info.get('value', 'N/A')} "
                        f"{sweep_info.get('unit', '')}"
                    ).strip(),
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
                            "Provided x-axis ignored for CSV "
                            "(shape mismatch or not ndarray)."
                            " Using index.",
                            stacklevel=2,
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
        warnings.warn(f"Could not write CSV file {csv_file}: {e}", stacklevel=2)
    except Exception as e:
        warnings.warn(
            f"An unexpected error occurred while writing CSV {csv_file}: {e}",
            stacklevel=2,
        )


def save_to_csv_json(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Save data to CSV and structured metadata to JSON.

    This function is maintained for backward compatibility and calls
    save_to_csv() and save_to_json() separately.

    Args:
        output_basename: Base path for output files (without extension)
        x: Abscissa data array(s) or None
        y: Intensity data array
        pars: Raw parameters dictionary
        original_file_path: Path to original data file
    """
    save_to_json(output_basename, pars, original_file_path)
    save_to_csv(output_basename, x, y, pars, original_file_path)


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
            f"directly in HDF5 attributes. Converting to string.",
            stacklevel=2,
        )
        h5_object.attrs[key] = str(value)
    except Exception as e:
        warnings.warn(
            f"Unexpected error storing attribute '{key}': "
            f"{type(e).__name__} - {e}. Skipping.",
            stacklevel=2,
        )


def save_to_hdf5(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Write data and metadata to ``<basename>.h5``.

    Datasets:
        ``/intensity``  : signal array (``y``)
        ``/abscissa``   : abscissa array (1D) or group ``/axis_0``,
                          ``/axis_1`` (2D)

    Metadata is written as attributes on the root group: the FAIR-mapped
    parameters, unmapped raw parameters, and the source-file path.

    Parameters
    ----------
    output_basename : pathlib.Path
        Base path; ``.h5`` is appended.
    x : np.ndarray, list of np.ndarray, or None
        Abscissa from :func:`epyr.eprload`.
    y : np.ndarray
        Signal array. Complex data is stored as two real datasets,
        ``/intensity_real`` and ``/intensity_imag``.
    pars : dict
        Raw Bruker parameters.
    original_file_path : str
        Source path kept as an HDF5 attribute.

    Returns
    -------
    None
        A ``UserWarning`` is emitted (no exception) if ``h5py`` is not
        installed; the file is then not written.

    Examples
    --------
    >>> from pathlib import Path
    >>> from epyr import eprload
    >>> from epyr.fair import save_to_hdf5
    >>> x, y, params, fp = eprload("examples/data/Rabi2D_GdCaWO4_13dB_3057G.DSC")
    >>> save_to_hdf5(Path("/tmp/demo"), x, y, params, fp)  # doctest: +SKIP
    """
    if not HAS_H5PY:
        warnings.warn(
            "h5py library not found. Skipping HDF5 output. "
            "Install with 'pip install h5py'",
            stacklevel=2,
        )
        return

    h5_file = output_basename.with_suffix(".h5")
    fair_meta, unmapped_meta = process_parameters(pars)

    logger.info(f"  Saving structured data and metadata to: {h5_file}")

    try:
        with h5py.File(h5_file, "w") as f:
            # Store global metadata
            f.attrs["original_file"] = original_file_path
            f.attrs["description"] = (
                "FAIR representation of EPR data converted from Bruker format."
            )
            f.attrs["conversion_timestamp"] = datetime.now().isoformat()
            f.attrs["converter_script_version"] = "epyr_fair_converter_v1.0"

            # Store structured FAIR metadata
            param_grp = f.create_group("metadata/parameters_fair")
            param_grp.attrs["description"] = (
                "Mapped parameters with units and descriptions."
            )

            for fair_key, info in fair_meta.items():
                item_grp = param_grp.create_group(fair_key)
                _try_set_h5_attr(item_grp, "value", info["value"])
                _try_set_h5_attr(item_grp, "unit", info["unit"])
                _try_set_h5_attr(item_grp, "description", info["description"])

            # Store unmapped parameters
            if unmapped_meta:
                unmap_grp = f.create_group("metadata/parameters_original")
                unmap_grp.attrs["description"] = (
                    "Parameters from the original file not found in the FAIR mapping."
                )
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
                _try_set_h5_attr(ds_x, "description", "X axis")
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
                            f"Error linking X dimension scale: "
                            f"{type(e).__name__} - {e}",
                            stacklevel=2,
                        )

                # Link Y dimension (second to last dimension)
                if current_ndim >= 2 and "y" in axis_datasets:
                    y_dim_index = current_ndim - 2
                    try:
                        dims[y_dim_index].label = "y"
                        dims[y_dim_index].attach_scale(axis_datasets["y"])
                    except Exception as e:
                        warnings.warn(
                            f"Error linking Y dimension scale: "
                            f"{type(e).__name__} - {e}",
                            stacklevel=2,
                        )

                # Link Z dimension (third to last dimension)
                if current_ndim >= 3 and "z" in axis_datasets:
                    z_dim_index = current_ndim - 3
                    try:
                        dims[z_dim_index].label = "z"
                        dims[z_dim_index].attach_scale(axis_datasets["z"])
                    except Exception as e:
                        warnings.warn(
                            f"Error linking Z dimension scale: "
                            f"{type(e).__name__} - {e}",
                            stacklevel=2,
                        )

    except IOError as e:
        warnings.warn(f"Could not write HDF5 file {h5_file}: {e}", stacklevel=2)
    except Exception as e:
        warnings.warn(
            f"An unexpected error occurred while writing HDF5 file {h5_file}: "
            f"{type(e).__name__} - {e}",
            stacklevel=2,
        )


def save_to_jpg(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
) -> None:
    """Write a preview figure to ``<basename>.jpg``.

    For 1D data, a single ``plot_1d`` figure. For 2D data, two files:
    ``<basename>_map.jpg`` and ``<basename>_waterfall.jpg``.

    Parameters
    ----------
    output_basename : pathlib.Path
        Base path; the JPG suffix is appended.
    x : np.ndarray, list of np.ndarray, or None
        Abscissa from :func:`epyr.eprload`.
    y : np.ndarray
        Signal array (1D or 2D).
    pars : dict
        Raw Bruker parameters, used for axis labels.
    original_file_path : str
        Source path, used as the figure title.

    Returns
    -------
    None
        Uses the ``Agg`` non-interactive backend; safe in scripts.

    Examples
    --------
    >>> from pathlib import Path
    >>> from epyr import eprload
    >>> from epyr.fair import save_to_jpg
    >>> x, y, params, fp = eprload("examples/data/130406SB_CaWO4_Er_CW_5K_20.DSC")
    >>> save_to_jpg(Path("/tmp/demo"), x, y, params, fp)  # doctest: +SKIP
    """
    try:
        import matplotlib

        matplotlib.use("Agg")  # Non-interactive backend
        import matplotlib.pyplot as plt

        from ..eprplot import plot_1d, plot_2d_map, plot_2d_waterfall
    except ImportError as e:
        warnings.warn(
            f"Could not import plotting modules: {e}. Skipping JPG export.",
            stacklevel=2,
        )
        return

    try:
        file_name = Path(original_file_path).name

        if y.ndim == 1:
            # 1D data - single plot
            logger.info(f"  Saving 1D plot to: {output_basename}.jpg")
            fig, ax = plot_1d(x, y, pars, title=file_name)
            plt.savefig(
                output_basename.with_suffix(".jpg"),
                dpi=200,
                format="jpg",
                bbox_inches="tight",
            )
            plt.close(fig)

        elif y.ndim == 2:
            # Symmetric color scale based on 98th percentile
            vlim = np.percentile(np.abs(np.real(y)), 98)

            # 2D data - both map and waterfall plots
            logger.info(f"  Saving 2D map plot to: {output_basename}_map.jpg")
            fig_map, ax_map = plot_2d_map(
                x,
                y,
                pars,
                title=f"{file_name} - Map",
                cmap="RdBu_r",
                vmin=-vlim,
                vmax=vlim,
            )
            plt.savefig(
                str(output_basename) + "_map.jpg",
                dpi=200,
                format="jpg",
                bbox_inches="tight",
            )
            plt.close(fig_map)

            logger.info(
                f"  Saving 2D waterfall plot to: {output_basename}_waterfall.jpg"
            )
            fig_waterfall, ax_waterfall = plot_2d_waterfall(
                x, y, pars, title=f"{file_name} - Waterfall"
            )
            plt.savefig(
                str(output_basename) + "_waterfall.jpg",
                dpi=200,
                format="jpg",
                bbox_inches="tight",
            )
            plt.close(fig_waterfall)

        else:
            warnings.warn(
                f"Cannot create JPG for {y.ndim}D data. Only 1D and 2D supported.",
                stacklevel=2,
            )

    except Exception as e:
        warnings.warn(
            f"Failed to create JPG for {original_file_path}: {e}", stacklevel=2
        )


def save_fair(
    output_basename: Path,
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    pars: Dict[str, Any],
    original_file_path: str,
    formats: Optional[List[str]] = None,
) -> None:
    """Save EPR data in specified FAIR formats.

    Args:
        output_basename: Base path for output files (without extension)
        x: Abscissa data array(s) or None
        y: Intensity data array
        pars: Raw parameters dictionary
        original_file_path: Path to original data file
        formats: List of output formats.
            Options: 'csv', 'json', 'hdf5', 'jpg', 'csv_json'
            - 'csv': Save data to CSV file only
            - 'json': Save metadata to JSON file only
            - 'hdf5': Save data and metadata to HDF5 file
            - 'jpg': Save visualization plots (1D: single plot, 2D: map + waterfall)
            - 'csv_json': Save both CSV and JSON (backward compatibility)
            Defaults to ['csv', 'json'].
    """
    if formats is None:
        formats = ["csv", "json"]
    # Handle individual formats
    if "csv" in formats:
        save_to_csv(output_basename, x, y, pars, original_file_path)

    if "json" in formats:
        save_to_json(output_basename, pars, original_file_path)

    # Handle backward compatibility format
    if "csv_json" in formats:
        save_to_csv_json(output_basename, x, y, pars, original_file_path)

    if "hdf5" in formats:
        save_to_hdf5(output_basename, x, y, pars, original_file_path)

    if "jpg" in formats:
        save_to_jpg(output_basename, x, y, pars, original_file_path)
