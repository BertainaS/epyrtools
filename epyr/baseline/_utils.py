# epyrtools/baseline/_utils.py

"""
Shared Utility Functions for Baseline Correction
==============================================

This module contains helper functions that are shared across different
baseline correction modules (e.g., _1d.py and _2d.py) within the
`epyrtools.baseline` sub-package.

Its main purpose is to handle the creation of boolean masks from
user-defined regions of interest and exclusion zones, for both 1D and 2D data.
"""

import warnings
from typing import List, Optional, Tuple

import numpy as np

# --- 1D Masking Utility ---


def _create_fit_mask(
    reference_axis: np.ndarray, 
    fit_window: Optional[Tuple[float, float]] = None, 
    exclude_regions: Optional[List[Tuple[float, float]]] = None
) -> np.ndarray:
    """Creates a boolean mask for data points to be used in 1D fitting.

    This helper function constructs a 1D boolean array where `True` values
    indicate which data points should be included in a fitting procedure.

    Args:
        reference_axis (np.ndarray): The axis (e.g., x-values or indices)
            against which `fit_window` and `exclude_regions` are defined.
        fit_window (tuple, optional): A `(start, end)` tuple specifying the
            primary region for fitting. Points outside this window are
            excluded. If None, the full range of `reference_axis` is
            initially considered. Defaults to None.
        exclude_regions (list of tuples, optional): A list of `(start, end)`
            tuples specifying regions to *exclude* from the fitting mask.
            These are applied after the `fit_window`. Defaults to None.

    Returns:
        np.ndarray:
            A boolean mask with the same size as `reference_axis`, where `True`
            indicates a point to include in the fit.
    """
    n_points = len(reference_axis)
    mask = np.ones(n_points, dtype=bool)

    # Apply the primary fit window (ROI) first
    if fit_window is not None:
        start_fw, end_fw = fit_window
        if start_fw >= end_fw:
            warnings.warn(
                f"Fit window start ({start_fw}) is not less than its end "
                f"({end_fw}). The window will be ignored, and the full range "
                f"considered (before exclusions).",
                UserWarning,
            )
        else:
            mask &= (reference_axis >= start_fw) & (reference_axis <= end_fw)

    # Apply exclusions on top of the current mask
    if exclude_regions is not None:
        for region_idx, region in enumerate(exclude_regions):
            start_ex, end_ex = region
            if start_ex >= end_ex:
                warnings.warn(
                    f"Exclusion region #{region_idx} start ({start_ex}) is not "
                    f"less than its end ({end_ex}). This region will be ignored.",
                    UserWarning,
                )
                continue
            # The ~ operator inverts the boolean mask for the exclusion region
            mask &= ~((reference_axis >= start_ex) & (reference_axis <= end_ex))
    return mask


# --- 2D Masking Utilities ---


def _get_slice_from_region(
    axis_data: np.ndarray, start_val: float, end_val: float
) -> slice:
    """Converts a value range in `axis_data` units to an index slice.

    Args:
        axis_data (np.ndarray): A sorted 1D array of axis coordinates.
        start_val (float): The starting value of the region.
        end_val (float): The ending value of the region.

    Returns:
        slice: An index slice `slice(start_index, end_index)` corresponding
               to the value range.
    """
    if start_val > end_val:  # Allow reversed specification for convenience
        start_val, end_val = end_val, start_val

    # searchsorted finds where elements should be inserted to maintain order.
    # 'left' includes the start value, 'right' includes the end value.
    idx_start = np.searchsorted(axis_data, start_val, side="left")
    idx_end = np.searchsorted(axis_data, end_val, side="right")

    # Clip to ensure indices are within array bounds for slicing
    idx_start = max(0, idx_start)
    idx_end = min(len(axis_data), idx_end)

    return slice(idx_start, idx_end)


def _create_fit_mask_2d(
    shape: tuple[int, int],
    x_axis_coords: np.ndarray = None,
    y_axis_coords: np.ndarray = None,
    fit_window_roi: tuple = None,
    exclude_regions: list = None,
) -> np.ndarray:
    """Creates a 2D boolean mask for fitting based on ROI and exclusions.

    Args:
        shape (tuple[int, int]): The (rows, cols) shape of the 2D data.
        x_axis_coords (np.ndarray, optional): 1D array of x-coordinates for
            all columns. If None, column indices are used.
        y_axis_coords (np.ndarray, optional): 1D array of y-coordinates for
            all rows. If None, row indices are used.
        fit_window_roi (tuple, optional): A region `((x_start, x_end), (y_start, y_end))`
            where fitting *should* occur. If None, the entire area is considered.
        exclude_regions (list, optional): List of regions `[((x1s,x1e),(y1s,y1e)), ...]`
            to *exclude* from fitting.

    Returns:
        np.ndarray: A 2D boolean mask of the given shape. `True` where fitting
                    should occur.
    """
    rows, cols = shape
    _x_axis = x_axis_coords if x_axis_coords is not None else np.arange(cols)
    _y_axis = y_axis_coords if y_axis_coords is not None else np.arange(rows)

    # Initialize mask based on the fit window (ROI)
    if fit_window_roi:
        (x_roi_start, x_roi_end), (y_roi_start, y_roi_end) = fit_window_roi
        x_slice_roi = _get_slice_from_region(_x_axis, x_roi_start, x_roi_end)
        y_slice_roi = _get_slice_from_region(_y_axis, y_roi_start, y_roi_end)
        mask = np.zeros(shape, dtype=bool)
        mask[y_slice_roi, x_slice_roi] = True
    else:
        mask = np.ones(shape, dtype=bool)

    # Set excluded regions in the mask to False
    if exclude_regions:
        for region in exclude_regions:
            (x_excl_start, x_excl_end), (y_excl_start, y_excl_end) = region
            x_slice_excl = _get_slice_from_region(_x_axis, x_excl_start, x_excl_end)
            y_slice_excl = _get_slice_from_region(_y_axis, y_excl_start, y_excl_end)
            mask[y_slice_excl, x_slice_excl] = False

    return mask
