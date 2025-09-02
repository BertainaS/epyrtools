# epyrtools/baseline/_1d.py

"""
1D Baseline Correction Algorithms
=================================

This module contains the core functions for performing baseline correction on
1D data arrays. These functions are exposed to the user through the
`epyrtools.baseline` sub-package.

It includes polynomial, constant offset, and exponential decay models.
"""

import warnings

import numpy as np
from scipy.optimize import curve_fit

# Import shared utility functions from within the same sub-package
from ._utils import _create_fit_mask

# --- Model Functions for Curve Fitting ---


def _stretched_exponential_decay_model(t, y0, amplitude, tau, beta):
    """Model: y0 + amplitude * exp(-(t / tau)**beta)."""
    if tau <= 1e-12:  # Prevent division by zero or extremely small tau
        return np.inf
    if not (1e-3 < beta <= 5.0001):  # Beta is typically 0 < beta <= 3
        return np.inf
    decay_arg = t / tau
    # Ensure argument to exponentiation is non-negative
    decay_arg_safe = np.maximum(decay_arg, 0)
    return y0 + amplitude * np.exp(-((decay_arg_safe) ** beta))


def _mono_exponential_decay_model(t, y0, amplitude, tau):
    """Model: y0 + amplitude * exp(-t / tau)."""
    if tau <= 1e-12:  # Prevent division by zero or extremely small tau
        return np.inf
    return y0 + amplitude * np.exp(-t / tau)


# --- Heuristic Initial Guess Generators for Exponential Models ---


def _heuristic_stretched_exp_guess(x_masked, y_masked):
    """Generates a heuristic initial guess for a stretched exponential fit."""
    guess_y0 = np.mean(
        y_masked[-max(1, len(y_masked) // 20) :]
    )  # Offset is mean of tail
    guess_A = y_masked[0] - guess_y0 if len(y_masked) > 0 else 0  # Amplitude
    if len(x_masked) > 1 and x_masked[-1] > x_masked[0]:
        guess_tau = (x_masked[-1] - x_masked[0]) / 3.0  # Decay time is ~1/3 of range
    else:
        guess_tau = 1.0
    guess_beta = 0.8  # Common value for stretched exponential
    return [guess_y0, guess_A, max(guess_tau, 1e-9), guess_beta]


def _heuristic_mono_exp_guess(x_masked, y_masked):
    """Generates a heuristic initial guess for a mono-exponential fit."""
    guess_y0 = np.mean(y_masked[-max(1, len(y_masked) // 20) :])
    guess_A = y_masked[0] - guess_y0 if len(y_masked) > 0 else 0
    if len(x_masked) > 1 and x_masked[-1] > x_masked[0]:
        guess_tau = (x_masked[-1] - x_masked[0]) / 3.0
    else:
        guess_tau = 1.0
    return [guess_y0, guess_A, max(guess_tau, 1e-9)]


# --- Exponential Baseline Core Logic ---


def _fit_exponential_baseline(
    y_data: np.ndarray,
    x_data: np.ndarray,
    model_func,
    param_names: list,
    default_initial_guess_generator,
    default_bounds: tuple,
    user_initial_guess: list = None,
    user_bounds: tuple = None,
    fit_region: tuple = None,
    exclude_regions: list = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """A generic helper to fit an exponential-like baseline model."""
    n_params = len(param_names)

    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")
    if not isinstance(x_data, np.ndarray) or x_data.ndim != 1:
        raise ValueError(f"{model_func.__name__} requires x_data as a 1D NumPy array.")
    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have the same length.")

    # Validate region formats
    if fit_region is not None and not (
        isinstance(fit_region, tuple) and len(fit_region) == 2
    ):
        raise ValueError("fit_region must be a (start_x, end_x) tuple.")
    if exclude_regions is not None:
        if not isinstance(exclude_regions, list):
            raise ValueError(
                "exclude_regions must be a list of (start_x, end_x) tuples."
            )
        for i, region in enumerate(exclude_regions):
            if not (isinstance(region, tuple) and len(region) == 2):
                raise ValueError(
                    f"Each item in exclude_regions (item {i}) must be a (start_x, end_x) tuple."
                )

    fit_mask = _create_fit_mask(
        x_data, fit_window=fit_region, exclude_regions=exclude_regions
    )
    x_masked, y_masked = x_data[fit_mask], y_data[fit_mask]

    if len(x_masked) < n_params:
        warnings.warn(
            f"Not enough points ({len(x_masked)}) to fit {model_func.__name__} "
            f"with {n_params} parameters after applying regions. "
            "Returning original data.",
            UserWarning,
        )
        return y_data, np.zeros_like(y_data), None

    # Handle initial guess
    if user_initial_guess is None:
        current_initial_guess = default_initial_guess_generator(x_masked, y_masked)
        param_order_str = ", ".join(param_names)
        warnings.warn(
            f"Using heuristic initial guess for {model_func.__name__} "
            f"([{param_order_str}]): {current_initial_guess}",
            UserWarning,
        )
    else:
        if not (
            isinstance(user_initial_guess, (list, np.ndarray))
            and len(user_initial_guess) == n_params
        ):
            raise ValueError(
                f"initial_guess must be a list/array of {n_params} values."
            )
        current_initial_guess = user_initial_guess

    # Handle parameter bounds
    current_bounds = default_bounds
    if user_bounds is not None:
        if not (
            isinstance(user_bounds, tuple)
            and len(user_bounds) == 2
            and len(user_bounds[0]) == n_params
            and len(user_bounds[1]) == n_params
        ):
            raise ValueError(
                f"bounds must be a tuple of two lists/arrays of length {n_params}: "
                "([lowers], [uppers])"
            )
        current_bounds = user_bounds

    try:
        popt, pcov = curve_fit(
            model_func,
            x_masked,
            y_masked,
            p0=current_initial_guess,
            bounds=current_bounds,
            maxfev=10000,
            method="trf",
        )

        baseline = model_func(x_data, *popt)
        y_corrected = y_data - baseline

        # Calculate parameter standard errors from covariance matrix
        if (
            pcov is not None
            and not np.any(np.isinf(pcov))
            and np.all(np.diag(pcov) >= 0)
        ):
            perr = np.sqrt(np.diag(pcov))
        else:
            perr = [np.nan] * n_params
            if pcov is not None:
                warnings.warn(
                    "Covariance matrix from fit is problematic. Std errors set to NaN.",
                    UserWarning,
                )

        # Compile results into a dictionary
        fit_parameters = {param_names[i]: popt[i] for i in range(n_params)}
        fit_parameters.update(
            {f"std_{param_names[i]}": perr[i] for i in range(n_params)}
        )

        return y_corrected, baseline, fit_parameters

    except RuntimeError as e:
        warnings.warn(
            f"{model_func.__name__} fit did not converge: {e}. Returning original data.",
            UserWarning,
        )
    except ValueError as e:
        warnings.warn(
            f"ValueError during {model_func.__name__} fit: {e}. Returning original data.",
            UserWarning,
        )
    except Exception as e:
        warnings.warn(
            f"Unexpected error in {model_func.__name__} fit: {e}. Returning original data.",
            UserWarning,
        )

    return y_data, np.zeros_like(y_data), None


# --- Public 1D Baseline Correction Functions ---


def baseline_polynomial(
    y_data: np.ndarray,
    x_data: np.ndarray = None,
    poly_order: int = 1,
    exclude_regions: list = None,
    roi: tuple = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Performs baseline correction by subtracting a polynomial fit.

    The polynomial is fitted to specific regions of the data, which are
    presumed to contain only the baseline. These regions are defined by
    specifying an overall Region of Interest (ROI) and/or by excluding regions
    that contain signals (e.g., peaks).

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        x_data (np.ndarray, optional): The x-axis data corresponding to y_data.
            If provided, `roi` and `exclude_regions` are interpreted in the
            units of `x_data`. If None, array indices are used. Defaults to None.
        poly_order (int, optional): Order of the polynomial to fit.
            Defaults to 1 (a linear fit).
        exclude_regions (list of tuples, optional): A list of `(start, end)`
            tuples specifying regions to *exclude* from the baseline fit.
            Interpreted in `x_data` units or indices. Defaults to None.
        roi (tuple, optional): A single `(start, end)` tuple specifying the
            Region of Interest where baseline fitting should *occur*. Data outside
            this ROI (and within `exclude_regions`) will not be used. If None,
            the entire dataset (minus `exclude_regions`) is used.
            Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            A tuple containing:
            - **y_corrected** (np.ndarray): The baseline-corrected `y_data`.
            - **baseline** (np.ndarray): The calculated polynomial baseline.

    Raises:
        ValueError: If inputs are invalid.
    """
    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")
    if poly_order < 0:
        raise ValueError("poly_order must be a non-negative integer.")

    n_points = len(y_data)

    # Determine the reference axis for fitting and evaluation
    if x_data is not None:
        if not isinstance(x_data, np.ndarray) or x_data.ndim != 1:
            raise ValueError("x_data must be a 1D NumPy array if provided.")
        if len(x_data) != n_points:
            raise ValueError("x_data and y_data must have the same length.")
        reference_axis = x_data
    else:
        reference_axis = np.arange(n_points)

    # Validate region formats
    if roi is not None and not (isinstance(roi, tuple) and len(roi) == 2):
        raise ValueError("roi must be a tuple of (start, end).")
    if exclude_regions is not None:
        if not isinstance(exclude_regions, list):
            raise ValueError("exclude_regions must be a list of (start, end) tuples.")
        for i, region in enumerate(exclude_regions):
            if not (isinstance(region, tuple) and len(region) == 2):
                raise ValueError(
                    f"Each item in exclude_regions (item {i}) must be a (start, end) tuple."
                )

    # Create the mask for points to include in the fit using the shared utility
    fit_mask = _create_fit_mask(
        reference_axis, fit_window=roi, exclude_regions=exclude_regions
    )

    num_points_for_fit = np.sum(fit_mask)
    if num_points_for_fit <= poly_order:
        warnings.warn(
            f"Not enough points ({num_points_for_fit}) for polynomial fit of "
            f"order {poly_order} after applying ROI/exclusions. "
            "Returning original data.",
            UserWarning,
        )
        return y_data, np.zeros_like(y_data)

    try:
        # Fit polynomial using only the masked (True) points
        coeffs = np.polyfit(reference_axis[fit_mask], y_data[fit_mask], poly_order)
        # Evaluate the polynomial over the entire original axis
        baseline = np.polyval(coeffs, reference_axis)
        y_corrected = y_data - baseline
        return y_corrected, baseline
    except (np.linalg.LinAlgError, ValueError) as e:
        warnings.warn(
            f"Polynomial baseline fit failed: {e}. Returning original data.",
            UserWarning,
        )
        return y_data, np.zeros_like(y_data)


def baseline_constant_offset(
    y_data: np.ndarray, offset_region_indices: tuple = None, method: str = "mean"
) -> tuple[np.ndarray, np.ndarray]:
    """Performs baseline correction by subtracting a constant offset.

    The offset is calculated from a specified region of `y_data` (using
    array indices) via its mean or median.

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        offset_region_indices (tuple, optional): A `(start_index, end_index)`
            tuple specifying the slice of `y_data` (exclusive of `end_index`)
            to calculate the offset from. If None, the entire array is used.
            Defaults to None.
        method (str, optional): The method to calculate the offset. Can be
            'mean' or 'median'. Defaults to 'mean'.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            A tuple containing:
            - **y_corrected** (np.ndarray): The baseline-corrected `y_data`.
            - **baseline** (np.ndarray): The calculated constant baseline array.

    Raises:
        ValueError: If inputs are invalid.
    """
    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")
    if method.lower() not in ["mean", "median"]:
        raise ValueError("Method must be 'mean' or 'median'.")

    # Determine the slice of data to use for offset calculation
    if offset_region_indices:
        if not (
            isinstance(offset_region_indices, tuple) and len(offset_region_indices) == 2
        ):
            raise ValueError(
                "offset_region_indices must be a (start_index, end_index) tuple."
            )
        try:
            start_idx, end_idx = map(int, offset_region_indices)
        except (ValueError, TypeError):
            raise ValueError(
                "offset_region_indices components must be convertible to int."
            )

        start_idx = max(0, start_idx)
        end_idx = min(len(y_data), end_idx)

        if start_idx >= end_idx:
            warnings.warn(
                f"Invalid offset_region_indices ({offset_region_indices}). "
                "Start index must be less than end index. Using whole array.",
                UserWarning,
            )
            region_for_offset = y_data
        else:
            region_for_offset = y_data[start_idx:end_idx]
    else:
        warnings.warn(
            "offset_region_indices not provided. Using whole array for offset.",
            UserWarning,
        )
        region_for_offset = y_data

    if region_for_offset.size == 0:
        warnings.warn(
            "Offset calculation region is empty. Returning original data.", UserWarning
        )
        return y_data, np.zeros_like(y_data)

    # Calculate the offset value
    offset_value = (
        np.mean(region_for_offset)
        if method.lower() == "mean"
        else np.median(region_for_offset)
    )

    baseline = np.full_like(y_data, offset_value)
    y_corrected = y_data - offset_value
    return y_corrected, baseline


def baseline_stretched_exponential(
    y_data: np.ndarray,
    x_data: np.ndarray,
    initial_guess: list = None,
    bounds: tuple = None,
    fit_region: tuple = None,
    exclude_regions: list = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Fits and subtracts a stretched exponential decay baseline.

    The model fitted is: ``y(x) = y0 + A * exp(-((x / τ)**β))``

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        x_data (np.ndarray): The x-axis data corresponding to y_data.
        initial_guess (list, optional): A list of initial guess values for the
            parameters `[y0, A, tau, beta]`. If None, a heuristic guess is
            generated. Defaults to None.
        bounds (tuple, optional): A tuple `([lowers], [uppers])` for parameter
            bounds. Defaults to `([-inf, -inf, 1e-12, 1e-3], [inf, inf, inf, 5.0])`.
        fit_region (tuple, optional): A `(start, end)` tuple specifying the region
            in `x_data` units to use for fitting. Defaults to the whole range.
        exclude_regions (list of tuples, optional): A list of `(start, end)`
            tuples to *exclude* from the fit. Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray, dict | None]:
            A tuple containing the corrected data, the baseline, and a dictionary
            of fit parameters. Returns `None` for parameters if the fit fails.
    """
    param_names = ["y0", "A", "tau", "beta"]
    default_bounds = ([-np.inf, -np.inf, 1e-12, 1e-3], [np.inf, np.inf, np.inf, 5.0])
    return _fit_exponential_baseline(
        y_data,
        x_data,
        _stretched_exponential_decay_model,
        param_names,
        _heuristic_stretched_exp_guess,
        default_bounds,
        user_initial_guess=initial_guess,
        user_bounds=bounds,
        fit_region=fit_region,
        exclude_regions=exclude_regions,
    )


def baseline_mono_exponential(
    y_data: np.ndarray,
    x_data: np.ndarray,
    initial_guess: list = None,
    bounds: tuple = None,
    fit_region: tuple = None,
    exclude_regions: list = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Fits and subtracts a mono-exponential decay baseline.

    The model fitted is: ``y(x) = y0 + A * exp(-(x / τ))``

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        x_data (np.ndarray): The x-axis data corresponding to y_data.
        initial_guess (list, optional): A list of initial guess values for the
            parameters `[y0, A, tau]`. If None, a heuristic guess is generated.
            Defaults to None.
        bounds (tuple, optional): A tuple `([lowers], [uppers])` for parameter
            bounds. Defaults to `([-inf, -inf, 1e-12], [inf, inf, inf])`.
        fit_region (tuple, optional): A `(start, end)` tuple specifying the region
            in `x_data` units to use for fitting. Defaults to the whole range.
        exclude_regions (list of tuples, optional): A list of `(start, end)`
            tuples to *exclude* from the fit. Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray, dict | None]:
            A tuple containing the corrected data, the baseline, and a dictionary
            of fit parameters. Returns `None` for parameters if the fit fails.
    """
    param_names = ["y0", "A", "tau"]
    default_bounds = ([-np.inf, -np.inf, 1e-12], [np.inf, np.inf, np.inf])
    return _fit_exponential_baseline(
        y_data,
        x_data,
        _mono_exponential_decay_model,
        param_names,
        _heuristic_mono_exp_guess,
        default_bounds,
        user_initial_guess=initial_guess,
        user_bounds=bounds,
        fit_region=fit_region,
        exclude_regions=exclude_regions,
    )
