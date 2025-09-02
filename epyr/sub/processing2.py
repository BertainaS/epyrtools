#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  6 15:57:27 2025

@author: sylvainbertaina
"""

import numpy as np
import warnings
# from scipy.signal import savgol_filter # If needed for other methods
from scipy.optimize import curve_fit # np.polyfit/polyval don't need this directly


# --- Baseline Correction Functions ---

def baseline_polynomial(
    y_data: np.ndarray,
    x_data: np.ndarray = None,
    poly_order: int = 1,
    exclude_regions: list = None,
    roi: tuple = None
) -> tuple[np.ndarray, np.ndarray]:
    """
    Performs baseline correction on 1D data by subtracting a polynomial fit.

    The polynomial is fitted to regions of the data presumed to be baseline.
    These regions can be specified by either defining an overall region of
    interest (ROI) for fitting, or by excluding signal regions from a
    broader fit.

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        x_data (np.ndarray, optional): The x-axis data corresponding to y_data.
            If provided, `exclude_regions` and `roi` are interpreted in
            the units of `x_data`. If None, indices are used. Defaults to None.
        poly_order (int, optional): Order of the polynomial to fit to the baseline.
            Defaults to 1 (linear).
        exclude_regions (list of tuples, optional): A list of (start, end) tuples
            specifying regions to *exclude* from the baseline fitting.
            Interpreted in `x_data` units if `x_data` is provided, otherwise as indices.
            Example: `[(x1, x2), (x3, x4)]`. Defaults to None (no regions explicitly excluded).
        roi (tuple, optional): A single (start, end) tuple specifying the region of
            interest (ROI) where the baseline fitting should *occur*. Data outside
            this ROI (and within any `exclude_regions` that might overlap) will not
            be used for the fit. Interpreted in `x_data` units if `x_data` is provided,
            otherwise as indices. Example: `(x_start_baseline, x_end_baseline)`.
            If None, the entire dataset (minus `exclude_regions`) is considered.
            Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            - y_corrected (np.ndarray): The baseline-corrected y_data.
            - baseline (np.ndarray): The calculated polynomial baseline.

    Raises:
        ValueError: If y_data is not 1D or if x_data and y_data lengths mismatch.
    """
    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")

    if x_data is not None:
        if not isinstance(x_data, np.ndarray) or x_data.ndim != 1:
            raise ValueError("x_data must be a 1D NumPy array if provided.")
        if len(x_data) != len(y_data):
            raise ValueError("x_data and y_data must have the same length.")
        x_axis_for_fitting = np.copy(x_data)
    else:
        x_axis_for_fitting = np.arange(len(y_data)) # Use indices if x_data is None

    y_to_fit = np.copy(y_data)

    # Create a mask for points to include in the baseline fit
    fit_mask = np.ones_like(y_to_fit, dtype=bool)

    if roi is not None:
        if len(roi) != 2:
            raise ValueError("roi must be a tuple of (start, end).")
        roi_start, roi_end = roi
        if roi_start >= roi_end:
            warnings.warn(f"ROI start ({roi_start}) is not less than ROI end ({roi_end}). ROI ignored.")
        else:
            # Initially, mask everything outside this ROI
            fit_mask[:] = False # Mask all
            # Unmask the ROI
            if x_data is not None: # ROI is in x_data units
                roi_active_mask = (x_axis_for_fitting >= roi_start) & (x_axis_for_fitting <= roi_end)
            else: # ROI is in index units
                roi_active_mask = (np.arange(len(y_to_fit)) >= roi_start) & \
                                  (np.arange(len(y_to_fit)) <= roi_end)
            fit_mask[roi_active_mask] = True

    if exclude_regions is not None:
        if not isinstance(exclude_regions, list):
            raise ValueError("exclude_regions must be a list of (start, end) tuples.")
        for region in exclude_regions:
            if not (isinstance(region, tuple) and len(region) == 2):
                raise ValueError("Each item in exclude_regions must be a (start, end) tuple.")
            start, end = region
            if start >= end:
                warnings.warn(f"Exclusion region start ({start}) is not less than end ({end}). Region ignored.")
                continue

            if x_data is not None: # Exclusion is in x_data units
                current_exclusion_mask = (x_axis_for_fitting >= start) & (x_axis_for_fitting <= end)
            else: # Exclusion is in index units
                current_exclusion_mask = (np.arange(len(y_to_fit)) >= start) & \
                                         (np.arange(len(y_to_fit)) <= end)
            fit_mask[current_exclusion_mask] = False # Exclude these points

    # Ensure there are enough points left for fitting
    num_points_for_fit = np.sum(fit_mask)
    if num_points_for_fit <= poly_order:
        warnings.warn(
            f"Not enough points ({num_points_for_fit}) for polynomial fit of order {poly_order} "
            "after applying ROI/exclusions. Returning original data."
        )
        return y_data, np.zeros_like(y_data)

    try:
        # Fit polynomial using only the masked (True) points
        coeffs = np.polyfit(x_axis_for_fitting[fit_mask], y_to_fit[fit_mask], poly_order)

        # Evaluate the polynomial over the entire original x_range (or indices)
        baseline = np.polyval(coeffs, x_axis_for_fitting)

        y_corrected = y_data - baseline
        return y_corrected, baseline
    except (np.linalg.LinAlgError, ValueError) as e:
        warnings.warn(f"Polynomial baseline fit failed: {e}. Returning original data.")
        return y_data, np.zeros_like(y_data)


def baseline_constant_offset(
    y_data: np.ndarray,
    offset_region_indices: tuple = None,
    method: str = 'mean'
) -> tuple[np.ndarray, np.ndarray]:
    """
    Performs baseline correction by subtracting a constant offset.

    The offset is calculated from a specified region of y_data (using indices)
    via its mean or median.

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        offset_region_indices (tuple, optional): A (start_index, end_index) tuple
            specifying the slice of y_data (exclusive of end_index, like Python slicing)
            to calculate the offset from. If None, the entire array is used,
            which is generally not recommended for data with signals.
            Defaults to None.
        method (str, optional): Method to calculate the offset within the region.
            Options are 'mean' or 'median'. Defaults to 'mean'.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            - y_corrected (np.ndarray): The baseline-corrected y_data.
            - baseline (np.ndarray): The calculated constant baseline (an array
                                     filled with the calculated offset value).

    Raises:
        ValueError: If y_data is not 1D or if the method is invalid.
    """
    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")

    if method.lower() not in ['mean', 'median']:
        raise ValueError("Method must be 'mean' or 'median'.")

    if offset_region_indices:
        if not (isinstance(offset_region_indices, tuple) and len(offset_region_indices) == 2):
            raise ValueError("offset_region_indices must be a (start_index, end_index) tuple.")
        start_idx, end_idx = offset_region_indices
        # Ensure indices are valid and form a valid slice
        start_idx = max(0, int(start_idx))
        end_idx = min(len(y_data), int(end_idx))

        if start_idx >= end_idx:
            warnings.warn(
                f"Invalid offset_region_indices ({offset_region_indices}). "
                "Start index must be less than end index. Using whole array for offset."
            )
            region_data_for_offset = y_data
        else:
            region_data_for_offset = y_data[start_idx:end_idx]
    else:
        warnings.warn("offset_region_indices not provided. Using whole array for offset calculation.")
        region_data_for_offset = y_data

    if region_data_for_offset.size == 0:
        warnings.warn("Offset calculation region is empty after slicing. Returning original data.")
        return y_data, np.zeros_like(y_data)

    if method.lower() == 'mean':
        offset_value = np.mean(region_data_for_offset)
    else: # 'median'
        offset_value = np.median(region_data_for_offset)

    baseline = np.full_like(y_data, offset_value)
    y_corrected = y_data - offset_value # Broadcasting subtracts scalar from array
    return y_corrected, baseline

# --- Stretched Exponential Decay Function ---
def stretched_exponential_decay(t, y0, amplitude, tau, beta):
    """
    Stretched exponential decay function (Kohlrausch–Williams–Watts type).
    f(t) = y0 + amplitude * exp(-(t / tau)**beta)
    """
    # Ensure t/tau is non-negative before raising to power beta, esp. if t can be negative.
    # However, for typical decays, t is expected to be non-negative.
    # We might need to handle negative t if the x-axis allows it and adjust the model.
    # For simplicity, assume t >= 0 and tau > 0.
    if tau <= 0: # tau must be positive
        return np.inf # Return a large number to penalize invalid tau in fitting
    if beta <= 0 or beta > 1.0001 : # Beta typically 0 < beta <= 1 (allow slight >1 for fitting robustness)
        return np.inf
    
    # Handle potential issues with (t/tau) being negative if t or x_offset not handled properly.
    # If t is always positive and starts from 0 (or t-t0), this is usually fine.
    decay_term = (t / tau)**beta
    return y0 + amplitude * np.exp(-decay_term)


def baseline_stretched_exponential(
    y_data: np.ndarray,
    x_data: np.ndarray,
    initial_guess: list = None,
    bounds: tuple = None,
    fit_region: tuple = None,
    exclude_regions: list = None
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Performs baseline correction by fitting and subtracting a stretched exponential decay.

    The function y(x) = y0 + A * exp(-((x - x_offset) / τ)**β) is fitted.
    This version assumes x_data is provided and that the decay starts effectively
    at or after x_data.min().

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        x_data (np.ndarray): The x-axis data corresponding to y_data. Must be provided
                             and should represent the variable against which the decay occurs
                             (e.g., time, field).
        initial_guess (list, optional): Initial guess for the parameters [y0, A, τ, β].
            - y0: baseline offset at infinity
            - A: amplitude of decay (can be negative if it's an increasing baseline)
            - τ: characteristic time
            - β: stretching exponent (0 < β <= 1)
            If None, some heuristics will be used to estimate, but providing
            good initial guesses is highly recommended for robust fitting.
        bounds (tuple, optional): Bounds for the parameters ([y0_min, A_min, τ_min, β_min],
                                 [y0_max, A_max, τ_max, β_max]).
            Example: `bounds=([0, -np.inf, 1e-9, 0.01], [np.inf, np.inf, np.inf, 1.0])`
            Defaults to `(-np.inf, np.inf)`.
        fit_region (tuple, optional): A (start_x, end_x) tuple specifying the
            region of x_data to use for fitting the decay. If None, the entire
            dataset (respecting exclude_regions) is used. Defaults to None.
        exclude_regions (list of tuples, optional): List of (start_x, end_x) tuples
            to exclude from fitting (e.g., signal peaks on top of the decay).
            Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray, dict]:
            - y_corrected (np.ndarray): Baseline-corrected y_data.
            - baseline (np.ndarray): The calculated stretched exponential baseline.
            - fit_params (dict): Dictionary of the fitted parameters (y0, A, τ, β)
                                 and their standard deviations (if calculable).
                                 Returns None for params if fit fails.
    Raises:
        ValueError: If y_data is not 1D or x_data is not provided or lengths mismatch.
    """
    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")
    if not isinstance(x_data, np.ndarray) or x_data.ndim != 1:
        raise ValueError("x_data must be a 1D NumPy array for stretched exponential fit.")
    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have the same length.")

    x_fit = np.copy(x_data)
    y_fit_data = np.copy(y_data)

    # Create a mask for points to include in the fit
    mask = np.ones_like(x_fit, dtype=bool)

    if fit_region:
        if len(fit_region) != 2: raise ValueError("fit_region must be a tuple (start_x, end_x).")
        start_x_fit, end_x_fit = fit_region
        if start_x_fit >= end_x_fit:
            warnings.warn(f"Fit region start ({start_x_fit}) not less than end ({end_x_fit}). Region ignored.")
        else:
            mask &= (x_fit >= start_x_fit) & (x_fit <= end_x_fit)

    if exclude_regions:
        if not isinstance(exclude_regions, list):
            raise ValueError("exclude_regions must be a list of (start_x, end_x) tuples.")
        for region in exclude_regions:
            if not (isinstance(region, tuple) and len(region) == 2):
                raise ValueError("Each item in exclude_regions must be a (start_x, end_x) tuple.")
            start_ex, end_ex = region
            if start_ex >= end_ex:
                warnings.warn(f"Exclusion region start ({start_ex}) not less than end ({end_ex}). Region ignored.")
                continue
            mask &= ~((x_fit >= start_ex) & (x_fit <= end_ex)) # Exclude these points

    # Apply mask
    x_masked = x_fit[mask]
    y_masked = y_fit_data[mask]

    if len(x_masked) < 4: # Need at least as many points as parameters
        warnings.warn("Not enough points to fit stretched exponential after applying regions. Returning original data.")
        return y_data, np.zeros_like(y_data), None

    # Heuristic initial guesses if not provided
    if initial_guess is None:
        # y0: Guess as the mean of the last few points of the *masked* data (or overall mean if fit_region is small)
        if len(y_masked) > 10:
            guess_y0 = np.mean(y_masked[-min(10, len(y_masked)//10):])
        else:
            guess_y0 = np.mean(y_masked)

        # A: Difference between first point (or max) and guess_y0
        # This assumes a decaying trend from the start of the masked data.
        first_masked_y = y_masked[0]
        guess_A = first_masked_y - guess_y0 # If decaying from high to low
        # If data might increase initially and then decay, or always increase (A < 0)
        # max_y = np.max(y_masked)
        # min_y = np.min(y_masked)
        # if abs(max_y - guess_y0) > abs(min_y - guess_y0):
        #     guess_A = max_y - guess_y0
        # else:
        #     guess_A = min_y - guess_y0

        # τ: Guess as roughly 1/3 of the range of x_masked where the decay happens.
        # This is a very rough guess and highly dependent on the data.
        if len(x_masked) > 1:
            guess_tau = (x_masked[-1] - x_masked[0]) / 3.0
            if guess_tau <= 0: guess_tau = (x_masked[-1] - x_masked[0]) if (x_masked[-1] - x_masked[0]) > 0 else 1.0
        else:
            guess_tau = 1.0
        guess_tau = max(guess_tau, 1e-9) # Ensure tau is positive

        # β: Start with 1 (simple exponential) or 0.5
        guess_beta = 0.8
        initial_guess = [guess_y0, guess_A, guess_tau, guess_beta]
        warnings.warn(f"Using heuristic initial guess for stretched exponential: {initial_guess}")

    # Default bounds if not provided
    if bounds is None:
        # y0, A, tau, beta
        # Tau must be > 0. Beta usually 0 < beta <= 1.
        lower_bounds = [-np.inf, -np.inf, 1e-12, 1e-3]  # Smallest tau and beta
        upper_bounds = [np.inf, np.inf, np.inf, 1.0] # Max beta = 1
        bounds = (lower_bounds, upper_bounds)
    else:
        # Ensure bounds are correctly structured
        if not (isinstance(bounds, tuple) and len(bounds) == 2 and
                isinstance(bounds[0], (list, np.ndarray)) and len(bounds[0]) == 4 and
                isinstance(bounds[1], (list, np.ndarray)) and len(bounds[1]) == 4):
            raise ValueError("Bounds must be a tuple of two lists/arrays of length 4: ([lowers], [uppers])")


    try:
        # --- Perform the fit on the masked data ---
        # It's crucial that x_masked corresponds to the part of the decay we want to model.
        # If the x_data doesn't start near the "beginning" of the decay (t=0 in model),
        # we might need to introduce an x_offset parameter in `stretched_exponential_decay`
        # and fit for it: e.g., ((x - x_offset) / tau).
        # For now, assuming x_masked[0] is close enough to the effective start.
        
        popt, pcov = curve_fit(
            stretched_exponential_decay,
            x_masked,  # Fit using the selected region's x values
            y_masked,
            p0=initial_guess,
            bounds=bounds,
            maxfev=5000 # Increase max iterations if convergence is an issue
        )

        # Calculate the baseline over the *entire original* x_data range
        baseline = stretched_exponential_decay(x_data, *popt)
        y_corrected = y_data - baseline

        # Extract parameter errors (standard deviations)
        perr = np.sqrt(np.diag(pcov)) if pcov is not None and not np.any(np.isinf(pcov)) else [np.nan]*len(popt)

        fit_parameters = {
            'y0': popt[0], 'A': popt[1], 'tau': popt[2], 'beta': popt[3],
            'std_y0': perr[0], 'std_A': perr[1], 'std_tau': perr[2], 'std_beta': perr[3]
        }
        return y_corrected, baseline, fit_parameters

    except RuntimeError as e:
        warnings.warn(f"Stretched exponential fit did not converge: {e}. Returning original data.")
        return y_data, np.zeros_like(y_data), None
    except ValueError as e: # Can happen from bounds or initial guess issues
        warnings.warn(f"ValueError during stretched exponential fit (check bounds/guesses): {e}. Returning original data.")
        return y_data, np.zeros_like(y_data), None
    except Exception as e: # Catch any other fitting errors
        warnings.warn(f"An unexpected error occurred during stretched exponential fit: {e}. Returning original data.")
        return y_data, np.zeros_like(y_data), None
    
    
# --- Mono-Exponential Decay Function ---
def mono_exponential_decay(t, y0, amplitude, tau):
    """
    Mono-exponential decay function.
    f(t) = y0 + amplitude * exp(-t / tau)
    """
    if tau <= 0: # tau must be positive
        return np.inf # Penalize invalid tau during fitting
    return y0 + amplitude * np.exp(-t / tau)


def baseline_mono_exponential(
    y_data: np.ndarray,
    x_data: np.ndarray,
    initial_guess: list = None,
    bounds: tuple = None,
    fit_region: tuple = None,
    exclude_regions: list = None
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Performs baseline correction by fitting and subtracting a mono-exponential decay.

    The function y(x) = y0 + A * exp(-(x - x_offset) / τ) is fitted.
    This version assumes x_data is provided and that the decay starts effectively
    at or after x_data.min().

    Args:
        y_data (np.ndarray): The 1D spectral data array.
        x_data (np.ndarray): The x-axis data corresponding to y_data. Must be provided.
        initial_guess (list, optional): Initial guess for parameters [y0, A, τ].
            - y0: baseline offset at infinity
            - A: amplitude of decay
            - τ: characteristic time
            If None, heuristics are used (providing good guesses is recommended).
        bounds (tuple, optional): Bounds for parameters ([y0_min, A_min, τ_min],
                                 [y0_max, A_max, τ_max]).
            Example: `bounds=([0, -np.inf, 1e-9], [np.inf, np.inf, np.inf])`
            Defaults to `(-np.inf, np.inf)`.
        fit_region (tuple, optional): (start_x, end_x) tuple specifying the x_data
            region for fitting. If None, the entire dataset (respecting
            exclude_regions) is used. Defaults to None.
        exclude_regions (list of tuples, optional): List of (start_x, end_x) tuples
            to exclude from fitting (e.g., signal peaks). Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray, dict]:
            - y_corrected (np.ndarray): Baseline-corrected y_data.
            - baseline (np.ndarray): The calculated mono-exponential baseline.
            - fit_params (dict): Dictionary of fitted parameters (y0, A, τ) and
                                 their standard deviations. None if fit fails.
    Raises:
        ValueError: If inputs are invalid.
    """
    if not isinstance(y_data, np.ndarray) or y_data.ndim != 1:
        raise ValueError("y_data must be a 1D NumPy array.")
    if not isinstance(x_data, np.ndarray) or x_data.ndim != 1:
        raise ValueError("x_data must be a 1D NumPy array for mono-exponential fit.")
    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have the same length.")

    x_fit = np.copy(x_data)
    y_fit_data = np.copy(y_data)

    # Create a mask for points to include in the fit (similar to stretched_exponential)
    mask = np.ones_like(x_fit, dtype=bool)
    if fit_region:
        if len(fit_region) != 2: raise ValueError("fit_region must be a tuple (start_x, end_x).")
        start_x_fit, end_x_fit = fit_region
        if start_x_fit >= end_x_fit: warnings.warn(f"Fit region start not less than end. Region ignored.")
        else: mask &= (x_fit >= start_x_fit) & (x_fit <= end_x_fit)

    if exclude_regions:
        if not isinstance(exclude_regions, list): raise ValueError("exclude_regions must be a list.")
        for region in exclude_regions:
            if not (isinstance(region, tuple) and len(region)==2): raise ValueError("Each exclude region must be (start,end).")
            start_ex, end_ex = region
            if start_ex >= end_ex: warnings.warn(f"Exclusion region start not less than end. Region ignored.")
            else: mask &= ~((x_fit >= start_ex) & (x_fit <= end_ex))

    x_masked = x_fit[mask]
    y_masked = y_fit_data[mask]

    if len(x_masked) < 3: # Need at least as many points as parameters (y0, A, τ)
        warnings.warn("Not enough points for mono-exponential fit after applying regions. Returning original data.")
        return y_data, np.zeros_like(y_data), None

    # Heuristic initial guesses if not provided
    if initial_guess is None:
        if len(y_masked) > 10:
            guess_y0 = np.mean(y_masked[-min(10, len(y_masked)//10):]) # Offset at long times
        else:
            guess_y0 = np.mean(y_masked)

        first_masked_y = y_masked[0]
        guess_A = first_masked_y - guess_y0 # Amplitude (assumes decay from start)

        if len(x_masked) > 1:
            guess_tau = (x_masked[-1] - x_masked[0]) / 3.0 # Characteristic time (very rough)
            if guess_tau <= 0 : guess_tau = (x_masked[-1] - x_masked[0]) if (x_masked[-1] - x_masked[0]) > 0 else 1.0
        else:
            guess_tau = 1.0
        guess_tau = max(guess_tau, 1e-9) # Ensure tau is positive

        initial_guess = [guess_y0, guess_A, guess_tau]
        warnings.warn(f"Using heuristic initial guess for mono-exponential: {initial_guess}")

    if bounds is None:
        # y0, A, tau
        lower_bounds = [-np.inf, -np.inf, 1e-12]  # Smallest tau
        upper_bounds = [np.inf, np.inf, np.inf]
        bounds = (lower_bounds, upper_bounds)
    else:
        if not (isinstance(bounds, tuple) and len(bounds) == 2 and
                isinstance(bounds[0], (list, np.ndarray)) and len(bounds[0]) == 3 and
                isinstance(bounds[1], (list, np.ndarray)) and len(bounds[1]) == 3):
            raise ValueError("Bounds must be a tuple of two lists/arrays of length 3 for mono-exponential.")


    try:
        popt, pcov = curve_fit(
            mono_exponential_decay,
            x_masked,
            y_masked,
            p0=initial_guess,
            bounds=bounds,
            maxfev=5000
        )

        baseline = mono_exponential_decay(x_data, *popt) # Calculate over full original x_data
        y_corrected = y_data - baseline

        perr = np.sqrt(np.diag(pcov)) if pcov is not None and not np.any(np.isinf(pcov)) else [np.nan]*len(popt)

        fit_parameters = {
            'y0': popt[0], 'A': popt[1], 'tau': popt[2],
            'std_y0': perr[0], 'std_A': perr[1], 'std_tau': perr[2]
        }
        return y_corrected, baseline, fit_parameters

    except RuntimeError as e:
        warnings.warn(f"Mono-exponential fit did not converge: {e}. Returning original data.")
        return y_data, np.zeros_like(y_data), None
    except ValueError as e:
        warnings.warn(f"ValueError during mono-exponential fit (check bounds/guesses): {e}. Returning original data.")
        return y_data, np.zeros_like(y_data), None
    except Exception as e:
        warnings.warn(f"An unexpected error occurred during mono-exponential fit: {e}. Returning original data.")
        return y_data, np.zeros_like(y_data), None