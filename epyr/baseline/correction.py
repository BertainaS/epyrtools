#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core baseline correction algorithms for EPR data.

This module contains the main baseline correction functions that work
directly with data from epyr.eprload().
"""

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from scipy.optimize import curve_fit

from ..logging_config import get_logger

logger = get_logger(__name__)

from .interactive import (
    RegionSelector,
    interactive_select_regions_1d,
    interactive_select_regions_2d,
    is_interactive_available,
)

# Import from our new modules
from .models import (
    bi_exponential_1d,
    exponential_1d,
    polynomial_1d,
    polynomial_2d,
    stretched_exponential_1d,
)
from .selection import (
    create_region_mask_1d,
    create_region_mask_2d,
    get_baseline_regions_1d,
    get_baseline_regions_2d,
)


def baseline_polynomial_1d(
    x: Union[np.ndarray, None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    order: int = 2,
    exclude_center: bool = True,
    center_fraction: float = 0.3,
    manual_regions: Optional[List[Tuple[float, float]]] = None,
    region_mode: str = "exclude",
    interactive: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Polynomial baseline correction for 1D EPR data.

    Fits a polynomial baseline to specified regions and subtracts it from
    the entire dataset. Ideal for CW EPR spectra with smooth baseline drifts.

    Args:
        x: X-axis data from eprload (can be None)
        y: 1D spectral data array from eprload
        params: Parameter dictionary from eprload (optional)
        order: Polynomial order (typically 1-4)
        exclude_center: If True, exclude center region from fitting
        center_fraction: Fraction of data to exclude from center
        manual_regions: List of manually selected regions as [(x1, x2), ...]
        region_mode: 'exclude' to exclude manual_regions, 'include' to use only manual_regions
        interactive: If True, open interactive region selector

    Returns:
        tuple: (corrected_data, baseline)
    """
    if y is None or y.ndim != 1:
        raise ValueError("y must be a 1D array")

    n_points = len(y)

    # Create x-coordinates
    if x is None or len(x) != n_points:
        x_coords = np.arange(n_points)
    else:
        x_coords = x

    # Interactive region selection if requested
    selected_regions = manual_regions if manual_regions is not None else []

    if interactive:
        if not is_interactive_available():
            logger.warning("⚠️  Interactive selection may not work in this environment.")
            logger.warning("   Consider using manual_regions parameter instead.")

        logger.info("🖱️ Interactive region selection enabled...")
        selected_regions = interactive_select_regions_1d(
            x_coords,
            y,
            f"Select regions to {region_mode.upper()} from baseline fitting",
        )
        logger.info(f"✅ Selected {len(selected_regions)} regions")

    # Create baseline mask
    mask = get_baseline_regions_1d(
        x_coords,
        y,
        exclude_center=exclude_center and not manual_regions,
        center_fraction=center_fraction,
        manual_regions=selected_regions,
        region_mode=region_mode,
    )

    # Fit polynomial to baseline regions
    x_fit = x_coords[mask]
    y_fit = y[mask]

    # Remove NaN values
    valid = np.isfinite(y_fit)
    x_fit = x_fit[valid]
    y_fit = y_fit[valid]

    if len(y_fit) < order + 1:
        warnings.warn(f"Not enough points ({len(y_fit)}) for polynomial order {order}")
        return y, np.zeros_like(y)

    try:
        # Fit polynomial
        coeffs = np.polyfit(x_fit, y_fit, order)

        # Evaluate baseline over full range
        baseline = np.polyval(coeffs, x_coords)

        # Subtract baseline
        corrected_data = y - baseline

        return corrected_data, baseline

    except np.linalg.LinAlgError as e:
        warnings.warn(f"Polynomial fitting failed: {e}")
        return y, np.zeros_like(y)


def baseline_polynomial_2d(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    order: Union[int, Tuple[int, int]] = 1,
    exclude_center: bool = True,
    center_fraction: float = 0.3,
    manual_regions: Optional[
        List[Tuple[Tuple[float, float], Tuple[float, float]]]
    ] = None,
    region_mode: str = "exclude",
    interactive: bool = False,
    use_real_part: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Polynomial baseline correction for 2D EPR data.

    Fits a 2D polynomial surface to specified regions and subtracts it from
    the entire dataset. Useful for 2D EPR measurements like DEER, Rabi oscillations, etc.

    Args:
        x: X-axis coordinates (None, 1D arrays, or meshgrids)
        y: 2D spectral data array
        params: Parameter dictionary from eprload (optional)
        order: Polynomial order (int for same order in both directions, or (order_x, order_y))
        exclude_center: If True, exclude center region from fitting
        center_fraction: Fraction of data to exclude from center
        manual_regions: List of manually selected regions as [((x1,x2), (y1,y2)), ...]
        region_mode: 'exclude' to exclude manual_regions, 'include' to use only manual_regions
        interactive: If True, open interactive region selector
        use_real_part: If True, use real part of complex data

    Returns:
        tuple: (corrected_data, baseline)
    """
    if y is None or y.ndim != 2:
        raise ValueError("y must be a 2D array")

    ny, nx = y.shape

    # Handle polynomial order
    if isinstance(order, int):
        order_x = order_y = order
    else:
        order_x, order_y = order

    # Create coordinate meshgrids
    if x is None:
        # No coordinates provided, use indices
        x_1d = np.arange(nx)
        y_1d = np.arange(ny)
        X, Y = np.meshgrid(x_1d, y_1d)
    elif isinstance(x, list) and len(x) == 2:
        # Two 1D coordinate arrays provided
        X, Y = np.meshgrid(x[1], x[0])  # Note: meshgrid swaps order
    elif isinstance(x, np.ndarray):
        if x.ndim == 1:
            # Single 1D array, assume it's x-coordinates
            y_1d = np.arange(ny)
            X, Y = np.meshgrid(x, y_1d)
        elif x.ndim == 3 and x.shape[0] == 2:
            # Two 2D meshgrids provided
            X, Y = x[1], x[0]  # Note: convention difference
        else:
            raise ValueError("Invalid x coordinate format for 2D data")
    else:
        raise ValueError("x must be None, list of 1D arrays, or ndarray")

    # Handle complex data
    if np.iscomplexobj(y):
        if use_real_part:
            data_for_fitting = np.real(y)
            logger.info("ℹ Using real part of complex 2D data for fitting")
        else:
            data_for_fitting = np.abs(y)
            logger.info("ℹ Using magnitude of complex 2D data for fitting")
    else:
        data_for_fitting = y

    # Interactive region selection if requested
    selected_regions = manual_regions if manual_regions is not None else []

    if interactive:
        if not is_interactive_available():
            logger.warning("⚠️  Interactive selection may not work in this environment.")

        logger.info("🖱️ Interactive 2D region selection enabled...")
        selected_regions = interactive_select_regions_2d(
            X,
            Y,
            data_for_fitting,
            f"Select regions to {region_mode.upper()} from baseline fitting",
        )
        logger.info(f"✅ Selected {len(selected_regions)} regions")

    # Create baseline mask
    mask = get_baseline_regions_2d(
        X,
        Y,
        data_for_fitting,
        exclude_center=exclude_center and not manual_regions,
        center_fraction=center_fraction,
        manual_regions=selected_regions,
        region_mode=region_mode,
    )

    # Prepare data for fitting
    X_flat = X[mask].flatten()
    Y_flat = Y[mask].flatten()
    Z_flat = data_for_fitting[mask].flatten()

    # Remove NaN values
    valid = np.isfinite(Z_flat)
    X_flat = X_flat[valid]
    Y_flat = Y_flat[valid]
    Z_flat = Z_flat[valid]

    min_points = (order_x + 1) * (order_y + 1)
    if len(Z_flat) < min_points:
        warnings.warn(
            f"Not enough points ({len(Z_flat)}) for polynomial order ({order_x}, {order_y})"
        )
        return y, np.zeros_like(data_for_fitting)

    try:
        # Prepare initial guess for polynomial coefficients
        initial_guess = np.zeros(min_points)
        initial_guess[0] = np.mean(Z_flat)  # Constant term

        # Fit 2D polynomial
        popt, _ = curve_fit(
            polynomial_2d, (X_flat, Y_flat), Z_flat, p0=initial_guess, maxfev=5000
        )

        # Evaluate baseline over full grid
        baseline = polynomial_2d((X.flatten(), Y.flatten()), *popt).reshape(y.shape)

        # Subtract baseline - preserve data type (real or complex)
        if np.iscomplexobj(y):
            if use_real_part:
                # Subtract from real part only
                corrected_data = y.copy()
                corrected_data.real -= baseline
            else:
                # Subtract from magnitude (this is tricky, use phase-preserving approach)
                magnitude = np.abs(y)
                phase = np.angle(y)
                corrected_magnitude = magnitude - baseline
                corrected_data = corrected_magnitude * np.exp(1j * phase)
        else:
            corrected_data = y - baseline

        return corrected_data, baseline

    except Exception as e:
        warnings.warn(f"2D polynomial fitting failed: {e}")
        return y, np.zeros_like(data_for_fitting)


def _prepare_data_for_exponential_fitting(
    x,
    y,
    use_real_part=True,
    exclude_initial=0,
    exclude_final=0,
    manual_regions=None,
    region_mode="exclude",
):
    """
    Prepare data for exponential baseline fitting.

    This helper function handles data preprocessing common to all exponential
    baseline correction functions.
    """
    if y is None or y.ndim != 1:
        raise ValueError("y must be a 1D array")

    n_points = len(y)

    # Create x-coordinates if needed
    if x is None or len(x) != n_points:
        x_coords = np.arange(n_points)
    else:
        x_coords = x

    # Handle complex data
    if np.iscomplexobj(y):
        if use_real_part:
            data_for_fitting = np.real(y)
            logger.info("ℹ Using real part of complex data for fitting")
        else:
            data_for_fitting = np.abs(y)
            logger.info("ℹ Using magnitude of complex data for fitting")
    else:
        data_for_fitting = y

    # Create baseline mask
    mask = get_baseline_regions_1d(
        x_coords,
        data_for_fitting,
        exclude_center=False,  # Don't exclude center for time-domain data
        exclude_initial=exclude_initial,
        exclude_final=exclude_final,
        manual_regions=manual_regions,
        region_mode=region_mode,
    )

    # Apply mask and remove invalid points
    x_fit = x_coords[mask]
    y_fit = data_for_fitting[mask]

    valid = np.isfinite(y_fit) & (y_fit > 0)  # Exponentials need positive values
    x_fit = x_fit[valid]
    y_fit = y_fit[valid]

    return x_coords, data_for_fitting, x_fit, y_fit


def _smart_exponential_initial_guess(x_fit, y_fit, model_type="stretched"):
    """
    Generate smart initial parameter guesses for exponential models.
    """
    if len(x_fit) == 0 or len(y_fit) == 0:
        raise ValueError("No valid data points for fitting")

    x_min, x_max = x_fit.min(), x_fit.max()
    y_min, y_max = y_fit.min(), y_fit.max()

    # Common parameters
    amplitude = y_max - y_min
    offset = y_min
    tau = (x_max - x_min) / 3  # Rough time constant

    if model_type == "stretched":
        # Stretched exponential parameters
        beta = 1.0  # Start with simple exponential
        return [amplitude, tau, beta, offset]

    elif model_type == "bi_exponential":
        # Bi-exponential parameters
        A1 = amplitude * 0.6
        A2 = amplitude * 0.4
        tau1 = tau * 0.3  # Fast component
        tau2 = tau * 3.0  # Slow component
        return [A1, tau1, A2, tau2, offset]

    else:  # simple exponential
        return [amplitude, tau, offset]


def baseline_stretched_exponential_1d(
    x: Union[np.ndarray, None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    use_real_part: bool = True,
    exclude_initial: int = 0,
    exclude_final: int = 0,
    manual_regions: Optional[List[Tuple[float, float]]] = None,
    region_mode: str = "exclude",
    interactive: bool = False,
    beta_range: Tuple[float, float] = (0.01, 5.0),
    initial_guess: Optional[Dict[str, float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stretched exponential baseline correction for 1D EPR data.

    Fits and removes a stretched exponential baseline of the form:
    baseline = offset + A * exp(-(x/tau)^beta)

    Particularly useful for EPR relaxation data (T1, T2 measurements) where
    the signal decays with a stretched exponential characteristic.

    Args:
        x: X-axis data from eprload (can be None)
        y: 1D spectral data array from eprload
        params: Parameter dictionary from eprload (optional)
        use_real_part: If True, fit only real part of complex data
        exclude_initial: Number of initial points to exclude from fitting
        exclude_final: Number of final points to exclude from fitting
        manual_regions: List of manually selected regions as [(x1, x2), ...]
        region_mode: 'exclude' to exclude manual_regions, 'include' to use only manual_regions
        interactive: If True, open interactive region selector
        beta_range: Range for stretching exponent (min, max)
        initial_guess: Dictionary with initial parameter guesses {'A': ..., 'tau': ..., 'beta': ..., 'offset': ...}

    Returns:
        tuple: (corrected_data, baseline)
    """
    try:
        # Prepare data
        x_coords, data_for_fitting, x_fit, y_fit = (
            _prepare_data_for_exponential_fitting(
                x,
                y,
                use_real_part,
                exclude_initial,
                exclude_final,
                manual_regions,
                region_mode,
            )
        )

        if interactive:
            if not is_interactive_available():
                logger.warning(
                    "⚠️  Interactive selection may not work in this environment."
                )

            logger.info(
                "🖱️ Interactive region selection for stretched exponential fitting..."
            )
            selected_regions = interactive_select_regions_1d(
                x_coords,
                data_for_fitting,
                "Select regions to include in stretched exponential fitting",
            )

            # Re-prepare data with interactive regions
            x_coords, data_for_fitting, x_fit, y_fit = (
                _prepare_data_for_exponential_fitting(
                    x,
                    y,
                    use_real_part,
                    exclude_initial,
                    exclude_final,
                    selected_regions,
                    "include",
                )
            )

        if len(y_fit) < 4:  # Need at least 4 points for 4 parameters
            warnings.warn(
                f"Not enough points ({len(y_fit)}) for stretched exponential fitting"
            )
            return y, np.zeros_like(data_for_fitting)

        # Initial parameter guess
        if initial_guess:
            p0 = [
                initial_guess.get("A", 1000),
                initial_guess.get("tau", 1000),
                initial_guess.get("beta", 1.0),
                initial_guess.get("offset", 0),
            ]
        else:
            p0 = _smart_exponential_initial_guess(x_fit, y_fit, "stretched")

        logger.debug(
            f"🔧 Initial guesses: A={p0[0]:.2e}, tau={p0[1]:.2e}, beta={p0[2]:.2f}, offset={p0[3]:.2e}"
        )

        # Parameter bounds
        bounds = (
            [0, 0, beta_range[0], -np.inf],  # Lower bounds
            [np.inf, np.inf, beta_range[1], np.inf],  # Upper bounds
        )

        # Fit stretched exponential
        popt, pcov = curve_fit(
            stretched_exponential_1d, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=5000
        )

        A_fit, tau_fit, beta_fit, offset_fit = popt

        # Calculate parameter uncertainties
        try:
            param_errors = np.sqrt(np.diag(pcov))
            logger.info(
                f"✅ Fit successful: A={A_fit:.2e}, tau={tau_fit:.2e}, beta={beta_fit:.2f}, offset={offset_fit:.2e}"
            )
            logger.info(
                f"📊 Parameter uncertainties: ΔA={param_errors[0]:.2e}, Δτ={param_errors[1]:.2e}, Δβ={param_errors[2]:.3f}, Δoffset={param_errors[3]:.2e}"
            )
        except:
            logger.info(
                f"✅ Fit successful: A={A_fit:.2e}, tau={tau_fit:.2e}, beta={beta_fit:.2f}, offset={offset_fit:.2e}"
            )

        # Evaluate baseline over full range
        baseline = stretched_exponential_1d(x_coords, *popt)

        # Subtract baseline from original data
        corrected_data = y - baseline

        return corrected_data, baseline

    except Exception as e:
        warnings.warn(f"Stretched exponential fitting failed: {e}")
        return y, np.zeros_like(y if not np.iscomplexobj(y) else np.real(y))


def baseline_bi_exponential_1d(
    x: Union[np.ndarray, None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    use_real_part: bool = True,
    exclude_initial: int = 0,
    exclude_final: int = 0,
    manual_regions: Optional[List[Tuple[float, float]]] = None,
    region_mode: str = "exclude",
    interactive: bool = False,
    tau_ratio_min: float = 2.5,
    initial_guess: Optional[Dict[str, float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Bi-exponential baseline correction for 1D EPR data.

    Fits and removes a bi-exponential baseline of the form:
    baseline = offset + A1*exp(-x/tau1) + A2*exp(-x/tau2)

    Useful for EPR data with multiple decay components or relaxation pathways.

    Args:
        x: X-axis data from eprload (can be None)
        y: 1D spectral data array from eprload
        params: Parameter dictionary from eprload (optional)
        use_real_part: If True, fit only real part of complex data
        exclude_initial: Number of initial points to exclude from fitting
        exclude_final: Number of final points to exclude from fitting
        manual_regions: List of manually selected regions as [(x1, x2), ...]
        region_mode: 'exclude' to exclude manual_regions, 'include' to use only manual_regions
        interactive: If True, open interactive region selector
        tau_ratio_min: Minimum ratio between tau2 and tau1 to ensure component separation
        initial_guess: Dictionary with initial parameter guesses

    Returns:
        tuple: (corrected_data, baseline)
    """
    try:
        # Prepare data
        x_coords, data_for_fitting, x_fit, y_fit = (
            _prepare_data_for_exponential_fitting(
                x,
                y,
                use_real_part,
                exclude_initial,
                exclude_final,
                manual_regions,
                region_mode,
            )
        )

        if interactive:
            if not is_interactive_available():
                logger.warning(
                    "⚠️  Interactive selection may not work in this environment."
                )

            logger.info("🖱️ Interactive region selection for bi-exponential fitting...")
            selected_regions = interactive_select_regions_1d(
                x_coords,
                data_for_fitting,
                "Select regions to include in bi-exponential fitting",
            )

            # Re-prepare data with interactive regions
            x_coords, data_for_fitting, x_fit, y_fit = (
                _prepare_data_for_exponential_fitting(
                    x,
                    y,
                    use_real_part,
                    exclude_initial,
                    exclude_final,
                    selected_regions,
                    "include",
                )
            )

        if len(y_fit) < 5:  # Need at least 5 points for 5 parameters
            warnings.warn(
                f"Not enough points ({len(y_fit)}) for bi-exponential fitting"
            )
            return y, np.zeros_like(data_for_fitting)

        # Initial parameter guess
        if initial_guess:
            p0 = [
                initial_guess.get("A1", 500),
                initial_guess.get("tau1", 100),
                initial_guess.get("A2", 500),
                initial_guess.get("tau2", 1000),
                initial_guess.get("offset", 0),
            ]
        else:
            p0 = _smart_exponential_initial_guess(x_fit, y_fit, "bi_exponential")

        logger.debug(
            f"🔧 Initial guesses: A1={p0[0]:.2e}, τ1={p0[1]:.2e}, A2={p0[2]:.2e}, τ2={p0[3]:.2e}, offset={p0[4]:.2e}"
        )

        # Custom fitting function with tau ratio constraint
        def constrained_bi_exponential(x, A1, tau1, A2, tau2, offset):
            # Enforce tau2 > tau1 * tau_ratio_min
            if tau2 < tau1 * tau_ratio_min:
                tau2 = tau1 * tau_ratio_min
            return bi_exponential_1d(x, A1, tau1, A2, tau2, offset)

        # Parameter bounds
        x_range = x_fit.max() - x_fit.min()
        bounds = (
            [0, 0, 0, 0, -np.inf],  # Lower bounds
            [np.inf, x_range * 2, np.inf, x_range * 10, np.inf],  # Upper bounds
        )

        # Fit bi-exponential
        popt, pcov = curve_fit(
            constrained_bi_exponential, x_fit, y_fit, p0=p0, bounds=bounds, maxfev=10000
        )

        A1_fit, tau1_fit, A2_fit, tau2_fit, offset_fit = popt

        # Ensure tau ordering
        if tau2_fit < tau1_fit * tau_ratio_min:
            tau2_fit = tau1_fit * tau_ratio_min

        logger.info(
            f"✅ Fit successful: A1={A1_fit:.2e}, τ1={tau1_fit:.2e}, A2={A2_fit:.2e}, τ2={tau2_fit:.2e}, offset={offset_fit:.2e}"
        )
        logger.info(f"📊 Time constant ratio: τ2/τ1 = {tau2_fit/tau1_fit:.2f}")

        # Evaluate baseline over full range
        baseline = bi_exponential_1d(
            x_coords, A1_fit, tau1_fit, A2_fit, tau2_fit, offset_fit
        )

        # Subtract baseline from original data
        corrected_data = y - baseline

        return corrected_data, baseline

    except Exception as e:
        warnings.warn(f"Bi-exponential fitting failed: {e}")
        return y, np.zeros_like(y if not np.iscomplexobj(y) else np.real(y))
