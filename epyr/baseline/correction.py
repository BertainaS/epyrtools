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
    interactive_select_regions_1d,
    interactive_select_regions_2d,
    is_interactive_available,
)

# Import from our new modules
from .models import (
    bi_exponential_1d,
    polynomial_2d,
    stretched_exponential_1d,
)
from .selection import (
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
    """Polynomial baseline correction for 1D EPR data.

    Fits a polynomial to selected baseline regions and subtracts it from
    the full dataset. Suited to CW EPR spectra with smooth drift.

    Parameters
    ----------
    x : np.ndarray or None
        Field axis from :func:`epyr.eprload`. If None or length mismatch,
        falls back to point indices.
    y : np.ndarray
        1D spectrum from :func:`epyr.eprload`.
    params : dict, optional
        Parameter dictionary from :func:`epyr.eprload` (unused for the
        fit itself; reserved for future axis-aware extensions).
    order : int, optional
        Polynomial order. Typical range 1-4. Default 2.
    exclude_center : bool, optional
        Exclude a centred fraction of the data from fitting (where the
        signal is expected). Default True.
    center_fraction : float, optional
        Width of the excluded central window, as a fraction of the data.
        Default 0.3.
    manual_regions : list of (float, float), optional
        Explicit (x_low, x_high) regions. Combined with ``region_mode``.
    region_mode : {'exclude', 'include'}, optional
        ``'exclude'`` removes ``manual_regions`` from the fit;
        ``'include'`` keeps only those regions. Default ``'exclude'``.
    interactive : bool, optional
        Launch the matplotlib region selector. Requires an interactive
        backend. Default False.

    Returns
    -------
    corrected : np.ndarray
        Baseline-subtracted spectrum.
    baseline : np.ndarray
        Fitted polynomial evaluated on the full axis.

    Examples
    --------
    >>> from epyr import eprload, baseline_polynomial_1d
    >>> x, y, params, _ = eprload("examples/data/130406SB_CaWO4_Er_CW_5K_20.DSC")
    >>> corrected, baseline = baseline_polynomial_1d(x, y, params, order=3)
    >>> corrected.shape == y.shape
    True

    Fit only on user-defined wings (exclude resonance region):

    >>> wings = [(3300, 3400), (3550, 3700)]
    >>> corrected, _ = baseline_polynomial_1d(
    ...     x, y, params, order=2,
    ...     manual_regions=wings, region_mode="include",
    ... )
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
            logger.warning("Interactive selection may not" " work in this environment.")
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
        warnings.warn(
            f"Not enough points ({len(y_fit)})" f" for polynomial order {order}",
            stacklevel=2,
        )
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
        warnings.warn(f"Polynomial fitting failed: {e}", stacklevel=2)
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
    """Polynomial baseline correction for 2D EPR data.

    Fits a 2D polynomial surface to selected regions and subtracts it from
    the full dataset. Useful for 2D experiments (DEER, Rabi 2D, etc.).

    Parameters
    ----------
    x : np.ndarray, list of np.ndarray, or None
        Axis data from :func:`epyr.eprload`. A two-element list ``[xa, ya]``
        sets both axes. None falls back to indices.
    y : np.ndarray
        2D spectrum, shape ``(ny, nx)``.
    params : dict, optional
        Parameter dictionary from :func:`epyr.eprload`.
    order : int or (int, int), optional
        Polynomial order. ``int`` applies the same order to both axes;
        a tuple sets ``(order_x, order_y)`` independently. Default 1.
    exclude_center : bool, optional
        Exclude a centred rectangular block of the data from the fit.
        Default True.
    center_fraction : float, optional
        Width of that block along each axis, as a fraction. Default 0.3.
    manual_regions : list of ((x1, x2), (y1, y2)), optional
        Explicit rectangular regions.
    region_mode : {'exclude', 'include'}, optional
        How to combine ``manual_regions`` with the fit mask. Default
        ``'exclude'``.
    interactive : bool, optional
        Launch the matplotlib 2D region selector. Default False.
    use_real_part : bool, optional
        For complex data, fit on ``y.real``. Default False.

    Returns
    -------
    corrected : np.ndarray
        Baseline-subtracted surface.
    baseline : np.ndarray
        Fitted polynomial evaluated on the full grid.

    Examples
    --------
    >>> from epyr import eprload, baseline_polynomial_2d
    >>> x, y, params, _ = eprload("examples/data/Rabi2D_GdCaWO4_13dB_3057G.DSC")
    >>> corrected, baseline = baseline_polynomial_2d(x, y, params, order=(2, 1))
    >>> corrected.shape == y.shape
    True
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
        # Two 1D coordinate arrays: x[0]=field axis (nx), x[1]=second axis (ny)
        # meshgrid(a, b) returns shape (len(b), len(a))
        # We need shape (ny, nx) to match y.shape
        X, Y = np.meshgrid(x[0], x[1])
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
            logger.warning("Interactive selection may not" " work in this environment.")

        logger.info("Interactive 2D region selection" " enabled...")
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
            f"Not enough points ({len(Z_flat)})"
            f" for polynomial order"
            f" ({order_x}, {order_y})",
            stacklevel=2,
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
                # Subtract from magnitude (tricky,
                # use phase-preserving approach)
                magnitude = np.abs(y)
                phase = np.angle(y)
                corrected_magnitude = magnitude - baseline
                corrected_data = corrected_magnitude * np.exp(1j * phase)
        else:
            corrected_data = y - baseline

        return corrected_data, baseline

    except Exception as e:
        warnings.warn(f"2D polynomial fitting failed: {e}", stacklevel=2)
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
    r"""Stretched-exponential baseline correction for 1D EPR data.

    Fits and removes a baseline of the form

    .. math::

        b(x) = \mathrm{offset} + A \exp\!\left[-(x/\tau)^{\beta}\right]

    Typical use: T2 echo decay envelope removal before peak analysis.

    Parameters
    ----------
    x : np.ndarray or None
        Time axis from :func:`epyr.eprload`. Falls back to indices if None
        or length mismatch.
    y : np.ndarray
        1D signal from :func:`epyr.eprload`.
    params : dict, optional
        Parameter dictionary from :func:`epyr.eprload`.
    use_real_part : bool, optional
        For complex ``y``, fit on ``y.real``. Default True.
    exclude_initial, exclude_final : int, optional
        Drop the first / last N points from the fit. Default 0.
    manual_regions : list of (float, float), optional
        Explicit (x_low, x_high) windows.
    region_mode : {'exclude', 'include'}, optional
        Combine ``manual_regions`` with the fit mask. Default ``'exclude'``.
    interactive : bool, optional
        Launch the matplotlib region selector. Default False.
    beta_range : (float, float), optional
        Bounds on the stretching exponent. Default ``(0.01, 5.0)``.
    initial_guess : dict, optional
        Seed values, any subset of ``{'A', 'tau', 'beta', 'offset'}``.

    Returns
    -------
    corrected : np.ndarray
        Baseline-subtracted signal.
    baseline : np.ndarray
        Fitted baseline on the full axis.

    Examples
    --------
    >>> from epyr import eprload, baseline_stretched_exponential_1d
    >>> path = "examples/data/ESEdecay_2D_rotgon_035_07.3K_h80_9.73687GHz_B3.DSC"
    >>> x, y, params, _ = eprload(path)
    >>> # take a single trace from the 2D dataset
    >>> trace = y[0]
    >>> corrected, baseline = baseline_stretched_exponential_1d(
    ...     x[0], trace, params, exclude_initial=5,
    ... )
    >>> corrected.shape == trace.shape
    True
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
                f"Not enough points ({len(y_fit)}) for stretched exponential fitting",
                stacklevel=2,
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
            f"Initial guesses: A={p0[0]:.2e},"
            f" tau={p0[1]:.2e},"
            f" beta={p0[2]:.2f},"
            f" offset={p0[3]:.2e}"
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
                f"Fit successful: A={A_fit:.2e},"
                f" tau={tau_fit:.2e},"
                f" beta={beta_fit:.2f},"
                f" offset={offset_fit:.2e}"
            )
            logger.info(
                f"Parameter uncertainties:"
                f" dA={param_errors[0]:.2e},"
                f" dtau={param_errors[1]:.2e},"
                f" dbeta={param_errors[2]:.3f},"
                f" doffset={param_errors[3]:.2e}"
            )
        except Exception:
            logger.info(
                f"Fit successful: A={A_fit:.2e},"
                f" tau={tau_fit:.2e},"
                f" beta={beta_fit:.2f},"
                f" offset={offset_fit:.2e}"
            )

        # Evaluate baseline over full range
        baseline = stretched_exponential_1d(x_coords, *popt)

        # Subtract baseline from original data
        corrected_data = y - baseline

        return corrected_data, baseline

    except Exception as e:
        warnings.warn(f"Stretched exponential fitting failed: {e}", stacklevel=2)
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
    r"""Bi-exponential baseline correction for 1D EPR data.

    Fits and removes a baseline of the form

    .. math::

        b(x) = \mathrm{offset} + A_1 e^{-x/\tau_1} + A_2 e^{-x/\tau_2}

    Suited to decays with two well-separated relaxation channels.

    Parameters
    ----------
    x : np.ndarray or None
        Time axis from :func:`epyr.eprload`. Falls back to indices if None.
    y : np.ndarray
        1D signal from :func:`epyr.eprload`.
    params : dict, optional
        Parameter dictionary from :func:`epyr.eprload`.
    use_real_part : bool, optional
        For complex ``y``, fit on ``y.real``. Default True.
    exclude_initial, exclude_final : int, optional
        Drop the first / last N points from the fit. Default 0.
    manual_regions : list of (float, float), optional
        Explicit (x_low, x_high) windows.
    region_mode : {'exclude', 'include'}, optional
        How to combine ``manual_regions`` with the fit mask.
    interactive : bool, optional
        Launch the matplotlib region selector. Default False.
    tau_ratio_min : float, optional
        Minimum ratio ``tau2 / tau1`` enforced during fitting to keep
        the two components separable. Default 2.5.
    initial_guess : dict, optional
        Seed values, any subset of ``{'A1', 'tau1', 'A2', 'tau2', 'offset'}``.

    Returns
    -------
    corrected : np.ndarray
        Baseline-subtracted signal.
    baseline : np.ndarray
        Fitted baseline on the full axis.

    Examples
    --------
    >>> from epyr import eprload, baseline_bi_exponential_1d
    >>> path = "examples/data/ESEdecay_2D_rotgon_035_07.3K_h80_9.73687GHz_B3.DSC"
    >>> x, y, params, _ = eprload(path)
    >>> trace = y[0]
    >>> corrected, baseline = baseline_bi_exponential_1d(
    ...     x[0], trace, params, tau_ratio_min=3.0,
    ... )
    >>> corrected.shape == trace.shape
    True
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
                f"Not enough points ({len(y_fit)}) for bi-exponential fitting",
                stacklevel=2,
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
            f"Initial guesses: A1={p0[0]:.2e},"
            f" t1={p0[1]:.2e},"
            f" A2={p0[2]:.2e},"
            f" t2={p0[3]:.2e},"
            f" offset={p0[4]:.2e}"
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
            f"Fit successful: A1={A1_fit:.2e},"
            f" t1={tau1_fit:.2e},"
            f" A2={A2_fit:.2e},"
            f" t2={tau2_fit:.2e},"
            f" offset={offset_fit:.2e}"
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
        warnings.warn(f"Bi-exponential fitting failed: {e}", stacklevel=2)
        return y, np.zeros_like(y if not np.iscomplexobj(y) else np.real(y))
