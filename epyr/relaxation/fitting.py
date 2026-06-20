"""
T1/T2 relaxation fitting module.

Fit time-domain EPR relaxation data (T1 recovery, T2 echo decay) with
mono-exponential, stretched-exponential, bi-exponential, inversion- or
saturation-recovery, or combined homogeneous/spectral-diffusion models.
"""

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from .._plot_utils import get_or_create_subplots
from ..logging_config import get_logger
from .models import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)

logger = get_logger(__name__)

SUPPORTED_MODELS = [
    "mono_exponential",
    "stretched_exponential",
    "biexponential",
    "inversion_recovery",
    "saturation_recovery",
    "gamma_gaussian_decay",
]


@dataclass
class RelaxationFitResult:
    """
    Container for relaxation fit results.

    Attributes
    ----------
    model : str
        Name of the fitted relaxation model.
    parameters : dict
        Fitted parameter values keyed by parameter name.
    parameter_errors : dict
        Standard errors of fitted parameters (square root of covariance diagonal).
    fitted_curve : np.ndarray
        Model evaluated at the fitted points (t_fit).
    residuals : np.ndarray
        Data minus model at the fitted points.
    r_squared : float
        Coefficient of determination R-squared.
    chi_squared : float
        Reduced chi-squared: sum of squared residuals divided by degrees of freedom.
    success : bool
        True if curve_fit converged.
    message : str
        Convergence message or error description.
    covariance_matrix : np.ndarray or None
        Full parameter covariance matrix returned by curve_fit.
    t_fit : np.ndarray or None
        Time values used for fitting, after NaN removal and masking.
    """

    model: str
    parameters: Dict[str, float]
    parameter_errors: Dict[str, float]
    fitted_curve: np.ndarray
    residuals: np.ndarray
    r_squared: float
    chi_squared: float
    success: bool
    message: str
    covariance_matrix: Optional[np.ndarray] = None
    t_fit: Optional[np.ndarray] = None

    def summary(self) -> str:
        """Return a formatted string summarizing fit quality and parameters."""
        lines = [f"=== Relaxation Fit Results - {self.model} ==="]
        lines.append(f"Success: {self.success}")
        if not self.success:
            lines.append(f"Error: {self.message}")
            return "\n".join(lines)

        lines.append(f"R2 = {self.r_squared:.6f}")
        lines.append(f"chi2 = {self.chi_squared:.6f}")
        lines.append("\nParameters:")

        for param, value in self.parameters.items():
            if param in self.parameter_errors:
                error = self.parameter_errors[param]
                lines.append(f"  {param}: {value:.6g} +/- {error:.6g}")
            else:
                lines.append(f"  {param}: {value:.6g}")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()


def fit_relaxation(
    t_data: np.ndarray,
    y_data: np.ndarray,
    model: str = "mono_exponential",
    initial_params: Optional[Dict[str, float]] = None,
    bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
    time_unit: str = "",
    **fit_kwargs,
) -> RelaxationFitResult:
    """
    Fit time-domain relaxation data with the specified decay/recovery model.

    Parameters
    ----------
    t_data : np.ndarray
        Time axis, in the unit chosen by the caller (e.g. ns or us).
    y_data : np.ndarray
        Relaxation signal (real-valued; take np.abs() of a complex echo
        signal before calling this function).
    model : str, optional
        Relaxation model name (default: 'mono_exponential'). See
        SUPPORTED_MODELS for the full list as models are added.
    initial_params : dict, optional
        Initial parameter guesses. Auto-estimated from data if None.
    bounds : dict, optional
        Parameter bounds as {name: (lower, upper)}, overriding data-derived
        defaults.
    mask : np.ndarray of bool, optional
        Boolean array of the same length as t_data. True selects a point for
        fitting; False excludes it. If None, all non-NaN points are used.
    plot : bool, optional
        Display a fit plot with residuals panel (default: True).
    time_unit : str, optional
        Cosmetic time-unit label for the plot axis and summary (e.g. 'ns').
        Does not affect fitting (default: '').
    **fit_kwargs
        Additional keyword arguments passed to scipy.optimize.curve_fit.

    Returns
    -------
    RelaxationFitResult
        Fit parameters, errors, statistics, fitted curve, and residuals.

    Examples
    --------
    >>> from epyr.relaxation import fit_relaxation
    >>> import numpy as np
    >>> t = np.linspace(0, 100, 200)
    >>> y = 5.0 * np.exp(-t / 15.0) + 1.0
    >>> result = fit_relaxation(t, y, 'mono_exponential', plot=False)
    >>> print(result.summary())
    """

    t_data = np.asarray(t_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    if len(t_data) != len(y_data):
        raise ValueError("t_data and y_data must have same length")

    if len(t_data) < 4:
        raise ValueError("Need at least 4 data points for fitting")

    if model not in SUPPORTED_MODELS:
        raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")

    valid_mask = ~(np.isnan(t_data) | np.isnan(y_data))
    if not np.any(valid_mask):
        raise ValueError("No valid data points (all NaN)")

    t_clean = t_data[valid_mask]
    y_clean = y_data[valid_mask]

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != t_data.shape:
            raise ValueError(
                f"mask shape {mask.shape} does not match t_data shape {t_data.shape}"
            )
        fit_mask = mask[valid_mask]
        if not np.any(fit_mask):
            raise ValueError("No valid data points after applying mask")
        t_fit = t_clean[fit_mask]
        y_fit = y_clean[fit_mask]
    else:
        t_fit = t_clean
        y_fit = y_clean

    fit_func, param_names, param_bounds = _get_fit_function(model)

    if initial_params is None:
        initial_params = _estimate_initial_params(t_fit, y_fit, model)

    initial_params = _validate_initial_params(initial_params, param_names, t_fit, y_fit)

    lower_bounds, upper_bounds = _setup_bounds(
        param_names, bounds, param_bounds, initial_params, t_fit, y_fit
    )

    p0 = [initial_params[name] for name in param_names]
    bounds_tuple = (lower_bounds, upper_bounds)

    default_kwargs = {"maxfev": 10000, "method": "trf"}
    default_kwargs.update(fit_kwargs)

    try:
        popt, pcov = curve_fit(
            fit_func, t_fit, y_fit, p0=p0, bounds=bounds_tuple, **default_kwargs
        )

        y_fitted = fit_func(t_fit, *popt)
        residuals = y_fit - y_fitted

        ss_res: float = float(np.sum(residuals**2))
        ss_tot: float = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        chi_squared = ss_res / (len(y_fit) - len(popt))

        param_dict = {name: value for name, value in zip(param_names, popt)}
        param_errors = {}

        if pcov is not None:
            param_std_errors = np.sqrt(np.diag(pcov))
            param_errors = {
                name: error for name, error in zip(param_names, param_std_errors)
            }

        result = RelaxationFitResult(
            model=model,
            parameters=param_dict,
            parameter_errors=param_errors,
            fitted_curve=y_fitted,
            residuals=residuals,
            r_squared=r_squared,
            chi_squared=chi_squared,
            success=True,
            message="Fit converged successfully",
            covariance_matrix=pcov,
            t_fit=t_fit,
        )

        if plot:
            _plot_fit_results(t_fit, y_fit, result, model, time_unit)

        return result

    except Exception as e:
        return RelaxationFitResult(
            model=model,
            parameters={},
            parameter_errors={},
            fitted_curve=np.array([]),
            residuals=np.array([]),
            r_squared=0.0,
            chi_squared=np.inf,
            success=False,
            message=str(e),
        )


def _get_fit_function(
    model: str,
) -> Tuple[Callable[..., np.ndarray], List[str], Dict[str, Tuple[float, float]]]:
    """
    Build the fit function, parameter names, and default bounds for a model.

    Parameters
    ----------
    model : str
        Relaxation model name.

    Returns
    -------
    fit_func : callable
        Function compatible with scipy.optimize.curve_fit, signature
        fit_func(t, *params).
    param_names : list of str
        Ordered parameter names matching fit_func's signature.
    param_bounds : dict
        Default bounds for each parameter as {name: (lower, upper)}.
    """

    if model == "mono_exponential":

        def fit_func(t, amplitude, T, offset):
            return mono_exponential(t, amplitude, T, offset)

        param_names = ["amplitude", "T", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "stretched_exponential":

        def fit_func(t, amplitude, T, beta, offset):
            return stretched_exponential(t, amplitude, T, beta, offset)

        param_names = ["amplitude", "T", "beta", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "beta": (0.05, 5.0),
            "offset": (-np.inf, np.inf),
        }

    elif model == "gamma_gaussian_decay":

        def fit_func(t, amplitude, Gamma0, GammaG, offset):
            return gamma_gaussian_decay(t, amplitude, Gamma0, GammaG, offset)

        param_names = ["amplitude", "Gamma0", "GammaG", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "Gamma0": (0.0, np.inf),
            "GammaG": (0.0, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "biexponential":

        def fit_func(t, amplitude1, tau1, amplitude2, tau2, offset):
            return biexponential(t, amplitude1, tau1, amplitude2, tau2, offset)

        param_names = ["amplitude1", "tau1", "amplitude2", "tau2", "offset"]
        param_bounds = {
            "amplitude1": (-np.inf, np.inf),
            "tau1": (1e-9, np.inf),
            "amplitude2": (-np.inf, np.inf),
            "tau2": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "inversion_recovery":

        def fit_func(t, amplitude, T1, offset):
            return inversion_recovery(t, amplitude, T1, offset)

        param_names = ["amplitude", "T1", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T1": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "saturation_recovery":

        def fit_func(t, amplitude, T1, offset):
            return saturation_recovery(t, amplitude, T1, offset)

        param_names = ["amplitude", "T1", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T1": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    else:
        raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")

    return fit_func, param_names, param_bounds


def _estimate_decay_rate(
    t: np.ndarray, y: np.ndarray, offset_guess: float, amplitude_guess: float
) -> float:
    """
    Estimate a decay rate via log-linear regression on the part of the curve
    that lies well above the asymptotic offset.

    Parameters
    ----------
    t : np.ndarray
        Time axis.
    y : np.ndarray
        Relaxation signal.
    offset_guess : float
        Estimated asymptotic offset.
    amplitude_guess : float
        Estimated amplitude at t = 0, relative to offset.

    Returns
    -------
    float
        Estimated decay rate, in inverse time units. Falls back to the
        inverse of the time range if too few points qualify for the
        regression.
    """

    deviation = np.abs(y - offset_guess)
    threshold = 0.05 * abs(amplitude_guess) if amplitude_guess != 0 else 0.0
    usable = deviation > threshold

    if np.sum(usable) >= 2:
        log_dev = np.log(deviation[usable])
        slope, _ = np.polyfit(t[usable], log_dev, 1)
        if slope < 0:
            return float(-slope)

    time_span = float(t[-1] - t[0]) or 1.0
    return 1.0 / time_span


def _estimate_initial_params(
    t: np.ndarray, y: np.ndarray, model: str
) -> Dict[str, float]:
    """
    Estimate initial fit parameters from the data.

    Parameters
    ----------
    t : np.ndarray
        Time axis.
    y : np.ndarray
        Relaxation signal.
    model : str
        Relaxation model name.

    Returns
    -------
    dict
        Initial parameter estimates keyed by parameter name.
    """

    offset_guess = float(y[-1])
    span = float(y[0] - y[-1])
    amplitude_guess = span if span != 0 else float(y.max() - y.min()) or 1.0

    if model in ("mono_exponential", "stretched_exponential", "gamma_gaussian_decay"):
        rate_guess = _estimate_decay_rate(t, y, offset_guess, amplitude_guess)
        T_guess = 1.0 / rate_guess if rate_guess > 0 else (t[-1] - t[0]) / 2

        if model == "mono_exponential":
            return {"amplitude": amplitude_guess, "T": T_guess, "offset": offset_guess}

        if model == "stretched_exponential":
            return {
                "amplitude": amplitude_guess,
                "T": T_guess,
                "beta": 1.0,
                "offset": offset_guess,
            }

        return {
            "amplitude": amplitude_guess,
            "Gamma0": rate_guess,
            "GammaG": rate_guess,
            "offset": offset_guess,
        }

    if model == "biexponential":
        time_span = float(t[-1] - t[0]) or 1.0
        return {
            "amplitude1": amplitude_guess / 2,
            "tau1": max(time_span * 0.1, 1e-6),
            "amplitude2": amplitude_guess / 2,
            "tau2": max(time_span, 1e-6),
            "offset": offset_guess,
        }

    if model in ("inversion_recovery", "saturation_recovery"):
        midpoint = y[0] + 0.5 * (y[-1] - y[0])
        crossing_idx = int(np.argmin(np.abs(y - midpoint)))
        half_life = max(float(t[crossing_idx] - t[0]), float(t[1] - t[0]))
        T1_guess = half_life / np.log(2)
        if model == "inversion_recovery":
            amplitude_guess = abs(float(y[-1] - y[0]) / 3.0) or 1.0
            offset_guess = float(y[0] + 2.0 * amplitude_guess)
        else:
            amplitude_guess = abs(float(y[-1] - y[0])) or 1.0
            offset_guess = float(y[0])
        return {
            "amplitude": amplitude_guess,
            "T1": T1_guess,
            "offset": offset_guess,
        }

    raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")


def _validate_initial_params(
    initial_params: Dict[str, float],
    param_names: List[str],
    t: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    """
    Fill in missing initial parameters with data-derived defaults.

    Parameters
    ----------
    initial_params : dict
        Partial or complete parameter dict.
    param_names : list of str
        Full list of parameters required by the fit function.
    t : np.ndarray
        Time axis, used to derive fallback values.
    y : np.ndarray
        Relaxation signal, used to derive fallback values.

    Returns
    -------
    dict
        Complete parameter dict with all entries in param_names present.
    """

    validated = initial_params.copy()
    time_span = float(t[-1] - t[0]) or 1.0

    for name in param_names:
        if name in validated:
            continue
        if name.startswith("amplitude"):
            validated[name] = float(y.max() - y.min()) or 1.0
        elif name in ("T", "T1") or name.startswith("tau"):
            validated[name] = time_span / 2
        elif name == "beta":
            validated[name] = 1.0
        elif name in ("Gamma0", "GammaG"):
            validated[name] = 1.0 / time_span
        elif name == "offset":
            validated[name] = float(y[-1])
        else:
            validated[name] = 1.0

    return validated


def _setup_bounds(
    param_names: List[str],
    user_bounds: Optional[Dict[str, Tuple[float, float]]],
    default_bounds: Dict[str, Tuple[float, float]],
    initial_params: Dict[str, float],
    t: np.ndarray,
    y: np.ndarray,
) -> Tuple[List[float], List[float]]:
    """
    Build ordered lower/upper bound lists for scipy.optimize.curve_fit.

    Parameters
    ----------
    param_names : list of str
        Ordered parameter names matching the fit function signature.
    user_bounds : dict or None
        Caller-supplied bounds as {name: (lower, upper)}.
    default_bounds : dict
        Default bounds per parameter as {name: (lower, upper)}.
    initial_params : dict
        Current initial values; used to ensure each value lies inside its bounds.
    t : np.ndarray
        Time axis, used to derive data-range bounds.
    y : np.ndarray
        Relaxation signal, used to derive data-range bounds.

    Returns
    -------
    lower_bounds : list of float
    upper_bounds : list of float
    """

    lower_bounds = []
    upper_bounds = []
    time_span = float(t[-1] - t[0]) or 1.0
    y_range = float(y.max() - y.min()) or 1.0

    for name in param_names:
        default_lower, default_upper = default_bounds[name]

        if user_bounds and name in user_bounds:
            lower, upper = user_bounds[name]
        else:
            lower, upper = default_lower, default_upper

        if name.startswith("amplitude"):
            if lower == -np.inf:
                lower = -10 * y_range
            if upper == np.inf:
                upper = 10 * y_range
        elif name in ("T", "T1") or name.startswith("tau"):
            if upper == np.inf:
                upper = 100 * time_span
        elif name in ("Gamma0", "GammaG"):
            dt = float(np.diff(t).mean()) or (time_span / len(t))
            if upper == np.inf:
                upper = 10.0 / dt
        elif name == "offset":
            if lower == -np.inf:
                lower = float(y.min()) - 10 * y_range
            if upper == np.inf:
                upper = float(y.max()) + 10 * y_range

        init_val = initial_params[name]
        if init_val <= lower:
            lower = (
                init_val * 0.1
                if init_val > 0
                else (init_val * 10 if init_val < 0 else -1.0)
            )
        if init_val >= upper:
            upper = (
                init_val * 10
                if init_val > 0
                else (init_val * 0.1 if init_val < 0 else 1.0)
            )

        lower_bounds.append(lower)
        upper_bounds.append(upper)

    return lower_bounds, upper_bounds


def _plot_fit_results(
    t: np.ndarray,
    y: np.ndarray,
    result: RelaxationFitResult,
    model: str,
    time_unit: str = "",
) -> None:
    """
    Plot data, fitted curve, and residuals for a single relaxation fit.

    Parameters
    ----------
    t : np.ndarray
        Time values used for fitting.
    y : np.ndarray
        Relaxation signal values used for fitting.
    result : RelaxationFitResult
        Fit result containing parameters, statistics, and curves.
    model : str
        Relaxation model name, shown in the plot title.
    time_unit : str, optional
        Unit label appended to the time axis (default: '').
    """

    fig, (ax1, ax2) = get_or_create_subplots(
        2, 1, gridspec_kw={"height_ratios": [3, 1]}
    )

    time_label = f"Time ({time_unit})" if time_unit else "Time"

    ax1.plot(t, y, "o", alpha=0.7, label="Data", color="#1f77b4")
    ax1.plot(
        t,
        result.fitted_curve,
        "-",
        label=f"{model} fit",
        color="#d62728",
    )
    ax1.set_xlabel(time_label)
    ax1.set_ylabel("Intensity")
    ax1.set_title(f"Relaxation Fitting - {model}")
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    results_lines = [
        f"R2 = {result.r_squared:.4f}",
        f"chi2 = {result.chi_squared:.2e}",
    ]
    for param, value in result.parameters.items():
        value_str = f"{value:.4g}"
        if param in result.parameter_errors:
            results_lines.append(
                f"{param}: {value_str} +/- {result.parameter_errors[param]:.2g}"
            )
        else:
            results_lines.append(f"{param}: {value_str}")

    ax1.text(
        0.98,
        0.98,
        "\n".join(results_lines),
        transform=ax1.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.8),
    )

    ax2.plot(t, result.residuals, "o-", alpha=0.7, color="#ff7f0e")
    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.set_xlabel(time_label)
    ax2.set_ylabel("Residuals")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


class RelaxationFitComparison(dict):
    """
    Dict of model name to RelaxationFitResult, with a tabular ``__repr__``.

    Behaves exactly like a plain ``dict`` (subscripting, iteration,
    ``.items()``, ...). Printing it shows R-squared, chi-squared, and every
    fitted parameter side by side for all models, so the full comparison is
    visible without looping over individual ``summary()`` calls.
    """

    def __repr__(self) -> str:
        if not self:
            return "RelaxationFitComparison({})"

        param_names = list(
            dict.fromkeys(
                name for result in self.values() for name in result.parameters
            )
        )

        headers = ["model", "success", "R2", "chi2"] + param_names
        rows = []
        for model, result in self.items():
            row = [model, str(result.success)]
            row.append(f"{result.r_squared:.6f}" if result.success else "-")
            row.append(f"{result.chi_squared:.4g}" if result.success else "-")
            for name in param_names:
                if name not in result.parameters:
                    row.append("-")
                    continue
                value_str = f"{result.parameters[name]:.4g}"
                if name in result.parameter_errors:
                    value_str += f" +/- {result.parameter_errors[name]:.2g}"
                row.append(value_str)
            rows.append(row)

        widths = [max(len(cell) for cell in column) for column in zip(headers, *rows)]

        def format_row(cells: List[str]) -> str:
            return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths))

        lines = [format_row(headers), format_row(["-" * w for w in widths])]
        lines.extend(format_row(row) for row in rows)

        return "\n".join(lines)


def fit_multiple_decays(
    t_data: np.ndarray,
    y_data: np.ndarray,
    models: Optional[List[str]] = None,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
) -> RelaxationFitComparison:
    """
    Fit relaxation data with multiple models and compare by reduced chi-squared.

    Parameters
    ----------
    t_data : np.ndarray
        Time axis, in the unit chosen by the caller.
    y_data : np.ndarray
        Relaxation signal.
    models : list of str, optional
        Models to fit. Default: ['mono_exponential', 'stretched_exponential',
        'biexponential']. The recovery models and 'gamma_gaussian_decay' are
        excluded by default since they assume a specific data shape rather
        than being interchangeable candidates for a generic decay.
    mask : np.ndarray of bool, optional
        Boolean array selecting points to include (True = include). Passed
        unchanged to each fit_relaxation call.
    plot : bool, optional
        Display a side-by-side comparison plot (default: True).

    Returns
    -------
    RelaxationFitComparison
        Dict subclass mapping model name to RelaxationFitResult for all
        attempted fits. Printing it shows a comparison table of R-squared,
        chi-squared, and every fitted parameter across all models.

    Notes
    -----
    Models are ranked by reduced chi-squared, not R-squared: R-squared is
    biased toward models with more free parameters, since extra parameters
    can only reduce the residual sum of squares on the same data.
    """

    if models is None:
        models = ["mono_exponential", "stretched_exponential", "biexponential"]

    results = RelaxationFitComparison()

    for model in models:
        try:
            result = fit_relaxation(t_data, y_data, model=model, mask=mask, plot=False)
            results[model] = result
        except Exception as e:
            logger.warning(f"Failed to fit {model}: {e}")
            results[model] = RelaxationFitResult(
                model=model,
                parameters={},
                parameter_errors={},
                fitted_curve=np.array([]),
                residuals=np.array([]),
                r_squared=0.0,
                chi_squared=np.inf,
                success=False,
                message=str(e),
            )

    successful_fits = {k: v for k, v in results.items() if v.success}

    if successful_fits and plot:
        _plot_comparison(t_data, y_data, successful_fits)

    logger.info("=== Relaxation Model Comparison ===")
    for model, result in results.items():
        if result.success:
            logger.info(
                f"{model:22s}: chi2_red = {result.chi_squared:.4g},"
                f" R2 = {result.r_squared:.6f}"
            )
        else:
            logger.info(f"{model:22s}: FAILED - {result.message}")

    if successful_fits:
        best_model = min(
            successful_fits.keys(), key=lambda k: successful_fits[k].chi_squared
        )
        logger.info(f"\nBest fit (lowest chi2_red): {best_model}")

    return results


def _plot_comparison(
    t: np.ndarray, y: np.ndarray, results: Dict[str, RelaxationFitResult]
) -> None:
    """
    Plot fitted curves and residuals for multiple relaxation models.

    Parameters
    ----------
    t : np.ndarray
        Full time data passed to fit_multiple_decays.
    y : np.ndarray
        Full relaxation signal passed to fit_multiple_decays.
    results : dict
        Mapping of model name to RelaxationFitResult; only successful fits
        are drawn.
    """

    fig, (ax1, ax2) = get_or_create_subplots(
        2, 1, gridspec_kw={"height_ratios": [3, 1]}
    )

    colors = ["#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

    ax1.plot(t, y, "o", alpha=0.7, label="Data", color="#1f77b4")

    for i, (model, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            t_plot = result.t_fit if result.t_fit is not None else t
            ax1.plot(
                t_plot,
                result.fitted_curve,
                "-",
                label=f"{model} (chi2_red={result.chi_squared:.3g})",
                color=color,
            )

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Intensity")
    ax1.set_title("Relaxation Fitting - Model Comparison")
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    for i, (model, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            t_plot = result.t_fit if result.t_fit is not None else t
            ax2.plot(
                t_plot,
                result.residuals,
                "o-",
                alpha=0.7,
                label=model,
                color=color,
            )

    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Residuals")
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
