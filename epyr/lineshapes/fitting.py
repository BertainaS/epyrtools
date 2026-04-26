"""
EPR signal fitting module.

Fit EPR spectra with Gaussian, Lorentzian, Voigt, or pseudo-Voigt lineshapes.
Supports absorption and derivative signals (orders 0-2), optional phase mixing,
and affine baseline correction.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from ..logging_config import get_logger
from .gaussian import gaussian
from .lorentzian import lorentzian
from .lshape import pseudo_voigt
from .voigtian import voigtian

logger = get_logger(__name__)


@dataclass
class FitResult:
    """
    Container for lineshape fit results.

    Attributes
    ----------
    shape_type : str
        Name of the fitted lineshape model.
    parameters : dict
        Fitted parameter values keyed by parameter name.
    parameter_errors : dict
        Standard errors of fitted parameters (square root of covariance diagonal).
    fitted_curve : np.ndarray
        Model evaluated at the fitted points (x_fit).
    residuals : np.ndarray
        Data minus model at the fitted points.
    r_squared : float
        Coefficient of determination R².
    chi_squared : float
        Reduced chi-squared: sum of squared residuals divided by degrees of freedom.
    success : bool
        True if curve_fit converged.
    message : str
        Convergence message or error description.
    covariance_matrix : np.ndarray or None
        Full parameter covariance matrix returned by curve_fit.
    x_fit : np.ndarray or None
        X values used for fitting, after NaN removal and masking.
    """

    shape_type: str
    parameters: Dict[str, float]
    parameter_errors: Dict[str, float]
    fitted_curve: np.ndarray
    residuals: np.ndarray
    r_squared: float
    chi_squared: float
    success: bool
    message: str
    covariance_matrix: Optional[np.ndarray] = None
    x_fit: Optional[np.ndarray] = None

    def summary(self) -> str:
        """Return a formatted string summarizing fit quality and parameters."""
        lines = [f"=== Fit Results - {self.shape_type.title()} ==="]
        lines.append(f"Success: {self.success}")
        if not self.success:
            lines.append(f"Error: {self.message}")
            return "\n".join(lines)

        lines.append(f"R² = {self.r_squared:.6f}")
        lines.append(f"χ² = {self.chi_squared:.6f}")
        lines.append("\nParameters:")

        for param, value in self.parameters.items():
            if param in self.parameter_errors:
                error = self.parameter_errors[param]
                lines.append(f"  {param}: {value:.6f} ± {error:.6f}")
            else:
                lines.append(f"  {param}: {value:.6f}")

        return "\n".join(lines)


def fit_epr_signal(
    x_data: np.ndarray,
    y_data: np.ndarray,
    shape_type: str = "gaussian",
    initial_params: Optional[Dict[str, float]] = None,
    bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    derivative: int = 0,
    fit_phase: bool = False,
    fit_baseline: bool = False,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
    **fit_kwargs,
) -> FitResult:
    """
    Fit EPR signal with specified lineshape function.

    Parameters
    ----------
    x_data : np.ndarray
        Magnetic field axis, in Gauss.
    y_data : np.ndarray
        EPR signal intensity (arbitrary units).
    shape_type : str, optional
        Lineshape model: 'gaussian', 'lorentzian', 'voigt', or 'pseudo_voigt'
        (default: 'gaussian').
    initial_params : dict, optional
        Initial parameter guesses. Auto-estimated from data if None.
        Keys depend on shape_type: basic models use {'center', 'width', 'amplitude'};
        voigt uses {'center', 'gaussian_width', 'lorentzian_width', 'amplitude'};
        pseudo_voigt adds 'alpha'; phase fitting adds 'phase'; baseline fitting
        adds 'baseline_slope' and 'baseline_offset'.
    bounds : dict, optional
        Parameter bounds as {name: (lower, upper)}, overriding data-derived defaults.
    derivative : int, optional
        Derivative order: 0 (absorption), 1, or 2. Fixed, not fitted (default: 0).
    fit_phase : bool, optional
        Fit the phase angle controlling absorption/dispersion mixing (default: False).
    fit_baseline : bool, optional
        Include an affine baseline ``a*x + b`` in the model. Adds 'baseline_slope'
        and 'baseline_offset' as fitted parameters (default: False).
    mask : np.ndarray of bool, optional
        Boolean array of the same length as x_data. True selects a point for
        fitting; False excludes it. Useful to reject artefacts or solvent peaks.
        If None, all non-NaN points are used.
    plot : bool, optional
        Display a fit plot with residuals panel (default: True).
    **fit_kwargs
        Additional keyword arguments passed to scipy.optimize.curve_fit.

    Returns
    -------
    FitResult
        Fit parameters, errors, statistics, fitted curve, and residuals.
        fitted_curve and residuals are defined on x_fit (masked points only).

    Examples
    --------
    >>> from epyr import eprload
    >>> x, y, params, filepath = eprload('data.DTA')
    >>> result = fit_epr_signal(x, y, 'gaussian')
    >>> print(result.summary())
    >>>
    >>> # 1st derivative with phase
    >>> result = fit_epr_signal(x, y, 'gaussian', derivative=1, fit_phase=True)
    >>>
    >>> # Exclude a spectral artefact between 3480-3510 G
    >>> mask = ~((x >= 3480) & (x <= 3510))
    >>> result = fit_epr_signal(x, y, 'gaussian', mask=mask)
    >>>
    >>> # 2nd derivative with explicit bounds
    >>> bounds = {'center': (3400, 3600), 'width': (5, 50), 'phase': (0, 3.14)}
    >>> result = fit_epr_signal(x, y, 'gaussian', derivative=2,
    ...                         fit_phase=True, bounds=bounds)
    """

    # Input validation
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)

    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have same length")

    if len(x_data) < 4:
        raise ValueError("Need at least 4 data points for fitting")

    # Remove NaN values
    valid_mask = ~(np.isnan(x_data) | np.isnan(y_data))
    if not np.any(valid_mask):
        raise ValueError("No valid data points (all NaN)")

    x_clean = x_data[valid_mask]
    y_clean = y_data[valid_mask]

    # Apply user mask
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != x_data.shape:
            raise ValueError(
                f"mask shape {mask.shape} does not match x_data shape {x_data.shape}"
            )
        fit_mask = mask[valid_mask]
        if not np.any(fit_mask):
            raise ValueError("No valid data points after applying mask")
        x_fit = x_clean[fit_mask]
        y_fit = y_clean[fit_mask]
    else:
        fit_mask = np.ones(len(x_clean), dtype=bool)
        x_fit = x_clean
        y_fit = y_clean

    # Validate derivative parameter
    if not isinstance(derivative, int) or derivative < 0 or derivative > 2:
        raise ValueError("derivative must be 0, 1, or 2")

    # Get fitting function and parameter info
    fit_func, param_names, param_bounds = _get_fit_function(
        shape_type, derivative, fit_phase, fit_baseline
    )

    # Estimate initial parameters if not provided
    if initial_params is None:
        initial_params = _estimate_initial_params(
            x_fit, y_fit, shape_type, derivative, fit_phase, fit_baseline
        )

    # Validate and complete initial parameters
    initial_params = _validate_initial_params(
        initial_params, param_names, x_fit, y_fit
    )

    # Setup parameter bounds
    lower_bounds, upper_bounds = _setup_bounds(
        param_names, bounds, param_bounds, initial_params, x_fit, y_fit
    )

    # Prepare parameters for fitting
    p0 = [initial_params[name] for name in param_names]
    bounds_tuple = (lower_bounds, upper_bounds)

    # Set default fitting options
    default_kwargs = {"maxfev": 10000, "method": "trf"}
    default_kwargs.update(fit_kwargs)

    try:
        # Perform the fit on masked points
        popt, pcov = curve_fit(
            fit_func, x_fit, y_fit, p0=p0, bounds=bounds_tuple, **default_kwargs
        )

        # Calculate fitted curve and residuals on fit points
        y_fitted = fit_func(x_fit, *popt)
        residuals = y_fit - y_fitted

        # Calculate statistics
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        chi_squared = ss_res / (len(y_fit) - len(popt))

        # Extract parameter values and errors
        param_dict = {name: value for name, value in zip(param_names, popt)}
        param_errors = {}

        if pcov is not None:
            param_std_errors = np.sqrt(np.diag(pcov))
            param_errors = {
                name: error for name, error in zip(param_names, param_std_errors)
            }

        # Create result object
        result = FitResult(
            shape_type=shape_type,
            parameters=param_dict,
            parameter_errors=param_errors,
            fitted_curve=y_fitted,
            residuals=residuals,
            r_squared=r_squared,
            chi_squared=chi_squared,
            success=True,
            message="Fit converged successfully",
            covariance_matrix=pcov,
            x_fit=x_fit,
        )

        # Create plot if requested
        if plot:
            _plot_fit_results(
                x_fit,
                y_fit,
                result,
                shape_type,
                derivative,
                fit_phase,
                x_all=x_clean,
                y_all=y_clean,
                fit_func=fit_func,
                popt=popt,
            )

        return result

    except Exception as e:
        # Return failed result
        return FitResult(
            shape_type=shape_type,
            parameters={},
            parameter_errors={},
            fitted_curve=np.array([]),
            residuals=np.array([]),
            r_squared=0.0,
            chi_squared=np.inf,
            success=False,
            message=str(e),
        )


def _make_simple_lineshape(
    func, derivative: int, fit_phase: bool
) -> Tuple[callable, List[str], Dict[str, Tuple[float, float]]]:
    """Build wrapper and metadata for single-width lineshapes (gaussian, lorentzian)."""
    if fit_phase:
        def _lineshape(x, center, width, amplitude, phase):
            return amplitude * func(
                x, center, width, derivative=derivative, phase=phase
            )
        param_names = ["center", "width", "amplitude", "phase"]
    else:
        def _lineshape(x, center, width, amplitude):
            return amplitude * func(x, center, width, derivative=derivative)
        param_names = ["center", "width", "amplitude"]

    param_bounds = {
        "center": (-np.inf, np.inf),
        "width": (0.001, np.inf),
        "amplitude": (-np.inf, np.inf),
        "phase": (-np.pi, np.pi),
    }
    return _lineshape, param_names, param_bounds


def _get_fit_function(
    shape_type: str, derivative: int, fit_phase: bool, fit_baseline: bool = False
) -> Tuple[callable, List[str], Dict[str, Tuple[float, float]]]:
    """
    Build the fitting function and parameter metadata for a given lineshape.

    Parameters
    ----------
    shape_type : str
        Lineshape model ('gaussian', 'lorentzian', 'voigt', 'pseudo_voigt').
    derivative : int
        Derivative order (0, 1, or 2).
    fit_phase : bool
        Include phase as a fitted parameter.
    fit_baseline : bool
        Wrap the lineshape with an affine baseline (a*x + b).

    Returns
    -------
    fit_func : callable
        Wrapped lineshape compatible with scipy.optimize.curve_fit.
    param_names : list of str
        Ordered parameter names matching the function signature.
    param_bounds : dict
        Default bounds for each parameter as {name: (lower, upper)}.
    """

    if shape_type == "gaussian":
        _lineshape, param_names, param_bounds = _make_simple_lineshape(
            gaussian, derivative, fit_phase
        )

    elif shape_type == "lorentzian":
        _lineshape, param_names, param_bounds = _make_simple_lineshape(
            lorentzian, derivative, fit_phase
        )

    elif shape_type == "voigt":
        if fit_phase:

            def _lineshape(
                x, center, gaussian_width, lorentzian_width, amplitude, phase
            ):
                return amplitude * voigtian(
                    x,
                    center,
                    (gaussian_width, lorentzian_width),
                    derivative=derivative,
                    phase=phase,
                )

            param_names = [
                "center",
                "gaussian_width",
                "lorentzian_width",
                "amplitude",
                "phase",
            ]
        else:

            def _lineshape(x, center, gaussian_width, lorentzian_width, amplitude):
                return amplitude * voigtian(
                    x, center, (gaussian_width, lorentzian_width), derivative=derivative
                )

            param_names = ["center", "gaussian_width", "lorentzian_width", "amplitude"]

        param_bounds = {
            "center": (-np.inf, np.inf),
            "gaussian_width": (0.001, np.inf),
            "lorentzian_width": (0.001, np.inf),
            "amplitude": (-np.inf, np.inf),
            "phase": (-np.pi, np.pi),
        }

    elif shape_type == "pseudo_voigt":
        if fit_phase:

            def _lineshape(x, center, width, amplitude, alpha, phase):
                return amplitude * pseudo_voigt(
                    x, center, width, eta=alpha, derivative=derivative, phase=phase
                )

            param_names = ["center", "width", "amplitude", "alpha", "phase"]
        else:

            def _lineshape(x, center, width, amplitude, alpha):
                return amplitude * pseudo_voigt(
                    x, center, width, eta=alpha, derivative=derivative
                )

            param_names = ["center", "width", "amplitude", "alpha"]

        param_bounds = {
            "center": (-np.inf, np.inf),
            "width": (0.001, np.inf),
            "amplitude": (-np.inf, np.inf),
            "alpha": (0.0, 1.0),
            "phase": (-np.pi, np.pi),
        }

    else:
        raise ValueError(
            f"Unsupported shape_type: {shape_type}. "
            "Choose from: 'gaussian', 'lorentzian', 'voigt', 'pseudo_voigt'"
        )

    # Wrap with affine baseline if requested
    if fit_baseline:
        n_lineshape_params = len(param_names)
        param_names = param_names + ["baseline_slope", "baseline_offset"]
        param_bounds["baseline_slope"] = (-np.inf, np.inf)
        param_bounds["baseline_offset"] = (-np.inf, np.inf)

        _inner = _lineshape

        def fit_func(x, *args):
            lineshape_args = args[:n_lineshape_params]
            baseline_slope = args[n_lineshape_params]
            baseline_offset = args[n_lineshape_params + 1]
            return _inner(x, *lineshape_args) + baseline_slope * x + baseline_offset

    else:
        fit_func = _lineshape

    return fit_func, param_names, param_bounds


def _estimate_initial_params(
    x: np.ndarray,
    y: np.ndarray,
    shape_type: str,
    derivative: int = 0,
    fit_phase: bool = False,
    fit_baseline: bool = False,
) -> Dict[str, float]:
    """
    Estimate initial fit parameters from the data.

    Dispatches to derivative-specific estimators for center, width, and amplitude,
    then appends shape-specific and optional parameters (phase, baseline).

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss.
    y : np.ndarray
        EPR signal values.
    shape_type : str
        Lineshape model ('gaussian', 'lorentzian', 'voigt', 'pseudo_voigt').
    derivative : int, optional
        Signal derivative order (0, 1, or 2). Default 0.
    fit_phase : bool, optional
        Include a phase parameter estimate. Default False.
    fit_baseline : bool, optional
        Include affine baseline parameter estimates. Default False.

    Returns
    -------
    dict
        Initial parameter estimates keyed by parameter name.
    """

    if derivative == 0:
        # Absorption signal: peak at center
        center, width, amplitude = _estimate_absorption_params(x, y)
    elif derivative == 1:
        # First derivative: bipolar signal, center at zero-crossing
        center, width, amplitude = _estimate_first_derivative_params(x, y, shape_type)
    elif derivative == 2:
        # Second derivative: central extremum with side lobes
        center, width, amplitude = _estimate_second_derivative_params(x, y, shape_type)
    else:
        # Fallback to absorption-like estimation
        center, width, amplitude = _estimate_absorption_params(x, y)

    # Shape-specific parameters
    initial_params = {
        "center": center,
        "amplitude": amplitude,
    }

    if shape_type in ["gaussian", "lorentzian", "pseudo_voigt"]:
        initial_params["width"] = max(width, np.diff(x).mean() * 2)

    if shape_type == "voigt":
        initial_params["gaussian_width"] = max(width * 0.6, np.diff(x).mean() * 2)
        initial_params["lorentzian_width"] = max(width * 0.6, np.diff(x).mean() * 2)

    if shape_type == "pseudo_voigt":
        initial_params["alpha"] = 0.5  # 50/50 mix

    if fit_phase:
        peak_region = np.abs(x - center) < width
        if np.sum(peak_region) > 2:
            y_peak, x_peak = y[peak_region], x[peak_region]
            avg_grad = np.mean(np.gradient(y_peak, x_peak))
            # Use π/4 when a dispersive component is detectable in the peak region
            dispersive = np.abs(avg_grad) > np.abs(np.mean(y_peak)) * 0.1
            initial_params["phase"] = np.pi / 4 if dispersive else 0.0
        else:
            initial_params["phase"] = 0.0

    # Estimate affine baseline parameters from data edges
    if fit_baseline:
        n_edge = max(2, len(x) // 10)
        x_edges = np.concatenate([x[:n_edge], x[-n_edge:]])
        y_edges = np.concatenate([y[:n_edge], y[-n_edge:]])
        if len(x_edges) >= 2:
            coeffs = np.polyfit(x_edges, y_edges, 1)
            initial_params["baseline_slope"] = coeffs[0]
            initial_params["baseline_offset"] = coeffs[1]
        else:
            initial_params["baseline_slope"] = 0.0
            initial_params["baseline_offset"] = 0.0

    return initial_params


def _estimate_absorption_params(
    x: np.ndarray, y: np.ndarray
) -> Tuple[float, float, float]:
    """
    Estimate center, width, and amplitude for an absorption signal (derivative=0).

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss.
    y : np.ndarray
        EPR absorption signal.

    Returns
    -------
    center : float
        Field position of the peak, in Gauss.
    width : float
        Estimated FWHM, in Gauss.
    amplitude : float
        Signed peak-to-peak amplitude.
    """

    amplitude = np.max(y) - np.min(y)
    if amplitude == 0:
        amplitude = 1.0

    if np.abs(np.min(y)) > np.abs(np.max(y)):
        amplitude = -amplitude
        peak_idx = np.argmin(y)
    else:
        peak_idx = np.argmax(y)

    center = x[peak_idx]

    # Estimate FWHM
    if amplitude > 0:
        half_max = np.min(y) + amplitude / 2
        above_half = y >= half_max
    else:
        half_max = np.max(y) + amplitude / 2
        above_half = y <= half_max

    if np.sum(above_half) > 1:
        indices = np.where(above_half)[0]
        width = x[indices[-1]] - x[indices[0]]
        if width <= 0:
            width = (x[-1] - x[0]) / 10
    else:
        width = (x[-1] - x[0]) / 10

    return center, width, amplitude


def _estimate_first_derivative_params(
    x: np.ndarray, y: np.ndarray, shape_type: str
) -> Tuple[float, float, float]:
    """
    Estimate center, width, and amplitude for a first-derivative EPR signal.

    The center is at the zero-crossing between the positive and negative extrema.
    Width is estimated from the extrema separation; amplitude is scaled by comparing
    the data peak-to-peak with a unit-amplitude model derivative.

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss.
    y : np.ndarray
        First-derivative EPR signal.
    shape_type : str
        Lineshape model, used to evaluate the unit-amplitude reference derivative.

    Returns
    -------
    center : float
        Estimated line center, in Gauss.
    width : float
        Estimated linewidth, in Gauss.
    amplitude : float
        Estimated amplitude scaling factor.
    """

    idx_max = np.argmax(y)
    idx_min = np.argmin(y)

    # Center is midpoint between extrema
    center = (x[idx_max] + x[idx_min]) / 2.0

    # Peak separation gives an estimate of the linewidth
    peak_sep = abs(x[idx_max] - x[idx_min])
    # For Gaussian: extrema at ±σ, so sep = 2σ = FWHM/√(2ln2) ≈ FWHM/1.18
    # For Lorentzian: extrema at ±γ/√3, so sep = 2γ/√3 = FWHM/√3 ≈ FWHM/1.73
    # Use ~1.4 as a general compromise
    width = max(peak_sep * 1.4, np.diff(x).mean() * 3)

    # Estimate amplitude by comparing data with unit-amplitude model derivative
    peak_to_peak = np.max(y) - np.min(y)
    amplitude = _estimate_amplitude_from_model(
        x, center, width, peak_to_peak, shape_type, derivative=1
    )

    return center, width, amplitude


def _estimate_second_derivative_params(
    x: np.ndarray, y: np.ndarray, shape_type: str
) -> Tuple[float, float, float]:
    """
    Estimate center, width, and amplitude for a second-derivative EPR signal.

    The central extremum (largest |y|) gives the line center. Width is estimated
    from the separation between the central peak and the flanking side lobes.

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss.
    y : np.ndarray
        Second-derivative EPR signal.
    shape_type : str
        Lineshape model, used to evaluate the unit-amplitude reference.

    Returns
    -------
    center : float
        Estimated line center, in Gauss.
    width : float
        Estimated linewidth, in Gauss.
    amplitude : float
        Estimated amplitude scaling factor.
    """

    # Central extremum (largest absolute value) gives the center
    idx_center = np.argmax(np.abs(y))
    center = x[idx_center]

    # Find the two side lobes (opposite sign from center)
    center_sign = np.sign(y[idx_center])
    opposite = np.where(np.sign(y) == -center_sign)[0]

    if len(opposite) > 1:
        # Distance from center to side lobes
        left_lobes = opposite[opposite < idx_center]
        right_lobes = opposite[opposite > idx_center]

        if len(left_lobes) > 0 and len(right_lobes) > 0:
            left_peak = left_lobes[np.argmin(y[left_lobes] * center_sign)]
            right_peak = right_lobes[np.argmin(y[right_lobes] * center_sign)]
            lobe_sep = x[right_peak] - x[left_peak]
            width = max(lobe_sep * 0.8, np.diff(x).mean() * 3)
        else:
            width = (x[-1] - x[0]) / 10
    else:
        width = (x[-1] - x[0]) / 10

    # Estimate amplitude from model comparison
    peak_to_peak = np.max(y) - np.min(y)
    amplitude = _estimate_amplitude_from_model(
        x, center, width, peak_to_peak, shape_type, derivative=2
    )

    return center, width, amplitude


def _estimate_amplitude_from_model(
    x: np.ndarray,
    center: float,
    width: float,
    data_peak_to_peak: float,
    shape_type: str,
    derivative: int,
) -> float:
    """
    Estimate the amplitude scaling factor from the model peak-to-peak.

    Evaluates the lineshape at unit amplitude with the estimated center and width,
    then scales to match the observed data peak-to-peak.

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss.
    center : float
        Estimated line center, in Gauss.
    width : float
        Estimated linewidth, in Gauss.
    data_peak_to_peak : float
        Observed peak-to-peak amplitude of the signal.
    shape_type : str
        Lineshape model used to compute the reference curve.
    derivative : int
        Derivative order of the reference curve (0, 1, or 2).

    Returns
    -------
    float
        Amplitude scaling factor. Returns data_peak_to_peak if the model
        peak-to-peak is zero or an exception occurs.
    """

    try:
        if shape_type == "voigt":
            y_model = voigtian(
                x, center, (width * 0.6, width * 0.6), derivative=derivative
            )
        elif shape_type == "gaussian":
            y_model = gaussian(x, center, width, derivative=derivative)
        elif shape_type == "lorentzian":
            y_model = lorentzian(x, center, width, derivative=derivative)
        elif shape_type == "pseudo_voigt":
            y_model = pseudo_voigt(x, center, width, eta=0.5)
            # pseudo_voigt doesn't support derivatives natively
            if derivative > 0:
                for _ in range(derivative):
                    y_model = np.gradient(y_model, x)
        else:
            return data_peak_to_peak

        model_ptp = np.max(y_model) - np.min(y_model)
        if model_ptp > 0:
            return data_peak_to_peak / model_ptp
        else:
            return data_peak_to_peak
    except Exception:
        return data_peak_to_peak


def _validate_initial_params(
    initial_params: Dict[str, float],
    param_names: List[str],
    x: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    """
    Fill in missing initial parameters with data-derived defaults.

    Parameters
    ----------
    initial_params : dict
        Partial or complete parameter dict provided by the caller or estimator.
    param_names : list of str
        Full list of parameters required by the fit function.
    x : np.ndarray
        Magnetic field axis, used to derive fallback values.
    y : np.ndarray
        EPR signal, used to derive fallback values.

    Returns
    -------
    dict
        Complete parameter dict with all entries in param_names present.
    """

    validated = initial_params.copy()

    # Ensure all required parameters are present
    for name in param_names:
        if name not in validated:
            if name == "center":
                validated[name] = x[np.argmax(np.abs(y))]
            elif name == "amplitude":
                validated[name] = np.max(y) - np.min(y)
            elif "width" in name:
                validated[name] = (x[-1] - x[0]) / 10
            elif name == "alpha":
                validated[name] = 0.5
            elif name == "phase":
                validated[name] = 0.0
            elif name in ("baseline_slope", "baseline_offset"):
                validated[name] = 0.0
            else:
                validated[name] = 1.0

    return validated


def _setup_bounds(
    param_names: List[str],
    user_bounds: Optional[Dict[str, Tuple[float, float]]],
    default_bounds: Dict[str, Tuple[float, float]],
    initial_params: Dict[str, float],
    x: np.ndarray,
    y: np.ndarray,
) -> Tuple[List[float], List[float]]:
    """
    Build ordered lower and upper bound lists for scipy.optimize.curve_fit.

    User-supplied bounds override defaults; remaining infinite bounds are
    replaced with data-range-derived values to keep the fit physically meaningful.

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
    x : np.ndarray
        Magnetic field axis, used to derive data-range bounds.
    y : np.ndarray
        EPR signal, used to derive data-range bounds.

    Returns
    -------
    lower_bounds : list of float
    upper_bounds : list of float
    """

    lower_bounds = []
    upper_bounds = []

    for name in param_names:
        # Get default bounds
        default_lower, default_upper = default_bounds[name]

        # Override with user bounds if provided
        if user_bounds and name in user_bounds:
            lower, upper = user_bounds[name]
        else:
            lower, upper = default_lower, default_upper

        # Make bounds reasonable based on data
        if name == "center":
            if lower == -np.inf:
                lower = x.min() - (x.max() - x.min())
            if upper == np.inf:
                upper = x.max() + (x.max() - x.min())
        elif name == "amplitude":
            if lower == -np.inf:
                lower = -10 * (np.max(y) - np.min(y))
            if upper == np.inf:
                upper = 10 * (np.max(y) - np.min(y))
        elif "width" in name:
            if upper == np.inf:
                upper = x.max() - x.min()
        elif name == "baseline_slope":
            # Slope bounded by data range ratio
            y_range = np.max(y) - np.min(y) if np.max(y) != np.min(y) else 1.0
            x_range = x.max() - x.min() if x.max() != x.min() else 1.0
            max_slope = 10 * y_range / x_range
            if lower == -np.inf:
                lower = -max_slope
            if upper == np.inf:
                upper = max_slope
        elif name == "baseline_offset":
            y_range = np.max(y) - np.min(y) if np.max(y) != np.min(y) else 1.0
            if lower == -np.inf:
                lower = np.min(y) - 10 * y_range
            if upper == np.inf:
                upper = np.max(y) + 10 * y_range

        # Ensure initial value is within bounds
        init_val = initial_params[name]
        if init_val <= lower:
            if init_val > 0:
                lower = init_val * 0.1
            elif init_val < 0:
                lower = init_val * 10
            else:
                lower = -1.0
        if init_val >= upper:
            if init_val > 0:
                upper = init_val * 10
            elif init_val < 0:
                upper = init_val * 0.1
            else:
                upper = 1.0

        lower_bounds.append(lower)
        upper_bounds.append(upper)

    return lower_bounds, upper_bounds


def _plot_fit_results(
    x: np.ndarray,
    y: np.ndarray,
    result: FitResult,
    shape_type: str,
    derivative: int = 0,
    fit_phase: bool = False,
    x_all: Optional[np.ndarray] = None,
    y_all: Optional[np.ndarray] = None,
    fit_func=None,
    popt=None,
):
    """
    Plot data, fitted curve, and residuals for a single fit.

    Parameters
    ----------
    x : np.ndarray
        X values used for fitting (masked subset).
    y : np.ndarray
        Y values used for fitting (masked subset).
    result : FitResult
        Fit result containing parameters, statistics, and curves.
    shape_type : str
        Lineshape model name, shown in the plot title.
    derivative : int, optional
        Derivative order, shown in the plot title (default: 0).
    fit_phase : bool, optional
        If True, the fitted phase value is shown in the title (default: False).
    x_all : np.ndarray, optional
        Full valid x array before masking. Excluded points are shown in gray
        and the model is drawn over this range.
    y_all : np.ndarray, optional
        Full valid y array before masking.
    fit_func : callable, optional
        Fit function used to evaluate the model over x_all.
    popt : np.ndarray, optional
        Optimal parameter vector from curve_fit.
    """

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [3, 1]}
    )

    has_excluded = x_all is not None and len(x_all) > len(x)

    # Show excluded points in gray background
    if has_excluded:
        ax1.plot(
            x_all,
            y_all,
            "o",
            markersize=4,
            alpha=0.3,
            color="gray",
            label="Excluded",
            zorder=1,
        )

    # Main plot - fitted data points
    ax1.plot(
        x, y, "o", markersize=4, alpha=0.7, label="Data (fitted)", color="#1f77b4",
        zorder=2,
    )

    # Model curve over full x range when available, otherwise only over fit points
    if fit_func is not None and popt is not None and x_all is not None:
        x_model = x_all
        y_model = fit_func(x_model, *popt)
    else:
        x_model = x
        y_model = result.fitted_curve

    ax1.plot(
        x_model,
        y_model,
        "-",
        linewidth=2,
        label=f"{shape_type.title()} fit",
        color="#d62728",
        zorder=3,
    )

    ax1.set_xlabel("Magnetic Field / G")
    ax1.set_ylabel("Intensity")

    # Build title with derivative and phase info
    title_parts = [f"EPR Signal Fitting - {shape_type.title()}"]
    if derivative > 0:
        title_parts.append(f"(d^{derivative})")
    if fit_phase:
        phase_val = result.parameters.get("phase", 0)
        title_parts.append(f"Phase: {phase_val:.3f} rad")

    ax1.set_title(" ".join(title_parts))
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Build fit results text with parameters and errors
    results_lines = [
        f"R² = {result.r_squared:.4f}",
        f"χ² = {result.chi_squared:.2e}",
    ]

    for param, value in result.parameters.items():
        value_str = f"{value:.4g}"
        if (
            param in result.parameter_errors
            and result.parameter_errors[param] is not None
        ):
            error_str = f"{result.parameter_errors[param]:.2g}"
            results_lines.append(f"{param}: {value_str} ± {error_str}")
        else:
            results_lines.append(f"{param}: {value_str}")

    textstr = "\n".join(results_lines)
    props = dict(boxstyle="round", facecolor="lightblue", alpha=0.8)
    ax1.text(
        0.02,
        0.98,
        textstr,
        transform=ax1.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="left",
        bbox=props,
    )

    # Residuals plot
    ax2.plot(x, result.residuals, "o-", markersize=3, alpha=0.7, color="#ff7f0e")
    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Magnetic Field / G")
    ax2.set_ylabel("Residuals")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def fit_multiple_shapes(
    x_data: np.ndarray,
    y_data: np.ndarray,
    shapes: List[str] = None,
    derivative: int = 0,
    fit_phase: bool = False,
    fit_baseline: bool = False,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
) -> Dict[str, FitResult]:
    """
    Fit EPR signal with multiple lineshape types and compare results.

    Parameters
    ----------
    x_data : np.ndarray
        Magnetic field axis, in Gauss.
    y_data : np.ndarray
        EPR signal intensity.
    shapes : list of str, optional
        Lineshape models to fit. Default: ['gaussian', 'lorentzian', 'pseudo_voigt'].
    derivative : int, optional
        Derivative order (0, 1, or 2). Fixed, not fitted (default: 0).
    fit_phase : bool, optional
        Fit the phase parameter in each model (default: False).
    fit_baseline : bool, optional
        Include an affine baseline in each model (default: False).
    mask : np.ndarray of bool, optional
        Boolean array selecting points to include (True = include).
        Passed unchanged to each fit_epr_signal call.
    plot : bool, optional
        Display a side-by-side comparison plot (default: True).

    Returns
    -------
    dict
        Mapping of shape name to FitResult for all attempted fits.
    """

    if shapes is None:
        shapes = ["gaussian", "lorentzian", "pseudo_voigt"]

    results = {}

    for shape in shapes:
        try:
            result = fit_epr_signal(
                x_data,
                y_data,
                shape,
                derivative=derivative,
                fit_phase=fit_phase,
                fit_baseline=fit_baseline,
                mask=mask,
                plot=False,
            )
            results[shape] = result
        except Exception as e:
            logger.warning(f"Failed to fit {shape}: {e}")
            results[shape] = FitResult(
                shape_type=shape,
                parameters={},
                parameter_errors={},
                fitted_curve=np.array([]),
                residuals=np.array([]),
                r_squared=0.0,
                chi_squared=np.inf,
                success=False,
                message=str(e),
            )

    # Find best fit based on R²
    successful_fits = {k: v for k, v in results.items() if v.success}

    if successful_fits and plot:
        _plot_comparison(x_data, y_data, successful_fits)

    # Print comparison
    logger.info("=== Fit Comparison ===")
    for shape, result in results.items():
        if result.success:
            logger.info(
                f"{shape:12s}: R² = {result.r_squared:.6f},"
                f" χ² = {result.chi_squared:.2e}"
            )
        else:
            logger.info(f"{shape:12s}: FAILED - {result.message}")

    if successful_fits:
        best_shape = max(
            successful_fits.keys(), key=lambda k: successful_fits[k].r_squared
        )
        logger.info(
            f"\nBest fit: {best_shape}"
            f" (R² = {successful_fits[best_shape].r_squared:.6f})"
        )

    return results


def _plot_comparison(x: np.ndarray, y: np.ndarray, results: Dict[str, FitResult]):
    """
    Plot fitted curves and residuals for multiple lineshape models.

    Parameters
    ----------
    x : np.ndarray
        Full X data passed to fit_multiple_shapes.
    y : np.ndarray
        Full Y data passed to fit_multiple_shapes.
    results : dict
        Mapping of shape name to FitResult; only successful fits are drawn.
    """

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [3, 1]}
    )

    colors = ["#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

    # Data
    ax1.plot(x, y, "o", markersize=4, alpha=0.7, label="Data", color="#1f77b4")

    # Fits
    for i, (shape, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            x_plot = result.x_fit if result.x_fit is not None else x
            ax1.plot(
                x_plot,
                result.fitted_curve,
                "-",
                linewidth=2,
                label=f"{shape.title()} (R²={result.r_squared:.4f})",
                color=color,
            )

    ax1.set_xlabel("Magnetic Field / G")
    ax1.set_ylabel("Intensity")
    ax1.set_title("EPR Signal Fitting - Shape Comparison")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Residuals
    for i, (shape, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            x_plot = result.x_fit if result.x_fit is not None else x
            ax2.plot(
                x_plot,
                result.residuals,
                "o-",
                markersize=2,
                alpha=0.7,
                label=f"{shape.title()}",
                color=color,
            )

    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Magnetic Field / G")
    ax2.set_ylabel("Residuals")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
