"""
EPR Signal Fitting Module

Provides comprehensive fitting capabilities for EPR signals using various lineshape functions.
Integrates with the eprload function and lineshapes module for complete analysis workflow.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from typing import Optional, Dict, Any, Tuple, Union, List
from dataclasses import dataclass
import warnings

from .gaussian import gaussian
from .lorentzian import lorentzian
from .voigtian import voigtian
from .lshape import pseudo_voigt
from .lineshape_class import Lineshape


@dataclass
class FitResult:
    """Container for fitting results"""
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

    def summary(self) -> str:
        """Generate a summary string of the fit results"""
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


def fit_epr_signal(x_data: np.ndarray,
                   y_data: np.ndarray,
                   shape_type: str = 'gaussian',
                   initial_params: Optional[Dict[str, float]] = None,
                   bounds: Optional[Dict[str, Tuple[float, float]]] = None,
                   derivative: int = 0,
                   fit_phase: bool = False,
                   plot: bool = True,
                   **fit_kwargs) -> FitResult:
    """
    Fit EPR signal with specified lineshape function.

    Parameters:
    -----------
    x_data : array
        X-axis data (magnetic field or frequency)
    y_data : array
        Y-axis data (EPR signal intensity)
    shape_type : str, default='gaussian'
        Type of lineshape to fit:
        - 'gaussian': Gaussian lineshape
        - 'lorentzian': Lorentzian lineshape
        - 'voigt': Voigt profile (requires 2 width parameters)
        - 'pseudo_voigt': Pseudo-Voigt profile
    initial_params : dict, optional
        Initial parameter guesses. If None, auto-estimated.
        Expected keys depend on shape_type and fit_phase:
        - Basic: {'center', 'width', 'amplitude'}
        - With phase: add 'phase'
        - Voigt: {'center', 'gaussian_width', 'lorentzian_width', 'amplitude'}
        - Pseudo-Voigt: {'center', 'width', 'amplitude', 'alpha'}
    bounds : dict, optional
        Parameter bounds as {param: (min, max)}
    derivative : int, default=0
        Derivative order to use (0, 1, 2). This is a FIXED parameter, not fitted.
    fit_phase : bool, default=False
        Whether to fit the phase parameter (absorption/dispersion mixing)
    plot : bool, default=True
        Whether to create a plot of the results
    **fit_kwargs : dict
        Additional keywords passed to scipy.optimize.curve_fit

    Returns:
    --------
    FitResult
        Object containing fit parameters, statistics, and curves

    Examples:
    ---------
    >>> # Load EPR data
    >>> from epyr import eprload
    >>> x, y, params, filepath = eprload('data.DTA')
    >>>
    >>> # Fit with Gaussian
    >>> result = fit_epr_signal(x, y, 'gaussian')
    >>> print(result.summary())
    >>>
    >>> # Fit 1st derivative with phase adjustment
    >>> result = fit_epr_signal(x, y, 'gaussian', derivative=1, fit_phase=True)
    >>>
    >>> # Fit with custom initial parameters including phase
    >>> initial = {'center': 3500, 'width': 10, 'amplitude': 1000, 'phase': 0.1}
    >>> result = fit_epr_signal(x, y, 'lorentzian', initial_params=initial, fit_phase=True)
    >>>
    >>> # Fit 2nd derivative with bounds
    >>> bounds = {'center': (3400, 3600), 'width': (5, 50), 'phase': (0, 3.14)}
    >>> result = fit_epr_signal(x, y, 'gaussian', derivative=2, fit_phase=True, bounds=bounds)
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

    # Validate derivative parameter
    if not isinstance(derivative, int) or derivative < 0 or derivative > 2:
        raise ValueError("derivative must be 0, 1, or 2")

    # Get fitting function and parameter info
    fit_func, param_names, param_bounds = _get_fit_function(shape_type, derivative, fit_phase)

    # Estimate initial parameters if not provided
    if initial_params is None:
        initial_params = _estimate_initial_params(x_clean, y_clean, shape_type, fit_phase)

    # Validate and complete initial parameters
    initial_params = _validate_initial_params(initial_params, param_names, x_clean, y_clean)

    # Setup parameter bounds
    lower_bounds, upper_bounds = _setup_bounds(param_names, bounds, param_bounds,
                                               initial_params, x_clean, y_clean)

    # Prepare parameters for fitting
    p0 = [initial_params[name] for name in param_names]
    bounds_tuple = (lower_bounds, upper_bounds)

    # Set default fitting options
    default_kwargs = {'maxfev': 10000, 'method': 'trf'}
    default_kwargs.update(fit_kwargs)

    try:
        # Perform the fit
        popt, pcov = curve_fit(fit_func, x_clean, y_clean,
                              p0=p0, bounds=bounds_tuple, **default_kwargs)

        # Calculate fitted curve and residuals
        y_fitted = fit_func(x_clean, *popt)
        residuals = y_clean - y_fitted

        # Calculate statistics
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_clean - np.mean(y_clean))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        chi_squared = ss_res / (len(y_clean) - len(popt))

        # Extract parameter values and errors
        param_dict = {name: value for name, value in zip(param_names, popt)}
        param_errors = {}

        if pcov is not None:
            param_std_errors = np.sqrt(np.diag(pcov))
            param_errors = {name: error for name, error in zip(param_names, param_std_errors)}

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
            covariance_matrix=pcov
        )

        # Create plot if requested
        if plot:
            _plot_fit_results(x_clean, y_clean, result, shape_type, derivative, fit_phase)

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
            message=str(e)
        )


def _get_fit_function(shape_type: str, derivative: int, fit_phase: bool) -> Tuple[callable, List[str], Dict[str, Tuple[float, float]]]:
    """Get the appropriate fitting function and parameter information"""

    if shape_type == 'gaussian':
        if fit_phase:
            def fit_func(x, center, width, amplitude, phase):
                return amplitude * gaussian(x, center, width, derivative=derivative, phase=phase)
            param_names = ['center', 'width', 'amplitude', 'phase']
        else:
            def fit_func(x, center, width, amplitude):
                return amplitude * gaussian(x, center, width, derivative=derivative)
            param_names = ['center', 'width', 'amplitude']

        param_bounds = {
            'center': (-np.inf, np.inf),
            'width': (0.001, np.inf),
            'amplitude': (-np.inf, np.inf),
            'phase': (-np.pi, np.pi)
        }

    elif shape_type == 'lorentzian':
        if fit_phase:
            def fit_func(x, center, width, amplitude, phase):
                return amplitude * lorentzian(x, center, width, derivative=derivative, phase=phase)
            param_names = ['center', 'width', 'amplitude', 'phase']
        else:
            def fit_func(x, center, width, amplitude):
                return amplitude * lorentzian(x, center, width, derivative=derivative)
            param_names = ['center', 'width', 'amplitude']

        param_bounds = {
            'center': (-np.inf, np.inf),
            'width': (0.001, np.inf),
            'amplitude': (-np.inf, np.inf),
            'phase': (-np.pi, np.pi)
        }

    elif shape_type == 'voigt':
        if fit_phase:
            def fit_func(x, center, gaussian_width, lorentzian_width, amplitude, phase):
                return amplitude * voigtian(x, center, (gaussian_width, lorentzian_width),
                                          derivative=derivative, phase=phase)
            param_names = ['center', 'gaussian_width', 'lorentzian_width', 'amplitude', 'phase']
        else:
            def fit_func(x, center, gaussian_width, lorentzian_width, amplitude):
                return amplitude * voigtian(x, center, (gaussian_width, lorentzian_width),
                                          derivative=derivative)
            param_names = ['center', 'gaussian_width', 'lorentzian_width', 'amplitude']

        param_bounds = {
            'center': (-np.inf, np.inf),
            'gaussian_width': (0.001, np.inf),
            'lorentzian_width': (0.001, np.inf),
            'amplitude': (-np.inf, np.inf),
            'phase': (-np.pi, np.pi)
        }

    elif shape_type == 'pseudo_voigt':
        if derivative != 0:
            warnings.warn("pseudo_voigt does not support derivatives. Using derivative=0.")

        if fit_phase:
            def fit_func(x, center, width, amplitude, alpha, phase):
                # Pseudo-voigt doesn't have native phase support, approximate with cos/sin mix
                abs_part = amplitude * pseudo_voigt(x, center, width, eta=alpha)
                disp_part = amplitude * np.gradient(pseudo_voigt(x, center, width, eta=alpha), x[1]-x[0])
                return np.cos(phase) * abs_part + np.sin(phase) * disp_part
            param_names = ['center', 'width', 'amplitude', 'alpha', 'phase']
        else:
            def fit_func(x, center, width, amplitude, alpha):
                return amplitude * pseudo_voigt(x, center, width, eta=alpha)
            param_names = ['center', 'width', 'amplitude', 'alpha']

        param_bounds = {
            'center': (-np.inf, np.inf),
            'width': (0.001, np.inf),
            'amplitude': (-np.inf, np.inf),
            'alpha': (0.0, 1.0),
            'phase': (-np.pi, np.pi)
        }

    else:
        raise ValueError(f"Unsupported shape_type: {shape_type}. "
                        "Choose from: 'gaussian', 'lorentzian', 'voigt', 'pseudo_voigt'")

    return fit_func, param_names, param_bounds


def _estimate_initial_params(x: np.ndarray, y: np.ndarray, shape_type: str, fit_phase: bool = False) -> Dict[str, float]:
    """Estimate initial parameters from data"""

    # Basic estimates
    amplitude = np.max(y) - np.min(y)
    if amplitude == 0:
        amplitude = 1.0

    # Handle negative peaks
    if np.abs(np.min(y)) > np.abs(np.max(y)):
        amplitude = -amplitude
        peak_idx = np.argmin(y)
    else:
        peak_idx = np.argmax(y)

    center = x[peak_idx]

    # Estimate width from full width at half maximum
    if amplitude > 0:
        half_max = np.min(y) + amplitude / 2
        above_half = y >= half_max
    else:
        half_max = np.max(y) + amplitude / 2  # amplitude is negative
        above_half = y <= half_max

    if np.sum(above_half) > 1:
        indices = np.where(above_half)[0]
        width = x[indices[-1]] - x[indices[0]]
        if width <= 0:
            width = (x[-1] - x[0]) / 10
    else:
        width = (x[-1] - x[0]) / 10

    # Shape-specific parameters
    initial_params = {
        'center': center,
        'amplitude': amplitude
    }

    if shape_type in ['gaussian', 'lorentzian', 'pseudo_voigt']:
        initial_params['width'] = max(width, np.diff(x).mean() * 2)

    if shape_type == 'voigt':
        initial_params['gaussian_width'] = max(width * 0.6, np.diff(x).mean() * 2)
        initial_params['lorentzian_width'] = max(width * 0.6, np.diff(x).mean() * 2)

    if shape_type == 'pseudo_voigt':
        initial_params['alpha'] = 0.5  # 50/50 mix

    # Add phase parameter if fitting phase
    if fit_phase:
        # Estimate phase from the sign and shape of the data
        # Simple heuristic: if first derivative of peak area is positive -> negative phase
        if len(y) > 2:
            peak_region = np.abs(x - center) < width
            if np.sum(peak_region) > 2:
                y_peak = y[peak_region]
                x_peak = x[peak_region]
                # Rough estimate based on asymmetry
                if len(y_peak) > 2:
                    gradient_estimate = np.gradient(y_peak, x_peak)
                    avg_gradient = np.mean(gradient_estimate)
                    # Phase estimation: if derivative-like -> π/2, if absorption-like -> 0
                    if np.abs(avg_gradient) > np.abs(np.mean(y_peak)) * 0.1:
                        initial_params['phase'] = np.pi/4  # Mixed
                    else:
                        initial_params['phase'] = 0.0  # Pure absorption
                else:
                    initial_params['phase'] = 0.0
            else:
                initial_params['phase'] = 0.0
        else:
            initial_params['phase'] = 0.0

    return initial_params


def _validate_initial_params(initial_params: Dict[str, float],
                           param_names: List[str],
                           x: np.ndarray,
                           y: np.ndarray) -> Dict[str, float]:
    """Validate and complete initial parameters"""

    validated = initial_params.copy()

    # Ensure all required parameters are present
    for name in param_names:
        if name not in validated:
            if name == 'center':
                validated[name] = x[np.argmax(np.abs(y))]
            elif name == 'amplitude':
                validated[name] = np.max(y) - np.min(y)
            elif 'width' in name:
                validated[name] = (x[-1] - x[0]) / 10
            elif name == 'alpha':
                validated[name] = 0.5
            elif name == 'phase':
                validated[name] = 0.0
            else:
                validated[name] = 1.0

    return validated


def _setup_bounds(param_names: List[str],
                 user_bounds: Optional[Dict[str, Tuple[float, float]]],
                 default_bounds: Dict[str, Tuple[float, float]],
                 initial_params: Dict[str, float],
                 x: np.ndarray,
                 y: np.ndarray) -> Tuple[List[float], List[float]]:
    """Setup parameter bounds for fitting"""

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
        if name == 'center':
            if lower == -np.inf:
                lower = x.min() - (x.max() - x.min())
            if upper == np.inf:
                upper = x.max() + (x.max() - x.min())
        elif name == 'amplitude':
            if lower == -np.inf:
                lower = -10 * (np.max(y) - np.min(y))
            if upper == np.inf:
                upper = 10 * (np.max(y) - np.min(y))
        elif 'width' in name:
            if upper == np.inf:
                upper = x.max() - x.min()

        # Ensure initial value is within bounds
        init_val = initial_params[name]
        if init_val <= lower:
            lower = init_val * 0.1 if init_val > 0 else init_val * 10
        if init_val >= upper:
            upper = init_val * 10 if init_val > 0 else init_val * 0.1

        lower_bounds.append(lower)
        upper_bounds.append(upper)

    return lower_bounds, upper_bounds


def _format_significant_figures(value: float, n_sig: int) -> str:
    """
    Format a number to n significant figures.

    Parameters:
    -----------
    value : float
        Number to format
    n_sig : int
        Number of significant figures

    Returns:
    --------
    str
        Formatted number string
    """
    if value == 0:
        return "0"

    # Handle negative values
    if value < 0:
        return "-" + _format_significant_figures(-value, n_sig)

    # Find the order of magnitude
    import math
    order = math.floor(math.log10(abs(value)))

    # Scale the number to have the first significant digit in the ones place
    scaled = value / (10 ** order)

    # Round to n_sig decimal places
    rounded = round(scaled, n_sig - 1)

    # Scale back
    result = rounded * (10 ** order)

    # Format based on the order of magnitude
    if order >= 4 or order < -3:
        # Use scientific notation for very large or very small numbers
        return f"{result:.{n_sig-1}e}"
    elif order >= 0:
        # Use regular notation, with appropriate decimal places
        decimal_places = max(0, n_sig - 1 - order)
        return f"{result:.{decimal_places}f}"
    else:
        # For numbers < 1, show enough decimal places
        decimal_places = n_sig - 1 - order
        return f"{result:.{decimal_places}f}"


def _plot_fit_results(x: np.ndarray, y: np.ndarray, result: FitResult, shape_type: str, derivative: int = 0, fit_phase: bool = False):
    """Create a plot showing the fit results"""

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.plot(x, y, 'o', markersize=4, alpha=0.7, label='Data', color='#1f77b4')
    ax1.plot(x, result.fitted_curve, '-', linewidth=2, label=f'{shape_type.title()} fit',
             color='#d62728')

    ax1.set_xlabel('Magnetic Field / G')
    ax1.set_ylabel('Intensity')

    # Build title with derivative and phase info
    title_parts = [f'EPR Signal Fitting - {shape_type.title()}']
    if derivative > 0:
        title_parts.append(f'(d^{derivative})')
    if fit_phase:
        phase_val = result.parameters.get('phase', 0)
        title_parts.append(f'Phase: {phase_val:.3f} rad')

    ax1.set_title(' '.join(title_parts))
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Build fit results text with parameters and errors
    results_lines = [f'R² = {result.r_squared:.4f}', f'χ² = {result.chi_squared:.2e}']

    # Add fitted parameters with 4 significant figures and errors with 2 significant figures
    for param, value in result.parameters.items():
        # Format value to 4 significant figures
        value_str = _format_significant_figures(value, 4)

        if param in result.parameter_errors and result.parameter_errors[param] is not None:
            error = result.parameter_errors[param]
            # Format error to 2 significant figures
            error_str = _format_significant_figures(error, 2)
            results_lines.append(f'{param}: {value_str} ± {error_str}')
        else:
            results_lines.append(f'{param}: {value_str}')

    textstr = '\n'.join(results_lines)
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='left', bbox=props)

    # Residuals plot
    ax2.plot(x, result.residuals, 'o-', markersize=3, alpha=0.7, color='#ff7f0e')
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Magnetic Field / G')
    ax2.set_ylabel('Residuals')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def fit_multiple_shapes(x_data: np.ndarray,
                       y_data: np.ndarray,
                       shapes: List[str] = None,
                       derivative: int = 0,
                       fit_phase: bool = False,
                       plot: bool = True) -> Dict[str, FitResult]:
    """
    Fit EPR signal with multiple lineshape types and compare results.

    Parameters:
    -----------
    x_data, y_data : arrays
        EPR signal data
    shapes : list, optional
        List of shapes to try. Default: ['gaussian', 'lorentzian', 'pseudo_voigt']
    derivative : int, default=0
        Derivative order to use (0, 1, 2). Fixed parameter.
    fit_phase : bool, default=False
        Whether to fit the phase parameter
    plot : bool
        Whether to create comparison plot

    Returns:
    --------
    dict
        Dictionary of {shape_type: FitResult} for all attempted fits
    """

    if shapes is None:
        shapes = ['gaussian', 'lorentzian', 'pseudo_voigt']

    results = {}

    for shape in shapes:
        try:
            result = fit_epr_signal(x_data, y_data, shape, derivative=derivative,
                                  fit_phase=fit_phase, plot=False)
            results[shape] = result
        except Exception as e:
            print(f"Failed to fit {shape}: {e}")
            results[shape] = FitResult(
                shape_type=shape,
                parameters={}, parameter_errors={},
                fitted_curve=np.array([]), residuals=np.array([]),
                r_squared=0.0, chi_squared=np.inf,
                success=False, message=str(e)
            )

    # Find best fit based on R²
    successful_fits = {k: v for k, v in results.items() if v.success}

    if successful_fits and plot:
        _plot_comparison(x_data, y_data, successful_fits)

    # Print comparison
    print("=== Fit Comparison ===")
    for shape, result in results.items():
        if result.success:
            print(f"{shape:12s}: R² = {result.r_squared:.6f}, χ² = {result.chi_squared:.2e}")
        else:
            print(f"{shape:12s}: FAILED - {result.message}")

    if successful_fits:
        best_shape = max(successful_fits.keys(),
                        key=lambda k: successful_fits[k].r_squared)
        print(f"\nBest fit: {best_shape} (R² = {successful_fits[best_shape].r_squared:.6f})")

    return results


def _plot_comparison(x: np.ndarray, y: np.ndarray, results: Dict[str, FitResult]):
    """Plot comparison of different fits"""

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10),
                                   gridspec_kw={'height_ratios': [3, 1]})

    colors = ['#d62728', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b']

    # Data
    ax1.plot(x, y, 'o', markersize=4, alpha=0.7, label='Data', color='#1f77b4')

    # Fits
    for i, (shape, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            ax1.plot(x, result.fitted_curve, '-', linewidth=2,
                    label=f'{shape.title()} (R²={result.r_squared:.4f})',
                    color=color)

    ax1.set_xlabel('Magnetic Field / G')
    ax1.set_ylabel('Intensity')
    ax1.set_title('EPR Signal Fitting - Shape Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Residuals
    for i, (shape, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            ax2.plot(x, result.residuals, 'o-', markersize=2, alpha=0.7,
                    label=f'{shape.title()}', color=color)

    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Magnetic Field / G')
    ax2.set_ylabel('Residuals')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()