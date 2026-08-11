"""
Lorentzian lineshape functions
Modern, optimized implementation for EPR spectroscopy
"""

from typing import Tuple, Union

import numpy as np

from ._validation import validate_abscissa


def lorentzian(
    x: np.ndarray,
    center: float,
    width: float,
    derivative: int = 0,
    phase: float = 0.0,
    return_both: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Area-normalized Lorentzian lineshape with derivatives and phase rotation.

    The Lorentzian profile is fundamental in magnetic resonance, representing
    homogeneous broadening from finite lifetimes and collision processes.

    Parameters:
    -----------
    x : array
        Abscissa points
    center : float
        Peak center position
    width : float
        Full width at half maximum (FWHM)
    derivative : int, default=0
        Derivative order:
        - 0: Standard lineshape
        - 1: First derivative
        - 2: Second derivative
        - -1: Integral from -∞
    phase : float, default=0.0
        Phase rotation in radians
        - 0: Pure absorption
        - π/2: Pure dispersion
    return_both : bool, default=False
        If True, return (absorption, dispersion) tuple

    Returns:
    --------
    array or tuple
        Lorentzian values, optionally with dispersion component

    Examples:
    ---------
    >>> x = np.linspace(-10, 10, 1000)
    >>> # Standard absorption Lorentzian
    >>> y = lorentzian(x, 0, 4)
    >>> # First derivative
    >>> dy = lorentzian(x, 0, 4, derivative=1)
    >>> # Dispersion mode
    >>> disp = lorentzian(x, 0, 4, phase=np.pi/2)
    >>> # Both absorption and dispersion
    >>> abs_part, disp_part = lorentzian(x, 0, 4, return_both=True)
    """

    x = validate_abscissa(x)

    # Input validation
    _validate_lorentzian_inputs(center, width, derivative, phase)

    # Normalized variable: u = (x - center) / gamma
    gamma = width / 2  # Half-width at half-maximum
    u = (x - center) / gamma

    # Compute absorption and dispersion components
    abs_part, disp_part = _compute_lorentzian_components(u, gamma, derivative)

    # Handle output based on phase and return options
    return _handle_lorentzian_output(abs_part, disp_part, phase, return_both)


def _validate_lorentzian_inputs(
    center: float, width: float, derivative: int, phase: float
) -> None:
    """Validate Lorentzian input parameters"""
    if not isinstance(center, (int, float)):
        raise ValueError("center must be a number")
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("width must be positive")
    if not isinstance(derivative, int) or derivative < -1:
        raise ValueError("derivative must be integer >= -1")
    if not isinstance(phase, (int, float)):
        raise ValueError("phase must be a real number")


def _compute_lorentzian_components(
    u: np.ndarray, gamma: float, derivative: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute absorption and dispersion components"""

    if derivative == -1:
        # Integral from -infinity
        abs_part = (1 / np.pi) * (np.arctan(u) + np.pi / 2)
        disp_part = (1 / np.pi) * np.log(1 + u**2) / 2

    elif derivative == 0:
        # Standard Lorentzian
        denominator = 1 + u**2
        abs_part = (1 / np.pi) / gamma / denominator
        disp_part = (1 / np.pi) / gamma * u / denominator

    elif derivative == 1:
        # First derivative
        denominator = (1 + u**2) ** 2
        abs_part = -(2 / np.pi) / gamma**2 * u / denominator
        disp_part = (1 / np.pi) / gamma**2 * (1 - u**2) / denominator

    elif derivative == 2:
        # Second derivative
        denominator = (1 + u**2) ** 3
        abs_part = (2 / np.pi) / gamma**3 * (3 * u**2 - 1) / denominator
        disp_part = -(4 / np.pi) / gamma**3 * u * (u**2 - 3) / denominator

    else:
        raise NotImplementedError(f"Derivative order {derivative} not implemented")

    return abs_part, disp_part


def _handle_lorentzian_output(
    abs_part: np.ndarray, disp_part: np.ndarray, phase: float, return_both: bool
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Handle output formatting based on phase and return options"""

    # Check if phase rotation is needed
    needs_phase_rotation = np.mod(phase, 2 * np.pi) != 0

    if needs_phase_rotation:
        # Apply phase rotation
        cos_p, sin_p = np.cos(phase), np.sin(phase)
        rotated_abs = cos_p * abs_part + sin_p * disp_part
        rotated_disp = -sin_p * abs_part + cos_p * disp_part

        if return_both:
            return rotated_abs, rotated_disp
        else:
            return rotated_abs
    else:
        # No phase rotation
        if return_both:
            return abs_part, disp_part
        else:
            return abs_part


# Convenience functions for common cases
def lorentzian_absorption(x: np.ndarray, center: float, width: float) -> np.ndarray:
    """Pure absorption Lorentzian"""
    return lorentzian(x, center, width)


def lorentzian_dispersion(x: np.ndarray, center: float, width: float) -> np.ndarray:
    """Pure dispersion Lorentzian"""
    return lorentzian(x, center, width, phase=np.pi / 2)


def lorentzian_derivative(
    x: np.ndarray, center: float, width: float, order: int = 1
) -> np.ndarray:
    """Lorentzian derivatives"""
    return lorentzian(x, center, width, derivative=order)
