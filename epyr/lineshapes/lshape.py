"""
General lineshape function - combines Gaussian and Lorentzian shapes
Modern implementation with pseudo-Voigt capability
"""

from typing import Tuple, Union

import numpy as np

from ._validation import validate_abscissa


def lshape(
    x: np.ndarray,
    center: float,
    width: Union[float, Tuple[float, float]],
    derivative: int = 0,
    alpha: float = 1.0,
    phase: float = 0.0,
) -> np.ndarray:
    """
    General normalized lineshape function.

    Computes a linear combination of Gaussian and Lorentzian lineshapes:
    alpha * Gaussian + (1-alpha) * Lorentzian

    This creates pseudo-Voigt profiles commonly used in spectroscopy.

    Parameters:
    -----------
    x : array
        Abscissa points
    center : float
        Peak center position
    width : float or (float, float)
        Full width at half maximum
        If single value: same width for both components
        If tuple: (gaussian_width, lorentzian_width)
    derivative : int, default=0
        Derivative order (0=function, 1=first derivative, 2=second, -1=integral)
    alpha : float, default=1.0
        Mixing parameter (0=pure Lorentzian, 1=pure Gaussian)
    phase : float, default=0.0
        Phase rotation (0=absorption, π/2=dispersion)

    Returns:
    --------
    array
        Lineshape values

    Examples:
    ---------
    >>> x = np.linspace(-10, 10, 1000)
    >>> # Pure Gaussian
    >>> gauss = lshape(x, 0, 5, alpha=1.0)
    >>> # Pure Lorentzian
    >>> lorentz = lshape(x, 0, 5, alpha=0.0)
    >>> # 50/50 mix (pseudo-Voigt)
    >>> mixed = lshape(x, 0, 5, alpha=0.5)
    >>> # Different widths for each component
    >>> mixed_widths = lshape(x, 0, (3, 7), alpha=0.3)
    """

    x = validate_abscissa(x)

    # Validate inputs
    if not isinstance(center, (int, float)):
        raise ValueError("center must be a number")
    if not isinstance(alpha, (int, float)) or not 0 <= alpha <= 1:
        raise ValueError("alpha must be between 0 and 1")
    if not isinstance(derivative, int) or derivative < -1:
        raise ValueError("derivative must be integer >= -1")
    if not isinstance(phase, (int, float)):
        raise ValueError("phase must be a number")

    # Handle width parameter
    if isinstance(width, (list, tuple)):
        if len(width) != 2:
            raise ValueError("width tuple must have exactly 2 values")
        width_gauss, width_lorentz = width
        if width_gauss <= 0 or width_lorentz <= 0:
            raise ValueError("all widths must be positive")
    else:
        if width <= 0:
            raise ValueError("width must be positive")
        width_gauss = width_lorentz = width

    from .gaussian import gaussian as _gaussian
    from .lorentzian import lorentzian as _lorentzian

    result = np.zeros_like(x, dtype=float)

    # Compute Gaussian component using the canonical implementation
    if alpha > 0:
        result += alpha * _gaussian(
            x, center, width_gauss, derivative=derivative, phase=phase
        )

    # Compute Lorentzian component using the canonical implementation
    if alpha < 1:
        result += (1 - alpha) * _lorentzian(
            x, center, width_lorentz, derivative=derivative, phase=phase
        )

    return result


# Convenience functions for common cases
def pseudo_voigt(
    x: np.ndarray,
    center: float,
    width: float,
    eta: float = 0.5,
    derivative: int = 0,
    phase: float = 0.0,
) -> np.ndarray:
    """
    Pseudo-Voigt profile: η*Lorentzian + (1-η)*Gaussian

    Parameters:
    -----------
    x : array
        Abscissa points
    center : float
        Peak center position
    width : float
        Full width at half maximum
    eta : float, default=0.5
        Mixing parameter (0=Gaussian, 1=Lorentzian)
    derivative : int, default=0
        Derivative order (0=function, 1=first derivative, 2=second)
    phase : float, default=0.0
        Phase rotation (0=absorption, π/2=dispersion)

    Returns:
    --------
    array
        Pseudo-Voigt profile values

    Examples:
    ---------
    >>> x = np.linspace(-10, 10, 1000)
    >>> # Standard pseudo-Voigt
    >>> y = pseudo_voigt(x, 0, 5, eta=0.5)
    >>> # First derivative
    >>> dy = pseudo_voigt(x, 0, 5, eta=0.5, derivative=1)
    """
    return lshape(x, center, width, derivative=derivative, alpha=1 - eta, phase=phase)
