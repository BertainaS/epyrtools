"""
Voigtian lineshape functions
Modern implementation of the Voigt profile (convolution of Gaussian and Lorentzian)
"""

from typing import Optional, Tuple, Union

import numpy as np
from scipy import special

from ._validation import validate_abscissa


def voigtian(
    x: np.ndarray,
    center: float,
    widths: Tuple[float, float],
    derivative: int = 0,
    phase: float = 0.0,
    return_both: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Area-normalized Voigt profile - convolution of Gaussian and Lorentzian.

    The Voigt profile models combined inhomogeneous (Gaussian) and
    homogeneous (Lorentzian) broadening mechanisms commonly found in
    magnetic resonance and optical spectroscopy.

    Parameters:
    -----------
    x : array
        Abscissa points
    center : float
        Peak center position
    widths : tuple of two floats
        (gaussian_fwhm, lorentzian_fwhm) - Full widths at half maximum
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
        Voigt profile values, optionally with dispersion component

    Examples:
    ---------
    >>> x = np.linspace(-10, 10, 1000)
    >>> # Equal Gaussian and Lorentzian widths
    >>> y = voigtian(x, 0, (4, 4))
    >>> # Gaussian-dominated profile
    >>> y_gauss = voigtian(x, 0, (6, 2))
    >>> # Lorentzian-dominated profile
    >>> y_lorentz = voigtian(x, 0, (2, 6))
    """

    x = validate_abscissa(x)

    # Input validation
    _validate_voigtian_inputs(center, widths, derivative, phase)

    gaussian_width, lorentzian_width = widths

    # Compute Voigt profile using Faddeeva function
    abs_part, disp_part = _voigt_faddeeva(
        x, center, gaussian_width, lorentzian_width, derivative
    )

    # Handle output based on phase and return options
    return _handle_voigt_output(abs_part, disp_part, phase, return_both)


def validate_voigt_widths(widths: Tuple[float, float]) -> None:
    """Check that widths is a valid (gaussian_width, lorentzian_width) pair.

    Both widths are FWHM values in the abscissa unit; they must be
    non-negative and at least one must be positive.

    Raises
    ------
    ValueError
        If widths is not a two-element tuple or violates the sign rules.
    """
    if not (isinstance(widths, (list, tuple)) and len(widths) == 2):
        raise ValueError(
            "widths must be a tuple of two values (gaussian_width, lorentzian_width)"
        )

    gaussian_width, lorentzian_width = widths
    if gaussian_width < 0 or lorentzian_width < 0:
        raise ValueError("both widths must be non-negative")

    if gaussian_width == 0 and lorentzian_width == 0:
        raise ValueError("at least one width must be positive")


def _validate_voigtian_inputs(
    center: float, widths: Tuple[float, float], derivative: int, phase: float
) -> None:
    """Validate Voigtian input parameters"""
    if not isinstance(center, (int, float)):
        raise ValueError("center must be a number")

    validate_voigt_widths(widths)

    if not isinstance(derivative, int) or derivative < -1:
        raise ValueError("derivative must be integer >= -1")

    if not isinstance(phase, (int, float)):
        raise ValueError("phase must be a real number")


def _voigt_faddeeva(
    x: np.ndarray,
    center: float,
    gauss_width: float,
    lorentz_width: float,
    derivative: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Voigt profile using Faddeeva function (complex error function)"""

    if gauss_width == 0:
        # Pure Lorentzian
        from .lorentzian import _compute_lorentzian_components

        gamma = lorentz_width / 2
        u = (x - center) / gamma
        return _compute_lorentzian_components(u, gamma, derivative)

    if lorentz_width == 0:
        # Pure Gaussian
        from .gaussian import _compute_gaussian_components

        sigma = gauss_width / (2 * np.sqrt(2 * np.log(2)))
        k = (x - center) / (sigma * np.sqrt(2))
        return _compute_gaussian_components(k, sigma, derivative)

    # True Voigt profile using Faddeeva function
    # Convert widths to standard parameters
    sigma = gauss_width / (2 * np.sqrt(2 * np.log(2)))  # Gaussian standard deviation
    gamma = lorentz_width / 2  # Lorentzian half-width

    # Normalized coordinates
    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))

    # Faddeeva function and normalization (used by all derivative orders except -1)
    w = special.wofz(z)
    norm = 1.0 / (sigma * np.sqrt(2 * np.pi))

    if derivative == 0:
        # Standard Voigt profile: V(x) = Re[w(z)] / (σ√(2π))
        abs_part = np.real(w) * norm
        disp_part = -np.imag(w) * norm

    elif derivative == -1:
        # Integral - approximate using individual components
        from .gaussian import gaussian
        from .lorentzian import lorentzian

        gauss_int = gaussian(x, center, gauss_width, derivative=-1)
        lorentz_int = lorentzian(x, center, lorentz_width, derivative=-1)
        weight = gauss_width / (gauss_width + lorentz_width)
        abs_part = weight * gauss_int + (1 - weight) * lorentz_int
        disp_part = np.zeros_like(x)

    elif derivative == 1:
        # Analytical first derivative using Faddeeva identity:
        # w'(z) = -2z·w(z) + 2i/√π
        dw = -2 * z * w + 2j / np.sqrt(np.pi)
        dzdx = 1.0 / (sigma * np.sqrt(2))
        result = dw * dzdx * norm
        abs_part = np.real(result)
        disp_part = -np.imag(result)

    elif derivative == 2:
        # Analytical second derivative using Faddeeva identities:
        # w'(z) = -2z·w(z) + 2i/√π
        # w''(z) = -2·w(z) + 4z²·w(z) - 4iz/√π
        d2w = -2 * w + 4 * z**2 * w - 4j * z / np.sqrt(np.pi)
        dzdx = 1.0 / (sigma * np.sqrt(2))
        result = d2w * dzdx**2 * norm
        abs_part = np.real(result)
        disp_part = -np.imag(result)

    else:
        raise NotImplementedError(
            f"Voigt derivative order {derivative} not implemented"
        )

    return abs_part, disp_part


def _handle_voigt_output(
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


# Convenience functions
def voigt_equal_widths(x: np.ndarray, center: float, width: float) -> np.ndarray:
    """Voigt profile with equal Gaussian and Lorentzian widths"""
    return voigtian(x, center, (width, width))


def voigt_gaussian_dominated(
    x: np.ndarray,
    center: float,
    gauss_width: float,
    lorentz_width: Optional[float] = None,
) -> np.ndarray:
    """Voigt profile dominated by Gaussian broadening"""
    if lorentz_width is None:
        lorentz_width = gauss_width * 0.3
    return voigtian(x, center, (gauss_width, lorentz_width))


def voigt_lorentzian_dominated(
    x: np.ndarray,
    center: float,
    lorentz_width: float,
    gauss_width: Optional[float] = None,
) -> np.ndarray:
    """Voigt profile dominated by Lorentzian broadening"""
    if gauss_width is None:
        gauss_width = lorentz_width * 0.3
    return voigtian(x, center, (gauss_width, lorentz_width))
