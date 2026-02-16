"""
Gaussian lineshape functions
Modern, optimized implementation for magnetic resonance spectroscopy
"""

from typing import Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from scipy import special


def gaussian(x, center, width, derivative=0, phase=0.0, return_both=False):
    """
    Area-normalized Gaussian lineshape with derivatives and phase rotation.

    The Gaussian profile represents inhomogeneous broadening from
    statistical distributions of local fields or magnetic environments.

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
        Gaussian values, optionally with dispersion component

    Examples:
    ---------
    >>> x = np.linspace(-10, 10, 1000)
    >>> # Standard absorption Gaussian
    >>> y = gaussian(x, 0, 4)
    >>> # First derivative
    >>> dy = gaussian(x, 0, 4, derivative=1)
    >>> # Dispersion mode
    >>> disp = gaussian(x, 0, 4, phase=np.pi/2)
    >>> # Both absorption and dispersion
    >>> abs_part, disp_part = gaussian(x, 0, 4, return_both=True)
    """

    x = np.asarray(x, dtype=float)

    # Input validation
    _validate_gaussian_inputs(center, width, derivative, phase)

    # Convert FWHM to standard deviation
    sigma = width / (2 * np.sqrt(2 * np.log(2)))
    k = (x - center) / (sigma * np.sqrt(2))

    # Compute absorption and dispersion components
    abs_part, disp_part = _compute_gaussian_components(k, sigma, derivative)

    # Handle output based on phase and return options
    return _handle_gaussian_output(abs_part, disp_part, phase, return_both)


def _validate_gaussian_inputs(center, width, derivative, phase):
    """Validate Gaussian input parameters"""
    if not isinstance(center, (int, float)):
        raise ValueError("center must be a number")
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("width must be positive")
    if not isinstance(derivative, int) or derivative < -1:
        raise ValueError("derivative must be integer >= -1")
    if not isinstance(phase, (int, float)):
        raise ValueError("phase must be a real number")


def _compute_gaussian_components(k, sigma, derivative):
    """Compute absorption and dispersion components"""

    if derivative == -1:
        # Integral from -infinity (cumulative distribution)
        abs_part = 0.5 * (1 + special.erf(k))
        # Dispersion integral - simplified approximation
        disp_part = np.sqrt(2 / np.pi) / (2 * sigma) * k * np.exp(-(k**2))

    elif derivative == 0:
        # Standard Gaussian
        abs_part = np.exp(-(k**2)) / (sigma * np.sqrt(2 * np.pi))
        # Dispersion uses Dawson function
        disp_part = 2 * special.dawsn(k) / (sigma * np.sqrt(np.pi))

    elif derivative >= 1:
        # Derivatives using Hermite polynomials
        # Absorption: H_n(k) * exp(-k²) / (σ√(2π)) * (-1/σ)^n * 2^(-n/2)
        prefactor = (
            1
            / (sigma * np.sqrt(2 * np.pi))
            * (-1 / sigma) ** derivative
            * 2 ** (-derivative / 2)
        )
        hermite_poly = _hermite_polynomial(k, derivative)
        abs_part = prefactor * hermite_poly * np.exp(-(k**2))

        # Dispersion derivatives using Dawson function
        dawson_prefactor = (
            2
            / (sigma * np.sqrt(np.pi))
            * (-1 / sigma) ** derivative
            * 2 ** (-derivative / 2)
        )
        dawson_deriv = _dawson_derivative(k, derivative)
        disp_part = dawson_prefactor * dawson_deriv

    else:
        raise NotImplementedError(f"Derivative order {derivative} not implemented")

    return abs_part, disp_part


def _handle_gaussian_output(abs_part, disp_part, phase, return_both):
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


def _hermite_polynomial(x, n):
    """
    Physicists' Hermite polynomials H_n(x) using recurrence.
    H_0(x) = 1, H_1(x) = 2x, H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x)
    """
    if n == 0:
        return np.ones_like(x)
    elif n == 1:
        return 2 * x
    elif n == 2:
        return 4 * x**2 - 2
    elif n == 3:
        return 8 * x**3 - 12 * x
    elif n == 4:
        return 16 * x**4 - 48 * x**2 + 12
    else:
        # Use recurrence for higher orders
        H_prev = np.ones_like(x)  # H_0
        H_curr = 2 * x  # H_1

        for k in range(2, n + 1):
            H_next = 2 * x * H_curr - 2 * (k - 1) * H_prev
            H_prev, H_curr = H_curr, H_next

        return H_curr


def _dawson_derivative(x, n):
    """Dawson function derivatives using recurrence relations"""
    F = special.dawsn(x)  # Dawson function F(x)

    if n == 0:
        return F
    elif n == 1:
        return 1 - 2 * x * F
    elif n == 2:
        return -4 * x - 2 * F + 4 * x**2 * F
    else:
        # For higher derivatives, use approximate recurrence
        # This is simplified - full implementation needs G_n polynomials
        prev_deriv = _dawson_derivative(x, n - 1)
        return (
            -2 * n * prev_deriv + 2 * x * _dawson_derivative(x, n - 1) if n > 0 else F
        )


# Convenience functions for common cases
def gaussian_absorption(x, center, width):
    """Pure absorption Gaussian"""
    return gaussian(x, center, width)


def gaussian_dispersion(x, center, width):
    """Pure dispersion Gaussian"""
    return gaussian(x, center, width, phase=np.pi / 2)


def gaussian_derivative(x, center, width, order=1):
    """Gaussian derivatives"""
    return gaussian(x, center, width, derivative=order)


def demo():
    """Interactive demonstration of Gaussian lineshapes"""

    x = np.linspace(-15, 15, 1000)

    # Modern colors
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Different widths
    ax = axes[0, 0]
    widths = [2, 4, 8]
    for i, width in enumerate(widths):
        y = gaussian(x, 0, width)
        ax.plot(x, y, color=colors[i], linewidth=2.5, label=f"FWHM = {width}")

    ax.set_title("Different Widths", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Derivatives
    ax = axes[0, 1]
    derivs = [0, 1, 2]
    labels = ["Function", "1st derivative", "2nd derivative"]

    for i, (deriv, label) in enumerate(zip(derivs, labels)):
        y = gaussian(x, 0, 6, derivative=deriv)
        ax.plot(x, y, color=colors[i], linewidth=2.5, label=label)

    ax.set_title("Derivatives", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Absorption vs Dispersion
    ax = axes[1, 0]
    abs_part, disp_part = gaussian(x, 0, 6, return_both=True)

    ax.plot(x, abs_part, color=colors[0], linewidth=2.5, label="Absorption")
    ax.plot(x, disp_part, color=colors[1], linewidth=2.5, label="Dispersion")

    ax.set_title("Absorption vs Dispersion", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Comparison with Lorentzian
    ax = axes[1, 1]
    gauss = gaussian(x, 0, 6)
    # Import lorentzian if available
    try:
        from .lorentzian import lorentzian

        lorentz = lorentzian(x, 0, 6)
        ax.plot(
            x,
            lorentz,
            color=colors[1],
            linewidth=2.5,
            label="Lorentzian",
            linestyle="--",
        )
    except ImportError:
        # Create simple Lorentzian for comparison
        gamma = 3  # half-width
        u = x / gamma
        lorentz = (1 / np.pi) / gamma / (1 + u**2)
        ax.plot(
            x,
            lorentz,
            color=colors[1],
            linewidth=2.5,
            label="Lorentzian",
            linestyle="--",
        )

    ax.plot(x, gauss, color=colors[0], linewidth=2.5, label="Gaussian")
    ax.set_title("Gaussian vs Lorentzian", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Style all subplots
    for ax in axes.flat:
        ax.set_xlabel("Position")
        ax.set_ylabel("Intensity")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    demo()
