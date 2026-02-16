"""
Lorentzian lineshape functions
Modern, optimized implementation for EPR spectroscopy
"""

from typing import Tuple, Union

import matplotlib.pyplot as plt
import numpy as np


def lorentzian(x, center, width, derivative=0, phase=0.0, return_both=False):
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

    x = np.asarray(x, dtype=float)

    # Input validation
    _validate_lorentzian_inputs(center, width, derivative, phase)

    # Normalized variable: u = (x - center) / gamma
    gamma = width / 2  # Half-width at half-maximum
    u = (x - center) / gamma

    # Compute absorption and dispersion components
    abs_part, disp_part = _compute_lorentzian_components(u, gamma, derivative)

    # Handle output based on phase and return options
    return _handle_lorentzian_output(abs_part, disp_part, phase, return_both)


def _validate_lorentzian_inputs(center, width, derivative, phase):
    """Validate Lorentzian input parameters"""
    if not isinstance(center, (int, float)):
        raise ValueError("center must be a number")
    if not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("width must be positive")
    if not isinstance(derivative, int) or derivative < -1:
        raise ValueError("derivative must be integer >= -1")
    if not isinstance(phase, (int, float)):
        raise ValueError("phase must be a real number")


def _compute_lorentzian_components(u, gamma, derivative):
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


def _handle_lorentzian_output(abs_part, disp_part, phase, return_both):
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
def lorentzian_absorption(x, center, width):
    """Pure absorption Lorentzian"""
    return lorentzian(x, center, width)


def lorentzian_dispersion(x, center, width):
    """Pure dispersion Lorentzian"""
    return lorentzian(x, center, width, phase=np.pi / 2)


def lorentzian_derivative(x, center, width, order=1):
    """Lorentzian derivatives"""
    return lorentzian(x, center, width, derivative=order)


def demo():
    """Interactive demonstration of Lorentzian lineshapes"""

    x = np.linspace(-15, 15, 1000)

    # Modern colors
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Different widths
    ax = axes[0, 0]
    widths = [2, 4, 8]
    for i, width in enumerate(widths):
        y = lorentzian(x, 0, width)
        ax.plot(x, y, color=colors[i], linewidth=2.5, label=f"FWHM = {width}")

    ax.set_title("Different Widths", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Derivatives
    ax = axes[0, 1]
    derivs = [0, 1, 2]
    labels = ["Function", "1st derivative", "2nd derivative"]

    for i, (deriv, label) in enumerate(zip(derivs, labels)):
        y = lorentzian(x, 0, 6, derivative=deriv)
        ax.plot(x, y, color=colors[i], linewidth=2.5, label=label)

    ax.set_title("Derivatives", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Absorption vs Dispersion
    ax = axes[1, 0]
    abs_part, disp_part = lorentzian(x, 0, 6, return_both=True)

    ax.plot(x, abs_part, color=colors[0], linewidth=2.5, label="Absorption")
    ax.plot(x, disp_part, color=colors[1], linewidth=2.5, label="Dispersion")

    ax.set_title("Absorption vs Dispersion", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Phase rotation
    ax = axes[1, 1]
    phases = [0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2]

    for i, phase in enumerate(phases):
        y = lorentzian(x, 0, 6, phase=phase)
        label = f"φ = {phase:.2f}" if i < 4 else "φ = π/2"
        ax.plot(x, y, color=colors[i], linewidth=2, label=label)

    ax.set_title("Phase Rotation", fontweight="bold")
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
