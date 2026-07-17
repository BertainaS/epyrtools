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
def pseudo_voigt(x, center, width, eta=0.5, derivative=0, phase=0.0):
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


def demo():
    """Demonstrate different lineshape combinations"""
    import matplotlib.pyplot as plt

    x = np.linspace(-15, 15, 1000)

    # Set1 colors
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Pure shapes
    ax = axes[0, 0]
    gauss = lshape(x, 0, 8, alpha=1.0)
    lorentz = lshape(x, 0, 8, alpha=0.0)

    ax.plot(x, gauss, color=colors[0], linewidth=2.5, label="Gaussian (α=1)")
    ax.plot(x, lorentz, color=colors[1], linewidth=2.5, label="Lorentzian (α=0)")
    ax.set_title("Pure Lineshapes", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Mixed shapes (pseudo-Voigt)
    ax = axes[0, 1]
    alphas = [0.25, 0.5, 0.75]
    for i, alpha in enumerate(alphas):
        mixed = lshape(x, 0, 8, alpha=alpha)
        ax.plot(x, mixed, color=colors[i], linewidth=2.5, label=f"α = {alpha}")

    ax.set_title("Pseudo-Voigt Profiles", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Different widths
    ax = axes[1, 0]
    narrow_gauss = lshape(x, -3, (4, 8), alpha=0.7)  # Narrow Gaussian + wide Lorentzian
    wide_gauss = lshape(x, 3, (12, 4), alpha=0.7)  # Wide Gaussian + narrow Lorentzian

    ax.plot(
        x,
        narrow_gauss,
        color=colors[0],
        linewidth=2.5,
        label="Narrow Gauss + Wide Lorentz",
    )
    ax.plot(
        x,
        wide_gauss,
        color=colors[1],
        linewidth=2.5,
        label="Wide Gauss + Narrow Lorentz",
    )
    ax.set_title("Different Component Widths", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    # Derivatives
    ax = axes[1, 1]
    center_func = lshape(x, 0, 8, derivative=0, alpha=0.5)
    first_deriv = lshape(x, 0, 8, derivative=1, alpha=0.5)
    second_deriv = lshape(x, 0, 8, derivative=2, alpha=0.5)

    ax.plot(x, center_func, color=colors[0], linewidth=2.5, label="Function")
    ax.plot(x, first_deriv, color=colors[1], linewidth=2.5, label="1st derivative")
    ax.plot(x, second_deriv, color=colors[2], linewidth=2.5, label="2nd derivative")
    ax.set_title("Derivatives (α=0.5)", fontweight="bold")
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
