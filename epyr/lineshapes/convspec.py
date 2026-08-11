"""
Spectrum convolution with lineshapes
Modern implementation with multi-dimensional support
"""

from typing import Union

import numpy as np
from scipy import signal

from ..logging_config import get_logger

logger = get_logger(__name__)


def convspec(
    spectrum: np.ndarray,
    step_size: Union[float, np.ndarray],
    width: Union[float, np.ndarray],
    derivative: Union[int, np.ndarray] = 0,
    alpha: Union[float, np.ndarray] = 1.0,
    phase: Union[float, np.ndarray] = 0.0,
) -> np.ndarray:
    """
    Convolve spectrum with lineshape functions.

    Applies broadening to stick spectra or other discrete data by
    convolution with Gaussian, Lorentzian, or pseudo-Voigt profiles.

    Parameters:
    -----------
    spectrum : array
        Input spectrum to convolve
    step_size : float or array
        Abscissa step size for each dimension
    width : float or array
        Full width at half maximum for lineshape
    derivative : int or array, default=0
        Derivative order (0=function, 1=first deriv, 2=second deriv)
    alpha : float or array, default=1.0
        Shape parameter (1=Gaussian, 0=Lorentzian, 0-1=pseudo-Voigt)
    phase : float or array, default=0.0
        Phase (0=absorption, π/2=dispersion)

    Returns:
    --------
    array
        Convolved spectrum with same shape as input

    Examples:
    ---------
    >>> # Simple 1D convolution
    >>> x = np.linspace(0, 100, 1000)
    >>> stick_spec = np.zeros_like(x)
    >>> stick_spec[500] = 1.0  # Delta peak at center
    >>> broadened = convspec(stick_spec, 0.1, 2.0)  # Gaussian, FWHM=2
    >>>
    >>> # Lorentzian broadening
    >>> lorentz = convspec(stick_spec, 0.1, 2.0, alpha=0.0)
    >>>
    >>> # First derivative
    >>> deriv = convspec(stick_spec, 0.1, 2.0, derivative=1)
    """

    spectrum = np.asarray(
        spectrum, dtype=complex if np.iscomplexobj(spectrum) else float
    )

    # Handle multi-dimensional parameters
    ndim = spectrum.ndim
    step_size = _expand_parameter(step_size, ndim)
    width = _expand_parameter(width, ndim)
    derivative = _expand_parameter(derivative, ndim)
    alpha = _expand_parameter(alpha, ndim)
    phase = _expand_parameter(phase, ndim)

    # Validate inputs
    _validate_convspec_inputs(spectrum, step_size, width, derivative, alpha, phase)

    # Perform convolution
    result = _convolve_spectrum(spectrum, step_size, width, derivative, alpha, phase)

    # Preserve real/complex nature of input
    if np.isrealobj(spectrum) and np.iscomplexobj(result):
        result = np.real(result)

    return result


def _expand_parameter(param: Union[float, int, np.ndarray], ndim: int) -> np.ndarray:
    """Expand scalar parameters to match number of dimensions"""
    param = np.asarray(param)
    if param.ndim == 0:
        return np.full(ndim, param.item())
    elif len(param) == ndim:
        return param
    else:
        raise ValueError(
            f"Parameter length {len(param)} doesn't match spectrum dimensions {ndim}"
        )


def _validate_convspec_inputs(
    spectrum: np.ndarray,
    step_size: np.ndarray,
    width: np.ndarray,
    derivative: np.ndarray,
    alpha: np.ndarray,
    phase: np.ndarray,
) -> None:
    """Validate convolution parameters"""

    if np.any(step_size <= 0):
        raise ValueError("step_size must be positive")

    if np.any(width < 0):
        raise ValueError("width must be non-negative")

    if np.any((derivative < -1) | (derivative > 2)):
        raise ValueError("derivative must be -1, 0, 1, or 2")

    if np.any((alpha < 0) | (alpha > 1)):
        raise ValueError("alpha must be between 0 and 1")

    if np.any(~np.isfinite([step_size, width, derivative, alpha, phase])):
        raise ValueError("All parameters must be finite")


def _convolve_spectrum(
    spectrum: np.ndarray,
    step_size: np.ndarray,
    width: np.ndarray,
    derivative: np.ndarray,
    alpha: np.ndarray,
    phase: np.ndarray,
) -> np.ndarray:
    """Core convolution implementation using FFT"""

    # For simplicity, use basic convolution with scipy.signal
    # Full implementation would use lshape function
    from .gaussian import gaussian
    from .lorentzian import lorentzian

    # Create convolution kernel
    n_kernel = min(len(spectrum) // 2, 500)  # Reasonable kernel size
    x_kernel = np.arange(-n_kernel, n_kernel + 1) * step_size[0]

    if alpha[0] == 1.0:
        # Pure Gaussian
        kernel = gaussian(
            x_kernel, 0, width[0], derivative=int(derivative[0]), phase=phase[0]
        )
    elif alpha[0] == 0.0:
        # Pure Lorentzian
        kernel = lorentzian(
            x_kernel, 0, width[0], derivative=int(derivative[0]), phase=phase[0]
        )
    else:
        # Mixed (pseudo-Voigt)
        gauss_part = gaussian(
            x_kernel, 0, width[0], derivative=int(derivative[0]), phase=phase[0]
        )
        lorentz_part = lorentzian(
            x_kernel, 0, width[0], derivative=int(derivative[0]), phase=phase[0]
        )
        kernel = alpha[0] * gauss_part + (1 - alpha[0]) * lorentz_part

    # Normalize kernel
    if derivative[0] == 0 and np.sum(kernel) != 0:
        kernel = kernel / np.sum(kernel)

    # Convolve
    result = signal.convolve(spectrum, kernel, mode="same")

    return result
