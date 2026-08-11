"""
Apodization windows for signal processing
Modern implementation with additional window types and features
"""

import numpy as np
from scipy import special

try:
    from ..logging_config import get_logger
except ImportError:
    import logging

    def get_logger(name):
        return logging.getLogger(name)


logger = get_logger(__name__)


def apowin(window_type, n_points, alpha=None, half_window=None):
    """
    Generate apodization windows for signal processing.

    Apodization windows are used to reduce spectral leakage and
    improve signal-to-noise ratio in Fourier transform spectroscopy.

    Parameters:
    -----------
    window_type : str
        Window type:
        - 'hamming' or 'ham': Hamming window
        - 'hann' or 'han': Hann (Hanning) window
        - 'blackman' or 'bla': Blackman window
        - 'bartlett' or 'bar': Bartlett (triangular) window
        - 'connes' or 'con': Connes window
        - 'cosine' or 'cos': Cosine window
        - 'welch' or 'wel': Welch window
        - 'kaiser' or 'kai': Kaiser window (needs alpha)
        - 'gaussian' or 'gau': Gaussian window (needs alpha)
        - 'exponential' or 'exp': Exponential window (needs alpha)
    n_points : int
        Number of points in the window
    alpha : float, optional
        Shape parameter for Kaiser, Gaussian, and Exponential windows
    half_window : str, optional
        Generate half window: 'left' (-1 to 0) or 'right' (0 to 1)

    Returns:
    --------
    array
        Normalized window values (peak = 1)

    Examples:
    ---------
    >>> # Hamming window
    >>> w = apowin('hamming', 256)
    >>> # Kaiser window with beta=6
    >>> w_kaiser = apowin('kaiser', 256, alpha=6)
    >>> # Half Hann window (right side)
    >>> w_half = apowin('hann', 128, half_window='right')
    """

    # Input validation
    if not isinstance(n_points, int) or n_points <= 0:
        raise ValueError("n_points must be a positive integer")

    # Normalize window type
    window_type = window_type.lower()

    # Handle abbreviated forms
    window_aliases = {
        "ham": "hamming",
        "han": "hann",
        "bla": "blackman",
        "bar": "bartlett",
        "con": "connes",
        "cos": "cosine",
        "wel": "welch",
        "kai": "kaiser",
        "gau": "gaussian",
        "exp": "exponential",
    }

    window_type = window_aliases.get(window_type, window_type)

    # Set coordinate range
    if half_window == "right":
        x = np.linspace(0, 1, n_points)
    elif half_window == "left":
        x = np.linspace(-1, 0, n_points)
    else:
        x = np.linspace(-1, 1, n_points)

    # Generate window
    window = _generate_window(window_type, x, alpha)

    # Normalize to peak value of 1
    if np.max(window) > 0:
        window = window / np.max(window)

    return window


def _generate_window(window_type, x, alpha):
    """Generate the window function values"""

    if window_type == "hamming":
        return 0.54 + 0.46 * np.cos(np.pi * x)

    elif window_type == "hann":
        return 0.5 + 0.5 * np.cos(np.pi * x)

    elif window_type == "blackman":
        return 0.42 + 0.5 * np.cos(np.pi * x) + 0.08 * np.cos(2 * np.pi * x)

    elif window_type == "bartlett":
        return 1 - np.abs(x)

    elif window_type == "connes":
        return (1 - x**2) ** 2

    elif window_type == "cosine":
        return np.cos(np.pi * x / 2)

    elif window_type == "welch":
        return 1 - x**2

    elif window_type == "kaiser":
        if alpha is None:
            raise ValueError("Kaiser window requires alpha parameter (typically 3-9)")
        return special.i0(alpha * np.sqrt(1 - x**2)) / special.i0(alpha)

    elif window_type == "gaussian":
        if alpha is None:
            raise ValueError(
                "Gaussian window requires alpha parameter (typically 0.6-1.2)"
            )
        return np.exp(-2 * x**2 / alpha**2)

    elif window_type == "exponential":
        if alpha is None:
            raise ValueError(
                "Exponential window requires alpha parameter (typically 2-6)"
            )
        return np.exp(-alpha * np.abs(x))

    else:
        raise ValueError(f"Unknown window type: {window_type}")
