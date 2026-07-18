"""
Pseudomodulation of EPR spectra.

Simulates the signal a lock-in amplifier would record under sinusoidal
field modulation, applied numerically to a spectrum recorded without field
modulation.

References:
    - Hyde, J.S., Pasenkiewicz-Gierula, M., Jesmanowicz, A., Antholine, W.E.,
      Appl. Magn. Reson. 1, 483-496 (1990).
"""

import numpy as np
from scipy import special

from ._validation import validate_abscissa


def pseudo_modulation(x, y, mod_amplitude, harmonic=1, pad=True):
    """
    Compute the pseudomodulated spectrum for simulated lock-in detection.

    Reproduces the signal a lock-in amplifier would record under sinusoidal
    field modulation of peak-to-peak amplitude `mod_amplitude`, applied to a
    spectrum recorded without field modulation.

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss (G) or millitesla (mT), uniformly
        spaced (ascending or descending).
    y : np.ndarray
        Spectrum values at each point in `x`, same length as `x`. Real or
        complex.
    mod_amplitude : float
        Peak-to-peak field modulation amplitude, in the same unit as `x`.
        Must be positive.
    harmonic : int, default=1
        Detection harmonic: 1 for first-harmonic (first-derivative-like)
        detection, 2 for second-harmonic detection.
    pad : bool, default=True
        Edge-pad `y` before the FFT convolution to suppress circular
        wraparound artifacts at the spectrum boundaries, then crop the
        result back to the original length.

    Returns
    -------
    np.ndarray
        Pseudomodulated spectrum, same length as `y`. Real if `y` is real,
        complex if `y` is complex.

    Raises
    ------
    ValueError
        If `x` and `y` have different lengths, `x` is not uniformly spaced,
        `mod_amplitude` is not positive, or `harmonic` is not 1 or 2.

    Notes
    -----
    For a sinusoidal modulation B(t) = B0 + (mod_amplitude/2) sin(theta),
    the signal detected at the n-th harmonic of the modulation frequency is

        S_n(B) = IFFT[ FFT(A)(k) * i^n * J_n(k * mod_amplitude / 2) ]

    where A is the input spectrum, k = 2*pi*fftfreq(N, dx) is the spatial
    frequency conjugate to the field axis, and J_n is the Bessel function of
    the first kind of order n. The i^n phase factor is required: J_n(-z) =
    (-1)^n J_n(z), so J_n(k * mod_amplitude/2) alone is odd in k for odd n,
    which breaks the Hermitian symmetry a real spectrum's FFT has and makes
    the IFFT purely imaginary. Multiplying by i^n restores Hermitian
    symmetry for every n, so the IFFT is real whenever A is real. In the
    small-amplitude limit, J_1(z) ~ z/2, so S_1(B) approaches
    (mod_amplitude/4) * dA/dB.

    Examples
    --------
    >>> x = np.linspace(-50, 50, 2000)
    >>> y = np.exp(-x**2 / (2 * 5.0**2))
    >>> y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)
    """
    x = validate_abscissa(x)
    y = np.asarray(y)
    _validate_pseudo_modulation_inputs(x, y, mod_amplitude, harmonic)

    dx = x[1] - x[0]
    is_complex_input = np.iscomplexobj(y)

    if pad:
        pad_width = len(y) // 2
        y_work = np.pad(y, pad_width, mode="edge")
    else:
        pad_width = 0
        y_work = y

    n = len(y_work)
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    transfer = (1j**harmonic) * special.jv(harmonic, k * mod_amplitude / 2)

    spectrum_fft = np.fft.fft(y_work)
    result = np.fft.ifft(spectrum_fft * transfer)

    if not is_complex_input:
        result = np.real(result)

    if pad_width:
        result = result[pad_width : pad_width + len(y)]

    return result


def _validate_pseudo_modulation_inputs(x, y, mod_amplitude, harmonic):
    """Validate pseudo_modulation input parameters."""
    if y.ndim != 1:
        raise ValueError("y must be a 1D array")
    if x.shape != y.shape:
        raise ValueError(f"Shape mismatch: x {x.shape} vs y {y.shape}")
    if len(x) < 2:
        raise ValueError("x must contain at least two points")

    steps = np.diff(x)
    if not np.allclose(steps, steps[0], rtol=1e-6):
        raise ValueError("x must be uniformly spaced")

    if mod_amplitude <= 0:
        raise ValueError("mod_amplitude must be positive")

    if harmonic not in (1, 2):
        raise ValueError("harmonic must be 1 or 2")
