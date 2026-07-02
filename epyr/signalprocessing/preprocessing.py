"""
Preprocessing utilities for time-domain EPR signal analysis.

Standalone, chainable functions for ESEEM, HYSCORE, and Rabi pipelines:
remove_baseline, apodize, zero_pad.
"""

from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

try:
    from ..logging_config import get_logger
except ImportError:
    import logging

    def get_logger(name):
        return logging.getLogger(name)


try:
    from .apowin import apowin
except ImportError:
    from apowin import apowin

logger = get_logger(__name__)


# =============================================================================
# remove_baseline
# =============================================================================


def remove_baseline(
    time: np.ndarray,
    signal: np.ndarray,
    method: str = "polynomial",
    order: int = 3,
    end_fraction: float = 0.15,
    axis: int = 1,
    plot: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit and subtract a slowly-varying background from a time-domain EPR signal.

    Parameters
    ----------
    time : np.ndarray
        Time axis, 1D array, any unit. For 2D data, must match the length
        of the signal along ``axis``.
    signal : np.ndarray
        Signal array, 1D or 2D.
    method : str
        ``'polynomial'`` or ``'exponential'``.
    order : int
        Polynomial degree; ignored when ``method='exponential'``.
    end_fraction : float
        Fraction of trailing points used to anchor the polynomial fit.
        For exponential fits, the full time axis is used.
    axis : int
        For 2D data: axis along which baseline is fitted and removed
        (0 = column-wise, 1 = row-wise). Default 1.
    plot : bool
        Show a before/after figure.

    Returns
    -------
    signal_corrected : np.ndarray
        Signal with baseline subtracted, same shape as ``signal``.
    baseline : np.ndarray
        Fitted baseline array, same shape as ``signal``.

    Raises
    ------
    ValueError
        If ``method`` is not recognised, if ``signal`` is not 1D or 2D,
        or if ``axis`` is invalid for 2D data.
    RuntimeError
        If the exponential fit fails to converge.

    Examples
    --------
    >>> time = np.linspace(0, 500, 256)
    >>> signal = np.exp(-time / 200) + 0.05 * np.random.randn(256)
    >>> corrected, baseline = remove_baseline(time, signal, method='exponential')
    """
    time = np.asarray(time, dtype=float)
    signal = np.asarray(signal, dtype=float)

    if signal.ndim == 1:
        corrected, baseline = _remove_baseline_1d(
            time, signal, method, order, end_fraction
        )
    elif signal.ndim == 2:
        corrected, baseline = _remove_baseline_2d(
            time, signal, method, order, end_fraction, axis
        )
    else:
        raise ValueError(f"signal must be 1D or 2D, got {signal.ndim}D")

    if plot:
        _plot_baseline(time, signal, corrected, baseline, method)

    return corrected, baseline


def _remove_baseline_1d(
    time: np.ndarray,
    signal: np.ndarray,
    method: str,
    order: int,
    end_fraction: float,
) -> Tuple[np.ndarray, np.ndarray]:
    n = len(signal)
    n_end = max(2, round(end_fraction * n))

    if method == "polynomial":
        t_end = time[-n_end:]
        s_end = signal[-n_end:]
        coeffs = np.polyfit(t_end, s_end, order)
        baseline = np.polyval(coeffs, time)

    elif method == "exponential":

        def _exp(t, A, tau, C):
            return A * np.exp(-(t - t[0]) / tau) + C

        A0 = float(signal[0] - signal[-1])
        tau0 = float((time[-1] - time[0]) / 3.0)
        C0 = float(np.mean(signal[-n_end:]))

        try:
            popt, _ = curve_fit(
                _exp,
                time,
                signal,
                p0=[A0, max(tau0, 1e-10), C0],
                bounds=([-np.inf, 1e-10, -np.inf], [np.inf, np.inf, np.inf]),
                maxfev=10000,
            )
            baseline = _exp(time, *popt)
        except RuntimeError as exc:
            raise RuntimeError(
                f"Exponential baseline fit failed: {exc}. "
                "Try method='polynomial' or adjust end_fraction."
            ) from exc

    else:
        raise ValueError(
            f"Unknown method: {method!r}. Use 'polynomial' or 'exponential'."
        )

    return signal - baseline, baseline


def _remove_baseline_2d(
    time: np.ndarray,
    signal: np.ndarray,
    method: str,
    order: int,
    end_fraction: float,
    axis: int,
) -> Tuple[np.ndarray, np.ndarray]:
    corrected = np.zeros_like(signal)
    baseline = np.zeros_like(signal)

    if axis == 1:
        for i in range(signal.shape[0]):
            corrected[i], baseline[i] = _remove_baseline_1d(
                time, signal[i], method, order, end_fraction
            )
    elif axis == 0:
        for j in range(signal.shape[1]):
            c, b = _remove_baseline_1d(time, signal[:, j], method, order, end_fraction)
            corrected[:, j] = c
            baseline[:, j] = b
    else:
        raise ValueError(f"axis must be 0 or 1 for 2D signals, got {axis}")

    return corrected, baseline


def _plot_baseline(
    time: np.ndarray,
    signal: np.ndarray,
    corrected: np.ndarray,
    baseline: np.ndarray,
    method: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    if signal.ndim == 1:
        axes[0].plot(time, signal, "b-", linewidth=1.5, label="Original")
        axes[0].plot(time, baseline, "r--", linewidth=2, label=f"{method} baseline")
        axes[0].set_xlabel("Time")
        axes[0].set_ylabel("Signal")
        axes[0].set_title("Signal and fitted baseline")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(time, corrected, "g-", linewidth=1.5)
        axes[1].set_xlabel("Time")
        axes[1].set_ylabel("Signal")
        axes[1].set_title("Baseline-corrected signal")
        axes[1].grid(True, alpha=0.3)
    else:
        im0 = axes[0].imshow(signal, aspect="auto", cmap="RdBu_r", origin="lower")
        axes[0].set_title("Original signal")
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].imshow(corrected, aspect="auto", cmap="RdBu_r", origin="lower")
        axes[1].set_title("Baseline-corrected signal")
        plt.colorbar(im1, ax=axes[1])

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()
