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

    if signal.ndim not in (1, 2):
        raise ValueError(f"signal must be 1D or 2D, got {signal.ndim}D")

    if signal.ndim == 2 and axis not in (0, 1):
        raise ValueError(f"axis must be 0 or 1 for 2D signals, got {axis}")

    expected_len = signal.shape[axis] if signal.ndim == 2 else signal.shape[0]
    if len(time) != expected_len:
        raise ValueError(
            f"time length {len(time)} does not match signal length"
            f" {expected_len} along axis {axis}."
        )

    logger.info(f"remove_baseline: {method}, signal shape {signal.shape}")

    if signal.ndim == 1:
        corrected, baseline = _remove_baseline_1d(
            time, signal, method, order, end_fraction
        )
    else:
        corrected, baseline = _remove_baseline_2d(
            time, signal, method, order, end_fraction, axis
        )

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
    corrected = np.zeros_like(signal, dtype=float)
    baseline = np.zeros_like(signal, dtype=float)

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


# =============================================================================
# apodize
# =============================================================================


def apodize(
    signal: np.ndarray,
    window: str = "hann",
    alpha: Optional[float] = None,
    half_window: Optional[str] = None,
    axis: Union[int, str] = "both",
    plot: bool = False,
) -> np.ndarray:
    """
    Apply an apodization window to reduce spectral leakage from signal truncation.

    Parameters
    ----------
    signal : np.ndarray
        Signal array, 1D or 2D.
    window : str
        Window type: any key accepted by ``apowin()``: ``'hann'``,
        ``'hamming'``, ``'blackman'``, ``'kaiser'``, ``'gaussian'``,
        ``'exponential'``, etc.
    alpha : float, optional
        Shape parameter for ``kaiser`` and ``gaussian`` windows.
    half_window : str, optional
        ``None`` (symmetric window), ``'left'`` (first half only),
        or ``'right'`` (second half only).
    axis : int or str
        For 2D data: ``0``, ``1``, or ``'both'``. When ``'both'``, a 2D
        window is built as the outer product of two 1D windows, one per
        axis. Ignored for 1D input.
    plot : bool
        Show a figure with the original signal, window shape, and
        windowed signal.

    Returns
    -------
    signal_windowed : np.ndarray
        Windowed signal, same shape as ``signal``.

    Raises
    ------
    ValueError
        If ``signal`` is not 1D or 2D, or if ``axis`` is invalid for 2D.

    Examples
    --------
    >>> signal = np.exp(-t / 120) * np.sin(2 * np.pi * 8.5e-3 * t)
    >>> windowed = apodize(signal, window='hann', half_window='right')
    """
    signal = np.asarray(signal, dtype=float)

    if signal.ndim == 1:
        windowed = _apodize_1d(signal, window, alpha, half_window)
    elif signal.ndim == 2:
        windowed = _apodize_2d(signal, window, alpha, half_window, axis)
    else:
        raise ValueError(f"signal must be 1D or 2D, got {signal.ndim}D")

    if plot:
        _plot_apodize(signal, windowed, window, half_window)

    return windowed


def _apodize_1d(
    signal: np.ndarray,
    window: str,
    alpha: Optional[float],
    half_window: Optional[str],
) -> np.ndarray:
    n = len(signal)
    w = apowin(window, n, alpha=alpha, half_window=half_window)
    return signal * w


def _apodize_2d(
    signal: np.ndarray,
    window: str,
    alpha: Optional[float],
    half_window: Optional[str],
    axis: Union[int, str],
) -> np.ndarray:
    n_rows, n_cols = signal.shape

    if axis == "both":
        w_rows = apowin(window, n_rows, alpha=alpha, half_window=half_window)
        w_cols = apowin(window, n_cols, alpha=alpha, half_window=half_window)
        return signal * np.outer(w_rows, w_cols)
    elif axis == 0:
        w = apowin(window, n_rows, alpha=alpha, half_window=half_window)
        return signal * w[:, np.newaxis]
    elif axis == 1:
        w = apowin(window, n_cols, alpha=alpha, half_window=half_window)
        return signal * w[np.newaxis, :]
    else:
        raise ValueError(f"axis must be 0, 1, or 'both' for 2D signals, got {axis!r}")


def _plot_apodize(
    signal: np.ndarray,
    windowed: np.ndarray,
    window: str,
    half_window: Optional[str],
) -> None:
    if signal.ndim == 1:
        n = len(signal)
        try:
            w = apowin(window, n, half_window=half_window)
        except Exception:
            w = np.ones(n)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        axes[0].plot(signal, "b-", linewidth=1.5)
        axes[0].set_title("Original signal")

        axes[1].plot(w, color="orange", linewidth=2)
        title = f"{window} window"
        if half_window:
            title += f" ({half_window} half)"
        axes[1].set_title(title)
        axes[1].set_ylim(-0.05, 1.1)

        axes[2].plot(windowed, "g-", linewidth=1.5)
        axes[2].set_title("Apodized signal")

        for ax in axes:
            ax.set_xlabel("Sample")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        im0 = axes[0].imshow(signal, aspect="auto", cmap="RdBu_r", origin="lower")
        axes[0].set_title("Original signal")
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].imshow(windowed, aspect="auto", cmap="RdBu_r", origin="lower")
        axes[1].set_title("Apodized signal")
        plt.colorbar(im1, ax=axes[1])

        for ax in axes:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()


# =============================================================================
# zero_pad
# =============================================================================


def zero_pad(
    signal: np.ndarray,
    factor: Optional[int] = None,
    n_points: Optional[int] = None,
    axis: Union[int, str] = -1,
    plot: bool = False,
) -> np.ndarray:
    """
    Pad a signal with trailing zeros to increase FFT frequency resolution.

    Exactly one of ``factor`` or ``n_points`` must be supplied.

    Parameters
    ----------
    signal : np.ndarray
        Signal array, 1D or 2D.
    factor : int, optional
        Multiplicative factor: output length = ``factor * N`` along ``axis``.
    n_points : int, optional
        Absolute output length along ``axis``.
    axis : int or str
        Axis to pad. For 1D input this parameter is ignored.
        For 2D input: ``0``, ``1``, ``-1`` (equivalent to ``1``), or
        ``'both'`` (pad each axis independently with the same factor or
        n_points). Default ``-1``.
    plot : bool
        Show a before/after figure.

    Returns
    -------
    signal_padded : np.ndarray
        Zero-padded signal, same dtype as ``signal``.

    Raises
    ------
    ValueError
        If both or neither of ``factor`` / ``n_points`` are given, if
        ``n_points`` is shorter than the signal length, if ``signal`` is
        not 1D or 2D, or if ``axis`` is invalid for 2D.

    Examples
    --------
    >>> signal = np.ones(256)
    >>> padded = zero_pad(signal, factor=4)   # 1024 points
    >>> padded = zero_pad(signal, n_points=1024)
    """
    signal = np.asarray(signal, dtype=float)

    if (factor is None) == (n_points is None):
        raise ValueError(
            "Provide exactly one of factor or n_points, not both or neither."
        )

    if signal.ndim == 1:
        padded = _zero_pad_1d(signal, factor, n_points)
    elif signal.ndim == 2:
        padded = _zero_pad_2d(signal, factor, n_points, axis)
    else:
        raise ValueError(f"signal must be 1D or 2D, got {signal.ndim}D")

    if plot:
        _plot_zero_pad(signal, padded)

    return padded


def _zero_pad_1d(
    signal: np.ndarray,
    factor: Optional[int],
    n_points: Optional[int],
) -> np.ndarray:
    n = len(signal)
    n_out = (factor * n) if factor is not None else n_points
    if n_out < n:
        raise ValueError(f"Target length {n_out} is shorter than signal length {n}.")
    padded = np.zeros(n_out, dtype=signal.dtype)
    padded[:n] = signal
    return padded


def _zero_pad_2d(
    signal: np.ndarray,
    factor: Optional[int],
    n_points: Optional[int],
    axis: Union[int, str],
) -> np.ndarray:
    n_rows, n_cols = signal.shape

    def _out_len(n: int) -> int:
        result = (factor * n) if factor is not None else n_points
        if result < n:
            raise ValueError(
                f"Target length {result} is shorter than signal length {n}."
            )
        return result

    if axis == "both":
        n_rows_out = _out_len(n_rows)
        n_cols_out = _out_len(n_cols)
        padded = np.zeros((n_rows_out, n_cols_out), dtype=signal.dtype)
        padded[:n_rows, :n_cols] = signal
    elif axis in (0,):
        n_rows_out = _out_len(n_rows)
        padded = np.zeros((n_rows_out, n_cols), dtype=signal.dtype)
        padded[:n_rows, :] = signal
    elif axis in (1, -1):
        n_cols_out = _out_len(n_cols)
        padded = np.zeros((n_rows, n_cols_out), dtype=signal.dtype)
        padded[:, :n_cols] = signal
    else:
        raise ValueError(
            f"axis must be 0, 1, -1, or 'both' for 2D signals, got {axis!r}"
        )

    return padded


def _plot_zero_pad(signal: np.ndarray, padded: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    if signal.ndim == 1:
        axes[0].plot(signal, "b-", linewidth=1.5)
        axes[0].set_title(f"Original signal ({len(signal)} points)")
        axes[0].set_xlabel("Sample")

        axes[1].plot(padded, "g-", linewidth=1.5)
        axes[1].axvline(
            len(signal) - 1,
            color="red",
            linestyle="--",
            alpha=0.7,
            label="Original data end",
        )
        axes[1].set_title(f"Zero-padded signal ({len(padded)} points)")
        axes[1].set_xlabel("Sample")
        axes[1].legend()
    else:
        im0 = axes[0].imshow(signal, aspect="auto", cmap="RdBu_r", origin="lower")
        axes[0].set_title(f"Original signal {signal.shape}")
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].imshow(padded, aspect="auto", cmap="RdBu_r", origin="lower")
        axes[1].set_title(f"Zero-padded signal {padded.shape}")
        plt.colorbar(im1, ax=axes[1])

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()
