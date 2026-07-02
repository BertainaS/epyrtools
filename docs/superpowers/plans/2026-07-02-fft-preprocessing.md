# FFT Preprocessing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three standalone, chainable preprocessing functions (`remove_baseline`, `apodize`, `zero_pad`) to `epyr/signalprocessing/` for ESEEM, HYSCORE, and Rabi time-domain FFT pipelines.

**Architecture:** New file `epyr/signalprocessing/preprocessing.py` holds all three functions. Each function auto-detects 1D vs 2D input via `signal.ndim`. Private helpers implement the core logic; public functions handle dispatch, validation, and plotting. `__init__.py` is updated to export the new names.

**Tech Stack:** NumPy, SciPy (`optimize.curve_fit`), Matplotlib, existing `apowin()` from `epyr/signalprocessing/apowin.py`.

## Global Constraints

- Python 3.8+ compatible type hints only (`Union`, `Optional`, `Tuple` from `typing`)
- `dtype=float` on all `np.asarray()` calls to avoid integer-array issues
- All public functions follow NumPy docstring format
- `plot=False` by default on all functions; plots use `Agg` backend in tests via `matplotlib.use("Agg")` at top of test file
- Pytest marks: `@pytest.mark.smoke`, `@pytest.mark.standard`, `@pytest.mark.deep`
- No print statements; use `logger = get_logger(__name__)`
- Commit message format: one imperative sentence, no AI attribution

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `epyr/signalprocessing/preprocessing.py` | Create | All three public functions + private helpers + plot helpers |
| `epyr/signalprocessing/__init__.py` | Modify | Export `remove_baseline`, `apodize`, `zero_pad` |
| `tests/test_preprocessing.py` | Create | Full test suite (smoke / standard / deep) |

---

## Task 1: `remove_baseline`

**Files:**
- Create: `epyr/signalprocessing/preprocessing.py`
- Create: `tests/test_preprocessing.py`

**Interfaces:**
- Produces: `remove_baseline(time, signal, method, order, end_fraction, axis, plot) -> Tuple[np.ndarray, np.ndarray]`

---

- [ ] **Step 1: Create `tests/test_preprocessing.py` with failing tests for `remove_baseline`**

```python
"""Tests for epyr.signalprocessing.preprocessing."""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from epyr.signalprocessing.preprocessing import remove_baseline


# ---------------------------------------------------------------------------
# remove_baseline — 1D polynomial
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_remove_baseline_polynomial_1d_shape():
    time = np.linspace(0, 500, 256)
    drift = 0.5 - 0.001 * time
    signal = np.sin(2 * np.pi * 8.5e-3 * time) * np.exp(-time / 120) + drift
    corrected, baseline = remove_baseline(time, signal, method="polynomial", order=1)
    assert corrected.shape == signal.shape
    assert baseline.shape == signal.shape


@pytest.mark.smoke
def test_remove_baseline_polynomial_1d_reduces_end_mean():
    time = np.linspace(0, 500, 256)
    drift = 0.5 * np.ones(256)
    signal = np.sin(2 * np.pi * 8.5e-3 * time) * np.exp(-time / 80) + drift
    corrected, _ = remove_baseline(time, signal, method="polynomial", order=1, end_fraction=0.2)
    # Trailing mean should be near zero after correction
    assert np.abs(np.mean(corrected[-30:])) < np.abs(np.mean(signal[-30:]))


# ---------------------------------------------------------------------------
# remove_baseline — 1D exponential
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_remove_baseline_exponential_1d_shape():
    time = np.linspace(0, 500, 256)
    decay = 1.5 * np.exp(-time / 200)
    signal = np.sin(2 * np.pi * 8.5e-3 * time) * 0.1 + decay
    corrected, baseline = remove_baseline(time, signal, method="exponential")
    assert corrected.shape == signal.shape
    assert baseline.shape == signal.shape


@pytest.mark.smoke
def test_remove_baseline_exponential_1d_subtracts_decay():
    time = np.linspace(0, 500, 256)
    decay = 2.0 * np.exp(-time / 150)
    signal = decay.copy()  # Pure exponential, no oscillation
    corrected, baseline = remove_baseline(time, signal, method="exponential")
    # After subtracting exponential baseline from itself, residual should be near zero
    assert np.max(np.abs(corrected)) < 0.1


@pytest.mark.smoke
def test_remove_baseline_unknown_method_raises():
    time = np.linspace(0, 100, 64)
    signal = np.ones(64)
    with pytest.raises(ValueError, match="Unknown method"):
        remove_baseline(time, signal, method="spline")


@pytest.mark.smoke
def test_remove_baseline_3d_raises():
    time = np.linspace(0, 100, 64)
    signal = np.ones((4, 8, 64))
    with pytest.raises(ValueError, match="1D or 2D"):
        remove_baseline(time, signal)


# ---------------------------------------------------------------------------
# remove_baseline — 2D
# ---------------------------------------------------------------------------

@pytest.mark.standard
def test_remove_baseline_2d_shape():
    time = np.linspace(0, 500, 128)
    drift = 0.5 * np.ones(128)
    n_traces = 16
    signal = np.tile(
        np.sin(2 * np.pi * 8.5e-3 * time) * np.exp(-time / 100) + drift,
        (n_traces, 1)
    )
    corrected, baseline = remove_baseline(time, signal, method="polynomial", order=2, axis=1)
    assert corrected.shape == signal.shape
    assert baseline.shape == signal.shape


@pytest.mark.standard
def test_remove_baseline_2d_axis0_shape():
    time = np.linspace(0, 200, 64)
    n_cols = 20
    signal = np.random.randn(64, n_cols) + np.linspace(1, 0, 64)[:, np.newaxis]
    corrected, baseline = remove_baseline(time, signal, method="polynomial", order=1, axis=0)
    assert corrected.shape == signal.shape


@pytest.mark.standard
def test_remove_baseline_2d_invalid_axis_raises():
    time = np.linspace(0, 100, 64)
    signal = np.ones((32, 64))
    with pytest.raises(ValueError, match="axis must be 0 or 1"):
        remove_baseline(time, signal, axis=2)


# ---------------------------------------------------------------------------
# remove_baseline — plot (smoke-tests that plot does not raise)
# ---------------------------------------------------------------------------

@pytest.mark.standard
def test_remove_baseline_plot_1d_no_error():
    time = np.linspace(0, 500, 64)
    signal = np.sin(2 * np.pi * 8e-3 * time) + 0.3
    remove_baseline(time, signal, method="polynomial", order=1, plot=True)
    plt.close("all")


@pytest.mark.standard
def test_remove_baseline_plot_2d_no_error():
    time = np.linspace(0, 500, 64)
    signal = np.tile(np.sin(2 * np.pi * 8e-3 * time) + 0.3, (8, 1))
    remove_baseline(time, signal, method="polynomial", order=1, plot=True)
    plt.close("all")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/sylvainbertaina/Documents/Cloud_CNRS/GitHub/epyrtools
pytest tests/test_preprocessing.py -v 2>&1 | head -30
```

Expected: `ModuleNotFoundError` or `ImportError` — `preprocessing` does not exist yet.

- [ ] **Step 3: Create `epyr/signalprocessing/preprocessing.py` with `remove_baseline`**

```python
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
        corrected, baseline = _remove_baseline_1d(time, signal, method, order, end_fraction)
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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_preprocessing.py -v -m "smoke or standard"
```

Expected: all green. If exponential tests are slow, check `maxfev`; 10000 is intentional for robustness.

- [ ] **Step 5: Commit**

```bash
git add epyr/signalprocessing/preprocessing.py tests/test_preprocessing.py
git commit -m "Add remove_baseline preprocessing function with tests"
```

---

## Task 2: `apodize`

**Files:**
- Modify: `epyr/signalprocessing/preprocessing.py` — append `apodize` and helpers
- Modify: `tests/test_preprocessing.py` — append `apodize` tests

**Interfaces:**
- Consumes: `apowin(window_type, n_points, alpha=None, half_window=None)` from `epyr/signalprocessing/apowin.py`
- Produces: `apodize(signal, window, alpha, half_window, axis, plot) -> np.ndarray`

---

- [ ] **Step 1: Append failing tests for `apodize` to `tests/test_preprocessing.py`**

Add after the existing `remove_baseline` tests:

```python
from epyr.signalprocessing.preprocessing import apodize


# ---------------------------------------------------------------------------
# apodize — 1D
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_apodize_1d_hann_shape():
    signal = np.sin(2 * np.pi * np.linspace(0, 1, 256))
    result = apodize(signal, window="hann")
    assert result.shape == signal.shape


@pytest.mark.smoke
def test_apodize_1d_hann_endpoints_near_zero():
    signal = np.ones(256)
    result = apodize(signal, window="hann")
    # Hann window is 0 at both ends
    assert result[0] < 0.01
    assert result[-1] < 0.01


@pytest.mark.smoke
def test_apodize_1d_hann_center_near_one():
    signal = np.ones(256)
    result = apodize(signal, window="hann")
    # Hann window peaks at 1 in the centre
    assert result[128] > 0.99


@pytest.mark.standard
def test_apodize_1d_half_window_right():
    signal = np.ones(128)
    result = apodize(signal, window="hann", half_window="right")
    # Right half: starts near 1, ends near 0
    assert result[0] > 0.95
    assert result[-1] < 0.05


@pytest.mark.standard
def test_apodize_1d_half_window_left():
    signal = np.ones(128)
    result = apodize(signal, window="hann", half_window="left")
    # Left half: starts near 0, ends near 1
    assert result[0] < 0.05
    assert result[-1] > 0.95


@pytest.mark.standard
def test_apodize_1d_kaiser_alpha():
    signal = np.ones(128)
    result = apodize(signal, window="kaiser", alpha=6)
    assert result.shape == signal.shape
    assert result[0] < result[64]  # Peak in the middle


@pytest.mark.smoke
def test_apodize_3d_raises():
    signal = np.ones((4, 8, 64))
    with pytest.raises(ValueError, match="1D or 2D"):
        apodize(signal)


# ---------------------------------------------------------------------------
# apodize — 2D
# ---------------------------------------------------------------------------

@pytest.mark.standard
def test_apodize_2d_both_shape():
    signal = np.ones((64, 128))
    result = apodize(signal, window="hann", axis="both")
    assert result.shape == signal.shape


@pytest.mark.standard
def test_apodize_2d_both_corners_near_zero():
    signal = np.ones((64, 128))
    result = apodize(signal, window="hann", axis="both")
    assert result[0, 0] < 0.01
    assert result[-1, -1] < 0.01


@pytest.mark.standard
def test_apodize_2d_axis0_shape():
    signal = np.ones((64, 128))
    result = apodize(signal, window="hann", axis=0)
    assert result.shape == signal.shape


@pytest.mark.standard
def test_apodize_2d_axis1_shape():
    signal = np.ones((64, 128))
    result = apodize(signal, window="hann", axis=1)
    assert result.shape == signal.shape


@pytest.mark.standard
def test_apodize_2d_invalid_axis_raises():
    signal = np.ones((64, 128))
    with pytest.raises(ValueError, match="axis must be"):
        apodize(signal, window="hann", axis=2)


# ---------------------------------------------------------------------------
# apodize — plot (smoke-tests that plot does not raise)
# ---------------------------------------------------------------------------

@pytest.mark.standard
def test_apodize_plot_1d_no_error():
    signal = np.ones(64)
    apodize(signal, window="hann", plot=True)
    plt.close("all")


@pytest.mark.standard
def test_apodize_plot_2d_no_error():
    signal = np.ones((32, 64))
    apodize(signal, window="hann", axis="both", plot=True)
    plt.close("all")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_preprocessing.py -k "apodize" -v 2>&1 | head -20
```

Expected: `ImportError` — `apodize` not yet defined.

- [ ] **Step 3: Append `apodize` implementation to `epyr/signalprocessing/preprocessing.py`**

Add after the `_plot_baseline` function:

```python
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
        Window type: any key accepted by ``apowin()`` — ``'hann'``,
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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_preprocessing.py -k "apodize" -v
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add epyr/signalprocessing/preprocessing.py tests/test_preprocessing.py
git commit -m "Add apodize preprocessing function with tests"
```

---

## Task 3: `zero_pad`

**Files:**
- Modify: `epyr/signalprocessing/preprocessing.py` — append `zero_pad` and helpers
- Modify: `tests/test_preprocessing.py` — append `zero_pad` tests

**Interfaces:**
- Produces: `zero_pad(signal, factor, n_points, axis, plot) -> np.ndarray`

---

- [ ] **Step 1: Append failing tests for `zero_pad` to `tests/test_preprocessing.py`**

Add after the existing `apodize` tests:

```python
from epyr.signalprocessing.preprocessing import zero_pad


# ---------------------------------------------------------------------------
# zero_pad — 1D
# ---------------------------------------------------------------------------

@pytest.mark.smoke
def test_zero_pad_factor_1d_length():
    signal = np.ones(128)
    result = zero_pad(signal, factor=4)
    assert len(result) == 512


@pytest.mark.smoke
def test_zero_pad_n_points_1d_length():
    signal = np.ones(128)
    result = zero_pad(signal, n_points=512)
    assert len(result) == 512


@pytest.mark.smoke
def test_zero_pad_1d_preserves_original():
    signal = np.arange(64, dtype=float)
    result = zero_pad(signal, factor=2)
    np.testing.assert_array_equal(result[:64], signal)


@pytest.mark.smoke
def test_zero_pad_1d_trailing_zeros():
    signal = np.ones(64)
    result = zero_pad(signal, factor=2)
    np.testing.assert_array_equal(result[64:], np.zeros(64))


@pytest.mark.smoke
def test_zero_pad_both_factor_and_n_points_raises():
    signal = np.ones(64)
    with pytest.raises(ValueError, match="exactly one"):
        zero_pad(signal, factor=2, n_points=256)


@pytest.mark.smoke
def test_zero_pad_neither_factor_nor_n_points_raises():
    signal = np.ones(64)
    with pytest.raises(ValueError, match="exactly one"):
        zero_pad(signal)


@pytest.mark.smoke
def test_zero_pad_n_points_shorter_than_signal_raises():
    signal = np.ones(128)
    with pytest.raises(ValueError, match="shorter than signal"):
        zero_pad(signal, n_points=64)


@pytest.mark.smoke
def test_zero_pad_3d_raises():
    signal = np.ones((4, 8, 64))
    with pytest.raises(ValueError, match="1D or 2D"):
        zero_pad(signal, factor=2)


# ---------------------------------------------------------------------------
# zero_pad — 2D
# ---------------------------------------------------------------------------

@pytest.mark.standard
def test_zero_pad_2d_axis1_shape():
    signal = np.ones((32, 64))
    result = zero_pad(signal, factor=2, axis=1)
    assert result.shape == (32, 128)


@pytest.mark.standard
def test_zero_pad_2d_axis0_shape():
    signal = np.ones((32, 64))
    result = zero_pad(signal, factor=2, axis=0)
    assert result.shape == (64, 64)


@pytest.mark.standard
def test_zero_pad_2d_both_shape():
    signal = np.ones((32, 64))
    result = zero_pad(signal, factor=2, axis="both")
    assert result.shape == (64, 128)


@pytest.mark.standard
def test_zero_pad_2d_both_n_points_shape():
    signal = np.ones((32, 64))
    result = zero_pad(signal, n_points=128, axis="both")
    assert result.shape == (128, 128)


@pytest.mark.standard
def test_zero_pad_2d_preserves_data():
    signal = np.random.randn(32, 64)
    result = zero_pad(signal, factor=2, axis=1)
    np.testing.assert_array_equal(result[:, :64], signal)
    np.testing.assert_array_equal(result[:, 64:], np.zeros((32, 64)))


@pytest.mark.standard
def test_zero_pad_2d_invalid_axis_raises():
    signal = np.ones((32, 64))
    with pytest.raises(ValueError, match="axis must be"):
        zero_pad(signal, factor=2, axis=3)


# ---------------------------------------------------------------------------
# zero_pad — plot
# ---------------------------------------------------------------------------

@pytest.mark.standard
def test_zero_pad_plot_1d_no_error():
    signal = np.ones(64)
    zero_pad(signal, factor=2, plot=True)
    plt.close("all")


@pytest.mark.standard
def test_zero_pad_plot_2d_no_error():
    signal = np.ones((32, 64))
    zero_pad(signal, factor=2, axis="both", plot=True)
    plt.close("all")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_preprocessing.py -k "zero_pad" -v 2>&1 | head -20
```

Expected: `ImportError` — `zero_pad` not yet defined.

- [ ] **Step 3: Append `zero_pad` implementation to `epyr/signalprocessing/preprocessing.py`**

Add after `_plot_apodize`:

```python
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
        raise ValueError(
            f"Target length {n_out} is shorter than signal length {n}."
        )
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
        return (factor * n) if factor is not None else n_points

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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_preprocessing.py -k "zero_pad" -v
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add epyr/signalprocessing/preprocessing.py tests/test_preprocessing.py
git commit -m "Add zero_pad preprocessing function with tests"
```

---

## Task 4: Export + deep integration test

**Files:**
- Modify: `epyr/signalprocessing/__init__.py` — add imports and `__all__` entries
- Modify: `tests/test_preprocessing.py` — append deep chaining test

**Interfaces:**
- Consumes: `remove_baseline`, `apodize`, `zero_pad` from Task 1–3
- Consumes: `analyze_frequencies` from `epyr.signalprocessing` (already exported)

---

- [ ] **Step 1: Append deep chaining test to `tests/test_preprocessing.py`**

Add at the end of the file:

```python
from epyr.signalprocessing import remove_baseline as rb, apodize as ap, zero_pad as zp
from epyr.signalprocessing import analyze_frequencies


# ---------------------------------------------------------------------------
# Deep: full ESEEM preprocessing chain
# ---------------------------------------------------------------------------

@pytest.mark.deep
def test_chain_eseem_pipeline():
    """Verify the full preprocessing chain produces a valid FFT result."""
    rng = np.random.default_rng(42)
    time = np.linspace(0, 500, 256)  # 500 ns, ns units
    rabi_freq = 8.5  # MHz
    decay = 1.5 * np.exp(-time / 200)
    oscillation = 0.3 * np.sin(2 * np.pi * rabi_freq * time * 1e-3)
    signal = oscillation + decay + 0.02 * rng.standard_normal(256)

    # Step 1: remove exponential baseline
    corrected, baseline = rb(time, signal, method="exponential")
    assert corrected.shape == (256,)
    # Baseline should track the slow decay
    assert np.mean(np.abs(baseline - decay)) < 0.2

    # Step 2: apodize with right half-Hann (signal starts high, decays)
    windowed = ap(corrected, window="hann", half_window="right")
    assert windowed.shape == (256,)
    assert windowed[-1] < windowed[0]  # right half: start > end

    # Step 3: zero-pad ×4
    padded = zp(windowed, factor=4)
    assert len(padded) == 1024
    np.testing.assert_array_equal(padded[256:], np.zeros(768))

    # Step 4: FFT — skip window and DC removal since already done
    result = analyze_frequencies(time, padded[:256], window=None, remove_dc=False, plot=False)
    # Dominant frequency should be near rabi_freq
    assert len(result["dominant_frequencies"]) > 0
    detected = result["dominant_frequencies"][0]
    assert abs(detected - rabi_freq) < 1.0  # within 1 MHz
```

- [ ] **Step 2: Run the deep test to verify it fails**

```bash
pytest tests/test_preprocessing.py -k "chain" -v 2>&1 | head -20
```

Expected: `ImportError` — `rb`, `ap`, `zp` not yet in `epyr.signalprocessing`.

- [ ] **Step 3: Update `epyr/signalprocessing/__init__.py`**

Read the current `__init__.py` first, then add the following imports after the existing `from .frequency_analysis import (...)` block:

```python
# Import preprocessing utilities
from .preprocessing import (
    apodize,
    remove_baseline,
    zero_pad,
)
```

And extend `__all__` to include:

```python
    # Preprocessing (standalone chainable functions)
    "remove_baseline",
    "apodize",
    "zero_pad",
```

- [ ] **Step 4: Run the full test suite**

```bash
pytest tests/test_preprocessing.py -v
```

Expected: all tests green across smoke / standard / deep.

Also run the smoke suite to check for regressions:

```bash
pytest -m smoke -v
```

Expected: all green, no regressions.

- [ ] **Step 5: Run code quality checks**

```bash
make format && make quality
```

Fix any issues reported before committing.

- [ ] **Step 6: Commit**

```bash
git add epyr/signalprocessing/__init__.py tests/test_preprocessing.py
git commit -m "Export remove_baseline, apodize, zero_pad from signalprocessing"
```

---

## Self-Review

**Spec coverage check:**

| Spec requirement | Covered by |
|---|---|
| `remove_baseline` polynomial 1D | Task 1 |
| `remove_baseline` exponential 1D | Task 1 |
| `remove_baseline` 2D (axis 0 and 1) | Task 1 |
| `remove_baseline` plot=True | Task 1 |
| `apodize` 1D full and half-window | Task 2 |
| `apodize` 2D both axes / axis=0 / axis=1 | Task 2 |
| `apodize` plot=True | Task 2 |
| `zero_pad` factor | Task 3 |
| `zero_pad` n_points | Task 3 |
| `zero_pad` both/neither raises | Task 3 |
| `zero_pad` 2D axis 0/1/both | Task 3 |
| `zero_pad` plot=True | Task 3 |
| Export from `__init__.py` | Task 4 |
| Chaining integration test | Task 4 |

All spec requirements are covered. No placeholders remain.
