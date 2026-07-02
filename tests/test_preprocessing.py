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
    corrected, _ = remove_baseline(
        time, signal, method="polynomial", order=1, end_fraction=0.2
    )
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
        np.sin(2 * np.pi * 8.5e-3 * time) * np.exp(-time / 100) + drift, (n_traces, 1)
    )
    corrected, baseline = remove_baseline(
        time, signal, method="polynomial", order=2, axis=1
    )
    assert corrected.shape == signal.shape
    assert baseline.shape == signal.shape


@pytest.mark.standard
def test_remove_baseline_2d_axis0_shape():
    time = np.linspace(0, 200, 64)
    n_cols = 20
    signal = np.random.randn(64, n_cols) + np.linspace(1, 0, 64)[:, np.newaxis]
    corrected, baseline = remove_baseline(
        time, signal, method="polynomial", order=1, axis=0
    )
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
