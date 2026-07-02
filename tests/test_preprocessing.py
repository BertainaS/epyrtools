"""Tests for epyr.signalprocessing.preprocessing."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from epyr.signalprocessing.preprocessing import apodize, remove_baseline

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


from epyr.signalprocessing.preprocessing import zero_pad  # noqa: E402

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
