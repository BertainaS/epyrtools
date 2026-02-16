"""
Apodization windows for signal processing
Modern implementation with additional window types and features
"""

from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from scipy import special

from ..logging_config import get_logger

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


def window_comparison(n_points=256):
    """
    Compare different window types side by side.

    Parameters:
    -----------
    n_points : int
        Number of points for each window
    """

    # Modern colors
    colors = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#a65628",
        "#f781bf",
        "#999999",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Basic windows (no parameters)
    ax = axes[0, 0]
    basic_windows = ["hamming", "hann", "blackman", "bartlett"]

    for i, win_type in enumerate(basic_windows):
        w = apowin(win_type, n_points)
        x = np.arange(n_points)
        ax.plot(x, w, color=colors[i], linewidth=2.5, label=win_type.capitalize())

    ax.set_title("Basic Windows", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.legend()
    ax.grid(alpha=0.3)

    # Specialized windows
    ax = axes[0, 1]
    specialized_windows = ["connes", "cosine", "welch"]

    for i, win_type in enumerate(specialized_windows):
        w = apowin(win_type, n_points)
        x = np.arange(n_points)
        ax.plot(x, w, color=colors[i], linewidth=2.5, label=win_type.capitalize())

    ax.set_title("Specialized Windows", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.legend()
    ax.grid(alpha=0.3)

    # Parametric windows
    ax = axes[1, 0]

    # Kaiser with different beta values
    betas = [3, 6, 9]
    for i, beta in enumerate(betas):
        w = apowin("kaiser", n_points, alpha=beta)
        x = np.arange(n_points)
        ax.plot(x, w, color=colors[i], linewidth=2.5, label=f"Kaiser (β={beta})")

    ax.set_title("Kaiser Windows", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.legend()
    ax.grid(alpha=0.3)

    # Half windows
    ax = axes[1, 1]

    # Show full vs half windows
    w_full = apowin("hann", n_points)
    w_left = apowin("hann", n_points // 2, half_window="left")
    w_right = apowin("hann", n_points // 2, half_window="right")

    x_full = np.arange(n_points)
    x_left = np.arange(n_points // 2)
    x_right = np.arange(n_points // 2, n_points)

    ax.plot(x_full, w_full, color=colors[0], linewidth=2.5, label="Full Hann")
    ax.plot(x_left, w_left, color=colors[1], linewidth=2.5, label="Left half")
    ax.plot(x_right, w_right, color=colors[2], linewidth=2.5, label="Right half")

    ax.set_title("Full vs Half Windows", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.legend()
    ax.grid(alpha=0.3)

    # Style all subplots
    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(0, n_points - 1)
        ax.set_ylim(0, 1.05)

    plt.tight_layout()
    plt.show()


def frequency_response_demo(n_points=256):
    """
    Show frequency response characteristics of different windows.
    """

    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    windows = ["hamming", "hann", "blackman", "kaiser"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Window shapes
    for i, win_type in enumerate(windows):
        if win_type == "kaiser":
            w = apowin(win_type, n_points, alpha=6)
        else:
            w = apowin(win_type, n_points)

        ax1.plot(w, color=colors[i], linewidth=2.5, label=win_type.capitalize())

    ax1.set_title("Window Functions", fontweight="bold")
    ax1.set_xlabel("Sample")
    ax1.set_ylabel("Amplitude")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Frequency response (magnitude)
    freqs = np.fft.fftfreq(n_points, 1.0)
    freqs_pos = freqs[: n_points // 2]

    for i, win_type in enumerate(windows):
        if win_type == "kaiser":
            w = apowin(win_type, n_points, alpha=6)
        else:
            w = apowin(win_type, n_points)

        # FFT and magnitude response
        W = np.fft.fft(w, n_points)
        magnitude_db = 20 * np.log10(np.abs(W[: n_points // 2]) + 1e-10)

        ax2.plot(
            freqs_pos,
            magnitude_db,
            color=colors[i],
            linewidth=2,
            label=win_type.capitalize(),
        )

    ax2.set_title("Frequency Response", fontweight="bold")
    ax2.set_xlabel("Normalized Frequency")
    ax2.set_ylabel("Magnitude (dB)")
    ax2.set_ylim(-120, 10)
    ax2.legend()
    ax2.grid(alpha=0.3)

    # Style
    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()


def apply_window_demo():
    """Demonstrate windowing effect on a test signal"""

    # Create test signal: sum of sinusoids
    n_points = 512
    t = np.arange(n_points)

    # Signal with two frequencies
    f1, f2 = 0.05, 0.12  # Normalized frequencies
    signal = (
        np.sin(2 * np.pi * f1 * t)
        + 0.5 * np.sin(2 * np.pi * f2 * t)
        + 0.1 * np.random.randn(n_points)
    )

    # Windows to compare
    windows = {
        "Rectangular (None)": np.ones(n_points),
        "Hamming": apowin("hamming", n_points),
        "Hann": apowin("hann", n_points),
        "Blackman": apowin("blackman", n_points),
    }

    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Time domain
    ax = axes[0, 0]
    ax.plot(t, signal, "k-", linewidth=1, alpha=0.7)
    ax.set_title("Original Signal", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.grid(alpha=0.3)

    # Windows
    ax = axes[0, 1]
    for i, (name, window) in enumerate(windows.items()):
        if name != "Rectangular (None)":
            ax.plot(t, window, color=colors[i - 1], linewidth=2, label=name)

    ax.set_title("Apodization Windows", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.legend()
    ax.grid(alpha=0.3)

    # Windowed signals
    ax = axes[1, 0]
    for i, (name, window) in enumerate(windows.items()):
        windowed_signal = signal * window
        if i < 2:  # Show only first two for clarity
            ax.plot(
                t,
                windowed_signal,
                color=colors[i],
                linewidth=1.5,
                alpha=0.8,
                label=f"{name} windowed",
            )

    ax.set_title("Windowed Signals", fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.legend()
    ax.grid(alpha=0.3)

    # Frequency spectra
    ax = axes[1, 1]
    freqs = np.fft.fftfreq(n_points, 1.0)[: n_points // 2]

    for i, (name, window) in enumerate(windows.items()):
        windowed_signal = signal * window
        spectrum = np.fft.fft(windowed_signal)
        magnitude_db = 20 * np.log10(np.abs(spectrum[: n_points // 2]) + 1e-10)

        ax.plot(
            freqs, magnitude_db, color=colors[i], linewidth=2, label=name, alpha=0.8
        )

    ax.set_title("Frequency Spectra", fontweight="bold")
    ax.set_xlabel("Normalized Frequency")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_ylim(-60, 40)
    ax.legend()
    ax.grid(alpha=0.3)

    # Style all subplots
    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()

    logger.info("\nWindowing Effects:")
    logger.info("- Rectangular: Sharp cutoff, high sidelobes")
    logger.info("- Hamming: Good sidelobe suppression, moderate resolution")
    logger.info("- Hann: Smooth transition, balanced performance")
    logger.info("- Blackman: Excellent sidelobe suppression, wider main lobe")


def demo():
    """Comprehensive demonstration of apodization windows"""

    logger.info("🪟 Apodization Windows Demo")
    logger.info("=" * 40)

    # Basic window examples
    logger.info("\n1. Basic window generation:")

    # Create some example windows
    n = 64

    hamming_win = apowin("hamming", n)
    logger.info(f"-> Hamming window: {n} points, peak = {np.max(hamming_win):.3f}")

    kaiser_win = apowin("kaiser", n, alpha=6)
    logger.info(f"-> Kaiser window (β=6): {n} points, peak = {np.max(kaiser_win):.3f}")

    # Half window
    half_win = apowin("hann", n // 2, half_window="right")
    logger.info(f"-> Half Hann window: {len(half_win)} points")

    # Show comparison plots
    logger.info("\n2. Showing window comparison...")
    window_comparison()

    logger.info("\n3. Showing frequency response characteristics...")
    frequency_response_demo()

    logger.info("\n4. Showing practical application to signals...")
    apply_window_demo()

    # Print window properties
    logger.info("\n>> Window Properties Summary:")
    logger.info("Window Type    | Main Lobe Width | Peak Sidelobe (dB) | Use Case")
    logger.info("-" * 70)
    logger.info(
        "Rectangular    |      2          |       -13          | Maximum resolution"
    )
    logger.info(
        "Hamming        |      4          |       -41          | General purpose"
    )
    logger.info(
        "Hann           |      4          |       -32          | Smooth spectrum"
    )
    logger.info("Blackman       |      6          |       -58          | Low sidelobes")
    logger.info(
        "Kaiser (β=6)   |      5          |       -50          | Adjustable tradeoff"
    )
    logger.info(
        "Gaussian       |   Variable      |       -55          | Smooth, parametric"
    )


if __name__ == "__main__":
    demo()
