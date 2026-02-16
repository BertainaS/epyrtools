"""
Frequency Analysis for Time-Domain EPR Signals

Comprehensive FFT-based frequency analysis with support for:
- 1D time-domain signals (Rabi, DEER, echo decay, etc.)
- 2D time-domain data with row-by-row 1D FFT processing
- 2D HYSCORE-type measurements with full 2D FFT

Includes DC offset removal and apodization windows for clean spectral analysis.
"""

from typing import Dict, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from scipy import fft
from scipy import signal as scipy_signal

from ..logging_config import get_logger

logger = get_logger(__name__)

try:
    from .apowin import apowin
except ImportError:
    # Handle direct execution
    from apowin import apowin


# ============================================================================
# Helper Functions for Common Operations
# ============================================================================


def _detect_time_units(time_data: np.ndarray) -> Tuple[str, str, float]:
    """
    Detect time unit and calculate time step in seconds.

    Parameters:
    -----------
    time_data : np.ndarray
        Time axis data

    Returns:
    --------
    tuple
        (time_unit, freq_unit, dt_seconds)
    """
    time_range = np.max(time_data) - np.min(time_data)
    dt_original = np.mean(np.diff(time_data))

    if time_range > 100:  # > 100 units, likely nanoseconds
        return "ns", "MHz", dt_original * 1e-9
    elif time_range > 1.0:  # 1-100 units, likely microseconds
        return "μs", "MHz", dt_original * 1e-6
    elif time_range > 0.01:  # 0.01-1 units, likely milliseconds
        return "ms", "kHz", dt_original * 1e-3
    elif time_range > 1e-6:  # 1e-6 to 0.01, likely seconds
        return "s", "Hz", dt_original
    else:  # Very small values, normalized time
        return "arb", "Hz", dt_original


def _convert_to_display_freq(frequencies_hz: np.ndarray, freq_unit: str) -> np.ndarray:
    """
    Convert frequencies from Hz to display units.

    Parameters:
    -----------
    frequencies_hz : np.ndarray
        Frequencies in Hz
    freq_unit : str
        Target unit ('MHz', 'kHz', or 'Hz')

    Returns:
    --------
    np.ndarray
        Frequencies in display units
    """
    if freq_unit == "MHz":
        return frequencies_hz / 1e6
    elif freq_unit == "kHz":
        return frequencies_hz / 1e3
    else:  # Hz
        return frequencies_hz


def _remove_dc_offset(
    signal: np.ndarray, axis: Optional[int] = None
) -> Tuple[np.ndarray, Union[float, np.ndarray]]:
    """
    Remove DC offset from signal.

    Parameters:
    -----------
    signal : np.ndarray
        Input signal
    axis : int, optional
        Axis along which to compute mean. None = mean over all data.

    Returns:
    --------
    tuple
        (processed_signal, dc_offset)
    """
    if axis is None:
        dc_offset = np.mean(signal)
        return signal - dc_offset, dc_offset
    else:
        dc_offset = np.mean(signal, axis=axis, keepdims=True)
        return signal - dc_offset, dc_offset


def _apply_window(
    signal: np.ndarray,
    window: Optional[str],
    window_alpha: Optional[float],
    axis: int = -1,
) -> np.ndarray:
    """
    Apply apodization window to signal.

    Parameters:
    -----------
    signal : np.ndarray
        Input signal
    window : str or None
        Window type ('hann', 'hamming', 'blackman', 'kaiser', None)
    window_alpha : float, optional
        Alpha parameter for Kaiser, Gaussian windows
    axis : int
        Axis along which to apply window (default: -1, last axis)

    Returns:
    --------
    np.ndarray
        Windowed signal
    """
    if window is None:
        return signal.copy()

    # Set default alpha for Kaiser and Gaussian
    if window in ["kaiser", "gaussian"] and window_alpha is None:
        window_alpha = 6.0

    # Get window size
    n_points = signal.shape[axis]

    # Generate window function
    if window_alpha is not None:
        window_func = apowin(window, n_points, alpha=window_alpha)
    else:
        window_func = apowin(window, n_points)

    # Apply window along specified axis
    if signal.ndim == 1:
        return signal * window_func
    elif signal.ndim == 2:
        if axis == -1 or axis == 1:
            return signal * window_func[np.newaxis, :]
        else:  # axis == 0
            return signal * window_func[:, np.newaxis]
    else:
        raise ValueError("Only 1D and 2D signals supported")


def analyze_frequencies(
    time_data: np.ndarray,
    signal_data: np.ndarray,
    window: Optional[str] = "hann",
    window_alpha: Optional[float] = None,
    zero_padding: int = 2,
    remove_dc: bool = True,
    plot: bool = True,
    freq_range: Optional[Tuple[float, float]] = None,
    **plot_kwargs,
) -> Dict:
    """
    FFT-based frequency analysis of time-domain EPR signals.

    This function performs clean FFT analysis to identify frequency components
    in time-dependent EPR signals, with proper DC offset removal.

    Parameters:
    -----------
    time_data : np.ndarray
        Time axis data (in ns, μs, or s)
    signal_data : np.ndarray
        EPR signal intensity vs time
    window : str or None, optional
        Apodization window type ('hann', 'hamming', 'blackman', 'kaiser', None)
        Default: 'hann'
    window_alpha : float, optional
        Alpha parameter for Kaiser, Gaussian windows (default: 6 for Kaiser)
    zero_padding : int, optional
        Zero padding factor (2 = double length, 4 = quadruple, etc.)
        Default: 2
    remove_dc : bool, optional
        Remove DC offset before analysis (recommended: True)
    plot : bool, optional
        Generate analysis plots. Default: True
    freq_range : tuple of float, optional
        Frequency range (min, max) to display in plots

    Returns:
    --------
    dict
        Analysis results containing:
        - 'frequencies': Frequency axis in appropriate units
        - 'power_spectrum': Power spectral density (normalized)
        - 'phase_spectrum': Phase spectrum
        - 'dominant_frequencies': List of peak frequencies
        - 'time_data': Original time data
        - 'processed_signal': Signal after DC removal and windowing
        - 'sampling_rate': Sampling rate in Hz
        - 'time_unit': Detected time unit
        - 'freq_unit': Frequency unit

    Examples:
    ---------
    >>> from epyr import eprload
    >>> from epyr.signalprocessing import analyze_frequencies
    >>>
    >>> # Load Rabi data
    >>> time, signal, params, _ = eprload('rabi_data.DTA')
    >>> result = analyze_frequencies(time, signal, window='hann', plot=True)
    >>> print(f"Dominant frequency: {result['dominant_frequencies'][0]:.3f} MHz")
    """

    # Input validation
    time_data = np.asarray(time_data)
    signal_data = np.asarray(signal_data)

    if time_data.shape != signal_data.shape:
        raise ValueError("Time and signal arrays must have the same shape")

    if len(time_data) < 4:
        raise ValueError("Need at least 4 data points for frequency analysis")

    logger.info(f"FFT Analysis of {len(signal_data)} data points")

    # Detect time units
    time_unit, freq_unit, dt_seconds = _detect_time_units(time_data)
    sampling_rate = 1.0 / dt_seconds
    logger.debug(f"Time unit: {time_unit}, Frequency unit: {freq_unit}")
    logger.debug(
        f"Sampling rate: {sampling_rate/{'MHz': 1e6, 'kHz': 1e3}.get(freq_unit, 1):.1f} {freq_unit}"
    )

    # Step 1: Remove DC offset (very important for EPR signals)
    if remove_dc:
        processed_signal, dc_offset = _remove_dc_offset(signal_data)
        logger.debug(f"Removed DC offset: {dc_offset:.6f}")
    else:
        processed_signal = signal_data.copy()
        logger.debug("DC offset not removed")

    # Step 2: Apply apodization window
    windowed_signal = _apply_window(processed_signal, window, window_alpha, axis=-1)
    if window is not None:
        if window_alpha is not None:
            logger.debug(f"Applied {window} window (alpha={window_alpha})")
        else:
            logger.debug(f"Applied {window} window")
    else:
        logger.debug("No window applied (rectangular)")

    # Step 3: Zero padding for better frequency resolution
    if zero_padding > 1:
        n_padded = len(windowed_signal) * zero_padding
        padded_signal = np.zeros(n_padded, dtype=windowed_signal.dtype)
        padded_signal[: len(windowed_signal)] = windowed_signal
        windowed_signal = padded_signal
        logger.debug(f"Zero padding: {len(processed_signal)} -> {n_padded} points")

    # Step 4: Perform FFT
    fft_result = fft.fft(windowed_signal)
    frequencies_hz = fft.fftfreq(len(windowed_signal), dt_seconds)

    # Take positive frequencies only
    n_pos = len(frequencies_hz) // 2
    frequencies_hz_pos = frequencies_hz[:n_pos]
    fft_positive = fft_result[:n_pos]

    # Convert frequencies to display units
    frequencies_display = _convert_to_display_freq(frequencies_hz_pos, freq_unit)

    # Step 5: Calculate power and phase spectra
    power_spectrum = np.abs(fft_positive) ** 2
    phase_spectrum = np.angle(fft_positive)

    # Normalize power spectrum
    if np.max(power_spectrum) > 0:
        power_spectrum = power_spectrum / np.max(power_spectrum)

    # Step 6: Find dominant frequencies (peaks above 10% of maximum)
    peak_threshold = 0.1
    peak_indices, _ = scipy_signal.find_peaks(power_spectrum, height=peak_threshold)
    dominant_frequencies_display = frequencies_display[peak_indices]

    # Sort by power (strongest first)
    if len(peak_indices) > 0:
        peak_powers = power_spectrum[peak_indices]
        sort_indices = np.argsort(peak_powers)[::-1]
        dominant_frequencies_display = dominant_frequencies_display[sort_indices]

    # Display results
    logger.info("Frequency Analysis Results:")
    logger.info(f"Frequency resolution: {frequencies_display[1]:.6f} {freq_unit}")
    logger.info(f"Maximum frequency: {frequencies_display[-1]:.3f} {freq_unit}")

    if len(dominant_frequencies_display) > 0:
        logger.info(f"Dominant frequencies ({freq_unit}):")
        for i, freq in enumerate(dominant_frequencies_display[:5]):  # Top 5
            if i < len(peak_indices):
                power_pct = power_spectrum[peak_indices[sort_indices[i]]] * 100
                logger.info(
                    f"  {i+1}. {freq:.6f} {freq_unit} (power: {power_pct:.1f}%)"
                )
    else:
        logger.info("No significant frequency peaks found")

    # Step 7: Create plots
    if plot:
        _plot_fft_analysis(
            time_data,
            signal_data,
            processed_signal,
            windowed_signal,
            frequencies_display,
            power_spectrum,
            phase_spectrum,
            dominant_frequencies_display,
            time_unit,
            freq_unit,
            freq_range,
            remove_dc,
            **plot_kwargs,
        )

    # Return results
    results = {
        "frequencies": frequencies_display,
        "power_spectrum": power_spectrum,
        "phase_spectrum": phase_spectrum,
        "dominant_frequencies": dominant_frequencies_display,
        "time_data": time_data,
        "processed_signal": processed_signal,
        "sampling_rate": sampling_rate,
        "time_unit": time_unit,
        "freq_unit": freq_unit,
        "dc_removed": remove_dc,
    }

    return results


def power_spectrum(
    time_data: np.ndarray,
    signal_data: np.ndarray,
    method: str = "welch",
    window: str = "hann",
    nperseg: Optional[int] = None,
    overlap: float = 0.5,
    remove_dc: bool = True,
    plot: bool = True,
) -> Dict:
    """
    Calculate power spectral density using Welch or periodogram methods.

    Parameters:
    -----------
    time_data : np.ndarray
        Time axis data
    signal_data : np.ndarray
        Signal data
    method : str
        Method: 'welch' or 'periodogram'
    window : str
        Window function for Welch method
    nperseg : int, optional
        Length of each segment for Welch method
    overlap : float
        Overlap fraction for Welch method (0-1)
    remove_dc : bool
        Remove DC offset before analysis
    plot : bool
        Generate plots

    Returns:
    --------
    dict
        Results with frequencies and power spectrum
    """

    # Remove DC offset if requested
    if remove_dc:
        signal_data, _ = _remove_dc_offset(signal_data)

    # Detect time units
    _, freq_unit, dt_seconds = _detect_time_units(time_data)
    sampling_rate = 1.0 / dt_seconds

    if method == "welch":
        if nperseg is None:
            nperseg = len(signal_data) // 4
        noverlap = int(nperseg * overlap)

        frequencies_hz, psd = scipy_signal.welch(
            signal_data,
            sampling_rate,
            window=window,
            nperseg=nperseg,
            noverlap=noverlap,
        )

    elif method == "periodogram":
        frequencies_hz, psd = scipy_signal.periodogram(
            signal_data, sampling_rate, window=window
        )
    else:
        raise ValueError(f"Unknown method: {method}")

    # Convert to display units
    frequencies = _convert_to_display_freq(frequencies_hz, freq_unit)

    # Normalize
    psd = psd / np.max(psd)

    if plot:
        plt.figure(figsize=(10, 6))
        plt.semilogy(frequencies, psd, linewidth=2)
        plt.xlabel(f"Frequency ({freq_unit})")
        plt.ylabel("Power Spectral Density")
        plt.title(f"Power Spectrum ({method.capitalize()} Method)")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return {
        "frequencies": frequencies,
        "psd": psd,
        "method": method,
        "freq_unit": freq_unit,
    }


def spectrogram_analysis(
    time_data: np.ndarray,
    signal_data: np.ndarray,
    window: str = "hann",
    nperseg: Optional[int] = None,
    overlap: float = 0.8,
    remove_dc: bool = True,
    plot: bool = True,
) -> Dict:
    """
    Time-frequency analysis using spectrogram.

    Parameters:
    -----------
    time_data : np.ndarray
        Time axis data
    signal_data : np.ndarray
        Signal data
    window : str
        Window function
    nperseg : int, optional
        Length of each segment
    overlap : float
        Overlap fraction (0-1)
    remove_dc : bool
        Remove DC offset
    plot : bool
        Generate spectrogram plot

    Returns:
    --------
    dict
        Results with time axis, frequencies, and spectrogram
    """

    # Remove DC offset if requested
    if remove_dc:
        signal_data, _ = _remove_dc_offset(signal_data)

    # Detect time units
    time_unit, freq_unit, dt_seconds = _detect_time_units(time_data)
    sampling_rate = 1.0 / dt_seconds

    if nperseg is None:
        nperseg = len(signal_data) // 8

    noverlap = int(nperseg * overlap)

    frequencies_hz, times_s, Sxx = scipy_signal.spectrogram(
        signal_data, sampling_rate, window=window, nperseg=nperseg, noverlap=noverlap
    )

    # Convert to display units
    frequencies = _convert_to_display_freq(frequencies_hz, freq_unit)

    # Convert time to original units
    time_offset = np.min(time_data)
    if time_unit == "ns":
        times = times_s / 1e-9 + time_offset
    elif time_unit == "μs":
        times = times_s / 1e-6 + time_offset
    elif time_unit == "ms":
        times = times_s / 1e-3 + time_offset
    else:
        times = times_s + time_offset

    if plot:
        plt.figure(figsize=(12, 8))
        plt.pcolormesh(
            times, frequencies, 10 * np.log10(Sxx + 1e-10), shading="gouraud"
        )
        plt.colorbar(label="Power (dB)")
        plt.xlabel(f"Time ({time_unit})")
        plt.ylabel(f"Frequency ({freq_unit})")
        plt.title("Spectrogram - Time-Frequency Analysis")
        plt.tight_layout()
        plt.show()

    return {
        "times": times,
        "frequencies": frequencies,
        "spectrogram": Sxx,
        "time_unit": time_unit,
        "freq_unit": freq_unit,
    }


def _plot_fft_analysis(
    time_data,
    signal_data,
    processed_signal,
    windowed_signal,
    frequencies,
    power_spectrum,
    phase_spectrum,
    dominant_frequencies,
    time_unit,
    freq_unit,
    freq_range,
    dc_removed,
    **plot_kwargs,
):
    """Helper function to create comprehensive FFT analysis plots."""

    figsize = plot_kwargs.get("figsize", (14, 10))
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Time domain - original and processed signal
    axes[0, 0].plot(time_data, signal_data, "b-", alpha=0.7, label="Original signal")
    if dc_removed:
        axes[0, 0].plot(
            time_data,
            processed_signal,
            "r-",
            linewidth=2,
            alpha=0.8,
            label="DC removed",
        )

    axes[0, 0].set_xlabel(f"Time ({time_unit})")
    axes[0, 0].set_ylabel("Signal Amplitude")
    axes[0, 0].set_title("Time Domain Signal")
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # Power spectrum (log scale)
    axes[0, 1].semilogy(frequencies, power_spectrum, "b-", linewidth=2)

    # Mark dominant frequencies
    for i, freq in enumerate(dominant_frequencies[:5]):
        if i < len(dominant_frequencies):
            axes[0, 1].axvline(
                freq,
                color="red",
                linestyle="--",
                alpha=0.7,
                label=f"Peak {i+1}: {freq:.3f}" if i < 3 else "",
            )

    axes[0, 1].set_xlabel(f"Frequency ({freq_unit})")
    axes[0, 1].set_ylabel("Normalized Power")
    axes[0, 1].set_title("Power Spectrum (Log Scale)")
    axes[0, 1].grid(True, alpha=0.3)

    if freq_range:
        axes[0, 1].set_xlim(freq_range)

    if len(dominant_frequencies) > 0:
        axes[0, 1].legend()

    # Processed signal ready for FFT (windowed + zero-padded)
    # Create time axis for the windowed signal (including zero padding)
    n_original = len(time_data)
    n_windowed = len(windowed_signal)

    # Time axis for windowed signal (extend original time range for zero padding)
    dt_original = np.mean(np.diff(time_data))
    time_start = time_data[0]
    time_windowed = time_start + np.arange(n_windowed) * dt_original

    axes[1, 0].plot(time_windowed, windowed_signal, "purple", linewidth=2)
    axes[1, 0].set_xlabel(f"Time ({time_unit})")
    axes[1, 0].set_ylabel("Signal Amplitude")
    axes[1, 0].set_title("Signal Sent to FFT (Windowed + Zero-Padded)")
    axes[1, 0].grid(True, alpha=0.3)

    # Add vertical line to show original data length
    if n_windowed > n_original:
        time_end_original = time_data[-1]
        axes[1, 0].axvline(
            time_end_original,
            color="red",
            linestyle="--",
            alpha=0.5,
            label=f"Original data end",
        )
        axes[1, 0].legend()

    # Power spectrum (linear scale)
    axes[1, 1].plot(frequencies, power_spectrum, "b-", linewidth=2)

    # Mark dominant frequencies
    for i, freq in enumerate(dominant_frequencies[:5]):
        if i < len(dominant_frequencies):
            axes[1, 1].axvline(freq, color="red", linestyle="--", alpha=0.7)

    axes[1, 1].set_xlabel(f"Frequency ({freq_unit})")
    axes[1, 1].set_ylabel("Normalized Power")
    axes[1, 1].set_title("Power Spectrum (Linear Scale)")
    axes[1, 1].grid(True, alpha=0.3)

    if freq_range:
        axes[1, 1].set_xlim(freq_range)

    # Style all subplots
    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()


def analyze_frequencies_2d(
    time_data: Union[np.ndarray, Tuple[np.ndarray, np.ndarray]],
    signal_data: np.ndarray,
    mode: Literal["row_by_row", "full_2d"] = "row_by_row",
    window: Optional[str] = "hann",
    window_alpha: Optional[float] = None,
    zero_padding: int = 2,
    remove_dc: bool = True,
    axis: int = 1,
    plot_result: bool = False,
    freq_range: Optional[Tuple[float, float]] = None,
    **plot_kwargs,
) -> Tuple:
    """
    FFT-based frequency analysis of 2D time-domain EPR signals.

    This function handles 2D EPR data with two processing modes:
    1. Row-by-row 1D FFT: Process each row/column independently (e.g., 2D Rabi)
    2. Full 2D FFT: Process both dimensions together (e.g., HYSCORE)

    Parameters:
    -----------
    time_data : np.ndarray or tuple of np.ndarray
        Time axis data. Can be:
        - Single 1D array: time axis for the FFT dimension
        - Tuple of two 1D arrays: (time_axis1, time_axis2) for 2D FFT
    signal_data : np.ndarray
        2D EPR signal intensity array (shape: n_traces x n_points)
    mode : str, optional
        Processing mode:
        - 'row_by_row': Apply 1D FFT to each row/column independently
        - 'full_2d': Apply 2D FFT to entire dataset (HYSCORE-type)
        Default: 'row_by_row'
    window : str or None, optional
        Apodization window type. Default: 'hann'
    window_alpha : float, optional
        Alpha parameter for Kaiser, Gaussian windows
    zero_padding : int, optional
        Zero padding factor. Default: 2
    remove_dc : bool, optional
        Remove DC offset before analysis. Default: True
    axis : int, optional
        Axis to process for row_by_row mode (0=columns, 1=rows). Default: 1
    plot_result : bool, optional
        Generate analysis plots. Default: False
    freq_range : tuple of float, optional
        Frequency range (min, max) to display in plots

    Returns:
    --------
    For mode='row_by_row':
        fq : np.ndarray
            Frequency axis (1D array)
        axis2 : np.ndarray
            Secondary axis (field, angle, trace index, etc.)
        spectrum : np.ndarray
            2D FFT spectrum magnitude (n_traces x n_frequencies)
        info : dict
            Analysis information (units, sampling_rate, mode, etc.)

    For mode='full_2d':
        fq1 : np.ndarray
            Frequency axis 1 (1D array)
        fq2 : np.ndarray
            Frequency axis 2 (1D array)
        spectrum : np.ndarray
            2D FFT spectrum magnitude (n_freq1 x n_freq2)
        info : dict
            Analysis information

    Examples:
    ---------
    >>> # Row-by-row 1D FFT (2D Rabi oscillations)
    >>> x_2d, y_2d, params, _ = eprload('rabi_2d.DTA')
    >>> fq, axis2, spectrum, info = analyze_frequencies_2d(
    ...     x_2d[0], y_2d, mode='row_by_row', plot_result=False)

    >>> # Full 2D FFT (HYSCORE)
    >>> x_hyscore, y_hyscore, params, _ = eprload('hyscore.DTA')
    >>> fq1, fq2, spectrum_2d, info = analyze_frequencies_2d(
    ...     (x_hyscore[0], x_hyscore[1]), y_hyscore,
    ...     mode='full_2d', plot_result=True)
    """

    # Input validation
    signal_data = np.asarray(signal_data)

    if signal_data.ndim != 2:
        raise ValueError(f"signal_data must be 2D array, got shape {signal_data.shape}")

    logger.info(f"2D FFT Analysis - Mode: {mode}")
    logger.info(f"Data shape: {signal_data.shape}")

    if mode == "row_by_row":
        return _analyze_2d_row_by_row(
            time_data,
            signal_data,
            window,
            window_alpha,
            zero_padding,
            remove_dc,
            axis,
            plot_result,
            freq_range,
            **plot_kwargs,
        )

    elif mode == "full_2d":
        return _analyze_2d_full(
            time_data,
            signal_data,
            window,
            window_alpha,
            zero_padding,
            remove_dc,
            plot_result,
            freq_range,
            **plot_kwargs,
        )

    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'row_by_row' or 'full_2d'")


def _analyze_2d_row_by_row(
    time_data,
    signal_data,
    window,
    window_alpha,
    zero_padding,
    remove_dc,
    axis,
    plot_result,
    freq_range,
    **plot_kwargs,
):
    """Row-by-row 1D FFT processing for 2D data"""

    # Extract time axis
    if isinstance(time_data, (tuple, list)):
        time_axis = time_data[axis]
        other_axis = time_data[1 - axis]
    else:
        time_axis = time_data
        other_axis = np.arange(signal_data.shape[1 - axis])

    time_axis = np.asarray(time_axis)
    other_axis = np.asarray(other_axis)

    # Determine units
    time_unit, freq_unit, dt_seconds = _detect_time_units(time_axis)
    sampling_rate = 1.0 / dt_seconds
    logger.debug(f"Time axis: {time_unit}, Frequency axis: {freq_unit}")
    logger.debug(f"Processing axis: {axis} (0=columns, 1=rows)")

    # Transpose if processing columns
    if axis == 0:
        signal_data = signal_data.T

    n_traces, n_points = signal_data.shape
    logger.info(f"Number of traces: {n_traces}, Points per trace: {n_points}")

    # Step 1: Remove DC offset
    if remove_dc:
        processed_signal, dc_offsets = _remove_dc_offset(signal_data, axis=1)
        logger.debug(
            f"Removed DC offset (mean across traces: {np.mean(dc_offsets):.6f})"
        )
    else:
        processed_signal = signal_data.copy()

    dc_removed_signal = processed_signal.copy()

    # Step 2: Apply window function
    windowed_signal = _apply_window(processed_signal, window, window_alpha, axis=1)
    if window is not None:
        if window_alpha is not None:
            logger.debug(f"Applied {window} window (alpha={window_alpha})")
        else:
            logger.debug(f"Applied {window} window")
    else:
        logger.debug("No window applied")

    # Step 3: Zero padding
    time_axis_extended = time_axis.copy()
    if zero_padding > 1:
        n_padded = n_points * zero_padding
        padded_signal = np.zeros((n_traces, n_padded), dtype=windowed_signal.dtype)
        padded_signal[:, :n_points] = windowed_signal

        # Extend time axis for zero-padded region
        dt = np.mean(np.diff(time_axis))
        time_extension = time_axis[-1] + dt * np.arange(1, n_padded - n_points + 1)
        time_axis_extended = np.concatenate([time_axis, time_extension])

        windowed_signal = padded_signal
        logger.debug(f"Zero padding: {n_points} -> {n_padded} points per trace")

    # Store the fully processed signal before FFT
    processed_signal_final = windowed_signal.copy()

    # Step 4: Perform FFT on each row
    fft_result = fft.fft(windowed_signal, axis=1)
    frequencies_hz = fft.fftfreq(windowed_signal.shape[1], dt_seconds)

    # Use fftshift to center zero frequency and get symmetric spectrum
    fft_result_shifted = fft.fftshift(fft_result, axes=1)
    frequencies_hz_shifted = fft.fftshift(frequencies_hz)

    # Convert to display units
    frequencies_display = _convert_to_display_freq(frequencies_hz_shifted, freq_unit)

    # Calculate spectrum magnitude and phase
    spectrum_magnitude = np.abs(fft_result_shifted)
    phase_spectrum = np.angle(fft_result_shifted)

    logger.info(f"Frequency resolution: {frequencies_display[1]:.6f} {freq_unit}")
    logger.info(f"Maximum frequency: {frequencies_display[-1]:.3f} {freq_unit}")

    # Create plots if requested
    if plot_result:
        _plot_2d_row_by_row(
            time_axis,
            time_axis_extended,
            other_axis,
            signal_data,
            dc_removed_signal,
            processed_signal_final,
            frequencies_display,
            spectrum_magnitude,
            time_unit,
            freq_unit,
            freq_range,
            axis,
            n_points,
            phase_spectrum,
            **plot_kwargs,
        )

    # Transpose back if needed
    if axis == 0:
        signal_data = signal_data.T
        spectrum_magnitude = spectrum_magnitude.T
        phase_spectrum = phase_spectrum.T

    # Create info dictionary
    info = {
        "mode": "row_by_row",
        "axis": axis,
        "time_unit": time_unit,
        "freq_unit": freq_unit,
        "sampling_rate": sampling_rate,
        "time_data": time_axis,
        "dc_removed": remove_dc,
        "window": window,
        "zero_padding": zero_padding,
        "phase_spectrum": phase_spectrum,
    }

    return frequencies_display, other_axis, spectrum_magnitude, info


def _analyze_2d_full(
    time_data,
    signal_data,
    window,
    window_alpha,
    zero_padding,
    remove_dc,
    plot_result,
    freq_range,
    **plot_kwargs,
):
    """Full 2D FFT processing for HYSCORE-type measurements"""

    # Extract time axes
    if not isinstance(time_data, (tuple, list)):
        raise ValueError(
            "For full 2D FFT, time_data must be tuple of (time_axis1, time_axis2)"
        )

    time_axis1, time_axis2 = time_data
    time_axis1 = np.asarray(time_axis1)
    time_axis2 = np.asarray(time_axis2)

    # Determine units for both axes
    time_unit1, freq_unit1, dt_seconds1 = _detect_time_units(time_axis1)
    time_unit2, freq_unit2, dt_seconds2 = _detect_time_units(time_axis2)

    sampling_rate1 = 1.0 / dt_seconds1
    sampling_rate2 = 1.0 / dt_seconds2

    logger.debug(
        f"Axis 1: {time_unit1} → {freq_unit1}, sampling rate: {sampling_rate1/{'MHz': 1e6, 'kHz': 1e3}.get(freq_unit1, 1):.1f} {freq_unit1}"
    )
    logger.debug(
        f"Axis 2: {time_unit2} → {freq_unit2}, sampling rate: {sampling_rate2/{'MHz': 1e6, 'kHz': 1e3}.get(freq_unit2, 1):.1f} {freq_unit2}"
    )

    n_points1, n_points2 = signal_data.shape
    logger.info(f"Data dimensions: {n_points1} x {n_points2}")

    # Remove DC offset
    if remove_dc:
        processed_signal, dc_offset = _remove_dc_offset(signal_data)
        logger.debug(f"Removed DC offset: {dc_offset:.6f}")
    else:
        processed_signal = signal_data.copy()

    # Apply 2D window function
    if window is not None:
        if window in ["kaiser", "gaussian"] and window_alpha is None:
            window_alpha = 6.0

        # Create 2D window as outer product of 1D windows
        if window_alpha is not None:
            window_func1 = apowin(window, n_points1, alpha=window_alpha)
            window_func2 = apowin(window, n_points2, alpha=window_alpha)
        else:
            window_func1 = apowin(window, n_points1)
            window_func2 = apowin(window, n_points2)

        window_2d = np.outer(window_func1, window_func2)
        windowed_signal = processed_signal * window_2d
        logger.debug(f"Applied 2D {window} window")
    else:
        windowed_signal = processed_signal.copy()
        logger.debug("No window applied")

    # Zero padding
    if zero_padding > 1:
        n_padded1 = n_points1 * zero_padding
        n_padded2 = n_points2 * zero_padding
        padded_signal = np.zeros((n_padded1, n_padded2), dtype=windowed_signal.dtype)
        padded_signal[:n_points1, :n_points2] = windowed_signal
        windowed_signal = padded_signal
        logger.debug(
            f"Zero padding: {n_points1}x{n_points2} -> {n_padded1}x{n_padded2}"
        )

    # Perform 2D FFT
    fft_result = fft.fft2(windowed_signal)
    frequencies_hz1 = fft.fftfreq(windowed_signal.shape[0], dt_seconds1)
    frequencies_hz2 = fft.fftfreq(windowed_signal.shape[1], dt_seconds2)

    # Use fftshift to center zero frequency and get symmetric spectrum
    fft_result_shifted = fft.fftshift(fft_result)
    frequencies_hz1_shifted = fft.fftshift(frequencies_hz1)
    frequencies_hz2_shifted = fft.fftshift(frequencies_hz2)

    # Convert to display units
    frequencies_display1 = _convert_to_display_freq(frequencies_hz1_shifted, freq_unit1)
    frequencies_display2 = _convert_to_display_freq(frequencies_hz2_shifted, freq_unit2)

    # Calculate spectrum magnitude and phase
    spectrum_magnitude = np.abs(fft_result_shifted)
    phase_spectrum = np.angle(fft_result_shifted)

    logger.info("Frequency resolution:")
    logger.info(f"  Axis 1: {frequencies_display1[1]:.6f} {freq_unit1}")
    logger.info(f"  Axis 2: {frequencies_display2[1]:.6f} {freq_unit2}")
    logger.info("Maximum frequencies:")
    logger.info(f"  Axis 1: {frequencies_display1[-1]:.3f} {freq_unit1}")
    logger.info(f"  Axis 2: {frequencies_display2[-1]:.3f} {freq_unit2}")

    # Create plots if requested
    if plot_result:
        _plot_2d_full(
            time_axis1,
            time_axis2,
            signal_data,
            processed_signal,
            frequencies_display1,
            frequencies_display2,
            spectrum_magnitude,
            time_unit1,
            time_unit2,
            freq_unit1,
            freq_unit2,
            freq_range,
            phase_spectrum,
            **plot_kwargs,
        )

    # Create info dictionary
    info = {
        "mode": "full_2d",
        "time_unit": (time_unit1, time_unit2),
        "freq_unit": (freq_unit1, freq_unit2),
        "sampling_rate": (sampling_rate1, sampling_rate2),
        "time_data": (time_axis1, time_axis2),
        "dc_removed": remove_dc,
        "window": window,
        "zero_padding": zero_padding,
        "phase_spectrum": phase_spectrum,
    }

    return frequencies_display1, frequencies_display2, spectrum_magnitude, info


def _plot_2d_row_by_row(
    time_axis,
    time_axis_extended,
    other_axis,
    signal_data,
    dc_removed_signal,
    processed_signal_final,
    frequencies,
    spectrum_magnitude,
    time_unit,
    freq_unit,
    freq_range,
    axis,
    n_points_original,
    phase_spectrum=None,
    **plot_kwargs,
):
    """4-panel plot: Original signal, FFT linear, FFT log, Phase"""

    figsize = plot_kwargs.get("figsize", (20, 10))
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Panel 1 (top-left): Original 2D signal
    im1 = axes[0, 0].imshow(
        signal_data,
        aspect="auto",
        cmap="RdBu_r",
        extent=[time_axis[0], time_axis[-1], other_axis[0], other_axis[-1]],
        origin="lower",
    )
    axes[0, 0].set_xlabel(f"Time ({time_unit})")
    axes[0, 0].set_ylabel("Trace Index")
    axes[0, 0].set_title("Original Signal")
    plt.colorbar(im1, ax=axes[0, 0], label="Amplitude")

    # Panel 2 (top-right): Magnitude spectrum (linear scale)
    im2 = axes[0, 1].imshow(
        spectrum_magnitude,
        aspect="auto",
        cmap="hot",
        extent=[frequencies[0], frequencies[-1], other_axis[0], other_axis[-1]],
        origin="lower",
    )
    axes[0, 1].set_xlabel(f"Frequency ({freq_unit})")
    axes[0, 1].set_ylabel("Trace Index")
    axes[0, 1].set_title("FFT Magnitude (linear scale)")
    if freq_range:
        axes[0, 1].set_xlim(freq_range)
    plt.colorbar(im2, ax=axes[0, 1], label="Magnitude")

    # Panel 3 (bottom-left): Magnitude spectrum (log scale)
    magnitude_log = np.log10(spectrum_magnitude + 1e-10)
    im3 = axes[1, 0].imshow(
        magnitude_log,
        aspect="auto",
        cmap="hot",
        extent=[frequencies[0], frequencies[-1], other_axis[0], other_axis[-1]],
        origin="lower",
    )
    axes[1, 0].set_xlabel(f"Frequency ({freq_unit})")
    axes[1, 0].set_ylabel("Trace Index")
    axes[1, 0].set_title("FFT Magnitude (log scale)")
    if freq_range:
        axes[1, 0].set_xlim(freq_range)
    plt.colorbar(im3, ax=axes[1, 0], label="log10(Magnitude)")

    # Panel 4 (bottom-right): Phase spectrum
    if phase_spectrum is not None:
        im4 = axes[1, 1].imshow(
            phase_spectrum,
            aspect="auto",
            cmap="twilight",
            extent=[frequencies[0], frequencies[-1], other_axis[0], other_axis[-1]],
            origin="lower",
            vmin=-np.pi,
            vmax=np.pi,
        )
        axes[1, 1].set_xlabel(f"Frequency ({freq_unit})")
        axes[1, 1].set_ylabel("Trace Index")
        axes[1, 1].set_title("Phase Spectrum")
        if freq_range:
            axes[1, 1].set_xlim(freq_range)
        plt.colorbar(im4, ax=axes[1, 1], label="Phase (rad)")
    else:
        axes[1, 1].text(
            0.5,
            0.5,
            "Phase spectrum\nnot available",
            ha="center",
            va="center",
            transform=axes[1, 1].transAxes,
        )
        axes[1, 1].set_title("Phase Spectrum")

    plt.tight_layout()
    plt.show()


def _plot_2d_full(
    time_axis1,
    time_axis2,
    signal_data,
    processed_signal,
    frequencies1,
    frequencies2,
    spectrum_magnitude,
    time_unit1,
    time_unit2,
    freq_unit1,
    freq_unit2,
    freq_range,
    phase_spectrum=None,
    **plot_kwargs,
):
    """4-panel plot: Original signal, FFT linear, FFT log, Phase"""

    figsize = plot_kwargs.get("figsize", (20, 10))
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Panel 1 (top-left): Original 2D time-domain signal
    im1 = axes[0, 0].imshow(
        signal_data,
        aspect="auto",
        cmap="RdBu_r",
        extent=[time_axis2[0], time_axis2[-1], time_axis1[0], time_axis1[-1]],
        origin="lower",
    )
    axes[0, 0].set_xlabel(f"Time 2 ({time_unit2})")
    axes[0, 0].set_ylabel(f"Time 1 ({time_unit1})")
    axes[0, 0].set_title("Original Signal")
    plt.colorbar(im1, ax=axes[0, 0], label="Amplitude")

    # Panel 2 (top-right): 2D Magnitude spectrum (linear scale)
    im2 = axes[0, 1].imshow(
        spectrum_magnitude,
        aspect="auto",
        cmap="hot",
        extent=[frequencies2[0], frequencies2[-1], frequencies1[0], frequencies1[-1]],
        origin="lower",
    )
    axes[0, 1].set_xlabel(f"Frequency 2 ({freq_unit2})")
    axes[0, 1].set_ylabel(f"Frequency 1 ({freq_unit1})")
    axes[0, 1].set_title("FFT Magnitude (linear scale)")
    if freq_range:
        axes[0, 1].set_xlim(freq_range)
        axes[0, 1].set_ylim(freq_range)
    plt.colorbar(im2, ax=axes[0, 1], label="Magnitude")

    # Panel 3 (bottom-left): 2D Magnitude spectrum (log scale)
    magnitude_log = np.log10(spectrum_magnitude + 1e-10)
    im3 = axes[1, 0].imshow(
        magnitude_log,
        aspect="auto",
        cmap="hot",
        extent=[frequencies2[0], frequencies2[-1], frequencies1[0], frequencies1[-1]],
        origin="lower",
    )
    axes[1, 0].set_xlabel(f"Frequency 2 ({freq_unit2})")
    axes[1, 0].set_ylabel(f"Frequency 1 ({freq_unit1})")
    axes[1, 0].set_title("FFT Magnitude (log scale)")
    if freq_range:
        axes[1, 0].set_xlim(freq_range)
        axes[1, 0].set_ylim(freq_range)
    plt.colorbar(im3, ax=axes[1, 0], label="log10(Magnitude)")

    # Panel 4 (bottom-right): 2D Phase spectrum
    if phase_spectrum is not None:
        im4 = axes[1, 1].imshow(
            phase_spectrum,
            aspect="auto",
            cmap="twilight",
            extent=[
                frequencies2[0],
                frequencies2[-1],
                frequencies1[0],
                frequencies1[-1],
            ],
            origin="lower",
            vmin=-np.pi,
            vmax=np.pi,
        )
        axes[1, 1].set_xlabel(f"Frequency 2 ({freq_unit2})")
        axes[1, 1].set_ylabel(f"Frequency 1 ({freq_unit1})")
        axes[1, 1].set_title("Phase Spectrum")
        if freq_range:
            axes[1, 1].set_xlim(freq_range)
            axes[1, 1].set_ylim(freq_range)
        plt.colorbar(im4, ax=axes[1, 1], label="Phase (rad)")
    else:
        axes[1, 1].text(
            0.5,
            0.5,
            "Phase spectrum\nnot available",
            ha="center",
            va="center",
            transform=axes[1, 1].transAxes,
        )
        axes[1, 1].set_title("Phase Spectrum")

    plt.tight_layout()
    plt.show()


def demo():
    """
    Simple demonstration of EPR FFT analysis.
    Shows clean frequency analysis with DC removal and windowing.
    """
    logger.info("EPR Signal Processing - Simplified FFT Analysis Demo")
    logger.info("=" * 60)
    logger.info("Focus on clean FFT analysis with proper DC removal")
    logger.info("")

    # Create synthetic Rabi oscillation
    t = np.linspace(0, 500, 256)  # 500 ns, 256 points
    rabi_freq = 8.5  # MHz
    decay_time = 120  # ns
    noise_level = 0.04
    dc_offset = 0.1  # Add DC offset to demonstrate removal

    # Clean Rabi signal with DC offset and noise
    clean_signal = np.sin(2 * np.pi * rabi_freq * t * 1e-3) * np.exp(-t / decay_time)
    noisy_signal = clean_signal + dc_offset + noise_level * np.random.randn(len(t))

    logger.info("Synthetic Rabi signal:")
    logger.info(f"  Target frequency: {rabi_freq} MHz")
    logger.info(f"  Decay time: {decay_time} ns")
    logger.info(f"  DC offset: {dc_offset}")
    logger.info(f"  Noise level: {noise_level:.1%}")
    logger.info(f"  Data points: {len(t)}")

    # Demo 1: Analysis with DC removal
    logger.info("")
    logger.info("=" * 50)
    logger.info("DEMO 1: FFT Analysis with DC Removal")
    logger.info("=" * 50)

    result_dc = analyze_frequencies(
        t,
        noisy_signal,
        window="hann",
        remove_dc=True,
        zero_padding=4,
        plot=True,
        freq_range=(0, 20),
    )

    if len(result_dc["dominant_frequencies"]) > 0:
        detected_freq = result_dc["dominant_frequencies"][0]
        error = abs(detected_freq - rabi_freq) / rabi_freq * 100
        logger.info("Results with DC removal:")
        logger.info(f"  Detected: {detected_freq:.3f} MHz")
        logger.info(f"  Error: {error:.2f}%")
        if error < 5:
            logger.info("  --> Excellent frequency detection!")

    # Demo 2: Comparison without DC removal
    logger.info("")
    logger.info("=" * 50)
    logger.info("DEMO 2: Comparison without DC Removal")
    logger.info("=" * 50)

    result_no_dc = analyze_frequencies(
        t,
        noisy_signal,
        window="hann",
        remove_dc=False,
        zero_padding=4,
        plot=True,
        freq_range=(0, 20),
    )

    # Demo 3: Window comparison
    logger.info("")
    logger.info("=" * 50)
    logger.info("DEMO 3: Window Function Effects")
    logger.info("=" * 50)

    windows = ["hann", "hamming", "blackman"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, window in enumerate(windows):
        result = analyze_frequencies(
            t, noisy_signal, window=window, remove_dc=True, plot=False
        )

        axes[i].semilogy(result["frequencies"], result["power_spectrum"])
        axes[i].set_title(f"{window.capitalize()} Window")
        axes[i].set_xlabel(f'Frequency ({result["freq_unit"]})')
        axes[i].set_ylabel("Power")
        axes[i].grid(True, alpha=0.3)
        axes[i].set_xlim(0, 20)

        # Mark dominant frequency
        if len(result["dominant_frequencies"]) > 0:
            peak_freq = result["dominant_frequencies"][0]
            axes[i].axvline(
                peak_freq,
                color="red",
                linestyle="--",
                alpha=0.7,
                label=f"{peak_freq:.2f} MHz",
            )
            axes[i].legend()

    plt.tight_layout()
    plt.show()

    # Demo 4: Power spectrum methods
    logger.info("")
    logger.info("=" * 50)
    logger.info("DEMO 4: Power Spectrum Methods")
    logger.info("=" * 50)

    psd_welch = power_spectrum(
        t, noisy_signal, method="welch", remove_dc=True, plot=True
    )
    logger.info("Welch method completed")

    psd_periodogram = power_spectrum(
        t, noisy_signal, method="periodogram", remove_dc=True, plot=True
    )
    logger.info("Periodogram method completed")

    logger.info("")
    logger.info("=" * 60)
    logger.info("DEMO COMPLETED!")
    logger.info("=" * 60)
    logger.info("Key Points Demonstrated:")
    logger.info("  * DC offset removal is crucial for clean spectra")
    logger.info("  * Window functions reduce spectral leakage")
    logger.info("  * Zero padding improves frequency resolution")
    logger.info("  * Multiple methods available for power spectra")
    logger.info("  * Automatic time unit detection (ns → MHz)")
    logger.info("Simplified module ready for EPR frequency analysis!")


if __name__ == "__main__":
    np.random.seed(42)  # For reproducible results
    demo()
