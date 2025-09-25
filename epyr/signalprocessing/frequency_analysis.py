"""
Simplified Frequency Analysis for Time-Domain EPR Signals

Focuses on FFT-based frequency analysis with proper DC offset removal
and apodization windows for clean spectral analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Union, Optional, Tuple, Dict
from scipy import signal as scipy_signal, fft

try:
    from .apowin import apowin
except ImportError:
    # Handle direct execution
    from apowin import apowin


def analyze_frequencies(time_data: np.ndarray, signal_data: np.ndarray,
                       window: Optional[str] = 'hann', window_alpha: Optional[float] = None,
                       zero_padding: int = 2, remove_dc: bool = True,
                       plot: bool = True, freq_range: Optional[Tuple[float, float]] = None,
                       **plot_kwargs) -> Dict:
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

    print(f"FFT Analysis of {len(signal_data)} data points")

    # Calculate time step
    dt_original = np.mean(np.diff(time_data))

    # Smart unit detection based on typical EPR time scales
    time_range = np.max(time_data) - np.min(time_data)

    if time_range > 100:  # > 100 units, likely nanoseconds (0-500 ns)
        time_unit = 'ns'
        freq_unit = 'MHz'
        dt_seconds = dt_original * 1e-9
    elif time_range > 1.0:  # 1-100 units, likely microseconds (0-50 μs)
        time_unit = 'μs'
        freq_unit = 'MHz'
        dt_seconds = dt_original * 1e-6
    elif time_range > 0.01:  # 0.01-1 units, likely milliseconds (0-10 ms)
        time_unit = 'ms'
        freq_unit = 'kHz'
        dt_seconds = dt_original * 1e-3
    elif time_range > 1e-6:  # 1e-6 to 0.01, likely seconds
        time_unit = 's'
        freq_unit = 'Hz'
        dt_seconds = dt_original
    else:  # Very small values, normalized time
        time_unit = 'arb'
        freq_unit = 'Hz'
        dt_seconds = dt_original

    sampling_rate = 1.0 / dt_seconds
    print(f"Time unit: {time_unit}, Frequency unit: {freq_unit}")
    print(f"Sampling rate: {sampling_rate/{'MHz': 1e6, 'kHz': 1e3}.get(freq_unit, 1):.1f} {freq_unit}")

    # Step 1: Remove DC offset (very important for EPR signals)
    processed_signal = signal_data.copy()
    if remove_dc:
        dc_offset = np.mean(signal_data)
        processed_signal = signal_data - dc_offset
        print(f"Removed DC offset: {dc_offset:.6f}")
    else:
        print("DC offset not removed")

    # Step 2: Apply apodization window
    if window is not None:
        if window in ['kaiser', 'gaussian'] and window_alpha is None:
            window_alpha = 6.0  # Default for Kaiser window

        if window_alpha is not None:
            window_func = apowin(window, len(processed_signal), alpha=window_alpha)
            print(f"Applied {window} window (alpha={window_alpha})")
        else:
            window_func = apowin(window, len(processed_signal))
            print(f"Applied {window} window")

        windowed_signal = processed_signal * window_func
    else:
        windowed_signal = processed_signal.copy()
        print("No window applied (rectangular)")

    # Step 3: Zero padding for better frequency resolution
    if zero_padding > 1:
        n_padded = len(windowed_signal) * zero_padding
        padded_signal = np.zeros(n_padded, dtype=windowed_signal.dtype)
        padded_signal[:len(windowed_signal)] = windowed_signal
        windowed_signal = padded_signal
        print(f"Zero padding: {len(processed_signal)} -> {n_padded} points")

    # Step 4: Perform FFT
    fft_result = fft.fft(windowed_signal)
    frequencies_hz = fft.fftfreq(len(windowed_signal), dt_seconds)

    # Take positive frequencies only
    n_pos = len(frequencies_hz) // 2
    frequencies_hz_pos = frequencies_hz[:n_pos]
    fft_positive = fft_result[:n_pos]

    # Convert frequencies to display units
    if freq_unit == 'MHz':
        frequencies_display = frequencies_hz_pos / 1e6
    elif freq_unit == 'kHz':
        frequencies_display = frequencies_hz_pos / 1e3
    else:  # Hz
        frequencies_display = frequencies_hz_pos

    # Step 5: Calculate power and phase spectra
    power_spectrum = np.abs(fft_positive)**2
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
    print(f"\nFrequency Analysis Results:")
    print(f"Frequency resolution: {frequencies_display[1]:.6f} {freq_unit}")
    print(f"Maximum frequency: {frequencies_display[-1]:.3f} {freq_unit}")

    if len(dominant_frequencies_display) > 0:
        print(f"\nDominant frequencies ({freq_unit}):")
        for i, freq in enumerate(dominant_frequencies_display[:5]):  # Top 5
            if i < len(peak_indices):
                power_pct = power_spectrum[peak_indices[sort_indices[i]]] * 100
                print(f"  {i+1}. {freq:.6f} {freq_unit} (power: {power_pct:.1f}%)")
    else:
        print("No significant frequency peaks found")

    # Step 7: Create plots
    if plot:
        _plot_fft_analysis(time_data, signal_data, processed_signal, windowed_signal,
                          frequencies_display, power_spectrum, phase_spectrum,
                          dominant_frequencies_display, time_unit, freq_unit,
                          freq_range, remove_dc, **plot_kwargs)

    # Return results
    results = {
        'frequencies': frequencies_display,
        'power_spectrum': power_spectrum,
        'phase_spectrum': phase_spectrum,
        'dominant_frequencies': dominant_frequencies_display,
        'time_data': time_data,
        'processed_signal': processed_signal,
        'sampling_rate': sampling_rate,
        'time_unit': time_unit,
        'freq_unit': freq_unit,
        'dc_removed': remove_dc
    }

    return results


def power_spectrum(time_data: np.ndarray, signal_data: np.ndarray,
                  method: str = 'welch', window: str = 'hann',
                  nperseg: Optional[int] = None, overlap: float = 0.5,
                  remove_dc: bool = True, plot: bool = True) -> Dict:
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
        signal_data = signal_data - np.mean(signal_data)

    # Get proper time step and sampling rate
    dt_original = np.mean(np.diff(time_data))
    time_range = np.max(time_data) - np.min(time_data)

    # Use same unit detection as analyze_frequencies
    if time_range > 100:
        dt_seconds = dt_original * 1e-9
        freq_unit = 'MHz'
    elif time_range > 0.1:
        dt_seconds = dt_original * 1e-6
        freq_unit = 'MHz'
    elif time_range > 0.001:
        dt_seconds = dt_original * 1e-3
        freq_unit = 'kHz'
    else:
        dt_seconds = dt_original
        freq_unit = 'Hz'

    sampling_rate = 1.0 / dt_seconds

    if method == 'welch':
        if nperseg is None:
            nperseg = len(signal_data) // 4
        noverlap = int(nperseg * overlap)

        frequencies_hz, psd = scipy_signal.welch(signal_data, sampling_rate,
                                       window=window, nperseg=nperseg,
                                       noverlap=noverlap)

    elif method == 'periodogram':
        frequencies_hz, psd = scipy_signal.periodogram(signal_data, sampling_rate,
                                             window=window)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Convert to display units
    if freq_unit == 'MHz':
        frequencies = frequencies_hz / 1e6
    elif freq_unit == 'kHz':
        frequencies = frequencies_hz / 1e3
    else:
        frequencies = frequencies_hz

    # Normalize
    psd = psd / np.max(psd)

    if plot:
        plt.figure(figsize=(10, 6))
        plt.semilogy(frequencies, psd, linewidth=2)
        plt.xlabel(f'Frequency ({freq_unit})')
        plt.ylabel('Power Spectral Density')
        plt.title(f'Power Spectrum ({method.capitalize()} Method)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return {'frequencies': frequencies, 'psd': psd, 'method': method, 'freq_unit': freq_unit}


def spectrogram_analysis(time_data: np.ndarray, signal_data: np.ndarray,
                        window: str = 'hann', nperseg: Optional[int] = None,
                        overlap: float = 0.8, remove_dc: bool = True,
                        plot: bool = True) -> Dict:
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
        signal_data = signal_data - np.mean(signal_data)

    dt_original = np.mean(np.diff(time_data))
    time_range = np.max(time_data) - np.min(time_data)

    # Unit detection
    if time_range > 100:
        dt_seconds = dt_original * 1e-9
        time_unit = 'ns'
        freq_unit = 'MHz'
    elif time_range > 0.1:
        dt_seconds = dt_original * 1e-6
        time_unit = 'μs'
        freq_unit = 'MHz'
    elif time_range > 0.001:
        dt_seconds = dt_original * 1e-3
        time_unit = 'ms'
        freq_unit = 'kHz'
    else:
        dt_seconds = dt_original
        time_unit = 's'
        freq_unit = 'Hz'

    sampling_rate = 1.0 / dt_seconds

    if nperseg is None:
        nperseg = len(signal_data) // 8

    noverlap = int(nperseg * overlap)

    frequencies_hz, times_s, Sxx = scipy_signal.spectrogram(signal_data, sampling_rate,
                                                window=window, nperseg=nperseg,
                                                noverlap=noverlap)

    # Convert to display units
    if freq_unit == 'MHz':
        frequencies = frequencies_hz / 1e6
    elif freq_unit == 'kHz':
        frequencies = frequencies_hz / 1e3
    else:
        frequencies = frequencies_hz

    # Convert time to original units
    time_offset = np.min(time_data)
    if time_unit == 'ns':
        times = times_s / 1e-9 + time_offset
    elif time_unit == 'μs':
        times = times_s / 1e-6 + time_offset
    elif time_unit == 'ms':
        times = times_s / 1e-3 + time_offset
    else:
        times = times_s + time_offset

    if plot:
        plt.figure(figsize=(12, 8))
        plt.pcolormesh(times, frequencies, 10*np.log10(Sxx + 1e-10), shading='gouraud')
        plt.colorbar(label='Power (dB)')
        plt.xlabel(f'Time ({time_unit})')
        plt.ylabel(f'Frequency ({freq_unit})')
        plt.title('Spectrogram - Time-Frequency Analysis')
        plt.tight_layout()
        plt.show()

    return {'times': times, 'frequencies': frequencies, 'spectrogram': Sxx,
            'time_unit': time_unit, 'freq_unit': freq_unit}


def _plot_fft_analysis(time_data, signal_data, processed_signal, windowed_signal,
                      frequencies, power_spectrum, phase_spectrum,
                      dominant_frequencies, time_unit, freq_unit,
                      freq_range, dc_removed, **plot_kwargs):
    """Helper function to create comprehensive FFT analysis plots."""

    figsize = plot_kwargs.get('figsize', (14, 10))
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Time domain - original and processed signal
    axes[0,0].plot(time_data, signal_data, 'b-', alpha=0.7, label='Original signal')
    if dc_removed:
        axes[0,0].plot(time_data, processed_signal, 'r-', linewidth=2,
                      alpha=0.8, label='DC removed')

    axes[0,0].set_xlabel(f'Time ({time_unit})')
    axes[0,0].set_ylabel('Signal Amplitude')
    axes[0,0].set_title('Time Domain Signal')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)

    # Power spectrum (log scale)
    axes[0,1].semilogy(frequencies, power_spectrum, 'b-', linewidth=2)

    # Mark dominant frequencies
    for i, freq in enumerate(dominant_frequencies[:5]):
        if i < len(dominant_frequencies):
            axes[0,1].axvline(freq, color='red', linestyle='--', alpha=0.7,
                            label=f'Peak {i+1}: {freq:.3f}' if i < 3 else '')

    axes[0,1].set_xlabel(f'Frequency ({freq_unit})')
    axes[0,1].set_ylabel('Normalized Power')
    axes[0,1].set_title('Power Spectrum (Log Scale)')
    axes[0,1].grid(True, alpha=0.3)

    if freq_range:
        axes[0,1].set_xlim(freq_range)

    if len(dominant_frequencies) > 0:
        axes[0,1].legend()

    # Processed signal ready for FFT (windowed + zero-padded)
    # Create time axis for the windowed signal (including zero padding)
    n_original = len(time_data)
    n_windowed = len(windowed_signal)

    # Time axis for windowed signal (extend original time range for zero padding)
    dt_original = np.mean(np.diff(time_data))
    time_start = time_data[0]
    time_windowed = time_start + np.arange(n_windowed) * dt_original

    axes[1,0].plot(time_windowed, windowed_signal, 'purple', linewidth=2)
    axes[1,0].set_xlabel(f'Time ({time_unit})')
    axes[1,0].set_ylabel('Signal Amplitude')
    axes[1,0].set_title('Signal Sent to FFT (Windowed + Zero-Padded)')
    axes[1,0].grid(True, alpha=0.3)

    # Add vertical line to show original data length
    if n_windowed > n_original:
        time_end_original = time_data[-1]
        axes[1,0].axvline(time_end_original, color='red', linestyle='--', alpha=0.5,
                         label=f'Original data end')
        axes[1,0].legend()

    # Power spectrum (linear scale)
    axes[1,1].plot(frequencies, power_spectrum, 'b-', linewidth=2)

    # Mark dominant frequencies
    for i, freq in enumerate(dominant_frequencies[:5]):
        if i < len(dominant_frequencies):
            axes[1,1].axvline(freq, color='red', linestyle='--', alpha=0.7)

    axes[1,1].set_xlabel(f'Frequency ({freq_unit})')
    axes[1,1].set_ylabel('Normalized Power')
    axes[1,1].set_title('Power Spectrum (Linear Scale)')
    axes[1,1].grid(True, alpha=0.3)

    if freq_range:
        axes[1,1].set_xlim(freq_range)

    # Style all subplots
    for ax in axes.flat:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.show()


def demo():
    """
    Simple demonstration of EPR FFT analysis.
    Shows clean frequency analysis with DC removal and windowing.
    """
    print("EPR Signal Processing - Simplified FFT Analysis Demo")
    print("=" * 60)
    print("Focus on clean FFT analysis with proper DC removal")
    print()

    # Create synthetic Rabi oscillation
    t = np.linspace(0, 500, 256)  # 500 ns, 256 points
    rabi_freq = 8.5  # MHz
    decay_time = 120  # ns
    noise_level = 0.04
    dc_offset = 0.1  # Add DC offset to demonstrate removal

    # Clean Rabi signal with DC offset and noise
    clean_signal = np.sin(2 * np.pi * rabi_freq * t * 1e-3) * np.exp(-t/decay_time)
    noisy_signal = clean_signal + dc_offset + noise_level * np.random.randn(len(t))

    print(f"Synthetic Rabi signal:")
    print(f"  Target frequency: {rabi_freq} MHz")
    print(f"  Decay time: {decay_time} ns")
    print(f"  DC offset: {dc_offset}")
    print(f"  Noise level: {noise_level:.1%}")
    print(f"  Data points: {len(t)}")

    # Demo 1: Analysis with DC removal
    print(f"\n" + "="*50)
    print("DEMO 1: FFT Analysis with DC Removal")
    print("="*50)

    result_dc = analyze_frequencies(t, noisy_signal, window='hann',
                                  remove_dc=True, zero_padding=4,
                                  plot=True, freq_range=(0, 20))

    if len(result_dc['dominant_frequencies']) > 0:
        detected_freq = result_dc['dominant_frequencies'][0]
        error = abs(detected_freq - rabi_freq) / rabi_freq * 100
        print(f"\nResults with DC removal:")
        print(f"  Detected: {detected_freq:.3f} MHz")
        print(f"  Error: {error:.2f}%")
        if error < 5:
            print("  --> Excellent frequency detection!")

    # Demo 2: Comparison without DC removal
    print(f"\n" + "="*50)
    print("DEMO 2: Comparison without DC Removal")
    print("="*50)

    result_no_dc = analyze_frequencies(t, noisy_signal, window='hann',
                                     remove_dc=False, zero_padding=4,
                                     plot=True, freq_range=(0, 20))

    # Demo 3: Window comparison
    print(f"\n" + "="*50)
    print("DEMO 3: Window Function Effects")
    print("="*50)

    windows = ['hann', 'hamming', 'blackman']
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, window in enumerate(windows):
        result = analyze_frequencies(t, noisy_signal, window=window,
                                   remove_dc=True, plot=False)

        axes[i].semilogy(result['frequencies'], result['power_spectrum'])
        axes[i].set_title(f'{window.capitalize()} Window')
        axes[i].set_xlabel(f'Frequency ({result["freq_unit"]})')
        axes[i].set_ylabel('Power')
        axes[i].grid(True, alpha=0.3)
        axes[i].set_xlim(0, 20)

        # Mark dominant frequency
        if len(result['dominant_frequencies']) > 0:
            peak_freq = result['dominant_frequencies'][0]
            axes[i].axvline(peak_freq, color='red', linestyle='--', alpha=0.7,
                           label=f'{peak_freq:.2f} MHz')
            axes[i].legend()

    plt.tight_layout()
    plt.show()

    # Demo 4: Power spectrum methods
    print(f"\n" + "="*50)
    print("DEMO 4: Power Spectrum Methods")
    print("="*50)

    psd_welch = power_spectrum(t, noisy_signal, method='welch', remove_dc=True, plot=True)
    print("Welch method completed")

    psd_periodogram = power_spectrum(t, noisy_signal, method='periodogram',
                                   remove_dc=True, plot=True)
    print("Periodogram method completed")

    print(f"\n" + "="*60)
    print("DEMO COMPLETED!")
    print("="*60)
    print("\nKey Points Demonstrated:")
    print("  * DC offset removal is crucial for clean spectra")
    print("  * Window functions reduce spectral leakage")
    print("  * Zero padding improves frequency resolution")
    print("  * Multiple methods available for power spectra")
    print("  * Automatic time unit detection (ns → MHz)")
    print(f"\nSimplified module ready for EPR frequency analysis!")


if __name__ == "__main__":
    np.random.seed(42)  # For reproducible results
    demo()