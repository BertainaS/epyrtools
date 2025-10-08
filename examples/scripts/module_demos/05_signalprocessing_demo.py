#!/usr/bin/env python3
"""
EPyR Tools Demo 05: Signal Processing and Frequency Analysis
============================================================

This script demonstrates the signal processing capabilities for time-domain EPR data.
The signalprocessing module provides tools for FFT analysis, apodization windows,
and power spectral density analysis.

Functions demonstrated:
- analyze_frequencies() - FFT-based frequency analysis with apodization
- power_spectrum() - Power spectral density using Welch and periodogram methods
- apowin() - Apodization window generation
- spectrogram_analysis() - Time-frequency analysis
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_apodization_windows():
    """Demonstrate apodization window functions."""
    print("=== EPyR Tools Signal Processing Demo - Apodization Windows ===")
    print()

    print("1. Apodization window generation:")
    print("-" * 35)

    # Generate different window types
    n_points = 256
    window_types = ['hamming', 'hann', 'blackman', 'bartlett', 'cosine', 'welch']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, window_type in enumerate(window_types):
        window = epyr.signalprocessing.apowin(window_type, n_points)

        axes[i].plot(window, linewidth=2)
        axes[i].set_title(f"{window_type.capitalize()} Window")
        axes[i].set_xlabel("Sample")
        axes[i].set_ylabel("Amplitude")
        axes[i].grid(True, alpha=0.3)
        axes[i].set_ylim(0, 1.1)

        print(f"  {window_type:<10}: Peak = {window.max():.3f}, Area = {window.sum():.1f}")

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "05_apodization_windows.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 05_apodization_windows.png")
    print()


def demo_kaiser_gaussian_windows():
    """Demonstrate parametric windows (Kaiser, Gaussian)."""
    print("2. Parametric windows (Kaiser and Gaussian):")
    print("-" * 45)

    n_points = 256

    # Kaiser windows with different beta values
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    beta_values = [2, 4, 6, 8, 10]
    for beta in beta_values:
        kaiser_window = epyr.signalprocessing.apowin('kaiser', n_points, alpha=beta)
        ax1.plot(kaiser_window, label=f'β = {beta}', linewidth=2)

    ax1.set_title("Kaiser Windows")
    ax1.set_xlabel("Sample")
    ax1.set_ylabel("Amplitude")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Gaussian windows with different alpha values
    alpha_values = [0.3, 0.4, 0.5, 0.6, 0.8]
    for alpha in alpha_values:
        gaussian_window = epyr.signalprocessing.apowin('gaussian', n_points, alpha=alpha)
        ax2.plot(gaussian_window, label=f'α = {alpha}', linewidth=2)

    ax2.set_title("Gaussian Windows")
    ax2.set_xlabel("Sample")
    ax2.set_ylabel("Amplitude")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "05_parametric_windows.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 05_parametric_windows.png")
    print()


def demo_time_domain_analysis():
    """Demonstrate frequency analysis of time-domain EPR data."""
    print("3. Time-domain EPR frequency analysis:")
    print("-" * 38)

    data_dir = Path(__file__).parent.parent.parent / "data"

    # Look for time-domain files
    time_files = [
        "2024_08_CaWO4171Yb_rabi_6K_6724G_18dB.DSC",
        "Rabi2D_GdCaWO4_6dB_3770G_2.DSC"
    ]

    for filename in time_files:
        file_path = data_dir / filename
        if file_path.exists():
            print(f"Loading real time-domain data: {filename}")
            try:
                time_data, signal_data, params, filepath = epyr.eprload(str(file_path), plot_if_possible=False)

                print(f"  Data shape: {signal_data.shape}")
                print(f"  Time range: {time_data.min():.1f} to {time_data.max():.1f} ns")
                print(f"  Signal is complex: {np.iscomplexobj(signal_data)}")

                # Use real part if complex
                if np.iscomplexobj(signal_data):
                    signal_for_analysis = signal_data.real
                else:
                    signal_for_analysis = signal_data

                # Frequency analysis with different windows
                print("  Performing frequency analysis...")
                result = epyr.signalprocessing.analyze_frequencies(
                    time_data, signal_for_analysis,
                    window='hann', remove_dc=True, plot=True
                )

                plt.savefig(Path(__file__).parent / f"05_freq_analysis_{Path(filename).stem}.png",
                           dpi=150, bbox_inches='tight')
                plt.show()
                print(f"  Saved as: 05_freq_analysis_{Path(filename).stem}.png")

                if 'dominant_frequencies' in result:
                    print(f"  Dominant frequencies: {result['dominant_frequencies'][:3]} MHz")
                break

            except Exception as e:
                print(f"  Error loading {filename}: {e}")
                continue
    else:
        print("No real time-domain data available, creating synthetic data...")
        create_synthetic_time_domain_analysis()

    print()


def create_synthetic_time_domain_analysis():
    """Create synthetic time-domain data for demonstration."""
    print("  Creating synthetic Rabi oscillation data...")

    # Generate synthetic Rabi oscillations with multiple frequencies
    t = np.linspace(0, 1000, 500)  # Time in ns
    dt = t[1] - t[0]

    # Multiple frequency components
    rabi_freq1 = 0.008  # MHz (8 kHz)
    rabi_freq2 = 0.025  # MHz (25 kHz)
    T2 = 300  # Coherence time in ns

    # Create signal with two Rabi frequencies and T2 decay
    signal = (0.7 * np.cos(2 * np.pi * rabi_freq1 * t) +
              0.3 * np.cos(2 * np.pi * rabi_freq2 * t)) * np.exp(-t/T2)

    # Add realistic noise
    signal += 0.05 * np.random.normal(size=len(t))

    print(f"  Time points: {len(t)}")
    print(f"  Time step: {dt:.1f} ns")
    print(f"  Expected frequencies: {rabi_freq1*1000:.1f} and {rabi_freq2*1000:.1f} kHz")

    # Frequency analysis
    result = epyr.signalprocessing.analyze_frequencies(
        t, signal, window='hann', remove_dc=True, plot=True,
        freq_range=(0, 0.1)  # Focus on low frequencies
    )

    plt.savefig(Path(__file__).parent / "05_synthetic_rabi_analysis.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 05_synthetic_rabi_analysis.png")

    if 'dominant_frequencies' in result:
        print(f"  Detected frequencies: {[f*1000 for f in result['dominant_frequencies'][:3]]} kHz")


def demo_power_spectrum():
    """Demonstrate power spectral density analysis."""
    print("4. Power spectral density analysis:")
    print("-" * 35)

    # Create noisy synthetic signal for PSD analysis
    t = np.linspace(0, 2000, 1000)  # 2 μs, 1000 points
    dt = t[1] - t[0]

    # Signal with multiple frequency components and noise
    freq1, freq2, freq3 = 0.015, 0.042, 0.087  # MHz
    signal = (np.sin(2 * np.pi * freq1 * t) +
              0.5 * np.sin(2 * np.pi * freq2 * t) +
              0.3 * np.sin(2 * np.pi * freq3 * t))

    # Add colored noise
    noise = np.random.normal(size=len(t))
    # Apply simple low-pass filter to noise for more realistic EPR noise
    from scipy import signal as scipy_signal
    b, a = scipy_signal.butter(4, 0.3)
    noise = scipy_signal.filtfilt(b, a, noise)
    signal += 0.4 * noise

    print(f"  Synthetic signal with frequencies: {freq1*1000:.1f}, {freq2*1000:.1f}, {freq3*1000:.1f} kHz")

    # Compare different PSD methods
    try:
        # Welch method
        psd_welch = epyr.signalprocessing.power_spectrum(
            t, signal, method='welch', window='hann', plot=True
        )

        plt.savefig(Path(__file__).parent / "05_power_spectrum_welch.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 05_power_spectrum_welch.png")

        # Periodogram method
        psd_period = epyr.signalprocessing.power_spectrum(
            t, signal, method='periodogram', window='hann', plot=True
        )

        plt.savefig(Path(__file__).parent / "05_power_spectrum_periodogram.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 05_power_spectrum_periodogram.png")

    except Exception as e:
        print(f"  Error in power spectrum analysis: {e}")
        print("  This feature may require the latest version of the signalprocessing module")

    print()


def demo_window_effects():
    """Demonstrate the effect of different windows on spectral analysis."""
    print("5. Window effects on frequency analysis:")
    print("-" * 40)

    # Create test signal with closely spaced frequencies
    t = np.linspace(0, 1000, 512)
    freq1, freq2 = 0.050, 0.052  # Very close frequencies (50 and 52 kHz)
    signal = np.sin(2 * np.pi * freq1 * t) + np.sin(2 * np.pi * freq2 * t)
    signal += 0.1 * np.random.normal(size=len(t))

    print(f"  Test signal with close frequencies: {freq1*1000:.0f} and {freq2*1000:.0f} kHz")

    # Compare different windows
    windows = ['none', 'hann', 'hamming', 'blackman']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, window in enumerate(windows):
        window_param = None if window == 'none' else window

        try:
            result = epyr.signalprocessing.analyze_frequencies(
                t, signal, window=window_param, remove_dc=True, plot=False
            )

            # Extract frequency and magnitude data for plotting
            if 'frequencies' in result and 'magnitude' in result:
                freq = result['frequencies']
                mag = result['magnitude']

                axes[i].plot(freq * 1000, mag)  # Convert to kHz
                axes[i].set_title(f"{window.capitalize() if window != 'none' else 'No'} Window")
                axes[i].set_xlabel("Frequency (kHz)")
                axes[i].set_ylabel("Magnitude")
                axes[i].set_xlim(40, 65)  # Focus on signal region
                axes[i].grid(True, alpha=0.3)

        except Exception as e:
            print(f"  Error with {window} window: {e}")
            axes[i].text(0.5, 0.5, f"Error: {e}", ha='center', va='center',
                        transform=axes[i].transAxes)
            axes[i].set_title(f"{window.capitalize()} Window - Error")

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "05_window_comparison.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 05_window_comparison.png")
    print("  Note: Different windows provide different frequency resolution vs. leakage tradeoffs")
    print()


def demo_spectrogram():
    """Demonstrate time-frequency analysis with spectrograms."""
    print("6. Time-frequency analysis (Spectrograms):")
    print("-" * 42)

    # Create signal with time-varying frequency (chirp)
    t = np.linspace(0, 1000, 1000)  # 1000 ns

    # Frequency chirp from 10 to 100 kHz
    f0, f1 = 0.01, 0.1  # Start and end frequencies in MHz
    # Linear chirp
    signal = np.sin(2 * np.pi * (f0 + (f1 - f0) * t / t[-1]) * t)

    # Add some amplitude modulation and noise
    signal *= (1 + 0.3 * np.sin(2 * np.pi * 0.005 * t))  # 5 kHz AM
    signal += 0.2 * np.random.normal(size=len(t))

    print(f"  Chirp signal: {f0*1000:.0f} to {f1*1000:.0f} kHz over {t[-1]:.0f} ns")

    try:
        # Spectrogram analysis
        result = epyr.signalprocessing.spectrogram_analysis(
            t, signal, window='hann', plot=True
        )

        plt.savefig(Path(__file__).parent / "05_spectrogram_analysis.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 05_spectrogram_analysis.png")

    except Exception as e:
        print(f"  Error in spectrogram analysis: {e}")
        print("  Creating basic spectrogram with scipy...")

        # Fallback: basic spectrogram with scipy
        from scipy import signal as scipy_signal

        fs = 1 / (t[1] - t[0]) * 1000  # Convert to MHz sampling rate
        f, t_spec, Sxx = scipy_signal.spectrogram(signal, fs, window='hann', nperseg=128)

        plt.figure(figsize=(10, 6))
        plt.pcolormesh(t_spec, f, 10 * np.log10(Sxx), shading='gouraud')
        plt.ylabel('Frequency (MHz)')
        plt.xlabel('Time (ns)')
        plt.title('Spectrogram - Time-Frequency Analysis')
        plt.colorbar(label='Power (dB)')
        plt.ylim(0, 0.15)

        plt.savefig(Path(__file__).parent / "05_basic_spectrogram.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 05_basic_spectrogram.png")

    print()


def main():
    """Run all signal processing demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_apodization_windows()
    demo_kaiser_gaussian_windows()
    demo_time_domain_analysis()
    demo_power_spectrum()
    demo_window_effects()
    demo_spectrogram()

    print("=== Signal Processing Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- apowin() generates various apodization windows for spectral leakage reduction")
    print("- analyze_frequencies() performs FFT analysis with proper DC removal and windowing")
    print("- power_spectrum() provides PSD analysis using Welch and periodogram methods")
    print("- spectrogram_analysis() reveals time-varying frequency content")
    print("- Window choice affects frequency resolution vs. spectral leakage tradeoffs")
    print("- Proper apodization is crucial for clean EPR frequency analysis")
    print()
    print("Generated plot files:")
    for plot_file in sorted(output_dir.glob("05_*.png")):
        print(f"  - {plot_file.name}")


if __name__ == "__main__":
    main()