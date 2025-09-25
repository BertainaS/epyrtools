"""
Signal Processing Module for EPR Time-Domain Data
================================================

This module provides comprehensive tools for frequency analysis of time-dependent EPR signals,
optimized for pulse EPR experiments including Rabi oscillations, DEER spectroscopy,
and other advanced time-domain techniques.

Key Features
------------
- FFT-based frequency analysis with automatic DC offset removal
- Advanced apodization windows for spectral leakage reduction
- Power spectral density analysis using Welch and periodogram methods
- Time-frequency analysis with spectrograms
- Automatic time unit detection (ns, μs, ms, s)
- Zero padding for enhanced frequency resolution

Main Functions
--------------
analyze_frequencies : FFT-based frequency analysis with comprehensive visualization
power_spectrum : Power spectral density using Welch or periodogram methods
spectrogram_analysis : Time-frequency analysis for evolving spectral content
apowin : Apodization window generation with multiple window types

Examples
--------
Basic frequency analysis of Rabi data::

    from epyr import eprload
    from epyr.signalprocessing import analyze_frequencies

    # Load time-domain EPR data
    time, signal, params, _ = eprload('rabi_data.DTA')

    # Analyze frequencies with DC removal and Hann window
    result = analyze_frequencies(time, signal, window='hann',
                               remove_dc=True, plot=True)

    print(f"Dominant frequency: {result['dominant_frequencies'][0]:.3f} MHz")

Advanced power spectrum analysis::

    from epyr.signalprocessing import power_spectrum

    # Welch method for noise reduction
    psd_result = power_spectrum(time, signal, method='welch',
                              window='hann', plot=True)

Notes
-----
This module is specifically designed for EPR time-domain signal analysis and includes
optimizations for typical EPR data characteristics including proper handling of complex
signals, automatic unit detection, and spectroscopic conventions.
"""

# Import apodization windows
from .apowin import apowin, window_comparison, frequency_response_demo, apply_window_demo

# Import frequency analysis tools
from .frequency_analysis import (
    analyze_frequencies, power_spectrum, spectrogram_analysis
)

__all__ = [
    # Apodization windows
    'apowin', 'window_comparison', 'frequency_response_demo', 'apply_window_demo',

    # Frequency analysis (simplified FFT-based functions)
    'analyze_frequencies', 'power_spectrum', 'spectrogram_analysis'
]

__version__ = "0.2.0"
__author__ = "EPyR Tools Development Team"