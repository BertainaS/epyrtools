"""
Signal Processing Module for EPR Time-Domain Data
================================================

This module provides comprehensive tools for frequency analysis of time-dependent EPR signals,
optimized for pulse EPR experiments including Rabi oscillations, DEER spectroscopy, HYSCORE,
and other advanced time-domain techniques.

Key Features
------------
- 1D FFT analysis with automatic DC offset removal and windowing
- 2D FFT analysis with two processing modes:
  - Row-by-row 1D FFT for 2D Rabi oscillations and similar experiments
  - Full 2D FFT for HYSCORE-type measurements with cross-peak analysis
- Advanced apodization windows for spectral leakage reduction
- Power spectral density analysis using Welch and periodogram methods
- Time-frequency analysis with spectrograms
- Automatic time unit detection (ns, μs, ms, s)
- Zero padding for enhanced frequency resolution

Main Functions
--------------
analyze_frequencies : 1D FFT-based frequency analysis with comprehensive visualization
analyze_frequencies_2d : 2D FFT analysis with row-by-row or full 2D processing modes
power_spectrum : Power spectral density using Welch or periodogram methods
spectrogram_analysis : Time-frequency analysis for evolving spectral content
apowin : Apodization window generation with multiple window types

Examples
--------
Basic 1D frequency analysis of Rabi data::

    from epyr import eprload
    from epyr.signalprocessing import analyze_frequencies

    # Load time-domain EPR data
    time, signal, params, _ = eprload('rabi_data.DTA')

    # Analyze frequencies with DC removal and Hann window
    result = analyze_frequencies(time, signal, window='hann',
                               remove_dc=True, plot=True)

    print(f"Dominant frequency: {result['dominant_frequencies'][0]:.3f} MHz")

2D frequency analysis - Row-by-row mode (2D Rabi oscillations)::

    from epyr.signalprocessing import analyze_frequencies_2d

    # Load 2D Rabi data
    x_2d, y_2d, params, _ = eprload('rabi_2d.DTA')

    # Process each trace independently with 1D FFT
    result = analyze_frequencies_2d(x_2d[0], y_2d, mode='row_by_row',
                                   window='hann', plot=True)

    print(f"Power spectrum shape: {result['power_spectrum'].shape}")

2D frequency analysis - Full 2D FFT (HYSCORE)::

    from epyr.signalprocessing import analyze_frequencies_2d

    # Load HYSCORE data
    x_hyscore, y_hyscore, params, _ = eprload('hyscore.DTA')

    # Apply full 2D FFT for cross-peak analysis
    result = analyze_frequencies_2d((x_hyscore[0], x_hyscore[1]), y_hyscore,
                                   mode='full_2d', window='hann', plot=True)

    print(f"2D power spectrum shape: {result['power_spectrum'].shape}")

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

The 2D FFT capabilities enable analysis of:
- 2D Rabi oscillations (row-by-row mode)
- HYSCORE spectra (full 2D mode)
- Other multi-dimensional time-domain EPR experiments
"""

# Import apodization windows
from .apowin import (
    apowin,
    apply_window_demo,
    frequency_response_demo,
    window_comparison,
)

# Import frequency analysis tools
from .frequency_analysis import (
    analyze_frequencies,
    analyze_frequencies_2d,
    power_spectrum,
    spectrogram_analysis,
)

__all__ = [
    # Apodization windows
    "apowin",
    "window_comparison",
    "frequency_response_demo",
    "apply_window_demo",
    # Frequency analysis (FFT-based functions)
    "analyze_frequencies",  # 1D FFT analysis
    "analyze_frequencies_2d",  # 2D FFT analysis (row-by-row or full 2D)
    "power_spectrum",  # Power spectral density
    "spectrogram_analysis",  # Time-frequency analysis
]

__version__ = "0.3.1"
__author__ = "EPyR Tools Development Team"
