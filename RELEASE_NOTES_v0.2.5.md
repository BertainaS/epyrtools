# EPyR Tools v0.2.5 Release Notes

**Release Date:** October 2025
**Release Type:** Feature Enhancement

## Overview

Version 0.2.5 introduces enhanced data handling capabilities and advanced signal processing features for time-domain EPR analysis. This release focuses on improving the flexibility of data loading and adding comprehensive 2D FFT analysis tools for Rabi oscillations and HYSCORE measurements.

## New Features

### Enhanced Data Loading (`eprload`)

The `eprload()` function now provides more control over signal component extraction:

#### New `return_type` Parameter
```python
# Return signal as-is (default, backward compatible)
x, y, params, filepath = eprload("data.dsc", return_type="default")

# Extract only real part (useful for CW-EPR)
x, y_real, params, filepath = eprload("data.dsc", return_type="real")

# Extract only imaginary part (useful for out-of-phase signals)
x, y_imag, params, filepath = eprload("data.dsc", return_type="imag")
```

#### Updated Default Behavior
- `plot_if_possible` parameter now defaults to `False` instead of `True`
- This change provides cleaner programmatic use without unexpected plot windows
- Users wanting automatic plotting should explicitly set `plot_if_possible=True`

**Migration Guide:**
```python
# Old code (v0.2.0)
x, y, params, filepath = eprload("data.dsc")  # Would plot automatically

# New code (v0.2.5) - equivalent behavior
x, y, params, filepath = eprload("data.dsc", plot_if_possible=True)

# Recommended for scripts (no plotting)
x, y, params, filepath = eprload("data.dsc")  # No automatic plot
```

### Advanced 2D FFT Analysis

Complete frequency analysis module for time-domain EPR measurements with two processing modes:

#### Row-by-Row 1D FFT Mode
For 2D datasets where each trace needs independent FFT processing (e.g., 2D Rabi oscillations):

```python
from epyr.signalprocessing import analyze_frequencies_2d

# Load 2D Rabi data
x_2d, y_2d, params, _ = eprload("rabi_2d.DTA")

# Row-by-row FFT processing
fq, field_axis, spectrum, info = analyze_frequencies_2d(
    x_2d[0],           # Time axis
    y_2d,              # 2D signal data
    mode='row_by_row',
    window='hann',
    zero_padding=2,
    remove_dc=True,
    plot_result=True   # Optional visualization
)

# Access results
frequencies = fq              # Frequency axis (MHz, kHz, or Hz)
spectrum_2d = spectrum        # 2D magnitude spectrum
freq_unit = info['freq_unit'] # Automatic unit detection
```

#### Full 2D FFT Mode
For true 2D time-domain measurements (e.g., HYSCORE):

```python
# Load HYSCORE data
x_hyscore, y_hyscore, params, _ = eprload("hyscore.DTA")

# Full 2D FFT
fq1, fq2, spectrum_2d, info = analyze_frequencies_2d(
    (x_hyscore[0], x_hyscore[1]),  # Both time axes
    y_hyscore,
    mode='full_2d',
    window='hann',
    zero_padding=2,
    plot_result=True
)

# Access results
freq_axis_1 = fq1
freq_axis_2 = fq2
spectrum_magnitude = spectrum_2d
phase_spectrum = info['phase_spectrum']
```

#### Key Features

**Automatic Time Unit Detection:**
- Automatically detects time units: nanoseconds (ns), microseconds (μs), milliseconds (ms), seconds (s)
- Converts to appropriate frequency units: MHz, kHz, or Hz
- No manual unit specification required

**Signal Processing Options:**
- DC offset removal (essential for clean spectra)
- Apodization windows: Hann, Hamming, Blackman, Kaiser, Gaussian
- Zero padding for improved frequency resolution
- Centered frequency spectrum (positive and negative frequencies via fftshift)

**Comprehensive Visualization:**
Four-panel display when `plot_result=True`:
1. Original time-domain signal
2. FFT magnitude (linear scale)
3. FFT magnitude (log scale for dynamic range)
4. Phase spectrum

**Return Values:**
- Tuple format: `(frequencies, axis2, spectrum, info_dict)`
- No automatic plotting by default (`plot_result=False`)
- Phase spectrum available in `info` dictionary

## Code Quality Improvements

### Refactored Signal Processing Module

The `frequency_analysis.py` module has been refactored to eliminate code duplication:

**New Helper Functions:**
- `_detect_time_units()`: Single source for time unit detection (eliminates 4 duplicate implementations)
- `_convert_to_display_freq()`: Unified frequency conversion
- `_remove_dc_offset()`: Standardized DC offset removal
- `_apply_window()`: Centralized apodization windowing

**Benefits:**
- Reduced code duplication by approximately 150 lines
- Single source of truth for common operations
- Easier maintenance and bug fixes
- Consistent behavior across all analysis functions

**Functions Updated:**
- `analyze_frequencies()` - 1D FFT analysis
- `power_spectrum()` - Power spectral density
- `spectrogram_analysis()` - Time-frequency analysis
- `_analyze_2d_row_by_row()` - 2D row-by-row FFT
- `_analyze_2d_full()` - Full 2D FFT

## Testing & Validation

### Test Coverage
- All new features tested with synthetic EPR data
- 2D FFT modes validated with Rabi oscillation and HYSCORE patterns
- Real/imaginary component extraction verified with complex signals
- Backward compatibility confirmed for existing workflows

### Verification Tests
```python
# Test 1D FFT
result = analyze_frequencies(time, signal, window='hann', plot=False)
assert 'frequencies' in result
assert 'power_spectrum' in result

# Test 2D row-by-row FFT
fq, axis2, spectrum, info = analyze_frequencies_2d(
    time_data, signal_2d, mode='row_by_row', plot_result=False)
assert fq.shape[0] > 0
assert spectrum.shape[0] == signal_2d.shape[0]

# Test return_type parameter
x, y_real, params, _ = eprload("data.dsc", return_type="real")
assert np.isrealobj(y_real)
```

## Backward Compatibility

### Breaking Changes
**None.** All changes are backward compatible.

### Default Behavior Changes
1. `eprload()` no longer plots automatically by default
   - **Impact**: Scripts will no longer show unexpected plot windows
   - **Migration**: Add `plot_if_possible=True` if automatic plotting is desired

2. `analyze_frequencies_2d()` returns tuple instead of dict
   - **Impact**: More Pythonic API, cleaner unpacking
   - **Migration**: None required for new code, function is new in v0.2.5

## Known Issues

None reported.

## Dependencies

No new dependencies added. All features use existing core dependencies:
- numpy >= 1.20.0
- scipy >= 1.7.0
- matplotlib >= 3.3.0

## Performance

- 2D FFT processing optimized for memory efficiency
- Row-by-row mode processes each trace independently
- Full 2D mode uses scipy's efficient fft2 implementation
- Helper function refactoring has no measurable performance impact

## Documentation Updates

- README updated with v0.2.5 features
- API documentation updated for new parameters
- Examples added for 2D FFT analysis
- Tutorial notebooks compatible with new features

## Contributors

- **Sylvain Bertaina** - Lead Developer & Maintainer
  - Email: sylvain.bertaina@cnrs.fr
  - Affiliation: Magnetism Group (MAG), IM2NP Laboratory, CNRS

## Installation

### Upgrade from v0.2.0
```bash
pip install --upgrade epyr-tools
```

### Fresh Installation
```bash
pip install epyr-tools==0.2.5
```

### Development Installation
```bash
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools
git checkout v0.2.5
pip install -e ".[dev,docs]"
```

## Next Steps

**Planned for v0.2.6:**
- Enhanced baseline correction with 2D support
- Additional lineshape fitting capabilities
- Improved documentation with more examples

## Feedback & Support

- **Issues**: [GitHub Issues](https://github.com/BertainaS/epyrtools/issues)
- **Discussions**: [GitHub Discussions](https://github.com/BertainaS/epyrtools/discussions)
- **Email**: sylvain.bertaina@cnrs.fr

---

**EPyR Tools v0.2.5** - Enhanced signal processing for modern EPR analysis
