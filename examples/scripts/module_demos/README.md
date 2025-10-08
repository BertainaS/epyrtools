# EPyR Tools Module Demonstrations

This directory contains comprehensive demonstration scripts showing the functionality of every module and function in EPyR Tools.

## Overview

Each script demonstrates a specific aspect of the EPyR Tools package with practical examples using real data from `examples/data/` or synthetic data when real data is not available.

## Demonstration Scripts

| Script | Module | Description |
|--------|--------|-------------|
| `00_master_demo.py` | **Master** | Runs all demonstrations with options for selection |
| `01_eprload_demo.py` | **Data Loading** | Bruker file loading, parameter extraction, scaling |
| `02_eprplot_demo.py` | **EPR Plotting** | 1D/2D visualization, complex data handling |
| `03_lineshapes_demo.py` | **Lineshapes** | Mathematical functions, fitting, derivatives |
| `04_baseline_demo.py` | **Baseline** | 1D/2D correction, exclusion regions, quality |
| `05_signalprocessing_demo.py` | **Signal Processing** | FFT analysis, apodization windows, spectrograms |
| `06_physics_demo.py` | **Physics** | Constants, unit conversions, EPR calculations |
| `07_fair_demo.py` | **FAIR Data** | Data conversion, validation, metadata |
| `08_performance_demo.py` | **Performance** | Memory monitoring, caching, optimization |
| `09_config_demo.py` | **Configuration** | Settings management, environment variables |
| `10_cli_demo.py` | **CLI Interface** | Command-line tools and workflows |

## Usage

### Run All Demonstrations
```bash
python 00_master_demo.py
```

### Run Specific Demonstrations
```bash
python 00_master_demo.py --modules 1,2,3
```

### Skip Plot Generation (Faster)
```bash
python 00_master_demo.py --no-plots
```

### Run Individual Scripts
```bash
python 01_eprload_demo.py
python 02_eprplot_demo.py
# ... etc
```

## What Each Demo Shows

### 01. Data Loading (eprload)
- Loading 1D and 2D Bruker EPR files (.DSC/.DTA, .PAR/.SPC)
- Parameter extraction and interpretation
- Data scaling options (n=scans, P=power, G=gain, T=temperature)
- Error handling and format detection
- Synthetic data generation for testing

### 02. EPR Plotting (eprplot)
- Specialized 1D EPR spectrum plotting
- 2D color maps and waterfall plots
- Complex data visualization (real/imaginary)
- Custom colormaps and styling
- Time-domain data plotting

### 03. Lineshape Analysis
- All lineshape functions: Gaussian, Lorentzian, Voigtian, Pseudo-Voigt
- Derivative lineshapes (1st and 2nd order)
- Phase effects and complex lineshapes
- Multi-component signal fitting
- Mathematical validation and HWHM relationships

### 04. Baseline Correction
- 1D polynomial baseline correction
- 2D baseline correction with exclusion regions
- Signal exclusion and quality assessment
- Different polynomial orders and validation
- Interactive baseline selection

### 05. Signal Processing
- Apodization windows (Hamming, Hann, Blackman, Kaiser, etc.)
- FFT-based frequency analysis with DC removal
- Power spectral density (Welch and periodogram methods)
- Time-frequency analysis (spectrograms)
- Window effects on spectral resolution

### 06. Physics Constants
- 2022 CODATA physical constants (SI and CGS)
- EPR frequency ↔ magnetic field conversions
- Energy unit conversions (eV, MHz, cm⁻¹, K)
- Temperature effects and thermal energy
- g-factor calculations and EPR band frequencies

### 07. FAIR Data Conversion
- Bruker to FAIR format conversion (CSV, JSON, HDF5)
- Metadata extraction and standardization
- Data validation and compliance checking
- Batch processing capabilities
- Parameter mapping and quality reporting

### 08. Performance Optimization
- Memory monitoring and usage optimization
- LRU caching for repeated file access
- Optimized data loading with chunking
- Performance benchmarking
- Configuration-based performance tuning

### 09. Configuration Management
- Hierarchical configuration system
- Setting and getting configuration values
- Environment variable integration (EPYR_*)
- Configuration file operations
- Practical configuration scenarios

### 10. CLI Interface
- All 9 CLI commands with examples
- Interactive plotting with measurement tools
- Batch processing workflows
- Configuration management from command line
- Help system and error handling

## Generated Output

Each demonstration generates:
- **Console output** showing function results and explanations
- **Plot files** (PNG format) demonstrating visualizations
- **Data files** (CSV, JSON, etc.) showing export capabilities
- **Configuration files** for testing persistence

## Key Features Demonstrated

✅ **Complete package coverage** - Every module and major function
✅ **Real data examples** - Uses actual EPR data when available
✅ **Synthetic data** - Creates realistic test data for all scenarios
✅ **Error handling** - Shows graceful degradation and error recovery
✅ **Best practices** - Demonstrates proper usage patterns
✅ **Integration** - Shows how modules work together
✅ **Performance** - Memory-efficient handling of large datasets
✅ **FAIR compliance** - Modern data management standards

## Requirements

- Python 3.8+
- NumPy, SciPy, Matplotlib (required)
- Pandas, H5py (optional, for full functionality)
- psutil (optional, for memory monitoring)
- tkinter (optional, for GUI components)

## Notes

- Demonstrations are designed to work without user interaction
- Plot generation can be disabled with `--no-plots` for faster execution
- Real EPR data files are used when available in `examples/data/`
- Synthetic data is generated when real files are not found
- All scripts include comprehensive error handling
- Console output provides educational explanations of each feature

This demonstration series provides a complete overview of EPyR Tools capabilities and serves as both tutorial and reference for users learning the package.