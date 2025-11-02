# EPyR Tools: Electron Paramagnetic Resonance (EPR) Tools in Python

<p align="center">
  <img src="Epyrtools_logo.jpg" alt="EPyR Tools Logo" height="240"/>
</p>

<p align="center">
  <img src="LogoL.png" alt="Institution Logo Left" height="120"/>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="LogoR.png" alt="Institution Logo Right" height="120"/>
</p>

| License | Tests | Documentation | Version |
|---------|-------|---------------|---------|
| [![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) | ![Tests Passing](https://img.shields.io/badge/tests-100%2B%20passed-brightgreen) | [![Documentation](https://img.shields.io/badge/docs-comprehensive-blue)](docs/) | ![Version](https://img.shields.io/badge/version-0.3.0-blue) |

**EPyR Tools** is a comprehensive Python package for Electron Paramagnetic Resonance (EPR) spectroscopy data analysis. It provides a complete toolkit for loading, processing, analyzing, and visualizing EPR data from Bruker spectrometers, with a focus on FAIR (Findable, Accessible, Interoperable, and Reusable) data principles.

From basic data loading to advanced quantitative analysis, EPyR Tools offers professional-grade capabilities for EPR researchers, with comprehensive documentation and interactive tutorials.

## Key Features

### **Data Loading & Formats**
- **Bruker File Support:** Load BES3T (.dta/.dsc) and ESP/WinEPR (.par/.spc) files seamlessly
- **Automatic Format Detection:** Smart file format recognition and parameter extraction
- **FAIR Data Conversion:** Export to CSV, JSON, and HDF5 formats with complete metadata
- **Batch Processing:** Handle multiple files efficiently with parallel processing
- **Plugin Architecture:** Extensible system for custom file formats and processing

### **Advanced Analysis**
- **Modern Baseline Correction:** Streamlined polynomial correction with interactive region selection
- **Peak Detection:** Automatic identification of EPR spectral features
- **g-Factor Calculations:** Precise electronic g-factor determination with field calibration
- **Hyperfine Analysis:** Pattern recognition and coupling constant extraction
- **Quantitative Integration:** Single and double integration for spin quantification
- **Lineshape Analysis:** Comprehensive suite of EPR lineshape functions (Gaussian, Lorentzian, Voigt, pseudo-Voigt)

### **Command Line Interface**
- **Complete CLI Suite:** 9 professional commands for all EPR workflows
- **Interactive Plotting:** `epyr-plot` with measurement tools for precise delta x/y analysis
- **Batch Processing:** `epyr-batch-convert` for high-throughput data processing  
- **Data Validation:** `epyr-validate` with FAIR compliance checking
- **Configuration Management:** `epyr-config` for system-wide settings
- **Interactive Tools:** `epyr-isotopes` GUI and system information display

### **Performance & Quality**
- **Memory Optimization:** Intelligent caching and memory management for large datasets
- **Parallel Processing:** Multi-core support for batch operations
- **Quality Assurance:** Comprehensive testing suite with 90+ tests
- **Code Standards:** Professional development workflow with pre-commit hooks
- **Security:** Automated vulnerability scanning and safe defaults

### **Visualization & Plotting**
- **Interactive Measurement Tools:** Click-to-measure delta x/y distances with real-time feedback
- **2D Spectral Maps:** Professional publication-quality EPR plots  
- **macOS Optimized:** Smooth interactive plotting with TkAgg backend
- **Customizable Styling:** Flexible plot configuration for different EPR experiments
- **Export Options:** High-resolution outputs for publications

### **Learning & Documentation**
- **Interactive Tutorials:** 4 comprehensive Jupyter notebooks (beginner → advanced)
- **Complete API Documentation:** Professional Sphinx-generated docs
- **Example Scripts:** Ready-to-use Python automation scripts including lineshape analysis
- **Best Practices Guide:** EPR analysis workflows and quality assessment

### **EPR-Specific Tools**
- **Physical Constants:** Comprehensive EPR-relevant constants library
- **Isotope Database:** Nuclear properties and magnetic parameters
- **Field-Frequency Conversion:** Precise EPR measurement calculations
- **Spectrometer Support:** Optimized for modern Bruker EPR systems

### **Professional Documentation**
- **[Complete Documentation](docs/)**: Comprehensive guides and API reference
- **[User Guide](docs/user_guide.md)**: Step-by-step tutorials and workflows
- **[CLI Reference](docs/cli_reference.md)**: Command-line interface documentation
- **[System Architecture](docs/README.md)**: Core modules and plugin system

## What's New in v0.3.0

**EPyR Tools v0.3.0** represents a **major code quality enhancement** with professional logging infrastructure:

### Professional Logging System
- **Centralized Logging**: Migrated 398 print statements to structured logging (95.7% coverage)
  - Professional logging infrastructure across 18 core modules
  - Hierarchical log levels: `DEBUG`, `INFO`, `WARNING`, `ERROR`
  - Timestamped messages with module names for better traceability
  - Consistent log format throughout the entire codebase

- **Enhanced Debugging**: Better production-ready monitoring capabilities
  - Configurable logging levels via standard Python logging
  - Log to files for permanent records
  - Integration with logging frameworks (Sentry, Logstash, etc.)
  - Full call chain tracing with module information

- **Backward Compatible**: 100% compatible with existing code
  - No breaking changes - all functionality preserved
  - Default logging configured to match previous print() behavior
  - All 100+ tests passing with no regressions
  - Optional custom logging configuration

### Modules Updated
Core modules refactored with professional logging:
- **Data Loading**: eprload, eprplot, cli
- **Signal Processing**: frequency_analysis, apowin
- **Physics**: conversions, units, constants
- **Baseline Correction**: baseline package (4 modules)
- **FAIR Conversion**: conversion, exporters, data_processing
- **Lineshape Analysis**: fitting, convspec
- **GUI**: isotope_gui

### Configuration Example
```python
import logging

# Configure EPyR logging (optional)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='epyr_analysis.log'
)

# Use EPyR Tools normally
from epyr import eprload
x, y, params, filepath = eprload("data.DTA")
# Logging messages written to epyr_analysis.log
```

### Benefits
- **Production Ready**: Professional logging standards for real-world deployments
- **Better Debugging**: Detailed trace information with timestamps and module names
- **Maintainability**: Single logging configuration point, consistent patterns
- **Flexibility**: Control output verbosity via standard Python logging configuration

**See Also:**
- [RELEASE_NOTES_v0.3.0.md](RELEASE_NOTES_v0.3.0.md) - Complete v0.3.0 details
- [RELEASE_NOTES_v0.2.5.md](RELEASE_NOTES_v0.2.5.md) - Enhanced data handling features
- [docs/changelog.rst](docs/changelog.rst) - Full version history

---

## Installation

### Prerequisites
- Python 3.8 or higher
- NumPy, matplotlib, pandas, h5py (automatically installed)

### Quick Install
```bash
pip install epyr-tools
```

### Development Installation
```bash
# Clone the repository
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools

# Install with development dependencies
pip install -e ".[dev,docs]"

# Set up pre-commit hooks
pre-commit install
```

### Verification
```bash
# Test installation
epyr --help
epyr-info

# Run test suite
make test
```

### Alternative: Manual Dependencies
```bash
# Clone the repository
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools

# Install dependencies manually
pip install -r requirements.txt

# Optional: Install development tools
pip install -r requirements-dev.txt
```

### Quick Test
```bash
# Verify installation
python -c "import epyr; print('EPyR Tools successfully installed!')"
```

## Getting Started

### Quick Demo: Loading EPR Data

![EPyR Tools eprload demo](eprload_simple_demo.gif)

### 1. Loading Data

The primary function for loading data is `epyr.eprload()`. It can be called with a file path or without arguments to open a file dialog.

```python
import epyr.eprload as eprload

# Open a file dialog to select a .dta, .dsc, .spc, or .par file
# The script will automatically plot the data.
x, y, params, file_path = eprload.eprload()

# Or, specify a file path directly:
# x, y, params, file_path = eprload.eprload('path/to/your/data.dsc')
```

### 2. Converting to FAIR Formats

Use the `epyr.fair` module to convert your Bruker files into more accessible formats (`.csv`, `.json`, `.h5`).

```python
from epyr.fair import convert_bruker_to_fair

# Opens a file dialog to select an input file and saves the converted
# files in the same directory.
convert_bruker_to_fair()

# You can also specify input and output paths:
# convert_bruker_to_fair('path/to/data.dsc', output_dir='path/to/output')
```

### 3. Modern Baseline Correction

The `epyr.baseline` module provides streamlined functions for correcting EPR baselines with interactive region selection. Compatible with `eprload()` data format.

```python
import epyr
import numpy as np

# Load EPR data using eprload
x, y, params, filepath = epyr.eprload("data.dsc")

# Simple automatic correction
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params)

# With manual region exclusion (e.g., exclude signal regions)
signal_regions = [(3340, 3360), (3380, 3400)]  # mT
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, 
    manual_regions=signal_regions,
    region_mode='exclude',
    order=2
)

# Interactive region selection for complex spectra
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params,
    interactive=True
)
```

### 4. Lineshape Analysis

The `epyr.lineshapes` module provides comprehensive EPR lineshape functions with advanced capabilities including derivatives, phase rotation, and convolution.

```python
from epyr.lineshapes import Lineshape, gaussian, lorentzian, voigtian
import numpy as np
import matplotlib.pyplot as plt

# Create magnetic field range
B = np.linspace(-10, 10, 1000)  # mT

# Method 1: Direct functions with advanced options
gauss = gaussian(B, center=0, width=4.0, derivative=0)  # Absorption
gauss_1st = gaussian(B, center=0, width=4.0, derivative=1)  # 1st derivative
lorentz = lorentzian(B, center=0, width=4.0, phase=0.0)  # Pure absorption

# True Voigt profiles (convolution of Gaussian and Lorentzian)
voigt = voigtian(B, center=0, sigma=2.0, gamma=2.0)

# Method 2: Unified Lineshape class for consistent API
shape = Lineshape('pseudo_voigt', width=4.0, alpha=0.5)  # 50/50 mix
mixed = shape(B, center=0)

# Plot comparison
plt.plot(B, gauss, label='Gaussian')
plt.plot(B, lorentz, label='Lorentzian')  
plt.plot(B, voigt, label='True Voigt')
plt.plot(B, mixed, label='Pseudo-Voigt')
plt.legend()
plt.show()
```

### 5. Simple EPR Plotting

The `epyr.eprplot` module provides simple, direct plotting functions for EPR data from `eprload()`.

```python
import epyr
import matplotlib.pyplot as plt

# Load EPR data
x, y, params, filepath = epyr.eprload("data.dsc")

# Plot 1D spectrum
epyr.eprplot.plot_1d(x, y, params, title="EPR Spectrum")

# Plot 2D data as color map
epyr.eprplot.plot_2d_map(x, y, params, title="2D EPR Map")

# Plot 2D data as waterfall plot
epyr.eprplot.plot_2d_waterfall(x, y, params, title="2D Waterfall", max_traces=30)

plt.show()
```

### 6. Interactive Command-Line Plotting

Experience the new interactive plotting with measurement tools:

```bash
# Interactive plotting with measurement tools
epyr-plot --interactive --measure

# Load specific file with measurements
epyr-plot spectrum.dsc --interactive --measure

# Quick analysis with scaling and save
epyr-plot data.dta --interactive --measure -s nG --save
```

**Measurement Features:**
- Click two points to measure Δx, Δy, and distance
- Real-time visual feedback with annotations
- macOS optimized for smooth performance
- Right-click to clear, 'c' key for quick clear

### 7. Isotope GUI

Run the interactive isotope GUI to explore nuclear data. Note that this requires the `pandas` library.

```python
from epyr.isotope_gui import run_gui

# This will launch the Tkinter GUI
run_gui()
```

## Interactive Tutorials

EPyR Tools includes comprehensive Jupyter notebook tutorials for hands-on learning:

### **Getting Started (Beginner)**
```bash
cd examples/notebooks
jupyter notebook 01_Getting_Started.ipynb
```
Learn EPR data loading, visualization, and FAIR data conversion.

### **Baseline Correction (Intermediate)**
```bash
jupyter notebook 06_Baseline_Correction_Functions_Complete.ipynb
```
Master modern polynomial baseline correction with interactive region selection.

### **Advanced Analysis (Expert)**
```bash
jupyter notebook 03_Advanced_Analysis.ipynb
```
Complete EPR analysis: g-factors, hyperfine structure, quantitative integration.

### **EPR Lineshape Analysis (Professional)**
```bash
jupyter notebook 04_EPR_Lineshape_Analysis.ipynb
```
Comprehensive lineshape analysis: Gaussian, Lorentzian, Voigt profiles, derivatives, and convolution techniques.

### **Example Scripts**
Ready-to-use automation scripts in `examples/scripts/`:
```bash
python examples/scripts/01_basic_loading.py
python examples/scripts/02_baseline_correction.py
python examples/scripts/04_lineshape_analysis.py
```

## Project Structure

```
epyrtools/
├── epyr/                           # Main package
│   ├── eprload.py                 # Core data loading (BES3T, ESP formats)
│   ├── eprplot.py                 # Simple EPR plotting functions
│   ├── baseline_correction.py     # Modern streamlined baseline correction
│   ├── baseline/                  # Compatibility layer for baseline functions
│   │   └── __init__.py           # Imports from baseline_correction
│   ├── fair/                     # FAIR data conversion
│   │   ├── conversion.py         # Format conversion tools
│   │   ├── exporters.py          # CSV, JSON, HDF5 export
│   │   └── parameter_mapping.py  # Metadata standardization
│   ├── lineshapes/               # EPR lineshape analysis
│   │   ├── gaussian.py           # Gaussian profiles with derivatives
│   │   ├── lorentzian.py         # Lorentzian profiles with phase rotation
│   │   ├── voigtian.py           # True Voigt convolution profiles
│   │   ├── lshape.py             # General lineshape functions
│   │   ├── convspec.py           # Spectrum convolution tools
│   │   └── lineshape_class.py    # Unified lineshape interface
│   ├── constants.py              # EPR physical constants
│   ├── plot.py                   # Advanced EPR plotting
│   ├── isotope_gui/             # Interactive isotope database
│   └── sub/                     # Utility modules
│       ├── loadBES3T.py         # BES3T format loader
│       ├── loadESP.py           # ESP format loader
│       └── utils.py             # File handling utilities
├── docs/                        # Sphinx API documentation
├── examples/                    # Comprehensive tutorial system
│   ├── notebooks/               # Interactive Jupyter tutorials
│   │   └── Getting_Started.ipynb    # Complete beginner tutorial
│   ├── scripts/                 # Python automation examples
│   └── data/                    # EPR measurement files
│       ├── *.DSC, *.DTA        # BES3T format files
│       ├── *.par, *.spc        # ESP format files
│       └── processed/          # Analysis results examples
├── tests/                       # Comprehensive test suite (44 tests)
├── requirements.txt             # Core dependencies
├── requirements-dev.txt         # Development tools
└── pyproject.toml              # Modern Python packaging
```

## Contributing & Support

### **Documentation**
- **API Reference:** [docs/](docs/) - Complete function documentation
- **Tutorials:** [examples/notebooks/](examples/notebooks/) - Interactive learning
- **Examples:** [examples/scripts/](examples/scripts/) - Ready-to-use code

### **Community**
- **Issues:** [GitHub Issues](https://github.com/BertainaS/epyrtools/issues)
- **Discussions:** Share EPR analysis workflows and tips
- **Contributing:** See contribution guidelines for code contributions

### **Quality Assurance**
- **44 passing tests** with pytest
- **Pre-commit hooks** for code quality
- **Type hints** and comprehensive docstrings
- **Professional packaging** with modern Python standards

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## Contributors

**Lead Developer & Maintainer:**
- **Sylvain Bertaina** - [sylvain.bertaina@cnrs.fr](mailto:sylvain.bertaina@cnrs.fr)

**Affiliation:**
- [Magnetism Group (MAG), IM2NP Laboratory](https://www.im2np.fr/fr/equipe-magnetisme-mag)
- CNRS (Centre National de la Recherche Scientifique)

---

**EPyR Tools** - *Professional EPR analysis for the Python ecosystem*
