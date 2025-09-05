# EPyR Tools: Electron Paramagnetic Resonance (EPR) Tools in Python

| License | Tests | Documentation | Version |
|---------|-------|---------------|---------|
| [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) | ![Tests Passing](https://img.shields.io/badge/tests-44%20passed-brightgreen) | [![Documentation](https://img.shields.io/badge/docs-available-blue)](docs/) | ![Version](https://img.shields.io/badge/version-0.1.1-blue) |

**EPyR Tools** is a comprehensive Python package for Electron Paramagnetic Resonance (EPR) spectroscopy data analysis. It provides a complete toolkit for loading, processing, analyzing, and visualizing EPR data from Bruker spectrometers, with a focus on FAIR (Findable, Accessible, Interoperable, and Reusable) data principles.

From basic data loading to advanced quantitative analysis, EPyR Tools offers professional-grade capabilities for EPR researchers, with comprehensive documentation and interactive tutorials.

## 🚀 Key Features

### **📁 Data Loading & Formats**
- **Bruker File Support:** Load BES3T (.dta/.dsc) and ESP/WinEPR (.par/.spc) files seamlessly
- **Automatic Format Detection:** Smart file format recognition and parameter extraction
- **FAIR Data Conversion:** Export to CSV, JSON, and HDF5 formats with complete metadata
- **Batch Processing:** Handle multiple files efficiently

### **📊 Advanced Analysis**
- **Baseline Correction:** Multiple algorithms (polynomial, exponential) with signal exclusion
- **Peak Detection:** Automatic identification of EPR spectral features
- **g-Factor Calculations:** Precise electronic g-factor determination with field calibration
- **Hyperfine Analysis:** Pattern recognition and coupling constant extraction
- **Quantitative Integration:** Single and double integration for spin quantification

### **📈 Visualization & Plotting**
- **2D Spectral Maps:** Professional publication-quality EPR plots
- **Interactive Plotting:** Real-time parameter adjustment and analysis
- **Customizable Styling:** Flexible plot configuration for different EPR experiments
- **Export Options:** High-resolution outputs for publications

### **🎓 Learning & Documentation**
- **Interactive Tutorials:** 3 comprehensive Jupyter notebooks (beginner → advanced)
- **Complete API Documentation:** Professional Sphinx-generated docs
- **Example Scripts:** Ready-to-use Python automation scripts
- **Best Practices Guide:** EPR analysis workflows and quality assessment

### **🧪 EPR-Specific Tools**
- **Physical Constants:** Comprehensive EPR-relevant constants library
- **Isotope Database:** Nuclear properties and magnetic parameters
- **Field-Frequency Conversion:** Precise EPR measurement calculations
- **Spectrometer Support:** Optimized for modern Bruker EPR systems

## 📦 Installation

### Prerequisites
- Python 3.8 or higher
- NumPy, matplotlib, pandas, h5py (automatically installed)

### Recommended: Install as Package
```bash
# Clone the repository
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools

# Install package with dependencies
pip install .

# For development (editable install)
pip install -e .

# With development tools
pip install -e .[dev]

# With documentation tools
pip install -e .[docs]

# Install everything
pip install -e .[all]
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

## 🚀 Getting Started

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

### 3. Baseline Correction

The `epyr.baseline` module provides functions for correcting distorted baselines. For more detailed examples, see the `notebooks/Demo_baseline_correction.ipynb` notebook.

```python
from epyr.baseline import baseline_polynomial
import numpy as np

# Assuming 'x' and 'y' are loaded spectral data
# For this example, let's create some sample data
x = np.linspace(0, 100, 500)
y = (0.1 * x + 5) + 10 * np.exp(-((x - 50)**2) / 10) # A peak on a sloped baseline

# Correct a linear baseline, excluding the peak region from the fit
y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1, exclude_regions=[(40, 60)])
```

### 4. Specialized Plotting

The `epyr.plot` module offers advanced plotting functions.

```python
from epyr.plot import plot_2d_map
import numpy as np
import matplotlib.pyplot as plt

# Create sample 2D data
x_axis = np.linspace(-10, 10, 100)
y_axis = np.linspace(-10, 10, 100)
XX, YY = np.meshgrid(x_axis, y_axis)
Z_data = np.exp(-(XX**2 + YY**2) / 8) # A 2D Gaussian peak

# Generate a 2D color map
fig, ax = plot_2d_map(x_axis, y_axis, Z_data, x_unit='mT', y_unit='GHz')
plt.show() # If running as a script
```

### 5. Isotope GUI

Run the interactive isotope GUI to explore nuclear data. Note that this requires the `pandas` library.

```python
from epyr.isotope_gui import run_gui

# This will launch the Tkinter GUI
run_gui()
```

## 📚 Interactive Tutorials

EPyR Tools includes comprehensive Jupyter notebook tutorials for hands-on learning:

### **Getting Started (Beginner)**
```bash
cd examples/notebooks
jupyter notebook 01_Getting_Started.ipynb
```
Learn EPR data loading, visualization, and FAIR data conversion.

### **Baseline Correction (Intermediate)**
```bash
jupyter notebook 02_Baseline_Correction.ipynb
```
Master polynomial and advanced baseline correction techniques.

### **Advanced Analysis (Expert)**
```bash
jupyter notebook 03_Advanced_Analysis.ipynb
```
Complete EPR analysis: g-factors, hyperfine structure, quantitative integration.

### **Example Scripts**
Ready-to-use automation scripts in `examples/scripts/`:
```bash
python examples/scripts/01_basic_loading.py
python examples/scripts/02_baseline_correction.py
```

## 📁 Project Structure

```
epyrtools/
├── epyr/                           # Main package
│   ├── eprload.py                 # Core data loading (BES3T, ESP formats)
│   ├── baseline/                  # Advanced baseline correction
│   │   ├── _1d.py                # 1D correction algorithms
│   │   ├── _2d.py                # 2D correction algorithms
│   │   └── _utils.py             # Correction utilities
│   ├── fair/                     # FAIR data conversion
│   │   ├── conversion.py         # Format conversion tools
│   │   ├── exporters.py          # CSV, JSON, HDF5 export
│   │   └── parameter_mapping.py  # Metadata standardization
│   ├── constants.py              # EPR physical constants
│   ├── plot.py                   # Advanced EPR plotting
│   ├── isotope_gui/             # Interactive isotope database
│   └── sub/                     # Low-level processing modules
├── docs/                        # Sphinx API documentation
├── examples/                    # Comprehensive tutorial system
│   ├── notebooks/               # Interactive Jupyter tutorials
│   │   ├── 01_Getting_Started.ipynb
│   │   ├── 02_Baseline_Correction.ipynb
│   │   └── 03_Advanced_Analysis.ipynb
│   ├── scripts/                 # Python automation examples
│   └── data/                    # Sample data organization
│       ├── BES3T/              # Modern Bruker format samples
│       ├── ESP/                # Legacy Bruker format samples
│       └── processed/          # Analysis results examples
├── tests/                       # Comprehensive test suite (44 tests)
├── requirements.txt             # Core dependencies
├── requirements-dev.txt         # Development tools
└── pyproject.toml              # Modern Python packaging
```

## 🤝 Contributing & Support

### **Documentation**
- **API Reference:** [docs/](docs/) - Complete function documentation
- **Tutorials:** [examples/notebooks/](examples/notebooks/) - Interactive learning
- **Examples:** [examples/scripts/](examples/scripts/) - Ready-to-use code

### **Community**
- **Issues:** [GitHub Issues](https://github.com/BertainaS/epyrtools/issues)
- **Discussions:** Share EPR analysis workflows and tips
- **Contributing:** See contribution guidelines for code contributions

### **Quality Assurance**
- ✅ **44 passing tests** with pytest
- ✅ **Pre-commit hooks** for code quality
- ✅ **Type hints** and comprehensive docstrings
- ✅ **Professional packaging** with modern Python standards

## 📄 License

This project is licensed under the **BSD 3-Clause License** - see the [LICENSE](LICENSE) file for details.

---

**EPyR Tools** - *Professional EPR analysis for the Python ecosystem* 🧲⚛️📊
