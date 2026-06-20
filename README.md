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
| [![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) | ![Tests Passing](https://img.shields.io/badge/tests-369%20passed-brightgreen) | [![Documentation](https://img.shields.io/badge/docs-Sphinx-blue)](docs/) | ![Version](https://img.shields.io/badge/version-0.4.0-blue) |

## What is EPyR Tools?

**EPyR Tools** is a Python package for analyzing Electron Paramagnetic Resonance (EPR) spectroscopy data from Bruker spectrometers. It covers the full path from raw instrument files to publication-ready results: loading BES3T and ESP/WinEPR binary formats, baseline correction, lineshape and T1/T2 relaxation fitting, FFT-based analysis of pulse-EPR time-domain data, and FAIR-compliant export to CSV, JSON, and HDF5.

The package targets EPR researchers who want a reproducible, scriptable Python workflow in place of proprietary vendor software, with a command-line interface for routine processing and a documented Python API for custom analysis.

## Key Features

### Data Loading & FAIR Conversion
- Read Bruker BES3T (`.DTA`/`.DSC`) and ESP/WinEPR (`.spc`/`.par`) files, 1D and 2D, with automatic format detection
- `epyr.fair`: export to CSV, JSON, and HDF5 with standardized, complete metadata
- Batch conversion of entire directories
- Plugin architecture for adding new file formats

### Baseline Correction (`epyr.baseline`)
- Polynomial, exponential, bi-exponential, and stretched-exponential models, for 1D and 2D data
- Automatic model selection ranked by AIC/BIC/R²
- Manual and interactive baseline-region selection

### Lineshape Analysis & Fitting (`epyr.lineshapes`)
- Gaussian, Lorentzian, true Voigt (Faddeeva-function convolution), and pseudo-Voigt lineshapes
- Absorption and 1st/2nd derivative forms, phase mixing, optional affine baseline term
- `fit_epr_signal()` and `fit_multiple_shapes()` for single-model and model-comparison fitting

### T1/T2 Relaxation Fitting (`epyr.relaxation`, new in v0.4.0)
- Six decay/recovery models: mono- and stretched-exponential, bi-exponential, inversion/saturation recovery, and combined homogeneous/spectral-diffusion (Gamma0/GammaG) echo decay
- `fit_relaxation()` for a single model, `fit_multiple_decays()` to rank candidates by reduced chi-squared
- Fit results print a readable parameter summary directly

### Signal Processing (`epyr.signalprocessing`)
- FFT-based frequency analysis for pulse-EPR time-domain data: Rabi oscillations, DEER, HYSCORE
- 1D and 2D FFT modes, apodization windows (Hann, Hamming, Blackman, Kaiser), automatic time-unit detection, zero-padding
- Power spectral density (Welch, periodogram) and spectrogram analysis

### Physics & Units (`epyr.physics`)
- CODATA 2022 physical constants in SI and CGS units (`GFREE`, `BMAGN`, `NMAGN`, `PLANCK`, ...)
- Field/frequency conversion (mT ↔ MHz ↔ cm⁻¹) via `unitconvert()` and dedicated helpers

### Command-Line Interface
- Nine commands covering the full workflow: `epyr-convert`, `epyr-baseline`, `epyr-batch-convert`, `epyr-config`, `epyr-info`, `epyr-isotopes`, `epyr-plot`, `epyr-validate`, `epyrview`
- `epyr-plot --interactive --measure`: click-to-measure delta x/y distance tool

### Performance & Isotope Database
- `OptimizedLoader` / `DataCache` for large files, with memory monitoring and streaming
- `epyr.isotopes` / `epyr-isotopes`: interactive periodic-table GUI with NMR frequency calculator and X/Q/W-band presets

## What's New in v0.4.0

Version 0.4.0 adds the **`epyr.relaxation`** package for T1/T2 relaxation fitting, complementing the existing field-domain lineshape fitting in `epyr.lineshapes.fitting`:

```python
from epyr.relaxation import fit_relaxation, fit_multiple_decays

# Single model
result = fit_relaxation(t, y, model="stretched_exponential")
print(result)
# === Relaxation Fit Results - stretched_exponential ===
# Success: True
# R2 = 0.998452
# ...

# Compare candidate models, ranked by reduced chi-squared (not R-squared,
# which is biased toward models with more free parameters)
results = fit_multiple_decays(t, y)
print(results)
# model                  success  R2        chi2       amplitude  T      ...
# mono_exponential       True     0.998391  0.0004842  2.005      1.310  ...
# stretched_exponential  True     0.998452  0.000476   1.976      1.309  ...
```

Fit plots in both `epyr.relaxation` and `epyr.lineshapes.fitting` now follow `matplotlib.rcParams` for figure size, marker size, line width, and font size, instead of a fixed layout.

See [docs/release_notes/v0.4.0.rst](docs/release_notes/v0.4.0.rst) for full details, or [docs/release_notes.rst](docs/release_notes.rst) for the complete version history.

## Installation

### Prerequisites
- Python 3.8 or higher
- NumPy, SciPy, matplotlib, pandas, h5py (installed automatically)

### Quick Install
```bash
pip install epyr-tools
```

### Development Installation
```bash
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools

# Install with development dependencies
pip install -e ".[dev,docs]"

# Set up pre-commit hooks
pre-commit install
```

### Verification
```bash
epyr --help
epyr-info
make test
```

## Getting Started

### 1. Loading Data

```python
import epyr

# Open a file dialog to select a .dta, .dsc, .spc, or .par file
x, y, params, filepath = epyr.eprload()

# Or specify a path directly:
# x, y, params, filepath = epyr.eprload('path/to/data.dsc')
```

### 2. Converting to FAIR Formats

```python
from epyr.fair import convert_bruker_to_fair

convert_bruker_to_fair('path/to/data.dsc', output_dir='path/to/output')
```

### 3. Baseline Correction

```python
import epyr

x, y, params, filepath = epyr.eprload("data.dsc")

# Automatic model selection
corrected, baseline, model_info = epyr.baseline.baseline_auto_1d(x, y, params)

# Or a specific model with manual region exclusion (e.g. signal regions, in mT)
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params,
    manual_regions=[(3340, 3360), (3380, 3400)],
    region_mode='exclude',
    order=2,
)
```

### 4. Lineshape Fitting

```python
from epyr.lineshapes import fit_epr_signal, fit_multiple_shapes

x, y, params, filepath = epyr.eprload('data.DTA')

# Single model
result = fit_epr_signal(x, y, 'gaussian')
print(result.summary())

# 1st-derivative signal with adjustable phase, comparing all lineshapes
results = fit_multiple_shapes(x, y, derivative=1, fit_phase=True)
```

### 5. T1/T2 Relaxation Fitting

```python
from epyr.relaxation import fit_relaxation, fit_multiple_decays

t, y, params, filepath = epyr.eprload('echo_decay.DTA')
y = abs(y)  # take the magnitude of a complex echo signal

result = fit_relaxation(t, y, model="stretched_exponential")
results = fit_multiple_decays(t, y)  # compare mono/stretched/bi-exponential
```

### 6. Time-Domain Signal Processing (FFT)

```python
from epyr.signalprocessing import analyze_frequencies

t, y, params, filepath = epyr.eprload('rabi_oscillation.DTA')
result = analyze_frequencies(t, y, window='hann', zero_padding=4)
print(f"Dominant frequency: {result['dominant_frequencies'][0]:.3f} MHz")
```

### 7. Plotting and the CLI

```python
import epyr

x, y, params, filepath = epyr.eprload("data.dsc")
epyr.plot_1d(x, y, params, title="EPR Spectrum")
```

```bash
# Interactive plot with click-to-measure delta x/y
epyr-plot spectrum.dsc --interactive --measure

# Batch FAIR conversion
epyr-batch-convert ./data --formats csv,json,hdf5
```

## Tutorials & Examples

### Jupyter Notebook Series
An eight-notebook tutorial series in `examples/notebooks/`, using real experimental data from `examples/data/`. Notebooks are committed without outputs; run them to generate figures.

```bash
cd examples/notebooks
jupyter lab 00_Tutorial_Series_Index.ipynb  # index and navigation
```

| Notebook | Topic |
|----------|-------|
| `01_Loading_and_Plotting.ipynb` | `eprload`, parameter inspection, 1D/2D plotting |
| `02_Baseline_Correction.ipynb` | Polynomial, automatic, and exponential baselines |
| `03_Lineshape_Analysis_and_Fitting.ipynb` | Gaussian/Lorentzian/Voigt, derivatives, fitting |
| `04_Relaxation_Fitting.ipynb` | T1/T2 decay/recovery models (new in v0.4.0) |
| `05_Signal_Processing_and_FFT.ipynb` | Frequency analysis of Rabi data, apodization |
| `06_FAIR_Conversion_and_Export.ipynb` | CSV/JSON/HDF5 export and validation |
| `07_Physics_Units_and_Constants.ipynb` | CODATA constants, field/frequency conversions |

### Standalone Example Scripts
Six short, self-contained scripts in `examples/clean/` exercising the public API end to end:

```bash
python examples/clean/01_basic_loading_and_plotting.py
python examples/clean/02_baseline_and_fitting.py
python examples/clean/03_advanced_fft_windows.py
python examples/clean/04_interactive_2d_slicer.py
python examples/clean/05_rabi_frequency_analysis.py
python examples/clean/06_relaxation_fitting.py
```

See [docs/tutorials/clean_examples.rst](docs/tutorials/clean_examples.rst) for a description of each script.

## Project Structure

```
epyrtools/
├── epyr/                          # Main package
│   ├── eprload.py                 # Core data loading (BES3T, ESP formats)
│   ├── eprplot.py                 # EPR plotting (1D, 2D map, waterfall, slicer)
│   ├── cli.py                     # Command-line interface (9 commands)
│   ├── config.py                  # Hierarchical configuration system
│   ├── performance.py             # OptimizedLoader, DataCache
│   ├── plugins.py                 # Plugin architecture
│   ├── logging_config.py          # Centralized logging
│   ├── isotope_gui.py             # Interactive isotope database GUI
│   ├── baseline/                  # Baseline correction (correction, selection, models, interactive)
│   ├── lineshapes/                # Gaussian, Lorentzian, Voigt, pseudo-Voigt, fitting
│   ├── relaxation/                # T1/T2 decay/recovery models and fitting
│   ├── signalprocessing/          # FFT frequency analysis, apodization windows
│   ├── physics/                   # CODATA constants and unit conversions
│   ├── fair/                      # FAIR conversion, exporters, validation
│   └── sub/                       # Bruker BES3T/ESP format loaders
├── docs/                          # Sphinx documentation and tutorials
├── examples/
│   ├── notebooks/                 # Jupyter tutorial series
│   ├── clean/                     # Six standalone end-to-end scripts
│   └── data/                      # Real EPR measurement files (CW, pulse, 2D)
├── tests/                         # Test suite (369 tests; smoke/standard/deep/scientific)
└── pyproject.toml                 # Packaging, dependencies, entry points
```

## Documentation

- [Full documentation](docs/): guides and API reference (Sphinx)
- [User guide](docs/user_guide.md): workflows and step-by-step tutorials
- [CLI reference](docs/cli_reference.md): command-line interface
- [API reference](docs/api_reference.md): public API
- [Release notes](docs/release_notes.rst): version history

## Testing & Quality

EPyR Tools follows a 4-level testing protocol (`pytest -m smoke|standard|deep|scientific`), with 369 tests covering basic functionality, broad feature coverage, edge cases, and scientific validation against NIST/CODATA values.

```bash
make test        # full suite
make test-cov    # with coverage report
make quality     # lint, type-check, security
```

## Contributing & Support

- **Issues:** [GitHub Issues](https://github.com/BertainaS/epyrtools/issues)
- **Contributing guide:** [docs/contributing.rst](docs/contributing.rst)

## License

This project is licensed under the **MIT License**, see [LICENSE](LICENSE) for details.

## Contributors

**Lead Developer & Maintainer:**
- **Sylvain Bertaina**, [sylvain.bertaina@cnrs.fr](mailto:sylvain.bertaina@cnrs.fr)

**Affiliation:**
- [Magnetism Group (MAG), IM2NP Laboratory](https://www.im2np.fr/fr/equipe-magnetisme-mag)
- CNRS (Centre National de la Recherche Scientifique)

---

**EPyR Tools**: EPR data analysis in Python.
