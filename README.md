# EPyR Tools: Electron Paramagnetic Resonance (EPR) Tools in Python

| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) |
|--------------|--------------------------------------------------------------------------------------------------------------------------|

**EPyR Tools** is a Python package designed for reading, processing, and visualizing Electron Paramagnetic Resonance (EPR) data. It provides a robust toolkit for handling proprietary data files from Bruker spectrometers (EMX, Elexis) and converting them into open, FAIR (Findable, Accessible, Interoperable, and Reusable) formats.

Beyond simple data loading, EPyR Tools offers advanced functionalities for baseline correction, specialized plotting, and interactive data exploration through a graphical user interface.

## Key Features

*   **Load Bruker Data:** Easily load data from Bruker BES3T (.dta, .dsc) and ESP/WinEPR (.par, .spc) files into Python.
*   **FAIR Data Conversion:** Convert proprietary Bruker files into accessible formats:
    *   **CSV & JSON:** Exports data to simple CSV files and creates detailed JSON files for metadata with human-readable parameter names.
    *   **HDF5:** Creates a self-contained, structured format that includes both data and metadata, ideal for long-term storage and analysis.
*   **Advanced Baseline Correction:** A suite of tools for 1D and 2D baseline removal, including polynomial, constant offset, and mono/stretched exponential decay models.
*   **Specialized Plotting:** Generate publication-quality plots tailored for EPR data, such as 2D spectral maps and waterfall plots for angular sweeps.
*   **Isotope GUI:** An interactive periodic table GUI to quickly look up nuclear isotope properties, such as spin, g-factor, and NMR frequency at a given magnetic field.

## Installation

EPyR Tools requires Python 3 and the following libraries. You can install them using the provided `requirements.txt` file:

```bash
pip install -r requirements.txt
```

## Quick Start

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

## Project Structure

*   `epyr/`: The main source directory for the package.
    *   `eprload.py`: Core module for loading Bruker data files.
    *   `fair.py`: Module for converting data to FAIR formats.
    *   `baseline/`: Sub-package containing 1D and 2D baseline correction algorithms.
    *   `plot.py`: Module for creating specialized 2D plots.
    *   `isotope_gui.py`: A Tkinter-based GUI for isotope data visualization.
    *   `sub/`: Contains low-level helper modules for data loading and processing.
*   `notebooks/`: Contains Jupyter notebooks demonstrating key features, such as an in-depth tutorial on baseline correction.
*   `requirements.txt`: A list of Python dependencies.

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.
