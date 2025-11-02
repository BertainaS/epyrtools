EPyR Tools Documentation
========================

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. image:: https://img.shields.io/badge/version-0.3.0-blue
   :alt: Version

.. image:: https://img.shields.io/badge/tests-100%2B%20passed-brightgreen
   :alt: Tests

**EPyR Tools** is a comprehensive Python package for Electron Paramagnetic Resonance (EPR) spectroscopy data analysis. It provides a complete toolkit for loading, processing, analyzing, and visualizing EPR data from Bruker spectrometers, with a focus on FAIR (Findable, Accessible, Interoperable, and Reusable) data principles.

From basic data loading to advanced quantitative analysis, EPyR Tools offers professional-grade capabilities for EPR researchers, with comprehensive documentation and interactive tutorials.

Key Features
---------------

**Data Loading & Formats**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Bruker File Support:** Load BES3T (.dta/.dsc) and ESP/WinEPR (.par/.spc) files seamlessly
* **Automatic Format Detection:** Smart file format recognition and parameter extraction
* **FAIR Data Conversion:** Export to CSV, JSON, and HDF5 formats with complete metadata
* **Batch Processing:** Handle multiple files efficiently

**Advanced Analysis**
~~~~~~~~~~~~~~~~~~~~~~~~
* **Baseline Correction:** Multiple algorithms (polynomial, exponential) with signal exclusion
* **Peak Detection:** Automatic identification of EPR spectral features
* **g-Factor Calculations:** Precise electronic g-factor determination with field calibration
* **Quantitative Integration:** Single and double integration for spin quantification
* **Lineshape Analysis:** Comprehensive EPR lineshape functions (Gaussian, Lorentzian, Voigt, pseudo-Voigt)

**Visualization & Plotting**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Interactive CLI Plotting:** Command-line plotting with measurement tools (``epyr-plot --interactive --measure``)
* **Delta X/Y Measurements:** Click two points to measure precise distances with visual feedback
* **2D Spectral Maps:** Professional publication-quality EPR plots
* **macOS Optimized:** Smooth interactive plotting with TkAgg backend
* **Customizable Styling:** Flexible plot configuration for different EPR experiments
* **Export Options:** High-resolution outputs for publications

**Learning & Documentation**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* **Interactive Tutorials:** Comprehensive Jupyter notebooks (beginner → advanced)
* **Complete API Documentation:** Professional Sphinx-generated docs
* **Example Scripts:** Ready-to-use Python automation scripts
* **Best Practices Guide:** EPR analysis workflows and quality assessment

**EPR-Specific Tools**
~~~~~~~~~~~~~~~~~~~~~~~~~
* **Physical Constants:** Comprehensive EPR-relevant constants library
* **Isotope Database:** Nuclear properties and magnetic parameters
* **Field-Frequency Conversion:** Precise EPR measurement calculations
* **Spectrometer Support:** Optimized for modern Bruker EPR systems

Quick Start
-----------

.. code-block:: python

   import epyr

   # Load EPR data (opens file dialog if no path given)
   x, y, params, filepath = epyr.eprload()

   # Apply baseline correction
   from epyr.baseline import baseline_polynomial
   y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)

   # Convert to FAIR formats (CSV, JSON, HDF5)
   from epyr.fair import convert_bruker_to_fair
   convert_bruker_to_fair(filepath, output_dir='./fair_data')

   # Generate EPR lineshapes for analysis
   from epyr.lineshapes import gaussian, lorentzian, Lineshape
   B = np.linspace(320, 340, 1000)  # mT
   gauss_line = gaussian(B, center=334.8, width=2.0)
   lorentz_line = lorentzian(B, center=334.8, width=2.0)
   
   # Or use the unified interface
   shape = Lineshape('pseudo_voigt', width=2.0, alpha=0.5)
   mixed_line = shape(B, center=334.8)

   # Create publication-quality plots
   import matplotlib.pyplot as plt
   plt.plot(x, y_corrected, 'b-', linewidth=1.5, label='Corrected Data')
   plt.plot(B, gauss_line, 'r--', linewidth=1.5, label='Gaussian Fit')
   plt.xlabel('Magnetic Field (G)')
   plt.ylabel('EPR Signal (a.u.)')
   plt.legend()
   plt.show()

For interactive plotting with measurement tools from command line:

.. code-block:: bash

   # Interactive plotting with measurement tools
   epyr-plot --interactive --measure

   # Load specific file with measurements
   epyr-plot spectrum.dsc --interactive --measure --save

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   tutorials/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/epyr
   api/modules
   cli_reference

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide

   contributing
   changelog
   release_notes

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
