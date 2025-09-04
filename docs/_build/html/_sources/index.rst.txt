EPyR Tools Documentation
========================

**EPyR Tools** is a Python package designed for reading, processing, and visualizing Electron Paramagnetic Resonance (EPR) data. It provides a robust toolkit for handling proprietary data files from Bruker spectrometers and converting them into open, FAIR (Findable, Accessible, Interoperable, and Reusable) formats.

.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
   :target: https://opensource.org/licenses/BSD-3-Clause
   :alt: License

Key Features
------------

* **Load Bruker Data:** Easily load data from Bruker BES3T (.dta, .dsc) and ESP/WinEPR (.par, .spc) files
* **FAIR Data Conversion:** Convert proprietary formats to CSV, JSON, and HDF5
* **Advanced Baseline Correction:** 1D and 2D baseline removal with multiple algorithms
* **Specialized Plotting:** Publication-quality EPR plots and 2D spectral maps  
* **Isotope GUI:** Interactive periodic table for nuclear isotope properties

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash

   pip install -e .

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   import epyr.eprload as eprload
   
   # Load EPR data
   x, y, params, file_path = eprload.eprload('data.dsc')
   
   # Apply baseline correction
   from epyr.baseline import baseline_polynomial
   y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)

   # Convert to FAIR format
   from epyr.fair import convert_bruker_to_fair
   convert_bruker_to_fair('data.dsc', output_dir='./fair_data')

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

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide
   
   contributing
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`