Quick Start Guide
=================

This guide will get you started with EPyR Tools in just a few minutes.

Basic Data Loading
------------------

EPyR Tools makes it easy to load EPR data from Bruker spectrometers:

.. code-block:: python

   import epyr

   # Load EPR data (opens file dialog if no path provided)
   x, y, params, filepath = epyr.eprload()

   # Or specify a file directly
   x, y, params, filepath = epyr.eprload('path/to/your/data.dsc')

   # Check what you loaded
   print(f"Data points: {len(y) if hasattr(y, '__len__') else y.shape}")
   print(f"Field range: {x[0].min():.1f} - {x[0].max():.1f} G" if isinstance(x, list) else f"{x.min():.1f} - {x.max():.1f} G")

Supported File Formats
~~~~~~~~~~~~~~~~~~~~~~~

* **BES3T format**: ``.dsc`` and ``.dta`` file pairs (modern Bruker)
* **ESP format**: ``.par`` and ``.spc`` file pairs (legacy Bruker)

The function automatically detects the format and returns:

* **x**: Field axis (or list of axes for 2D data)
* **y**: Intensity data (1D array or 2D matrix)
* **params**: Dictionary of experimental parameters
* **filepath**: Path to the loaded file

Baseline Correction
-------------------

Remove baseline drift and artifacts from your spectra:

.. code-block:: python

   from epyr.baseline import baseline_polynomial

   # Simple polynomial baseline correction
   y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)

   # Advanced correction with signal exclusion
   exclude_regions = [(3300, 3400)]  # Exclude peak regions
   y_corrected, baseline = baseline_polynomial(
       y, x_data=x, poly_order=2, exclude_regions=exclude_regions
   )

Available Methods
~~~~~~~~~~~~~~~~~

* **Polynomial**: Orders 0-5 (constant, linear, quadratic, etc.)
* **Exponential**: Single and stretched exponential decay models
* **Signal exclusion**: Exclude peaks and features from baseline fitting

FAIR Data Conversion
--------------------

Convert proprietary Bruker files to open, accessible formats:

.. code-block:: python

   from epyr.fair import convert_bruker_to_fair

   # Convert single file to all formats (CSV, JSON, HDF5)
   convert_bruker_to_fair('data.dsc', output_dir='./converted')

   # The conversion creates:
   # - data.csv: Tabulated data with metadata headers
   # - data.json: Complete parameters and measurement info
   # - data.h5: Efficient binary format with full metadata

FAIR Benefits
~~~~~~~~~~~~~

* **Findable**: Standardized metadata and naming
* **Accessible**: No proprietary software needed
* **Interoperable**: Works with Excel, R, Python, MATLAB
* **Reusable**: Complete experimental documentation

Visualization
-------------

Create publication-quality EPR plots:

.. code-block:: python

   import matplotlib.pyplot as plt

   # For 1D spectra
   if not isinstance(x, list):
       plt.figure(figsize=(10, 6))
       plt.plot(x, y, 'b-', linewidth=1.5)
       plt.xlabel('Magnetic Field (G)')
       plt.ylabel('EPR Signal (a.u.)')
       plt.grid(True, alpha=0.3)
       plt.show()

   # For 2D data, EPyR Tools provides specialized plotting
   from epyr.plot import plot_2d_map

   if isinstance(x, list) and len(x) > 1:
       # Create 2D color map
       fig, ax = plot_2d_map(x[0], x[1], y, x_unit='G', y_unit='ns')
       plt.show()

Advanced Features
~~~~~~~~~~~~~~~~~

* **2D spectral maps**: Color plots with customizable scales
* **Interactive plotting**: Real-time parameter adjustment
* **Export options**: High-resolution outputs for publications

Complete Example
----------------

Here's a complete workflow from loading to analysis:

.. code-block:: python

   import epyr
   from epyr.baseline import baseline_polynomial
   from epyr.fair import convert_bruker_to_fair
   import matplotlib.pyplot as plt
   import numpy as np

   # 1. Load EPR data
   print("Loading EPR data...")
   x, y, params, filepath = epyr.eprload('example.dsc')

   # 2. Display basic info
   if isinstance(x, list):
       print(f"2D data: {y.shape}")
       print(f"Complex data: {np.iscomplexobj(y)}")
   else:
       print(f"1D data: {len(y)} points")
       print(f"Field range: {x.min():.1f} - {x.max():.1f} G")

   # 3. Apply baseline correction (for 1D data)
   if not isinstance(x, list):
       print("Applying baseline correction...")
       y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)

       # Plot before/after
       plt.figure(figsize=(12, 5))
       plt.subplot(1, 2, 1)
       plt.plot(x, y, 'b-', label='Original')
       plt.plot(x, baseline, 'r--', label='Baseline')
       plt.legend()
       plt.title('Before Correction')

       plt.subplot(1, 2, 2)
       plt.plot(x, y_corrected, 'g-', label='Corrected')
       plt.legend()
       plt.title('After Correction')
       plt.tight_layout()
       plt.show()

   # 4. Convert to FAIR formats
   print("Converting to FAIR formats...")
   convert_bruker_to_fair(filepath, output_dir='./fair_data')
   print("✅ Conversion complete!")

   # 5. Display key parameters
   print("\n📋 Key Parameters:")
   key_params = ['MWFQ', 'MWPW', 'HCF', 'HSW', 'AVGS']
   for param in key_params:
       if param in params:
           print(f"  {param}: {params[param]}")

Next Steps
----------

* **Explore Examples**: Check ``examples/notebooks/`` for interactive tutorials
* **Run Scripts**: Try the automated scripts in ``examples/scripts/``
* **Read API Docs**: Browse the complete API reference
* **Join Community**: Report issues and contribute on GitHub

Common Workflows
~~~~~~~~~~~~~~~~

1. **Routine Analysis**: Load → Baseline Correct → Plot → Export
2. **Data Conversion**: Load → Convert to FAIR → Archive
3. **Batch Processing**: Script multiple files → Automated analysis
4. **Research**: Interactive notebooks → Custom analysis → Publication plots
