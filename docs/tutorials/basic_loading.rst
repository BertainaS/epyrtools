Basic Loading Tutorial
======================

This tutorial covers the fundamentals of loading EPR data with EPyR Tools.

Interactive Notebook
--------------------

**File**: ``examples/notebooks/01_Basic_Loading.ipynb``

The most comprehensive way to learn basic loading is through our interactive Jupyter notebook. It provides:

* Step-by-step guidance with explanations
* Real EPR data examples
* Interactive code cells you can modify
* Troubleshooting help and error handling

To run the notebook:

.. code-block:: bash

   cd examples/notebooks
   jupyter notebook 01_Basic_Loading.ipynb

Core Concepts
-------------

File Formats
~~~~~~~~~~~~

EPyR Tools supports two main Bruker file formats:

**BES3T Format (Modern Bruker)**
  * ``.dsc`` files contain parameters and metadata
  * ``.dta`` files contain the actual spectral data
  * Always come as pairs - both files needed

**ESP Format (Legacy Bruker)**
  * ``.par`` files contain parameters
  * ``.spc`` files contain spectral data
  * Also come as pairs

Data Types
~~~~~~~~~~

EPR data can be:

**1D Data**
  * Single spectrum vs magnetic field
  * Returns: x (field array), y (intensity array)
  * Example: CW-EPR spectrum

**2D Data**
  * Parameter-dependent spectra (time, power, angle, etc.)
  * Returns: x (list of axes), y (2D intensity matrix)
  * Example: Rabi oscillation, ENDOR

**Complex vs Real**
  * Pulsed EPR often produces complex data (I/Q detection)
  * CW-EPR typically produces real data
  * Complex data visualization uses magnitude: ``np.abs(data)``

Basic Usage Examples
--------------------

Simple Loading
~~~~~~~~~~~~~~

.. code-block:: python

   import epyr

   # Open file dialog to select data
   x, y, params, filepath = epyr.eprload()

   # Or specify file directly
   x, y, params, filepath = epyr.eprload('path/to/spectrum.dsc')

   # Check what you loaded
   print(f"Loaded: {filepath}")
   if isinstance(x, list):
       print(f"2D data: {y.shape}")
   else:
       print(f"1D data: {len(y)} points")

Handling Different Data Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt

   # Load data
   x, y, params, filepath = epyr.eprload('example.dsc')

   if isinstance(x, list) and len(x) > 1:
       # 2D data
       print(f"2D EPR data: {y.shape}")
       print(f"Complex data: {np.iscomplexobj(y)}")

       # For visualization, use magnitude if complex
       y_plot = np.abs(y) if np.iscomplexobj(y) else y

       # Plot as 2D color map
       plt.figure(figsize=(10, 6))
       plt.imshow(y_plot, aspect='auto', origin='lower', cmap='viridis')
       plt.colorbar(label='Signal (a.u.)')
       plt.xlabel('Field Points')
       plt.ylabel('Parameter Points')
       plt.title('2D EPR Data')
       plt.show()

   else:
       # 1D data
       x_array = x[0] if isinstance(x, list) else x
       print(f"1D EPR data: {len(x_array)} points")
       print(f"Field range: {x_array.min():.1f} - {x_array.max():.1f} G")

       # Simple 1D plot
       plt.figure(figsize=(10, 6))
       plt.plot(x_array, y, 'b-', linewidth=1.5)
       plt.xlabel('Magnetic Field (G)')
       plt.ylabel('EPR Signal (a.u.)')
       plt.title('EPR Spectrum')
       plt.grid(True, alpha=0.3)
       plt.show()

Parameter Extraction
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Common EPR parameters
   key_params = {
       'MWFQ': 'Microwave Frequency (Hz)',
       'MWPW': 'Microwave Power (dB)',
       'HCF': 'Center Field (G)',
       'HSW': 'Sweep Width (G)',
       'AVGS': 'Number of Averages',
       'TE': 'Temperature (K)',
       'MA': 'Modulation Amplitude (G)',
   }

   print("📋 Experimental Parameters:")
   for param, description in key_params.items():
       if param in params:
           value = params[param]
           print(f"  {description}: {value}")

   # Calculate derived quantities
   if 'MWFQ' in params and 'HCF' in params:
       freq_ghz = float(params['MWFQ']) / 1e9
       field_g = float(params['HCF'])

       # Approximate g-factor at center field
       h = 6.626e-34  # Planck constant
       mu_b = 9.274e-24  # Bohr magneton
       g_factor = (h * freq_ghz * 1e9) / (mu_b * field_g * 1e-4)
       print(f"  Center g-factor: {g_factor:.3f}")

Error Handling
--------------

Robust Loading
~~~~~~~~~~~~~~

.. code-block:: python

   def safe_load_epr(file_path):
       """Safely load EPR data with error handling."""
       try:
           x, y, params, filepath = epyr.eprload(file_path)

           if x is None or y is None:
               print(f"❌ Failed to load data from {file_path}")
               return None

           print(f"✅ Successfully loaded {file_path}")
           return x, y, params, filepath

       except FileNotFoundError:
           print(f"❌ File not found: {file_path}")
       except Exception as e:
           print(f"❌ Error loading {file_path}: {e}")

       return None

   # Usage
   result = safe_load_epr('spectrum.dsc')
   if result is not None:
       x, y, params, filepath = result
       # Continue with analysis...

Common Issues
~~~~~~~~~~~~~

**Missing .dta or .spc file**

.. code-block:: python

   from pathlib import Path

   def check_file_pairs(dsc_file):
       """Check if required data file exists."""
       dsc_path = Path(dsc_file)

       if dsc_path.suffix.lower() == '.dsc':
           dta_file = dsc_path.with_suffix('.dta')
           if not dta_file.exists():
               print(f"❌ Missing data file: {dta_file}")
               return False
       elif dsc_path.suffix.lower() == '.par':
           spc_file = dsc_path.with_suffix('.spc')
           if not spc_file.exists():
               print(f"❌ Missing data file: {spc_file}")
               return False

       return True

**Case sensitivity issues**

.. code-block:: python

   def find_epr_files(directory):
       """Find EPR files handling case variations."""
       from pathlib import Path

       data_dir = Path(directory)
       epr_files = []

       # Check both upper and lower case
       for pattern in ['*.dsc', '*.DSC', '*.par', '*.PAR']:
           epr_files.extend(data_dir.glob(pattern))

       return epr_files

Data Export
-----------

Simple Export
~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd

   # Export 1D data to CSV
   if not isinstance(x, list):
       df = pd.DataFrame({
           'Field_G': x,
           'Intensity': y
       })
       df.to_csv('spectrum.csv', index=False)
       print("✅ Data exported to spectrum.csv")

   # Export parameters to JSON
   import json
   with open('parameters.json', 'w') as f:
       json.dump(params, f, indent=2)

Batch Processing
~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path

   def process_directory(data_dir):
       """Process all EPR files in a directory."""
       results = {}

       for epr_file in Path(data_dir).glob('*.dsc'):
           print(f"Processing {epr_file.name}...")

           result = safe_load_epr(epr_file)
           if result is not None:
               x, y, params, filepath = result

               # Store basic info
               results[epr_file.stem] = {
                   'data_type': '2D' if isinstance(x, list) and len(x) > 1 else '1D',
                   'complex': np.iscomplexobj(y),
                   'shape': y.shape,
                   'frequency': params.get('MWFQ', 'Unknown'),
                   'temperature': params.get('TE', 'Unknown')
               }

       return results

   # Process all files
   results = process_directory('examples/data')
   for filename, info in results.items():
       print(f"{filename}: {info['data_type']} data, shape {info['shape']}")

Best Practices
--------------

1. **Always check data validity**

   .. code-block:: python

      if x is None or y is None:
          print("Failed to load data")
          return

2. **Handle both 1D and 2D data**

   .. code-block:: python

      if isinstance(x, list) and len(x) > 1:
          # 2D data processing
          pass
      else:
          # 1D data processing
          pass

3. **Use magnitude for complex data visualization**

   .. code-block:: python

      y_display = np.abs(y) if np.iscomplexobj(y) else y

4. **Preserve metadata**

   .. code-block:: python

      # Always keep parameters for reproducibility
      analysis_info = {
          'original_file': str(filepath),
          'parameters': params,
          'processing_date': datetime.now().isoformat()
      }

5. **Validate critical parameters**

   .. code-block:: python

      required_params = ['MWFQ', 'HCF', 'HSW']
      missing = [p for p in required_params if p not in params]
      if missing:
          print(f"Warning: Missing parameters {missing}")

Next Steps
----------

After mastering basic loading:

1. **Try Baseline Correction**: Learn to remove baseline drift
2. **Explore FAIR Conversion**: Convert data to open formats
3. **Advanced Visualization**: Create publication-quality plots
4. **Quantitative Analysis**: Extract g-factors and coupling constants

Additional Resources
--------------------

* **Example Scripts**: ``examples/scripts/01_basic_loading.py``
* **Sample Data**: ``examples/data/`` contains real EPR measurements
* **API Reference**: Complete function documentation
* **Community**: GitHub issues and discussions for help
