Example Scripts
===============

EPyR Tools includes several ready-to-use Python scripts that demonstrate common EPR analysis workflows. These scripts are located in ``examples/scripts/`` and can be run directly or adapted for your own analysis.

Running the Scripts
-------------------

All scripts are self-contained and can be run from the command line:

.. code-block:: bash

   cd epyrtools
   python examples/scripts/01_basic_loading.py
   python examples/scripts/02_baseline_correction.py
   python examples/scripts/03_fair_conversion.py

The scripts will automatically find EPR data files in ``examples/data/`` and process them.

Script 1: Basic Loading (01_basic_loading.py)
----------------------------------------------

**Purpose**: Demonstrates automatic EPR data loading and visualization

**Features**:
- Finds all EPR files in the data directory
- Handles both 1D and 2D data formats automatically
- Creates professional-quality plots
- Extracts and displays key experimental parameters
- Saves plots with descriptive filenames

**Output**:
- Console: Data statistics and parameter information
- Files: PNG plots saved in the scripts directory

**Key Code Sections**:

.. code-block:: python

   # Automatic file discovery
   for dsc_file in data_dir.glob("*.dsc"):
       dta_file = dsc_file.with_suffix(".dta")
       if dta_file.exists():
           sample_files.append(("BES3T", dsc_file))

   # Smart data type handling
   if isinstance(x, list) and len(x) > 1:
       # 2D data plotting with color maps and slices
       y_plot = np.abs(y) if np.iscomplexobj(y) else y
       plt.imshow(y_plot, aspect='auto', origin='lower')
   else:
       # 1D data plotting
       plt.plot(x, y, 'b-', linewidth=1.5)

**Sample Output**:

.. code-block:: text

   EPyR Tools - Basic Data Loading Example
   ========================================

   Loading BES3T file: 130406SB_CaWO4_Er_CW_5K_20.DSC
   ✅ Data type: 1D
   📏 Data points: 1024
   🧲 Field range: 100.0 to 6100.0 G
   📈 Signal range: -7.79e+04 to 8.75e+04
   📋 Key Parameters:
     • Microwave Frequency (Hz): 9704197000.0
     • Microwave Power (dB): 1.013e-05
   ✅ Plot saved: 130406SB_CaWO4_Er_CW_5K_20_plot.png

Script 2: Baseline Correction (02_baseline_correction.py)
----------------------------------------------------------

**Purpose**: Comprehensive baseline correction workflow with method comparison

**Features**:
- Automatically selects 1D spectra for baseline correction
- Compares multiple correction algorithms (constant, linear, quadratic)
- Demonstrates advanced exclusion region handling
- Shows visual comparison of correction methods
- Calculates quality metrics (RMS)

**Output**:
- Console: Processing information and quality metrics
- Files: Comparison plots showing before/after correction

**Key Code Sections**:

.. code-block:: python

   # Method comparison
   corrections = [
       ("Original", y, None),
       ("Constant Offset", *baseline_polynomial(y, x_data=x, poly_order=0)),
       ("Linear", *baseline_polynomial(y, x_data=x, poly_order=1)),
       ("Quadratic", *baseline_polynomial(y, x_data=x, poly_order=2)),
   ]

   # Advanced exclusion handling
   exclude_regions = [(center_field - width/2, center_field + width/2)
                     for peak_idx in peak_indices]
   y_corrected, baseline = baseline_polynomial(
       y, x_data=x, poly_order=2, exclude_regions=exclude_regions
   )

**Sample Output**:

.. code-block:: text

   EPyR Tools - Baseline Correction Example
   ==========================================
   Processing file: 130406SB_CaWO4_Er_CW_5K_20.DSC
   Using 1024 data points for baseline correction
   ✅ Comparison plot saved: 130406SB_CaWO4_Er_CW_5K_20_baseline_correction.png

   Demonstrating baseline correction with signal exclusion...
   Excluding 4 signal regions from baseline fit
   ✅ Exclusion example saved: 130406SB_CaWO4_Er_CW_5K_20_exclusion_correction.png

Script 3: FAIR Conversion (03_fair_conversion.py)
--------------------------------------------------

**Purpose**: Batch conversion of proprietary Bruker files to FAIR formats

**Features**:
- Processes entire data directories automatically
- Creates CSV, JSON, and HDF5 output formats
- Preserves complete experimental metadata
- Generates format comparison analysis
- Shows how to read converted data

**Output**:
- Console: Conversion progress and file information
- Files: FAIR format files (CSV, JSON, HDF5) and comparison plots

**Key Code Sections**:

.. code-block:: python

   # Batch processing
   for file_path in epr_files:
       x, y, params, filepath = epyr.eprload(str(file_path), plot_if_possible=False)

       # Use FAIR conversion module
       convert_bruker_to_fair(str(file_path), output_dir=str(output_dir))

       # Demonstrate reading converted formats
       demonstrate_converted_data(output_dir, base_name, is_2d, is_complex)

   # Format comparison analysis
   create_format_comparison(output_dir)

**Sample Output**:

.. code-block:: text

   EPyR Tools - FAIR Data Conversion Example
   ===========================================
   Converting Bruker EPR files to FAIR formats (CSV, JSON, HDF5)

   Found 4 EPR files for conversion:
     1. [BES3T] Rabi2D_GdCaWO4_13dB_3057G.DSC (403.9 KB)
     2. [BES3T] 130406SB_CaWO4_Er_CW_5K_20.DSC (2.7 KB)

   🔄 Processing: 130406SB_CaWO4_Er_CW_5K_20.DSC
   📊 Data type: 1D
   🔢 Complex data: No
   ✅ Converted to 3 formats:
     - 130406SB_CaWO4_Er_CW_5K_20.csv (26.7 KB)
     - 130406SB_CaWO4_Er_CW_5K_20.json (18.7 KB)
     - 130406SB_CaWO4_Er_CW_5K_20.h5 (142.8 KB)

Customizing Scripts
-------------------

The scripts are designed to be easily modified for your specific needs:

Changing Input Directory
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # In any script, modify this line:
   data_dir = examples_dir / "data"  # Original
   data_dir = Path("/path/to/your/data")  # Custom

Adding File Filters
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Filter by filename pattern
   for file_path in data_dir.glob("*CaWO4*.dsc"):
       # Process only files containing "CaWO4"

   # Filter by date
   from datetime import datetime, timedelta
   recent = datetime.now() - timedelta(days=30)
   for file_path in data_dir.glob("*.dsc"):
       if file_path.stat().st_mtime > recent.timestamp():
           # Process only recent files

Custom Parameter Extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Add custom parameters to display
   key_params = {
       "MWFQ": "Microwave Frequency (Hz)",
       "TE": "Temperature (K)",  # Add temperature
       "SAMP": "Sample Name",    # Add sample info
       # Add your parameters here
   }

   # Calculate derived quantities
   if 'HSW' in params and 'RES' in params:
       field_resolution = params['HSW'] / params['RES']
       print(f"    Field Resolution: {field_resolution:.2f} G/point")

Automated Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def process_epr_file(file_path):
       """Complete analysis pipeline for one file."""
       # 1. Load data
       x, y, params, filepath = epyr.eprload(str(file_path))

       # 2. Apply baseline correction (if 1D)
       if not isinstance(x, list):
           y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)
       else:
           y_corrected = y

       # 3. Extract key info
       analysis_result = {
           'file': file_path.name,
           'data_type': '2D' if isinstance(x, list) else '1D',
           'frequency': params.get('MWFQ'),
           'temperature': params.get('TE'),
           'snr': np.ptp(y) / np.std(y) if not isinstance(x, list) else None
       }

       # 4. Save processed data
       output_name = file_path.stem + "_processed"
       if not isinstance(x, list):
           np.savetxt(f"{output_name}.txt",
                     np.column_stack([x, y_corrected]),
                     header="Field(G) Intensity(a.u.)")

       return analysis_result

Script Integration
------------------

Using Scripts in Jupyter Notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Import script functions into notebooks
   import sys
   sys.path.append('examples/scripts')

   from script_01_basic_loading import load_and_plot_example
   from script_02_baseline_correction import baseline_correction_example

   # Run script functions interactively
   load_and_plot_example()

Creating Custom Scripts
~~~~~~~~~~~~~~~~~~~~~~~~

Use the existing scripts as templates:

.. code-block:: python

   #!/usr/bin/env python3
   """
   Custom EPR Analysis Script
   =========================

   Description of your analysis workflow
   """

   import sys
   from pathlib import Path

   # Add EPyR Tools to path
   sys.path.insert(0, str(Path(__file__).parent.parent.parent))

   import epyr
   from epyr.baseline import baseline_polynomial

   def your_analysis_function():
       """Your custom analysis workflow."""
       # Load data
       # Process data
       # Generate results
       pass

   if __name__ == "__main__":
       your_analysis_function()

Advanced Features
-----------------

The scripts demonstrate several advanced EPyR Tools features:

Complex Data Handling
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Automatic detection and handling of complex EPR data
   if np.iscomplexobj(y):
       y_plot = np.abs(y)  # Use magnitude for visualization
       data_info = "(Complex data - showing magnitude)"
   else:
       y_plot = y
       data_info = "(Real data)"

Multi-dimensional Data
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Handle 2D datasets (time series, field sweeps, etc.)
   if isinstance(x, list) and len(x) > 1:
       x_axis = x[0]  # First dimension (usually field)
       y_axis = x[1]  # Second dimension (time, power, etc.)

       # Plot 2D color map
       plt.imshow(y_plot, extent=[x_axis.min(), x_axis.max(),
                                 y_axis.min(), y_axis.max()])

Professional Plotting
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Publication-quality plot settings
   plt.figure(figsize=(12, 8), dpi=150)
   plt.plot(x, y, 'b-', linewidth=2)
   plt.xlabel('Magnetic Field (G)', fontsize=14)
   plt.ylabel('EPR Signal (a.u.)', fontsize=14)
   plt.title('EPR Spectrum', fontsize=16)
   plt.grid(True, alpha=0.3)
   plt.tight_layout()
   plt.savefig('spectrum.png', dpi=300, bbox_inches='tight')

Error Handling
~~~~~~~~~~~~~~

.. code-block:: python

   # Robust error handling for batch processing
   try:
       x, y, params, filepath = epyr.eprload(str(file_path))
       if x is None or y is None:
           print(f"❌ Failed to load {file_path.name}")
           continue
   except Exception as e:
       print(f"❌ Error processing {file_path.name}: {e}")
       continue

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Efficient processing for large datasets
   if y.size > 1e6:  # Large dataset
       print("Large dataset detected, optimizing processing...")
       # Use memory-efficient operations
       # Process data in chunks if needed

Next Steps
----------

After exploring the example scripts:

1. **Modify for Your Data**: Adapt scripts to your specific EPR experiments
2. **Create Analysis Pipelines**: Combine multiple scripts for complete workflows
3. **Automate Processing**: Set up batch processing for routine analysis
4. **Share and Collaborate**: Contribute improved scripts back to the community

The scripts provide a solid foundation for EPR data analysis and can be extended to meet virtually any analysis requirement.
