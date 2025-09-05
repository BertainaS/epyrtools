Tutorials and Examples
======================

Learn EPyR Tools through hands-on tutorials and practical examples.

Interactive Notebooks
----------------------

Our Jupyter notebooks provide interactive tutorials for learning EPyR Tools:

.. toctree::
   :maxdepth: 1

   basic_loading
   baseline_correction
   fair_conversion

**Location**: ``examples/notebooks/``

To run the notebooks:

.. code-block:: bash

   cd examples/notebooks
   jupyter notebook

Basic Loading Tutorial
~~~~~~~~~~~~~~~~~~~~~~

**File**: ``01_Basic_Loading.ipynb``

Learn the fundamentals of EPR data loading:

* Loading BES3T and ESP format files
* Understanding 1D vs 2D data structures
* Handling complex EPR data
* Basic parameter extraction
* Simple visualization techniques
* Data export for external analysis

**Prerequisites**: Basic Python knowledge
**Duration**: 20-30 minutes

Baseline Correction Tutorial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Coming Soon**: ``02_Baseline_Correction.ipynb``

Master baseline correction techniques:

* Understanding baseline artifacts
* Polynomial correction methods
* Signal region exclusion
* Advanced exponential models
* Quality assessment metrics
* Best practices for different data types

**Prerequisites**: Basic Loading Tutorial
**Duration**: 30-45 minutes

FAIR Data Conversion Tutorial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Coming Soon**: ``03_FAIR_Conversion.ipynb``

Convert proprietary data to open formats:

* FAIR data principles
* CSV export with metadata
* JSON parameter documentation
* HDF5 for large datasets
* Cross-platform compatibility
* Integration with other tools

**Prerequisites**: Basic Loading Tutorial
**Duration**: 25-35 minutes

Example Scripts
---------------

Ready-to-use Python scripts for common workflows:

.. toctree::
   :maxdepth: 1

   script_examples

**Location**: ``examples/scripts/``

Basic Loading Script
~~~~~~~~~~~~~~~~~~~~

**File**: ``01_basic_loading.py``

Automated EPR data loading and visualization:

.. code-block:: bash

   python examples/scripts/01_basic_loading.py

Features:

* Finds all EPR files in data directory
* Handles both 1D and 2D data automatically
* Creates publication-quality plots
* Extracts and displays key parameters
* Saves plots with descriptive names

Baseline Correction Script
~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``02_baseline_correction.py``

Comprehensive baseline correction workflow:

.. code-block:: bash

   python examples/scripts/02_baseline_correction.py

Features:

* Multiple correction algorithms
* Visual comparison of methods
* Advanced exclusion region handling
* Quality metrics calculation
* Automated method selection

FAIR Conversion Script
~~~~~~~~~~~~~~~~~~~~~~

**File**: ``03_fair_conversion.py``

Batch conversion to FAIR formats:

.. code-block:: bash

   python examples/scripts/03_fair_conversion.py

Features:

* Processes entire data directories
* Creates multiple output formats
* Preserves complete metadata
* Generates conversion reports
* Format comparison analysis

Sample Data
-----------

The ``examples/data/`` directory contains sample EPR data files:

**BES3T Format Files**:

* ``130406SB_CaWO4_Er_CW_5K_20.DSC/.DTA`` - 1D CW-EPR spectrum
* ``Rabi2D_GdCaWO4_13dB_3057G.DSC/.DTA`` - 2D pulsed EPR (Rabi oscillation)
* ``Rabi2D_GdCaWO4_6dB_3770G_2.DSC/.DTA`` - 2D pulsed EPR (different conditions)

**ESP Format Files**:

* ``2014_03_19_MgO_300K_111_fullrotation33dB.par/.spc`` - Angular rotation study

**Data Characteristics**:

* **1D Data**: Single EPR spectra vs magnetic field
* **2D Data**: Parameter-dependent EPR (time, power, angle, etc.)
* **Complex Data**: Pulsed EPR experiments with I/Q detection
* **Real Data**: CW-EPR and integrated signals

Common Workflows
----------------

Routine EPR Analysis
~~~~~~~~~~~~~~~~~~~~

1. **Load Data**:

   .. code-block:: python

      import epyr
      x, y, params, filepath = epyr.eprload('spectrum.dsc')

2. **Check Data Type**:

   .. code-block:: python

      if isinstance(x, list):
          print(f"2D data: {y.shape}")
      else:
          print(f"1D data: {len(y)} points")

3. **Apply Baseline Correction**:

   .. code-block:: python

      from epyr.baseline import baseline_polynomial
      y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)

4. **Visualize Results**:

   .. code-block:: python

      import matplotlib.pyplot as plt
      plt.plot(x, y_corrected)
      plt.xlabel('Magnetic Field (G)')
      plt.ylabel('EPR Signal (a.u.)')
      plt.show()

5. **Export Data**:

   .. code-block:: python

      from epyr.fair import convert_bruker_to_fair
      convert_bruker_to_fair(filepath, output_dir='./results')

Research Data Management
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Batch Processing**:

   .. code-block:: python

      from pathlib import Path
      from epyr.fair import convert_bruker_to_fair

      data_dir = Path('raw_data')
      for epr_file in data_dir.glob('*.dsc'):
          convert_bruker_to_fair(epr_file, output_dir='fair_data')

2. **Quality Control**:

   .. code-block:: python

      # Check data integrity
      if x is None or y is None:
          print(f"Failed to load {filepath}")
          continue

      # Validate parameters
      required_params = ['MWFQ', 'HCF', 'HSW']
      missing = [p for p in required_params if p not in params]
      if missing:
          print(f"Missing parameters: {missing}")

3. **Metadata Documentation**:

   .. code-block:: python

      # Document experimental conditions
      metadata = {
          'sample': params.get('SAMP', 'Unknown'),
          'temperature': params.get('TE', 'Not recorded'),
          'frequency': params.get('MWFQ', 'Not specified'),
          'power': params.get('MWPW', 'Not specified'),
      }

Publication Preparation
~~~~~~~~~~~~~~~~~~~~~~~

1. **High-Quality Plots**:

   .. code-block:: python

      plt.figure(figsize=(8, 6), dpi=300)
      plt.plot(x, y_corrected, 'k-', linewidth=2)
      plt.xlabel('Magnetic Field (G)', fontsize=14)
      plt.ylabel('EPR Signal (a.u.)', fontsize=14)
      plt.tick_params(labelsize=12)
      plt.tight_layout()
      plt.savefig('figure1.png', dpi=300, bbox_inches='tight')

2. **Data Archiving**:

   .. code-block:: python

      # Create comprehensive archive
      archive_dir = Path('publication_data')
      archive_dir.mkdir(exist_ok=True)

      # Save processed data
      np.savetxt(archive_dir / 'processed_spectrum.txt',
                 np.column_stack([x, y_corrected]),
                 header='Field(G) Intensity(a.u.)')

      # Save metadata
      import json
      with open(archive_dir / 'metadata.json', 'w') as f:
          json.dump(params, f, indent=2)

3. **Reproducible Analysis**:

   .. code-block:: python

      # Document analysis parameters
      analysis_log = {
          'epyr_version': epyr.__version__,
          'baseline_method': 'polynomial',
          'baseline_order': 1,
          'processing_date': datetime.now().isoformat(),
          'input_file': str(filepath),
      }

Troubleshooting
---------------

Common Issues and Solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**"No module named 'epyr'"**

* Check installation: ``pip list | grep epyr``
* Reinstall: ``pip install -e . --force-reinstall``
* Verify Python path: ``sys.path`` should include project directory

**"Failed to load EPR file"**

* Check file pairs: ``.dsc`` needs ``.dta``, ``.par`` needs ``.spc``
* Verify file permissions and path
* Try with different file to isolate issue

**"Baseline correction fails"**

* Check if data is 1D: ``y.ndim == 1``
* Verify field data: ``len(x) == len(y)``
* Try lower polynomial order first

**"Complex data visualization issues"**

* Use magnitude: ``np.abs(complex_data)``
* Check data shape for 2D: ``isinstance(x, list)``
* Try different plotting approaches for 2D data

Getting Help
~~~~~~~~~~~~

* **Documentation**: Check API reference and examples
* **GitHub Issues**: Report bugs and ask questions
* **Community**: Join discussions and share experiences
* **Email**: Contact maintainers for specific problems

Next Steps
----------

After completing these tutorials:

1. **Explore Advanced Features**: Try custom analysis scripts
2. **Join Community**: Contribute examples and improvements
3. **Apply to Research**: Use EPyR Tools in your scientific work
4. **Share Results**: Publish reproducible EPR data and analyses
