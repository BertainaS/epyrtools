Installation Guide
==================

EPyR Tools can be installed in several ways depending on your needs.

Requirements
------------

* Python 3.8 or higher
* NumPy >= 1.20.0
* Matplotlib >= 3.3.0
* SciPy >= 1.7.0
* h5py >= 3.1.0 (for HDF5 export)
* pandas >= 1.3.0 (for enhanced CSV handling)

Recommended Installation
------------------------

The recommended way to install EPyR Tools is using pip with an editable installation:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/BertainaS/epyrtools.git
   cd epyrtools

   # Install with all dependencies
   pip install -e .

Development Installation
------------------------

For development work, install with development dependencies:

.. code-block:: bash

   # Install with development tools
   pip install -e .[dev]

   # Or install with all optional dependencies
   pip install -e .[all]

This includes:

* **Development tools**: black, isort, flake8, pre-commit
* **Testing**: pytest, pytest-cov
* **Documentation**: sphinx, sphinx-rtd-theme, myst-parser

Manual Installation
-------------------

If you prefer to manage dependencies manually:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/BertainaS/epyrtools.git
   cd epyrtools

   # Install core dependencies
   pip install -r requirements.txt

   # Optional: Install development dependencies
   pip install -r requirements-dev.txt

Verification
------------

To verify your installation:

.. code-block:: python

   import epyr
   print(f"EPyR Tools version: {epyr.__version__}")

   # Test basic functionality
   from epyr.baseline import baseline_polynomial
   from epyr.fair import convert_bruker_to_fair
   print("✅ EPyR Tools installed successfully!")

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**ImportError: No module named 'epyr'**

* Ensure you're in the correct Python environment
* Try: ``pip install -e . --force-reinstall``
* Check that the installation completed without errors

**Missing dependencies**

* Install missing packages: ``pip install package_name``
* Or reinstall with all dependencies: ``pip install -e .[all]``

**h5py installation issues**

* On macOS with Apple Silicon: ``pip install --no-use-pep517 h5py``
* On Linux: ``sudo apt-get install libhdf5-dev`` first
* On Windows: Use conda instead: ``conda install h5py``

Platform-Specific Notes
~~~~~~~~~~~~~~~~~~~~~~~~

**Windows**

* Use Anaconda or Miniconda for easier dependency management
* Some binary packages may require Microsoft Visual C++ Build Tools

**macOS**

* Xcode Command Line Tools may be required: ``xcode-select --install``
* For Apple Silicon Macs, ensure you're using compatible package versions

**Linux**

* Install system dependencies: ``sudo apt-get install build-essential python3-dev``
* For HDF5 support: ``sudo apt-get install libhdf5-dev``

Getting Help
------------

If you encounter installation issues:

1. Check the `GitHub Issues <https://github.com/BertainaS/epyrtools/issues>`_
2. Create a new issue with your system details and error message
3. Include output of ``pip list`` and ``python --version``
