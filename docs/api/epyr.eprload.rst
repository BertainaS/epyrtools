epyr.eprload module
===================

The ``epyr.eprload`` module provides the main interface for loading EPR data files from Bruker spectrometers.

Main Function
-------------

.. automodule:: epyr.eprload
   :members: eprload
   :undoc-members:
   :show-inheritance:

The ``eprload()`` function is the primary entry point for loading EPR data. It automatically detects the file format and calls the appropriate loader function.

Usage Examples
--------------

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   import epyr.eprload as eprload
   
   # Open file dialog to select EPR data
   x, y, params, filepath = eprload.eprload()
   
   # Or specify file directly
   x, y, params, filepath = eprload.eprload('spectrum.dsc')

With Parameters
~~~~~~~~~~~~~~~

.. code-block:: python

   # Disable plotting
   x, y, params, filepath = eprload.eprload('data.dsc', plot_if_possible=False)
   
   # Apply scaling during load
   x, y, params, filepath = eprload.eprload('data.dsc', scaling='nPGT')
   
   # Save processed data
   x, y, params, filepath = eprload.eprload('data.dsc', save_if_possible=True)

Supported File Formats
-----------------------

BES3T Format
~~~~~~~~~~~~
- **.dsc files**: Descriptor files containing measurement parameters
- **.dta files**: Binary data files containing spectral data

ESP/WinEPR Format  
~~~~~~~~~~~~~~~~~
- **.par files**: Parameter files with measurement settings
- **.spc files**: Binary spectrum data files

Return Values
-------------

The ``eprload()`` function returns a tuple containing:

* **x** (*numpy.ndarray*): X-axis data (typically magnetic field in Gauss)
* **y** (*numpy.ndarray*): Y-axis data (EPR signal intensity) 
* **params** (*dict*): Dictionary of measurement parameters from the file
* **filepath** (*str*): Path to the loaded data file

All Functions
-------------

.. automodule:: epyr.eprload
   :members:
   :undoc-members:
   :show-inheritance: