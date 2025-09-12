epyr.fair module
================

The ``epyr.fair`` module provides functionality for converting proprietary Bruker EPR data into FAIR (Findable, Accessible, Interoperable, and Reusable) formats.

Overview
--------

FAIR data principles ensure that scientific data is:

* **Findable**: Data has metadata and identifiers
* **Accessible**: Data can be retrieved with standard protocols  
* **Interoperable**: Data uses standard vocabularies and formats
* **Reusable**: Data has clear licenses and provenance

This module converts Bruker EPR files to open formats: CSV, JSON, and HDF5.

Main Functions
--------------

Conversion Functions
~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.fair.conversion
   :members: convert_bruker_to_fair
   :undoc-members:
   :show-inheritance:

Data Processing
~~~~~~~~~~~~~~~

.. automodule:: epyr.fair.data_processing
   :members:
   :undoc-members:
   :show-inheritance:

Export Functions
~~~~~~~~~~~~~~~~

.. automodule:: epyr.fair.exporters
   :members: save_to_csv_json, save_to_hdf5
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

Basic Conversion
~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.fair import convert_bruker_to_fair
   
   # Convert single file (creates CSV, JSON, HDF5)
   convert_bruker_to_fair('spectrum.dsc')
   
   # Specify output directory
   convert_bruker_to_fair('spectrum.dsc', output_dir='./fair_data/')
   
   # Batch conversion
   convert_bruker_to_fair()  # Opens file dialog for selection

Output Formats
--------------

CSV Format
~~~~~~~~~~
Simple comma-separated values for data:

.. code-block:: text

   # EPR Spectrum Data
   # Original file: spectrum.dsc
   # Converted: 2025-09-02
   field_gauss,intensity_au
   3200.0,0.123
   3201.0,0.125
   ...

JSON Format  
~~~~~~~~~~~
Structured metadata with human-readable parameter names:

.. code-block:: json

   {
     "original_file": "spectrum.dsc",
     "conversion_info": {
       "timestamp": "2025-09-02T10:30:00",
       "epyr_version": "0.1.6"
     },
     "measurement_parameters": {
       "microwave_frequency": {
         "value": 9.4e9,
         "unit": "Hz",
         "description": "Microwave frequency"
       },
       "magnetic_field_center": {
         "value": 3350.0,
         "unit": "G", 
         "description": "Center magnetic field"
       }
     },
     "data": {
       "field_axis": [3200.0, 3201.0, ...],
       "intensity": [0.123, 0.125, ...]
     }
   }

HDF5 Format
~~~~~~~~~~~
Self-contained hierarchical format with full metadata:

.. code-block:: text

   spectrum.h5
   ├── data/
   │   ├── intensity          # EPR signal data
   │   ├── field_axis        # Magnetic field axis
   │   └── ...
   ├── metadata/
   │   ├── parameters_fair/   # FAIR-mapped parameters
   │   │   ├── microwave_frequency/
   │   │   ├── magnetic_field_center/
   │   │   └── ...
   │   └── parameters_original/  # Unmapped original parameters
   └── attributes             # Global metadata

Parameter Mapping
-----------------

The FAIR converter maps Bruker parameter names to standardized, human-readable names:

.. automodule:: epyr.fair.parameter_mapping
   :members:
   :undoc-members:
   :show-inheritance:

Common Parameter Mappings:

* ``MWFQ`` → ``microwave_frequency`` (Hz)
* ``MWPW`` → ``microwave_power`` (dB)  
* ``AVGS`` → ``number_of_averages``
* ``RCAG`` → ``receiver_gain`` (dB)
* ``STMP`` → ``sample_temperature`` (K)

Complete API
------------

.. automodule:: epyr.fair
   :members:
   :undoc-members:
   :show-inheritance: