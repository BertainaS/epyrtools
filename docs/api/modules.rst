Complete Module Reference
=========================

This page provides a complete overview of all EPyR Tools modules and submodules.

.. autosummary::
   :toctree: generated/
   :recursive:

   epyr

Main Package
------------

.. automodule:: epyr
   :members:
   :undoc-members:
   :show-inheritance:

Module Tree
-----------

epyr
├── eprload           # EPR data loading
├── baseline/         # Baseline correction algorithms  
│   ├── _1d          # 1D baseline correction
│   ├── _2d          # 2D baseline correction
│   └── _utils       # Baseline utilities
├── fair/            # FAIR data conversion
│   ├── conversion   # Main conversion functions
│   ├── data_processing  # Data processing utilities
│   ├── exporters    # Export to CSV/JSON/HDF5
│   └── parameter_mapping  # Bruker parameter mapping
├── plot             # EPR-specific plotting
├── constants        # Physical constants and conversions
├── isotope_gui/     # Interactive isotope database
│   ├── main_window  # Main GUI application
│   ├── gui_components  # GUI widgets and components
│   ├── isotope_data    # Nuclear isotope database
│   └── gui_helpers     # GUI utility functions
└── sub/             # Internal/legacy modules
    ├── baseline2    # Legacy baseline functions
    ├── loadBES3T    # BES3T format loader
    ├── loadESP      # ESP format loader
    ├── processing2  # Legacy processing functions
    └── utils        # Internal utilities