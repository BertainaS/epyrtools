epyr.plot module
================

``epyr.plot`` is a backward-compatible alias for :mod:`epyr.eprplot`, EPyR Tools'
plotting module for EPR spectroscopy data.

Overview
--------

* **1D plots**: absorption/derivative EPR traces
* **2D color maps**: field vs. angle/frequency/time
* **Waterfall plots**: stacked series of 1D spectra
* **Interactive slicer**: browse a 2D dataset one slice at a time

Main Functions
--------------

.. automodule:: epyr.eprplot
   :members: plot_1d, plot_2d_map, plot_2d_waterfall, plot_2d_slicer
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

1D Plotting
~~~~~~~~~~~

.. code-block:: python

   from epyr.plot import plot_1d
   from epyr import eprload

   x, y, params, _ = eprload('spectrum.DTA')
   fig, ax = plot_1d(x, y, params, title='EPR Spectrum')

2D Color Map
~~~~~~~~~~~~

.. code-block:: python

   from epyr.plot import plot_2d_map
   import numpy as np

   field_axis = np.linspace(3200, 3400, 200)  # Gauss
   angle_axis = np.linspace(0, 180, 37)       # Degrees
   epr_2d_data = load_your_2d_data()          # Shape: (37, 200)

   fig, ax = plot_2d_map(
       [field_axis, angle_axis], epr_2d_data,
       title='EPR Angular Dependence', cmap='magma',
   )

Waterfall Plot
~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.plot import plot_2d_waterfall

   fig, ax = plot_2d_waterfall(
       [field_axis, angle_axis], epr_2d_data,
       offset_factor=0.5, max_traces=20,
   )

Interactive Slicer
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.plot import plot_2d_slicer

   plot_2d_slicer([field_axis, angle_axis], epr_2d_data, slice_direction='horizontal')

Complete API
------------

.. automodule:: epyr.eprplot
   :members:
   :undoc-members:
   :show-inheritance:
