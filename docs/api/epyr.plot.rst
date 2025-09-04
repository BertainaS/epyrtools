epyr.plot module
================

The ``epyr.plot`` module provides specialized plotting functions for EPR spectroscopy data visualization.

Overview
--------

This module offers publication-quality plotting tools specifically designed for EPR data:

* **2D spectral maps**: Field vs. angle/frequency/time plots
* **Waterfall plots**: Stacked 1D spectra
* **EPR-specific styling**: Appropriate axis labels, colormaps, and formatting
* **Interactive plotting**: Integration with matplotlib for customization

Main Classes and Functions
--------------------------

Plot Configuration
~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.plot
   :members: EPRPlotConfig
   :undoc-members:
   :show-inheritance:

2D Plotting Functions
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.plot
   :members: plot_2d_spectral_map
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

2D Spectral Mapping
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.plot import plot_2d_spectral_map
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Create 2D EPR data (e.g., field vs. angle)
   field_axis = np.linspace(3200, 3400, 200)  # Gauss
   angle_axis = np.linspace(0, 180, 37)       # Degrees
   
   # Your 2D EPR data array (shape: angle x field)
   epr_2d_data = load_your_2d_data()  # Shape: (37, 200)
   
   # Create 2D plot
   fig, ax = plot_2d_spectral_map(
       field_axis, angle_axis, epr_2d_data,
       x_unit='G', y_unit='°',
       title='EPR Angular Dependence'
   )
   
   plt.show()

Custom Styling
~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.plot import EPRPlotConfig
   
   # Configure EPR-specific plot settings
   config = EPRPlotConfig()
   
   # Use configuration for consistent styling
   fig, ax = plot_2d_spectral_map(
       x_axis, y_axis, data_2d,
       figsize=config.DEFAULT_FIGSIZE_2D,
       cmap=config.DEFAULT_CMAP
   )

Advanced Plotting Options
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Customize colormap and ranges
   fig, ax = plot_2d_spectral_map(
       field_axis, angle_axis, epr_2d_data,
       x_unit='mT', y_unit='°',
       title='Single Crystal EPR',
       cmap='RdBu_r',  # Red-blue colormap
       vmin=-1000, vmax=1000,  # Set intensity range
       interpolation='bilinear'
   )
   
   # Add custom colorbar label
   cbar = ax.collections[0].colorbar
   cbar.set_label('EPR Signal (a.u.)', rotation=270, labelpad=20)
   
   # Save high-resolution figure
   fig.savefig('epr_2d.png', dpi=300, bbox_inches='tight')

Plot Types
----------

2D Spectral Maps
~~~~~~~~~~~~~~~~
For visualizing 2D EPR datasets:

* **Field vs. Angle**: Single crystal EPR measurements
* **Field vs. Time**: Time-resolved EPR studies  
* **Field vs. Temperature**: Variable temperature EPR
* **Field vs. Frequency**: Multi-frequency EPR

Waterfall Plots
~~~~~~~~~~~~~~~
For displaying series of 1D spectra:

* **Angular series**: Multiple orientations
* **Temperature series**: Variable temperature studies
* **Time series**: Kinetic measurements
* **Power series**: Saturation studies

Styling Options
---------------

Colormaps
~~~~~~~~~
EPR-appropriate colormaps:

* ``viridis``: Default, perceptually uniform
* ``plasma``: High contrast for weak signals
* ``RdBu_r``: Red-blue for absorption/emission
* ``seismic``: Blue-white-red for derivatives

Axis Labels
~~~~~~~~~~~
Automatic unit formatting:

* **Magnetic field**: G, mT, T
* **Frequency**: Hz, MHz, GHz  
* **Angle**: degrees, radians
* **Temperature**: K, °C
* **Time**: s, ms, μs

Plot Configuration
------------------

.. automodule:: epyr.plot
   :members: EPRPlotConfig
   :undoc-members:
   :show-inheritance:

The ``EPRPlotConfig`` class provides default settings optimized for EPR data visualization.

Complete API
------------

.. automodule:: epyr.plot
   :members:
   :undoc-members:
   :show-inheritance: