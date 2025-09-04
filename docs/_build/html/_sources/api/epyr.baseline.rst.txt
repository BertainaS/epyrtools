epyr.baseline module
====================

The ``epyr.baseline`` module provides functions for correcting baseline distortions in EPR spectra.

Overview
--------

Baseline correction is essential for EPR spectral analysis. This module provides several algorithms:

* **Polynomial correction**: Fits polynomial baselines (constant, linear, quadratic, etc.)
* **Exponential correction**: Fits mono- and stretched-exponential decay baselines
* **2D correction**: Specialized algorithms for 2D EPR datasets

Main Functions
--------------

1D Baseline Correction
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.baseline._1d
   :members: baseline_polynomial, baseline_constant_offset, baseline_mono_exponential, baseline_stretched_exponential
   :undoc-members:
   :show-inheritance:

2D Baseline Correction  
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.baseline._2d
   :members: baseline_polynomial_2d
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

Polynomial Baseline Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.baseline import baseline_polynomial
   import numpy as np
   
   # Load your EPR data (x, y arrays)
   x = np.linspace(3200, 3400, 1000)  # Magnetic field in G
   y = epr_spectrum  # Your EPR data
   
   # Linear baseline correction
   y_corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=1)
   
   # Exclude signal regions from baseline fit
   signal_regions = [(3340, 3360), (3380, 3400)]
   y_corrected, baseline = baseline_polynomial(
       y, x_data=x, poly_order=2, exclude_regions=signal_regions
   )

Exponential Baseline Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.baseline import baseline_mono_exponential
   
   # For spectra with exponential baseline drift
   y_corrected, baseline = baseline_mono_exponential(
       y, x_data=x, exclude_regions=[(3340, 3380)]
   )

2D Baseline Correction
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.baseline import baseline_polynomial_2d
   import numpy as np
   
   # 2D EPR data (e.g., field vs. angle)
   field_axis = np.linspace(3200, 3400, 200)
   angle_axis = np.linspace(0, 180, 37) 
   spectrum_2d = epr_data_2d  # Shape: (37, 200)
   
   # 2D polynomial baseline correction
   corrected_2d, baseline_2d = baseline_polynomial_2d(
       spectrum_2d, x_axis=field_axis, y_axis=angle_axis, 
       poly_order_x=1, poly_order_y=1
   )

Baseline Algorithms
-------------------

Polynomial Baseline
~~~~~~~~~~~~~~~~~~~
Fits polynomials of specified order to the baseline:

* **Order 0**: Constant offset correction
* **Order 1**: Linear drift correction  
* **Order 2**: Quadratic baseline correction
* **Higher orders**: For complex baseline shapes

Exponential Baseline
~~~~~~~~~~~~~~~~~~~~
For spectra with exponential decay components:

* **Mono-exponential**: ``y = A * exp(-x/τ) + C``
* **Stretched-exponential**: ``y = A * exp(-(x/τ)^β) + C``

Utilities
---------

.. automodule:: epyr.baseline._utils
   :members:
   :undoc-members:
   :show-inheritance:

Complete API
------------

.. automodule:: epyr.baseline
   :members:
   :undoc-members:
   :show-inheritance: