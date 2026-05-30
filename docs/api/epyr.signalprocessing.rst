epyr.signalprocessing module
============================

The ``epyr.signalprocessing`` package provides FFT-based frequency analysis
for time-dependent EPR signals: Rabi oscillations, DEER, HYSCORE, and other
pulse-EPR experiments.

Frequency analysis
------------------

.. automodule:: epyr.signalprocessing.frequency_analysis
   :members:
   :undoc-members:
   :show-inheritance:

Apodization windows
-------------------

.. automodule:: epyr.signalprocessing.apowin
   :members:
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

1D Rabi analysis
~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr import eprload
   from epyr.signalprocessing import analyze_frequencies

   # Load a Rabi oscillation trace (time-domain)
   t, y, params, _ = eprload("examples/data/Rabi2D_GdCaWO4_13dB_3057G.DSC")
   trace = y[0]  # first row of the 2D dataset

   # FFT with DC offset removal, Hann window, 4x zero padding
   result = analyze_frequencies(
       t[0], trace,
       window="hann",
       zero_pad=4,
       remove_dc=True,
   )
   print(f"Dominant Rabi frequency: {result['dominant_frequencies'][0]:.3f} MHz")

2D processing modes
~~~~~~~~~~~~~~~~~~~

The 2D analysis supports two strategies:

* **Row-by-row 1D FFT** -- for 2D Rabi datasets, each row is an oscillation.
* **Full 2D FFT** -- for HYSCORE, exposes cross-peaks in the (f1, f2) plane.

.. code-block:: python

   from epyr.signalprocessing import analyze_frequencies_2d

   result = analyze_frequencies_2d(t, y, params, mode="row_by_row")
   print(f"Power spectrum shape: {result['power_spectrum'].shape}")
