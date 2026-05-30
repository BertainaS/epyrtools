Clean End-to-End Examples
=========================

The ``examples/clean/`` directory contains five short, standalone Python
scripts that exercise the public EPyR Tools API on the bundled datasets.
Each script is self-contained: copy, run, edit, repeat.

Run all five from the repository root with:

.. code-block:: bash

   python examples/clean/01_basic_loading_and_plotting.py
   python examples/clean/02_baseline_and_fitting.py
   python examples/clean/03_advanced_fft_windows.py
   python examples/clean/04_interactive_2d_slicer.py
   python examples/clean/05_rabi_frequency_analysis.py

01 -- Basic loading and plotting
--------------------------------

**File:** ``examples/clean/01_basic_loading_and_plotting.py``

Loads a 1D CW EPR spectrum (CaWO4:Er, BES3T format) and a 2D rotation
pattern (MgO, ESP format), then plots them with :func:`epyr.plot_1d`,
:func:`epyr.plot_2d_map` and :func:`epyr.plot_2d_waterfall`.

Covers:

- BES3T (``.dsc``/``.dta``) and ESP (``.par``/``.spc``) loading
- 1D vs 2D abscissa handling (``x`` returned as array vs list of arrays)
- Reading common parameters from the ``params`` dict
- Three plotting modes side-by-side

Use this as the smoke test that your install is wired correctly.

02 -- Baseline correction and lineshape fitting
-----------------------------------------------

**File:** ``examples/clean/02_baseline_and_fitting.py``

Four mini-experiments in one script:

1. Polynomial baseline on a real CW spectrum (CaWO4:Er).
2. Polynomial baseline on a synthetic Lorentzian on a sloping background.
3. First-derivative Lorentzian fit on a real CuSO4 ESP file.
4. Model comparison (Gaussian vs Lorentzian vs Voigt vs pseudo-Voigt) on
   a synthetic two-line spectrum.

Exercises :func:`epyr.baseline.baseline_polynomial_1d`,
:func:`epyr.lineshapes.fit_epr_signal` and
:func:`epyr.lineshapes.fit_multiple_shapes`.

03 -- Advanced FFT windowing
----------------------------

**File:** ``examples/clean/03_advanced_fft_windows.py``

Compares apodization windows (rectangular, Hann, Hamming, Blackman, Kaiser)
on a synthetic Rabi-like trace, then on real Rabi data. Demonstrates the
spectral-leakage / frequency-resolution trade-off.

Touches :func:`epyr.signalprocessing.analyze_frequencies` and
:func:`epyr.signalprocessing.apowin`.

04 -- Interactive 2D slicer
---------------------------

**File:** ``examples/clean/04_interactive_2d_slicer.py``

Launches :func:`epyr.plot_2d_slicer` on the 2D Rabi rotation dataset and
on the MgO angular-rotation file. Use the slider (or arrow keys, depending
on backend) to scan slices in both directions.

Requires an interactive matplotlib backend. In Jupyter, run
``%matplotlib widget`` first.

05 -- Rabi frequency analysis
-----------------------------

**File:** ``examples/clean/05_rabi_frequency_analysis.py``

Full row-by-row FFT analysis on two Rabi 2D datasets (Gd:CaWO4 at 13 dB
and 6 dB microwave attenuation). Recovers the Rabi frequency and the
expected scaling with sqrt(power), and prints a comparison table.

Useful as a template for any 2D experiment where each row is a
time-domain oscillation.
