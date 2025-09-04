epyr.constants module
====================

The ``epyr.constants`` module provides physical constants, unit conversions, and nuclear data relevant to EPR spectroscopy.

Overview
--------

This module contains:

* **Fundamental constants**: Planck constant, Bohr magneton, etc.
* **EPR-specific constants**: g-factors, magnetogyric ratios
* **Unit conversions**: Gauss ↔ Tesla, frequency ↔ field
* **Nuclear isotope data**: Spins, abundances, NMR frequencies

Physical Constants
------------------

Fundamental Constants
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.constants
   :members: PLANCK_CONSTANT, BOHR_MAGNETON, NUCLEAR_MAGNETON, ELECTRON_G_FACTOR
   :undoc-members:
   :show-inheritance:

EPR-Specific Constants
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.constants
   :members: MHZ_PER_GAUSS_PER_G_FACTOR
   :undoc-members:
   :show-inheritance:

Unit Conversions
----------------

Magnetic Field Units
~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.constants
   :members: GAUSS_TO_TESLA, TESLA_TO_GAUSS
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

Basic Constants
~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr import constants
   
   # Fundamental constants
   h = constants.PLANCK_CONSTANT      # 6.626e-34 J⋅s
   μ_B = constants.BOHR_MAGNETON      # 9.274e-24 J/T
   g_e = constants.ELECTRON_G_FACTOR  # ~2.002
   
   # Calculate EPR frequency
   B_field = 3350  # Gauss
   B_tesla = B_field * constants.GAUSS_TO_TESLA
   frequency = g_e * μ_B * B_tesla / h
   print(f"EPR frequency: {frequency/1e9:.3f} GHz")

Field-Frequency Conversions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Convert between magnetic field and EPR frequency
   field_gauss = 3350  # Typical X-band EPR field
   g_factor = 2.0      # Free electron g-factor
   
   # Calculate corresponding frequency
   freq_hz = field_gauss * g_factor * constants.MHZ_PER_GAUSS_PER_G_FACTOR * 1e6
   print(f"Frequency: {freq_hz/1e9:.2f} GHz")
   
   # Unit conversions
   field_tesla = field_gauss * constants.GAUSS_TO_TESLA
   field_mt = field_tesla * 1000  # Convert to mT
   print(f"Field: {field_gauss} G = {field_tesla:.4f} T = {field_mt:.1f} mT")

g-Factor Calculations
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate g-factor from EPR parameters
   frequency_hz = 9.4e9  # 9.4 GHz X-band
   field_gauss = 3350    # Resonance field
   
   # g = hν / (μ_B × B)
   field_tesla = field_gauss * constants.GAUSS_TO_TESLA
   g_factor = (constants.PLANCK_CONSTANT * frequency_hz) / \
              (constants.BOHR_MAGNETON * field_tesla)
   
   print(f"g-factor: {g_factor:.6f}")

Nuclear Isotope Data
--------------------

Isotope Database
~~~~~~~~~~~~~~~~

.. automodule:: epyr.constants
   :members: ISOTOPE_DATA
   :undoc-members:
   :show-inheritance:

Common EPR-Relevant Isotopes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The isotope database includes nuclear properties for common isotopes encountered in EPR:

**Hydrogen isotopes:**
- ¹H: I=1/2, 99.98% abundance, γ=42.577 MHz/T
- ²H (D): I=1, 0.02% abundance, γ=6.536 MHz/T

**Carbon isotopes:**
- ¹²C: I=0, 98.9% abundance (EPR silent)
- ¹³C: I=1/2, 1.1% abundance, γ=10.705 MHz/T

**Nitrogen isotopes:**
- ¹⁴N: I=1, 99.63% abundance, γ=3.078 MHz/T
- ¹⁵N: I=1/2, 0.37% abundance, γ=-4.316 MHz/T

Usage Example:

.. code-block:: python

   # Access nuclear data
   if hasattr(constants, 'ISOTOPE_DATA'):
       isotopes = constants.ISOTOPE_DATA
       
       # Get proton data
       if 'H' in isotopes and '1' in isotopes['H']:
           proton = isotopes['H']['1']
           spin = proton.get('spin', 0.5)
           gamma = proton.get('gamma', 42.577)  # MHz/T
           abundance = proton.get('abundance', 99.98)  # %
           
           print(f"¹H: I={spin}, γ={gamma} MHz/T, {abundance}% abundant")

EPR Parameter Calculations
--------------------------

Hyperfine Coupling
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate hyperfine coupling in different units
   A_gauss = 10.0  # Hyperfine coupling in Gauss
   A_mhz = A_gauss * constants.MHZ_PER_GAUSS_PER_G_FACTOR
   A_tesla = A_gauss * constants.GAUSS_TO_TESLA
   
   print(f"Hyperfine coupling: {A_gauss} G = {A_mhz:.2f} MHz = {A_tesla*1e4:.1f} mT")

Nuclear Zeeman Effect
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Nuclear Zeeman splitting for ¹H in EPR field
   B_field = 3350 * constants.GAUSS_TO_TESLA  # Tesla
   gamma_proton = 42.577e6  # Hz/T (magnetogyric ratio)
   
   nuclear_freq = gamma_proton * B_field  # Nuclear Larmor frequency
   print(f"Proton Larmor frequency: {nuclear_freq/1e6:.2f} MHz")

Complete API
------------

.. automodule:: epyr.constants
   :members:
   :undoc-members:
   :show-inheritance: