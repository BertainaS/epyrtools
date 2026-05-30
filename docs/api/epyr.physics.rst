epyr.physics module
===================

The ``epyr.physics`` package collects CODATA 2022 physical constants, EPR/NMR
unit conversions, and a general-purpose ``unitconvert()`` interface.

It replaces the legacy flat ``epyr.constants`` module. For backwards
compatibility, ``epyr.constants`` is still importable and resolves to
``epyr.physics.constants``.

Overview
--------

* :mod:`epyr.physics.constants` -- fundamental constants in SI and CGS
  (``GFREE``, ``BMAGN``, ``PLANCK``, ``HBAR``, ``CLIGHT``, ``BOLTZM``,
  ``NMAGN``, ``ECHARGE``, ``EVOLT``, plus ``*_CGS`` variants and callable
  accessors for backwards compatibility).
* :mod:`epyr.physics.conversions` -- direct frequency/field conversions
  (``mhz_to_mt``, ``mt_to_mhz``, ``cm_inv_to_mhz``, ``mhz_to_cm_inv``)
  and pre-formatted lookup tables.
* :mod:`epyr.physics.units` -- a general ``unitconvert(value, from, to)``
  interface backed by the conversion functions.

Constants
---------

.. automodule:: epyr.physics.constants
   :members:
   :undoc-members:
   :show-inheritance:

Conversions
-----------

.. automodule:: epyr.physics.conversions
   :members:
   :undoc-members:
   :show-inheritance:

Units
-----

.. automodule:: epyr.physics.units
   :members:
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

Constants
~~~~~~~~~

.. code-block:: python

   from epyr.physics import PLANCK, BMAGN, GFREE

   # Compute the resonance frequency at B = 0.335 T with g = g_e
   B = 0.335              # Tesla
   nu = GFREE * BMAGN * B / PLANCK
   print(f"EPR frequency: {nu / 1e9:.3f} GHz")  # ~9.4 GHz

Field/frequency conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.physics import mhz_to_mt, mt_to_mhz

   # X-band
   B_mT = mhz_to_mt(9500, g_factor=2.0023)
   print(f"X-band field for g=g_e: {B_mT:.1f} mT")  # ~339 mT

   # Reverse direction
   nu_MHz = mt_to_mhz(339.0, g_factor=2.0023)
   print(f"{nu_MHz / 1000:.2f} GHz")

Energy unit cheat-sheet
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.physics import energy_conversion_table

   # Print cm^-1 / MHz / GHz / eV / K / meV mapping for typical ZFS values
   energy_conversion_table(energies_cm_inv=[1, 10, 100])

g-factor from EPR parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.physics import PLANCK, BMAGN

   nu_Hz = 9.4e9               # X-band microwave frequency
   B_T = 3350e-4               # resonance field, Tesla (3350 G)
   g = PLANCK * nu_Hz / (BMAGN * B_T)
   print(f"g = {g:.6f}")        # ~2.003

Nuclear isotope data
--------------------

Nuclear properties (spin, abundance, magnetogyric ratio) are not stored in
:mod:`epyr.physics`. They live in ``epyr/sub/isotopedata.txt`` (EasySpin
format) and are surfaced by the interactive GUI documented at
:doc:`epyr.isotope_gui`.
