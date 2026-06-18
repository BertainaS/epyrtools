Release Notes
=============

.. toctree::
   :maxdepth: 2
   :caption: Version History:

   release_notes/v0.4.0
   release_notes/v0.3.9
   release_notes/v0.3.8
   release_notes/v0.3.6
   release_notes/v0.3.5
   release_notes/v0.2.0

Version 0.4.0 (Latest)
-----------------------

**Release Date:** June 2026

New ``epyr.relaxation`` package for T1/T2 relaxation fitting: six
decay/recovery models (mono- and stretched-exponential, bi-exponential,
inversion/saturation recovery, Gamma0/GammaG), ``fit_relaxation()``, and
``fit_multiple_decays()`` (ranked by reduced chi-squared). See
:doc:`release_notes/v0.4.0`.

Version 0.3.9
--------------

**Release Date:** May 2026

Bruker loader refactor (cyclomatic complexity 60 and 89 split into small
private helpers), full type annotations on ``eprload()``, explicit
imports in ``epyr/__init__.py``, NumPy docstrings with worked examples
on 53% of the public API, and tooling fixes (``make format`` and
``make quality-fix`` both return 0). See :doc:`release_notes/v0.3.9`.

Version 0.3.8
-------------

**Release Date:** April 2026

Documentation-quality release with a fitting-module refactor and CLI
output fixes. ``epyr.lineshapes.fitting`` split into focused private
helpers; public API unchanged. See :doc:`release_notes/v0.3.8`.

Version 0.3.6
-------------

**Release Date:** February 2026

Maintenance release: dead-code removal, test-API alignment with the
lineshape refactor, ``plot_2d_waterfall`` ``clip_factor`` option, and
``.flake8`` tuning. See :doc:`release_notes/v0.3.6`.

Version 0.3.5
-------------

**Release Date:** February 2026

Analytical Voigt derivatives via Faddeeva-function identities, pseudo-Voigt
derivative and phase support, matplotlib 3.8+ compatibility for
``plot_2d_map`` and ``vmin``/``vmax`` color scale control. See
:doc:`release_notes/v0.3.5`.

Previous Versions
-----------------

Version 0.3.1
~~~~~~~~~~~~~

- Automatic file extension detection for EPR data loading
- Affine baseline fitting option in ``fit_epr_signal()``
- Fixed filename parsing for files with multiple dots

Version 0.3.0
~~~~~~~~~~~~~

- Migrated 398 print statements to centralized logging system
- Structured logging infrastructure across 18 core modules

Version 0.2.x
~~~~~~~~~~~~~

- Stable core functionality with validated test suite
- FAIR data conversion with 235+ parameter mappings
- 2D FFT analysis module
- ``return_type`` parameter in ``eprload()``

For complete details, see :doc:`release_notes/v0.2.0`.

Version 0.1.x
~~~~~~~~~~~~~

Early development versions focusing on core EPR data loading and basic
analysis capabilities.

- Initial Bruker file format support
- Basic plotting functionality
- Foundation plugin architecture
- Early FAIR data conversion
