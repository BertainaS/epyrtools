Release Notes
=============

.. toctree::
   :maxdepth: 2
   :caption: Version History:

   release_notes/v0.3.5
   release_notes/v0.2.0

Version 0.3.5 (Latest)
----------------------

**Release Date:** February 2026

**Status:** Stable release with analytical Voigt derivatives and plotting improvements

EPyR Tools v0.3.5 focuses on **lineshape fitting accuracy** and **visualization quality**.

Key Highlights
~~~~~~~~~~~~~~

**Lineshape Fitting**

- Analytical Voigt derivatives via Faddeeva function identities (replaces numerical approach)
- Pseudo-Voigt now fully supports derivatives and phase fitting
- Robust initial parameter estimation for 1st and 2nd derivative signals

**Visualization**

- Fixed ``plot_2d_map`` colorbar crash with matplotlib 3.8+
- Added ``vmin``/``vmax`` color scale control for 2D maps
- Thinner default line width (0.75) for cleaner plots
- JPG export: 200 DPI, ``RdBu_r`` colormap, symmetric 98th percentile scale

**Bug Fixes**

- Fixed ``baseline_polynomial_2d`` meshgrid shape mismatch
- Fixed ``lshape()`` Dawson function returning complex values with phase

For complete details, see :doc:`release_notes/v0.3.5`.

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

Early development versions focusing on core EPR data loading and basic analysis capabilities.

- Initial Bruker file format support
- Basic plotting functionality
- Foundation plugin architecture
- Early FAIR data conversion
