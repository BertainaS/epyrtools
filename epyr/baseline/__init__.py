# epyrtools/baseline/__init__.py
"""
Baseline Correction Sub-package
===============================

Provides a suite of functions for performing baseline correction on 1D and 2D data.
"""

from ._1d import (
    baseline_constant_offset,
    baseline_mono_exponential,
    baseline_polynomial,
    baseline_stretched_exponential,
)
from ._2d import baseline_polynomial_2d

# This controls what `from epyrtools.baseline import *` imports.
__all__ = [
    "baseline_polynomial",
    "baseline_constant_offset",
    "baseline_mono_exponential",
    "baseline_stretched_exponential",
    "baseline_polynomial_2d",
]
