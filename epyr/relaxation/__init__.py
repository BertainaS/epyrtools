"""
EPyR Tools - Relaxation Fitting Module

Time-domain decay and recovery models for T1/T2 EPR relaxation measurements,
with a curve-fitting layer mirroring epyr.lineshapes.fitting.
"""

from .fitting import (
    RelaxationFitComparison,
    RelaxationFitResult,
    fit_multiple_decays,
    fit_relaxation,
)
from .models import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)

__all__ = [
    "mono_exponential",
    "stretched_exponential",
    "biexponential",
    "inversion_recovery",
    "saturation_recovery",
    "gamma_gaussian_decay",
    "fit_relaxation",
    "fit_multiple_decays",
    "RelaxationFitResult",
    "RelaxationFitComparison",
]
