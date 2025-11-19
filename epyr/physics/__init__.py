"""
EPR Physics Module

Physical constants and unit conversion utilities for EPR/NMR spectroscopy.
All values from 2022 CODATA recommendations with proper units and uncertainties.
"""

# Import constants and functions
from .constants import (  # SI constants (direct values); CGS constants (direct values); Backward compatibility functions; EPR-specific functions
    AVOGADRO,
    AVOGADRO_CGS,
    BMAGN,
    BMAGN_CGS,
    BOLTZM,
    BOLTZM_CGS,
    CLIGHT,
    CLIGHT_CGS,
    ECHARGE,
    ECHARGE_CGS,
    EVOLT,
    EVOLT_CGS,
    GFREE,
    GFREE_CGS,
    HBAR,
    HBAR_CGS,
    NMAGN,
    NMAGN_CGS,
    PLANCK,
    PLANCK_CGS,
    avogadro,
    bmagn,
    boltzm,
    clight,
    constants_summary,
    echarge,
    evolt,
    frequency_to_magnetic_field,
    gamma_hz,
    gfree,
    hbar,
    magnetic_field_to_frequency,
    nmagn,
    planck,
    thermal_energy,
    wavelength_to_frequency,
)

# Import direct conversion functions
from .conversions import (
    cm_inv_to_mhz,
    energy_conversion_table,
    frequency_field_conversion_table,
    mhz_to_cm_inv,
    mhz_to_mt,
    mt_to_mhz,
)

# Import unit conversion utilities
from .units import demo_conversions, list_conversions, unitconvert

__all__ = [
    # SI constants (direct values - preferred)
    "GFREE",
    "BMAGN",
    "PLANCK",
    "HBAR",
    "CLIGHT",
    "BOLTZM",
    "AVOGADRO",
    "NMAGN",
    "ECHARGE",
    "EVOLT",
    # CGS constants (direct values)
    "GFREE_CGS",
    "BMAGN_CGS",
    "PLANCK_CGS",
    "HBAR_CGS",
    "CLIGHT_CGS",
    "BOLTZM_CGS",
    "AVOGADRO_CGS",
    "NMAGN_CGS",
    "ECHARGE_CGS",
    "EVOLT_CGS",
    # Backward compatibility functions
    "gfree",
    "bmagn",
    "planck",
    "hbar",
    "clight",
    "boltzm",
    "avogadro",
    "nmagn",
    "echarge",
    "evolt",
    # EPR-specific functions
    "gamma_hz",
    "magnetic_field_to_frequency",
    "frequency_to_magnetic_field",
    "thermal_energy",
    "wavelength_to_frequency",
    # Unit conversions
    "unitconvert",
    "list_conversions",
    "demo_conversions",
    # Direct conversion functions
    "mhz_to_mt",
    "mt_to_mhz",
    "cm_inv_to_mhz",
    "mhz_to_cm_inv",
    "frequency_field_conversion_table",
    "energy_conversion_table",
    # Utilities
    "constants_summary",
]

# Version info
__version__ = "0.3.1"
__author__ = "EPyR Tools Development Team"
