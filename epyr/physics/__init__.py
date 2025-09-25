"""
EPR Physics Module

Physical constants and unit conversion utilities for EPR/NMR spectroscopy.
All values from 2022 CODATA recommendations with proper units and uncertainties.
"""

# Import constants and functions
from .constants import (
    # SI constants (direct values)
    GFREE, BMAGN, PLANCK, HBAR, CLIGHT, BOLTZM, AVOGADRO, NMAGN, ECHARGE, EVOLT,
    # CGS constants (direct values)
    GFREE_CGS, BMAGN_CGS, PLANCK_CGS, HBAR_CGS, CLIGHT_CGS, BOLTZM_CGS, AVOGADRO_CGS, NMAGN_CGS, ECHARGE_CGS, EVOLT_CGS,
    # Backward compatibility functions
    gfree, bmagn, planck, hbar, clight, boltzm, avogadro, nmagn, echarge, evolt,
    # EPR-specific functions
    gamma_hz, magnetic_field_to_frequency, frequency_to_magnetic_field,
    thermal_energy, wavelength_to_frequency, constants_summary
)

# Import unit conversion utilities
from .units import (
    unitconvert, list_conversions, demo_conversions
)

# Import direct conversion functions
from .conversions import (
    mhz_to_mt, mt_to_mhz, cm_inv_to_mhz, mhz_to_cm_inv,
    frequency_field_conversion_table, energy_conversion_table
)

__all__ = [
    # SI constants (direct values - preferred)
    'GFREE', 'BMAGN', 'PLANCK', 'HBAR', 'CLIGHT', 'BOLTZM', 'AVOGADRO', 'NMAGN', 'ECHARGE', 'EVOLT',

    # CGS constants (direct values)
    'GFREE_CGS', 'BMAGN_CGS', 'PLANCK_CGS', 'HBAR_CGS', 'CLIGHT_CGS', 'BOLTZM_CGS', 'AVOGADRO_CGS', 'NMAGN_CGS', 'ECHARGE_CGS', 'EVOLT_CGS',

    # Backward compatibility functions
    'gfree', 'bmagn', 'planck', 'hbar', 'clight', 'boltzm', 'avogadro', 'nmagn', 'echarge', 'evolt',

    # EPR-specific functions
    'gamma_hz', 'magnetic_field_to_frequency', 'frequency_to_magnetic_field',
    'thermal_energy', 'wavelength_to_frequency',

    # Unit conversions
    'unitconvert', 'list_conversions', 'demo_conversions',

    # Direct conversion functions
    'mhz_to_mt', 'mt_to_mhz', 'cm_inv_to_mhz', 'mhz_to_cm_inv',
    'frequency_field_conversion_table', 'energy_conversion_table',

    # Utilities
    'constants_summary'
]

# Version info
__version__ = "0.2.0"
__author__ = "EPyR Tools Development Team"