"""
Unit conversion utilities for EPR/NMR spectroscopy

Converts between common spectroscopic units: cm⁻¹, eV, K, mT, MHz
All conversions use 2022 CODATA physical constants for accuracy.
"""

from typing import Optional, Union

import numpy as np

from ..logging_config import get_logger
from .constants import bmagn, boltzm, clight, evolt, gfree, planck

logger = get_logger(__name__)


def unitconvert(
    value: Union[float, np.ndarray],
    units: str,
    g_factor: Optional[Union[float, np.ndarray]] = None,
) -> Union[float, np.ndarray]:
    """
    Convert between spectroscopic units.

    Supported conversions:
    - cm⁻¹ ↔ eV, K, mT, MHz
    - eV ↔ cm⁻¹, K, mT, MHz
    - K ↔ cm⁻¹, eV, mT, MHz
    - mT ↔ cm⁻¹, eV, K, MHz
    - MHz ↔ cm⁻¹, eV, K, mT

    Parameters
    ----------
    value : float or array
        Input value(s) to convert
    units : str
        Conversion string in format 'unit_from->unit_to'
        e.g., 'cm^-1->MHz', 'eV->mT'
    g_factor : float or array, optional
        g-factor for magnetic field conversions
        Defaults to free electron g-factor (2.002319...)

    Returns
    -------
    float or array
        Converted value(s)

    Examples
    --------
    >>> # Convert wavenumbers to frequency
    >>> freq = unitconvert(1000, 'cm^-1->MHz')  # 1000 cm⁻¹ to MHz
    >>> print(f"{freq:.3f} MHz")

    >>> # Convert with custom g-factor
    >>> field = unitconvert(100, 'cm^-1->mT', g_factor=2.005)
    >>> print(f"{field:.2f} mT")

    >>> # Vector conversion
    >>> energies = np.array([100, 200, 300])  # cm⁻¹
    >>> temps = unitconvert(energies, 'cm^-1->K')
    >>> print(f"Temperatures: {temps}")
    """
    if g_factor is None:
        g_factor = gfree()

    # Dictionary of conversion functions
    conversions = {
        # From cm⁻¹
        "cm^-1->eV": lambda v: v * 100 * clight() * planck() / evolt(),
        "cm^-1->K": lambda v: v * 100 * clight() * planck() / boltzm(),
        "cm^-1->mT": lambda v: v
        / g_factor
        * (planck() / bmagn() / 1e-3)
        * 100
        * clight(),
        "cm^-1->MHz": lambda v: v * 100 * clight() / 1e6,
        # From eV
        "eV->cm^-1": lambda v: v * evolt() / (100 * clight() * planck()),
        "eV->K": lambda v: v * evolt() / boltzm(),
        "eV->mT": lambda v: v / g_factor / bmagn() / 1e-3 * evolt(),
        "eV->MHz": lambda v: v * evolt() / planck() / 1e6,
        # From K
        "K->cm^-1": lambda v: v * boltzm() / (100 * clight() * planck()),
        "K->eV": lambda v: v * boltzm() / evolt(),
        "K->mT": lambda v: v / g_factor / bmagn() / 1e-3 * boltzm(),
        "K->MHz": lambda v: v * boltzm() / planck() / 1e6,
        # From mT
        "mT->cm^-1": lambda v: v
        * g_factor
        / (planck() / bmagn() / 1e-3)
        / (100 * clight()),
        "mT->eV": lambda v: v * g_factor * bmagn() * 1e-3 / evolt(),
        "mT->K": lambda v: v * g_factor * bmagn() * 1e-3 / boltzm(),
        "mT->MHz": lambda v: v * g_factor * (1e-3 * bmagn() / planck() / 1e6),
        # From MHz
        "MHz->cm^-1": lambda v: v * 1e6 / (100 * clight()),
        "MHz->eV": lambda v: v * 1e6 * planck() / evolt(),
        "MHz->K": lambda v: v * 1e6 * planck() / boltzm(),
        "MHz->mT": lambda v: v / g_factor * (planck() / bmagn() / 1e-3) * 1e6,
    }

    # Check if conversion exists (case-sensitive first)
    if units in conversions:
        return conversions[units](value)

    # Check case-insensitive
    units_lower = units.lower()
    conversions_lower = {k.lower(): v for k, v in conversions.items()}

    if units_lower in conversions_lower:
        # Find the correct case
        correct_units = next(k for k in conversions.keys() if k.lower() == units_lower)
        raise ValueError(
            f"Case mismatch. You provided: '{units}'. Did you mean '{correct_units}'?"
        )

    raise ValueError(
        f"Unknown unit conversion: '{units}'. "
        f"Supported: {', '.join(sorted(conversions.keys()))}"
    )


def list_conversions():
    """List all supported unit conversions."""
    conversions = [
        "cm^-1 <-> eV, K, mT, MHz",
        "eV <-> cm^-1, K, mT, MHz",
        "K <-> cm^-1, eV, mT, MHz",
        "mT <-> cm^-1, eV, K, MHz",
        "MHz <-> cm^-1, eV, K, mT",
    ]

    logger.info("Supported Unit Conversions")
    logger.info("=" * 35)
    for conv in conversions:
        logger.info(f"  {conv}")

    logger.info("Physical Constants Used:")
    logger.info(f"  Speed of light: {clight():.0f} m⋅s⁻¹")
    logger.info(f"  Planck constant: {planck():.2e} J⋅s")
    logger.info(f"  Bohr magneton: {bmagn():.2e} J⋅T⁻¹")
    logger.info(f"  Boltzmann constant: {boltzm():.2e} J⋅K⁻¹")
    logger.info(f"  Electron volt: {evolt():.2e} J")
    logger.info(f"  Free electron g-factor: {gfree():.8f}")


def demo_conversions():
    """Demonstrate common unit conversions in EPR spectroscopy."""
    logger.info("EPR Unit Conversion Examples")
    logger.info("=" * 40)

    # Example 1: X-band EPR field calculation
    freq_ghz = 9.5  # GHz
    freq_mhz = freq_ghz * 1000
    field_mt = unitconvert(freq_mhz, "MHz->mT")
    logger.info(f"X-band EPR ({freq_ghz} GHz):")
    logger.info(f"  Resonant field: {field_mt:.1f} mT")

    # Example 2: Energy scale conversions
    energy_wn = 1000  # cm⁻¹
    energy_ev = unitconvert(energy_wn, "cm^-1->eV")
    energy_k = unitconvert(energy_wn, "cm^-1->K")
    energy_mhz = unitconvert(energy_wn, "cm^-1->MHz")

    logger.info(f"Energy scale comparisons for {energy_wn} cm⁻¹:")
    logger.info(f"  {energy_ev:.6f} eV")
    logger.info(f"  {energy_k:.1f} K")
    logger.info(f"  {energy_mhz/1000:.1f} GHz")

    # Example 3: Temperature to field conversion
    temp_k = 4.2  # Liquid helium temperature
    temp_wn = unitconvert(temp_k, "K->cm^-1")
    temp_field = unitconvert(temp_k, "K->mT")

    logger.info(f"Thermal energy at {temp_k} K:")
    logger.info(f"  {temp_wn:.3f} cm⁻¹")
    logger.info(f"  {temp_field:.3f} mT equivalent field")

    # Example 4: Vector conversion
    fields = np.array([100, 200, 300, 400])  # mT
    freqs = unitconvert(fields, "mT->MHz", g_factor=2.003)

    logger.info("Field to frequency conversion (g=2.003):")
    for b, f in zip(fields, freqs):
        logger.info(f"  {b} mT → {f/1000:.3f} GHz")

    # Example 5: Different g-factors
    field = 340  # mT
    g_factors = np.array([2.000, 2.003, 2.010, 2.100])
    frequencies = unitconvert(field, "mT->MHz", g_factor=g_factors)

    logger.info(f"Frequency at {field} mT for different g-factors:")
    for g, f in zip(g_factors, frequencies):
        logger.info(f"  g = {g:.3f} → {f/1000:.3f} GHz")


if __name__ == "__main__":
    list_conversions()
    logger.info("")
    demo_conversions()
