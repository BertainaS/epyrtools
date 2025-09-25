"""
Unit conversion utilities for EPR/NMR spectroscopy

Converts between common spectroscopic units: cm⁻¹, eV, K, mT, MHz
All conversions use 2022 CODATA physical constants for accuracy.
"""

import numpy as np
from typing import Union, Optional
from .constants import (
    clight, planck, bmagn, boltzm, evolt, gfree
)


def unitconvert(value: Union[float, np.ndarray], 
                units: str, 
                g_factor: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
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
        "cm^-1->mT": lambda v: v / g_factor * (planck() / bmagn() / 1e-3) * 100 * clight(),
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
        "mT->cm^-1": lambda v: v * g_factor / (planck() / bmagn() / 1e-3) / (100 * clight()),
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
        raise ValueError(f"Case mismatch. You provided: '{units}'. Did you mean '{correct_units}'?")
    
    raise ValueError(f"Unknown unit conversion: '{units}'. "
                    f"Supported: {', '.join(sorted(conversions.keys()))}")


def list_conversions():
    """List all supported unit conversions."""
    conversions = [
        "cm^-1 <-> eV, K, mT, MHz",
        "eV <-> cm^-1, K, mT, MHz",
        "K <-> cm^-1, eV, mT, MHz",
        "mT <-> cm^-1, eV, K, MHz",
        "MHz <-> cm^-1, eV, K, mT"
    ]

    print("Supported Unit Conversions")
    print("=" * 35)
    for conv in conversions:
        print(f"  {conv}")

    print(f"\nPhysical Constants Used:")
    print(f"  Speed of light: {clight():.0f} m⋅s⁻¹")
    print(f"  Planck constant: {planck():.2e} J⋅s")
    print(f"  Bohr magneton: {bmagn():.2e} J⋅T⁻¹")
    print(f"  Boltzmann constant: {boltzm():.2e} J⋅K⁻¹")
    print(f"  Electron volt: {evolt():.2e} J")
    print(f"  Free electron g-factor: {gfree():.8f}")


def demo_conversions():
    """Demonstrate common unit conversions in EPR spectroscopy."""
    print("EPR Unit Conversion Examples")
    print("=" * 40)
    
    # Example 1: X-band EPR field calculation
    freq_ghz = 9.5  # GHz
    freq_mhz = freq_ghz * 1000
    field_mt = unitconvert(freq_mhz, 'MHz->mT')
    print(f"X-band EPR ({freq_ghz} GHz):")
    print(f"  Resonant field: {field_mt:.1f} mT")
    
    # Example 2: Energy scale conversions
    energy_wn = 1000  # cm⁻¹
    energy_ev = unitconvert(energy_wn, 'cm^-1->eV')
    energy_k = unitconvert(energy_wn, 'cm^-1->K')
    energy_mhz = unitconvert(energy_wn, 'cm^-1->MHz')
    
    print(f"\nEnergy scale comparisons for {energy_wn} cm⁻¹:")
    print(f"  {energy_ev:.6f} eV")
    print(f"  {energy_k:.1f} K") 
    print(f"  {energy_mhz/1000:.1f} GHz")
    
    # Example 3: Temperature to field conversion
    temp_k = 4.2  # Liquid helium temperature
    temp_wn = unitconvert(temp_k, 'K->cm^-1')
    temp_field = unitconvert(temp_k, 'K->mT')
    
    print(f"\nThermal energy at {temp_k} K:")
    print(f"  {temp_wn:.3f} cm⁻¹")
    print(f"  {temp_field:.3f} mT equivalent field")
    
    # Example 4: Vector conversion
    fields = np.array([100, 200, 300, 400])  # mT
    freqs = unitconvert(fields, 'mT->MHz', g_factor=2.003)
    
    print(f"\nField to frequency conversion (g=2.003):")
    for b, f in zip(fields, freqs):
        print(f"  {b} mT → {f/1000:.3f} GHz")
    
    # Example 5: Different g-factors
    field = 340  # mT
    g_factors = np.array([2.000, 2.003, 2.010, 2.100])
    frequencies = unitconvert(field, 'mT->MHz', g_factor=g_factors)
    
    print(f"\nFrequency at {field} mT for different g-factors:")
    for g, f in zip(g_factors, frequencies):
        print(f"  g = {g:.3f} → {f/1000:.3f} GHz")


if __name__ == "__main__":
    list_conversions()
    print()
    demo_conversions()