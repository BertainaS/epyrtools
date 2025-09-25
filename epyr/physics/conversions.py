"""
Direct conversion functions for EPR spectroscopy

Provides simple, direct conversion functions between common EPR units:
MHz, mT, cm-1, and related energy/field conversions.
"""

import numpy as np
from typing import Union, Optional, List, Tuple
from .constants import (
    gfree, bmagn, planck, clight, boltzm, evolt
)


def mhz_to_mt(frequency_mhz: Union[float, np.ndarray],
              g_factor: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
    """
    Convert frequency in MHz to magnetic field in mT.

    Uses the fundamental EPR relation: B = h*nu / (g*mu_B)

    Parameters:
    -----------
    frequency_mhz : float or array
        Frequency in MHz
    g_factor : float or array, optional
        g-factor (defaults to free electron g-factor = 2.002319...)

    Returns:
    --------
    float or array
        Magnetic field in mT

    Examples:
    ---------
    >>> # X-band EPR frequency
    >>> field = mhz_to_mt(9500)  # 9.5 GHz
    >>> print(f"Field: {field:.1f} mT")

    >>> # Different g-factors
    >>> fields = mhz_to_mt(9500, g_factor=[2.000, 2.005, 2.010])
    >>> print(f"Fields: {fields}")
    """
    # Store original input types for return type determination
    freq_is_scalar = np.isscalar(frequency_mhz)
    g_is_scalar = g_factor is None or np.isscalar(g_factor)

    if g_factor is None:
        g_factor = gfree()

    # Convert MHz to Hz, then use fundamental relation
    frequency_hz = np.asarray(frequency_mhz) * 1e6
    g_factor_arr = np.asarray(g_factor)
    field_tesla = frequency_hz * planck() / (g_factor_arr * bmagn())
    field_mt = field_tesla * 1000  # Tesla to mT

    # Return scalar if both inputs were scalar
    if freq_is_scalar and g_is_scalar:
        return float(field_mt)
    return field_mt


def mt_to_mhz(field_mt: Union[float, np.ndarray],
              g_factor: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
    """
    Convert magnetic field in mT to frequency in MHz.

    Uses the fundamental EPR relation: nu = g*mu_B*B / h

    Parameters:
    -----------
    field_mt : float or array
        Magnetic field in mT
    g_factor : float or array, optional
        g-factor (defaults to free electron g-factor = 2.002319...)

    Returns:
    --------
    float or array
        Frequency in MHz

    Examples:
    ---------
    >>> # What frequency for 340 mT field?
    >>> freq = mt_to_mhz(340)
    >>> print(f"Frequency: {freq:.1f} MHz")

    >>> # Array of fields
    >>> freqs = mt_to_mhz([100, 200, 300, 400])
    >>> print(f"Frequencies: {freqs}")
    """
    # Store original input types for return type determination
    field_is_scalar = np.isscalar(field_mt)
    g_is_scalar = g_factor is None or np.isscalar(g_factor)

    if g_factor is None:
        g_factor = gfree()

    # Convert mT to Tesla, then use fundamental relation
    field_tesla = np.asarray(field_mt) * 1e-3
    g_factor_arr = np.asarray(g_factor)
    frequency_hz = g_factor_arr * bmagn() * field_tesla / planck()
    frequency_mhz = frequency_hz / 1e6  # Hz to MHz

    # Return scalar if both inputs were scalar
    if field_is_scalar and g_is_scalar:
        return float(frequency_mhz)
    return frequency_mhz


def cm_inv_to_mhz(wavenumber_cm_inv: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert wavenumber in cm^-1 to frequency in MHz.

    Uses the relation: nu = c * wavenumber
    where c is the speed of light.

    Parameters:
    -----------
    wavenumber_cm_inv : float or array
        Wavenumber in cm^-1

    Returns:
    --------
    float or array
        Frequency in MHz

    Examples:
    ---------
    >>> # Convert 1000 cm^-1 to MHz
    >>> freq = cm_inv_to_mhz(1000)
    >>> print(f"Frequency: {freq:.3e} MHz")

    >>> # Array conversion
    >>> freqs = cm_inv_to_mhz([100, 500, 1000, 2000])
    >>> print(f"Frequencies: {freqs}")
    """
    # wavenumber in cm^-1 * c in m/s * 100 (cm->m) / 1e6 (Hz->MHz)
    wavenumber_cm_inv = np.asarray(wavenumber_cm_inv)
    frequency_mhz = wavenumber_cm_inv * clight() * 100 / 1e6

    # Return scalar if input was scalar
    if np.isscalar(wavenumber_cm_inv):
        return float(frequency_mhz)
    return frequency_mhz


def mhz_to_cm_inv(frequency_mhz: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Convert frequency in MHz to wavenumber in cm^-1.

    Uses the relation: wavenumber = nu / c
    where c is the speed of light.

    Parameters:
    -----------
    frequency_mhz : float or array
        Frequency in MHz

    Returns:
    --------
    float or array
        Wavenumber in cm^-1

    Examples:
    ---------
    >>> # Convert 30000 MHz (30 GHz) to cm^-1
    >>> wn = mhz_to_cm_inv(30000)
    >>> print(f"Wavenumber: {wn:.3f} cm^-1")

    >>> # Array conversion
    >>> wns = mhz_to_cm_inv([1000, 5000, 10000, 30000])
    >>> print(f"Wavenumbers: {wns}")
    """
    # frequency in MHz * 1e6 (MHz->Hz) / (c in m/s * 100 (m->cm))
    frequency_mhz = np.asarray(frequency_mhz)
    wavenumber_cm_inv = frequency_mhz * 1e6 / (clight() * 100)

    # Return scalar if input was scalar
    if np.isscalar(frequency_mhz):
        return float(wavenumber_cm_inv)
    return wavenumber_cm_inv


def frequency_field_conversion_table(frequencies_ghz: Optional[List[float]] = None,
                                   g_factors: Optional[List[float]] = None) -> None:
    """
    Print a conversion table between frequency and magnetic field.

    Parameters:
    -----------
    frequencies_ghz : list of float, optional
        Frequencies in GHz to include in table
        Default: Common EPR frequencies [1, 3, 9.5, 34, 95, 263]
    g_factors : list of float, optional
        g-factors to include in table
        Default: [2.000, 2.002, 2.005, 2.010]
    """
    if frequencies_ghz is None:
        frequencies_ghz = [1, 3, 9.5, 34, 95, 263]  # Common EPR bands

    if g_factors is None:
        g_factors = [2.000, 2.002, 2.005, 2.010]

    print("EPR Frequency to Magnetic Field Conversion Table")
    print("=" * 60)
    print(f"{'Freq (GHz)':<12}", end="")
    for g in g_factors:
        print(f"{'g=' + str(g):<12}", end="")
    print()
    print("-" * 60)

    for freq_ghz in frequencies_ghz:
        freq_mhz = freq_ghz * 1000
        print(f"{freq_ghz:<12.1f}", end="")
        for g in g_factors:
            field_mt = mhz_to_mt(freq_mhz, g_factor=g)
            print(f"{field_mt:<12.1f}", end="")
        print()

    print("-" * 60)
    print("All fields in mT")


def energy_conversion_table(energies_cm_inv: Optional[List[float]] = None) -> None:
    """
    Print a conversion table between different energy units.

    Parameters:
    -----------
    energies_cm_inv : list of float, optional
        Energies in cm^-1 to include in table
        Default: [1, 10, 100, 1000, 5000, 10000]
    """
    if energies_cm_inv is None:
        energies_cm_inv = [1, 10, 100, 1000, 5000, 10000]

    print("Energy Unit Conversion Table")
    print("=" * 70)
    print(f"{'cm^-1':<10} {'MHz':<15} {'GHz':<10} {'eV':<12} {'K':<10} {'meV':<8}")
    print("-" * 70)

    for wn in energies_cm_inv:
        # Convert cm^-1 to other units
        freq_mhz = cm_inv_to_mhz(wn)
        freq_ghz = freq_mhz / 1000

        # Using physical constants for other conversions
        energy_j = wn * 100 * clight() * planck()  # cm^-1 to J
        energy_ev = energy_j / evolt()  # J to eV
        energy_mev = energy_ev * 1000  # eV to meV
        temp_k = energy_j / boltzm()  # J to K

        print(f"{wn:<10.0f} {freq_mhz:<15.3e} {freq_ghz:<10.3f} "
              f"{energy_ev:<12.6f} {temp_k:<10.3f} {energy_mev:<8.3f}")

    print("-" * 70)


def conversion_examples():
    """Show practical examples of EPR unit conversions."""
    print("EPR Unit Conversion Examples")
    print("=" * 40)

    # Example 1: X-band EPR
    print("1. X-band EPR (9.5 GHz)")
    freq_xband = 9500  # MHz
    field_xband = mhz_to_mt(freq_xband)
    print(f"   Frequency: {freq_xband} MHz")
    print(f"   Field: {field_xband:.1f} mT")
    print(f"   Check: {mt_to_mhz(field_xband):.0f} MHz")
    print()

    # Example 2: Q-band EPR
    print("2. Q-band EPR (~34 GHz)")
    freq_qband = 34000  # MHz
    field_qband = mhz_to_mt(freq_qband)
    print(f"   Frequency: {freq_qband} MHz")
    print(f"   Field: {field_qband:.0f} mT")
    print()

    # Example 3: Energy scale conversions
    print("3. Energy scale conversions")
    energy_wn = 1000  # cm^-1
    energy_freq = cm_inv_to_mhz(energy_wn)
    print(f"   {energy_wn} cm^-1 = {energy_freq:.3e} MHz")
    print(f"   Check: {mhz_to_cm_inv(energy_freq):.1f} cm^-1")
    print()

    # Example 4: Different g-factors
    print("4. g-factor effects at X-band")
    freq = 9500  # MHz
    g_values = [2.000, 2.002, 2.005, 2.010]
    fields = mhz_to_mt(freq, g_factor=g_values)

    print(f"   Frequency: {freq} MHz")
    for g, b in zip(g_values, fields):
        print(f"   g = {g:.3f}: {b:.1f} mT")
    print()

    # Example 5: Field sweep range
    print("5. Typical EPR field sweep (g = 2.002)")
    center_field = 340  # mT
    sweep_width = 20   # mT
    fields = np.array([center_field - sweep_width/2, center_field, center_field + sweep_width/2])
    freqs = mt_to_mhz(fields, g_factor=2.002)

    print(f"   Field range: {fields[0]:.1f} - {fields[2]:.1f} mT")
    print(f"   Frequency range: {freqs[0]:.1f} - {freqs[2]:.1f} MHz")
    print(f"   Frequency width: {freqs[2] - freqs[0]:.1f} MHz")


if __name__ == "__main__":
    conversion_examples()
    print()
    frequency_field_conversion_table()
    print()
    energy_conversion_table()