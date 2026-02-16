"""
Physical constants for EPR/NMR spectroscopy
All values from 2022 CODATA recommendations with proper units and uncertainties
Constants available in both SI and CGS units.
"""

from typing import Optional, Tuple, Union

import numpy as np

from ..logging_config import get_logger

logger = get_logger(__name__)


# ============================================================================
# FUNDAMENTAL PHYSICAL CONSTANTS (2022 CODATA) - SI Units
# ============================================================================

# Free electron g-factor (dimensionless)
# The g-factor relates the magnetic moment to the angular momentum
# for a free electron: μ = -g μ_B S/ℏ
# 2022 CODATA recommended values
GFREE = 2.00231930436092

# Bohr magneton in SI units (J T⁻¹)
# The fundamental quantum of magnetic moment for an electron.
# μ_B = eℏ/(2m_e)
# 2022 CODATA recommended values
BMAGN = 9.2740100657e-24  # J T⁻¹

# Planck constant in SI units (J⋅s = J⋅Hz⁻¹)
# Fundamental constant relating energy to frequency: E = hν
# Exact value since 2019 SI redefinition.
PLANCK = 6.62607015e-34  # J⋅s

# Reduced Planck constant ℏ = h/(2π) in SI units (J⋅s)
HBAR = PLANCK / (2 * np.pi)

# Speed of light in vacuum in SI units (m⋅s⁻¹)
# Exact value by definition since 1983.
CLIGHT = 299792458  # m⋅s⁻¹

# Boltzmann constant in SI units (J⋅K⁻¹)
# Relates temperature to energy: E = k_B T
# Exact value since 2019 SI redefinition.
BOLTZM = 1.380649e-23  # J⋅K⁻¹ (exact)

# Avogadro constant in SI units (mol⁻¹)
# Number of constituent particles per mole.
# Exact value since 2019 SI redefinition.
AVOGADRO = 6.02214076e23  # mol⁻¹ (exact)

# Nuclear magneton in SI units (J⋅T⁻¹)
# Fundamental quantum of magnetic moment for nucleons.
# μ_N = eℏ/(2m_p) where m_p is the proton mass.
# 2022 CODATA recommended values
NMAGN = 5.0507837393e-27  # J⋅T⁻¹

# Elementary charge in SI units (C)
# Exact value since 2019 SI redefinition
ECHARGE = 1.602176634e-19  # C (exact)

# Electron volt in SI units (J)
# Energy equivalent to one electron volt: 1 eV = e × 1 V
# Exact since charge is exact
EVOLT = ECHARGE  # J


# ============================================================================
# FUNDAMENTAL PHYSICAL CONSTANTS - CGS Units
# ============================================================================

# Free electron g-factor (dimensionless, same in all unit systems)
GFREE_CGS = GFREE

# Bohr magneton in CGS units (erg G⁻¹)
# Conversion: 1 J T⁻¹ = 1000 erg G⁻¹
BMAGN_CGS = BMAGN * 1000  # erg G⁻¹

# Planck constant in CGS units (erg⋅s)
# Conversion: 1 J = 1e7 erg
PLANCK_CGS = PLANCK * 1e7  # erg⋅s

# Reduced Planck constant ℏ = h/(2π) in CGS units (erg⋅s)
HBAR_CGS = PLANCK_CGS / (2 * np.pi)

# Speed of light in vacuum in CGS units (cm⋅s⁻¹)
# Conversion: 1 m = 100 cm
CLIGHT_CGS = CLIGHT * 100  # cm⋅s⁻¹

# Boltzmann constant in CGS units (erg⋅K⁻¹)
# Conversion: 1 J = 1e7 erg
BOLTZM_CGS = BOLTZM * 1e7  # erg⋅K⁻¹

# Avogadro constant in CGS units (mol⁻¹, same in all unit systems)
AVOGADRO_CGS = AVOGADRO

# Nuclear magneton in CGS units (erg⋅G⁻¹)
# Conversion: 1 J T⁻¹ = 1000 erg G⁻¹
NMAGN_CGS = NMAGN * 1000  # erg⋅G⁻¹

# Elementary charge in CGS units (esu, same as statcoulomb)
# Conversion: 1 C = 2.997924580e9 esu (from c in units)
ECHARGE_CGS = ECHARGE * 2.997924580e9  # esu

# Electron volt in CGS units (erg)
# Conversion: 1 J = 1e7 erg
EVOLT_CGS = EVOLT * 1e7  # erg


# ============================================================================
# BACKWARD COMPATIBILITY FUNCTIONS
# ============================================================================


def gfree(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Free electron g-factor (dimensionless).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Free electron g-factor (dimensionless)

    References:
    -----------
    2022 CODATA recommended values
    """
    uncertainty = 0.00000000000036
    if return_uncertainty:
        return GFREE, uncertainty
    return GFREE


def bmagn(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Bohr magneton in SI units (J T⁻¹).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Bohr magneton in J T⁻¹

    References:
    -----------
    2022 CODATA recommended values
    """
    uncertainty = 0.0000000029e-24
    if return_uncertainty:
        return BMAGN, uncertainty
    return BMAGN


def planck() -> float:
    """
    Planck constant in SI units (J⋅s = J⋅Hz⁻¹).

    Returns:
    --------
    float
        Planck constant in J⋅s

    References:
    -----------
    2019 SI redefinition, exact value
    """
    return PLANCK


def hbar() -> float:
    """
    Reduced Planck constant ℏ = h/(2π) in SI units (J⋅s).

    Returns:
    --------
    float
        Reduced Planck constant in J⋅s
    """
    return HBAR


def clight() -> float:
    """
    Speed of light in vacuum in SI units (m⋅s⁻¹).

    Returns:
    --------
    float
        Speed of light in m⋅s⁻¹

    References:
    -----------
    1983 SI redefinition, exact value
    """
    return CLIGHT


def boltzm(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Boltzmann constant in SI units (J⋅K⁻¹).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Boltzmann constant in J⋅K⁻¹

    References:
    -----------
    2019 SI redefinition, exact value
    """
    if return_uncertainty:
        return BOLTZM, 0.0  # Exact since 2019
    return BOLTZM


def avogadro(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Avogadro constant in SI units (mol⁻¹).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Avogadro constant in mol⁻¹

    References:
    -----------
    2019 SI redefinition, exact value
    """
    if return_uncertainty:
        return AVOGADRO, 0.0  # Exact since 2019
    return AVOGADRO


def nmagn(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Nuclear magneton in SI units (J⋅T⁻¹).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Nuclear magneton in J⋅T⁻¹

    References:
    -----------
    2022 CODATA recommended values
    """
    uncertainty = 0.0000000016e-27
    if return_uncertainty:
        return NMAGN, uncertainty
    return NMAGN


def echarge(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Elementary charge in SI units (C).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Elementary charge in C

    References:
    -----------
    2019 SI redefinition, exact value
    """
    if return_uncertainty:
        return ECHARGE, 0.0  # Exact since 2019
    return ECHARGE


def evolt(return_uncertainty: bool = False) -> Union[float, Tuple[float, float]]:
    """
    Electron volt in SI units (J).

    Parameters:
    -----------
    return_uncertainty : bool
        If True, return (value, standard_uncertainty)

    Returns:
    --------
    float or tuple
        Electron volt in J

    References:
    -----------
    2019 SI redefinition, exact value (e × 1V)
    """
    if return_uncertainty:
        return EVOLT, 0.0  # Exact since 2019
    return EVOLT


# ============================================================================
# EPR-SPECIFIC FUNCTIONS
# ============================================================================


def gamma_hz(g_factor: Optional[float] = None) -> float:
    """
    Calculate gyromagnetic ratio in Hz/T for any g-factor.

    The gyromagnetic ratio relates frequency to magnetic field: ν = γ B
    where γ = g μ_B / h

    Parameters:
    -----------
    g_factor : float, optional
        g-factor (defaults to free electron g-factor)

    Returns:
    --------
    float
        Gyromagnetic ratio in Hz/T

    Examples:
    ---------
    >>> # Free electron gyromagnetic ratio
    >>> gamma_e = gamma_hz()
    >>> print(f"Free electron: {gamma_e:.3e} Hz/T")

    >>> # Custom g-factor
    >>> gamma_custom = gamma_hz(2.005)
    >>> print(f"g=2.005: {gamma_custom:.3e} Hz/T")

    >>> # Calculate resonance frequency
    >>> B = 0.34  # Tesla (X-band field)
    >>> freq = gamma_hz() * B
    >>> print(f"X-band frequency: {freq/1e9:.2f} GHz")
    """
    if g_factor is None:
        g_factor = GFREE

    return g_factor * BMAGN / PLANCK


def magnetic_field_to_frequency(
    B_tesla: float, g_factor: Optional[float] = None
) -> float:
    """
    Convert magnetic field to resonance frequency.

    Uses the fundamental EPR/NMR relation: ν = γB = gμ_B B/h

    Parameters:
    -----------
    B_tesla : float
        Magnetic field in Tesla
    g_factor : float, optional
        g-factor (defaults to free electron g-factor)

    Returns:
    --------
    float
        Resonance frequency in Hz

    Examples:
    ---------
    >>> # X-band EPR at ~9.5 GHz
    >>> B = 0.34  # Tesla
    >>> freq = magnetic_field_to_frequency(B)  # Hz
    >>> print(f"Frequency: {freq/1e9:.2f} GHz")
    """
    return gamma_hz(g_factor) * B_tesla


def frequency_to_magnetic_field(
    freq_hz: float, g_factor: Optional[float] = None
) -> float:
    """
    Convert frequency to magnetic field.

    Parameters:
    -----------
    freq_hz : float
        Frequency in Hz
    g_factor : float, optional
        g-factor (defaults to free electron g-factor)

    Returns:
    --------
    float
        Magnetic field in Tesla

    Examples:
    ---------
    >>> # What field for 9.5 GHz EPR?
    >>> freq = 9.5e9  # Hz
    >>> B = frequency_to_magnetic_field(freq)
    >>> print(f"Magnetic field: {B*1000:.1f} mT")
    """
    return freq_hz / gamma_hz(g_factor)


def thermal_energy(temperature_k: float) -> float:
    """
    Thermal energy k_B T at given temperature.

    Parameters:
    -----------
    temperature_k : float
        Temperature in Kelvin

    Returns:
    --------
    float
        Thermal energy in J

    Examples:
    ---------
    >>> # Room temperature thermal energy
    >>> E_thermal = thermal_energy(295)  # K
    >>> print(f"kT = {E_thermal/(1.602176634e-19):.3f} meV")
    """
    return BOLTZM * temperature_k


def wavelength_to_frequency(wavelength_m: float) -> float:
    """
    Convert wavelength to frequency.

    Parameters:
    -----------
    wavelength_m : float
        Wavelength in meters

    Returns:
    --------
    float
        Frequency in Hz

    Examples:
    ---------
    >>> # 3 cm microwave wavelength
    >>> freq = wavelength_to_frequency(0.03)  # m
    >>> print(f"Frequency: {freq/1e9:.1f} GHz")
    """
    return CLIGHT / wavelength_m


# ============================================================================
# PHYSICAL CONSTANTS SUMMARY
# ============================================================================


def constants_summary():
    """Print summary of all physical constants with units and values."""

    logger.info("EPyR Tools Physical Constants")
    logger.info("=" * 50)
    logger.info("All values from 2022 CODATA recommendations")
    logger.info("")

    logger.info("SI Units:")
    si_constants = [
        ("Free electron g-factor", GFREE, "dimensionless"),
        ("Bohr magneton", BMAGN, "J⋅T⁻¹"),
        ("Nuclear magneton", NMAGN, "J⋅T⁻¹"),
        ("Planck constant", PLANCK, "J⋅s"),
        ("Reduced Planck constant", HBAR, "J⋅s"),
        ("Speed of light", CLIGHT, "m⋅s⁻¹"),
        ("Boltzmann constant", BOLTZM, "J⋅K⁻¹"),
        ("Avogadro constant", AVOGADRO, "mol⁻¹"),
        ("Elementary charge", ECHARGE, "C"),
        ("Electron volt", EVOLT, "J"),
    ]

    for name, value, unit in si_constants:
        logger.info(f"  {name:<30}: {value:.6e} {unit}")

    logger.info("")
    logger.info("CGS Units:")
    cgs_constants = [
        ("Free electron g-factor", GFREE_CGS, "dimensionless"),
        ("Bohr magneton", BMAGN_CGS, "erg⋅G⁻¹"),
        ("Nuclear magneton", NMAGN_CGS, "erg⋅G⁻¹"),
        ("Planck constant", PLANCK_CGS, "erg⋅s"),
        ("Reduced Planck constant", HBAR_CGS, "erg⋅s"),
        ("Speed of light", CLIGHT_CGS, "cm⋅s⁻¹"),
        ("Boltzmann constant", BOLTZM_CGS, "erg⋅K⁻¹"),
        ("Avogadro constant", AVOGADRO_CGS, "mol⁻¹"),
        ("Elementary charge", ECHARGE_CGS, "esu"),
        ("Electron volt", EVOLT_CGS, "erg"),
    ]

    for name, value, unit in cgs_constants:
        logger.info(f"  {name:<30}: {value:.6e} {unit}")

    logger.info("")
    logger.info("EPR Examples:")
    # Free electron gyromagnetic ratio
    gamma_free = gamma_hz()
    logger.info(f"  Free electron γ/2π: {gamma_free:.3e} Hz/T")

    # X-band EPR
    freq_xband = 9.5e9  # Hz
    B_xband = frequency_to_magnetic_field(freq_xband)
    logger.info(f"  X-band EPR (9.5 GHz): {B_xband*1000:.1f} mT")

    # Room temperature thermal energy
    E_thermal = thermal_energy(295)  # K
    logger.info(f"  Room temperature kT: {E_thermal/(1.602176634e-19)*1000:.2f} meV")


if __name__ == "__main__":
    constants_summary()
