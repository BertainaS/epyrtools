#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  6 09:52:16 2025

@author: sylvainbertaina
"""
"""
Module defining common physical constants relevant to physics,
magnetism, and electron spin resonance (ESR).

Values are based on CODATA 2018 recommended values.
Constants are provided in both SI (International System) units and
CGS (Gaussian) units where applicable and meaningful.

Note on CGS: Gaussian units are used. Formulas often change between
SI and CGS (e.g., factors of 4π, c, handling of μ₀/ε₀).
These constants represent the numerical values in the respective units.
"""

import math

# --- Conversion Factors (Internal Use) ---
_J_PER_ERG = 1.0e-7
_ERG_PER_J = 1.0e7
_M_PER_CM = 0.01
_CM_PER_M = 100.0
_KG_PER_G = 0.001
_G_PER_KG = 1000.0
_T_PER_GAUSS = 1.0e-4
_GAUSS_PER_T = 1.0e4
# Magnetic moment: J/T -> erg/G. Energy = μ⋅B. J = (J/T)⋅T; erg = (erg/G)⋅G
# (J/T) / (erg/G) = (J/erg) * (G/T) = 1e7 * 1e-4 = 1e3
_JT_PER_ERGG = 1.0e-3
_ERGG_PER_JT = 1.0e3
# Gyromagnetic ratio: rad/s/T -> rad/s/G
# (rad/s/T) / (rad/s/G) = G/T = 1e-4
_RADST_PER_RADSG = 1.0e-4
_RADSG_PER_RADST = 1.0e4
# Charge: C -> statC (esu) approx c/10 (where c is in cm/s)
_C_cm_s = 299792458.0 * _CM_PER_M # c in cm/s
_STATC_PER_C = _C_cm_s / 10.0
_C_PER_STATC = 1.0 / _STATC_PER_C

# --- Fundamental Constants (CODATA 2018) - SI ---

SPEED_OF_LIGHT_SI = 299792458.0  # Unit: m⋅s⁻¹
"""(c) Speed of light in vacuum (SI). Exact."""

VACUUM_PERMEABILITY_SI = 1.25663706212e-6 # Unit: N⋅A⁻² or H⋅m⁻¹
"""(μ₀) Magnetic constant / Vacuum permeability (SI). CODATA 2018."""
# Note: Before 2019 redefinition, μ₀ was exactly 4π×10⁻⁷ H/m. Now it has uncertainty.

VACUUM_PERMITTIVITY_SI = 8.8541878128e-12 # Unit: F⋅m⁻¹
"""(ε₀) Electric constant / Vacuum permittivity (SI). CODATA 2018. Calculated from μ₀ and c."""
# ε₀ = 1 / (μ₀ * c²)

PLANCK_CONSTANT_SI = 6.62607015e-34  # Unit: J⋅s
"""(h) Planck constant (SI). Exact."""

REDUCED_PLANCK_CONSTANT_SI = PLANCK_CONSTANT_SI / (2 * math.pi) # Unit: J⋅s
"""(ħ) Reduced Planck constant (h-bar = h/2π) (SI)."""

ELEMENTARY_CHARGE_SI = 1.602176634e-19  # Unit: C (Coulomb)
"""(e) Elementary charge (magnitude) (SI). Exact."""

ELECTRON_MASS_SI = 9.1093837015e-31  # Unit: kg
"""(mₑ) Electron mass (SI). CODATA 2018."""

PROTON_MASS_SI = 1.67262192369e-27 # Unit: kg
"""(mₚ) Proton mass (SI). CODATA 2018."""

FINE_STRUCTURE_CONSTANT = 7.2973525693e-3 # Unit: Dimensionless
"""(α) Fine-structure constant (dimensionless). CODATA 2018. Same in SI & CGS."""
# α = e² / (4πε₀ ħ c)  (SI) = e² / (ħ c) (Gaussian CGS)

# --- Fundamental Constants (CODATA 2018) - CGS (Gaussian) ---

SPEED_OF_LIGHT_CGS = SPEED_OF_LIGHT_SI * _CM_PER_M # Unit: cm⋅s⁻¹
"""(c) Speed of light in vacuum (CGS)."""

# μ₀ and ε₀ are often implicitly 1 or handled by factors of 4π/c in CGS formulas.
# We don't define direct CGS equivalents for them as standalone constants.

PLANCK_CONSTANT_CGS = PLANCK_CONSTANT_SI * _ERG_PER_J # Unit: erg⋅s
"""(h) Planck constant (CGS)."""

REDUCED_PLANCK_CONSTANT_CGS = REDUCED_PLANCK_CONSTANT_SI * _ERG_PER_J # Unit: erg⋅s
"""(ħ) Reduced Planck constant (h-bar) (CGS)."""

ELEMENTARY_CHARGE_CGS = ELEMENTARY_CHARGE_SI * _STATC_PER_C # Unit: statC or esu
"""(e) Elementary charge (magnitude) (CGS)."""

ELECTRON_MASS_CGS = ELECTRON_MASS_SI * _G_PER_KG # Unit: g
"""(mₑ) Electron mass (CGS)."""

PROTON_MASS_CGS = PROTON_MASS_SI * _G_PER_KG # Unit: g
"""(mₚ) Proton mass (CGS)."""


# --- Magnetic Moments & g-factors (CODATA 2018) ---

BOHR_MAGNETON_SI = 9.2740100783e-24  # Unit: J⋅T⁻¹
"""(μ<0xE2><0x82><0x8B>) Bohr magneton (SI). μ<0xE2><0x82><0x8B> = eħ / (2mₑ)."""

BOHR_MAGNETON_CGS = BOHR_MAGNETON_SI * _ERGG_PER_JT # Unit: erg⋅G⁻¹
"""(μ<0xE2><0x82><0x8B>) Bohr magneton (CGS)."""

NUCLEAR_MAGNETON_SI = 5.0507837461e-27  # Unit: J⋅T⁻¹
"""(μ<0xE2><0x82><0x99>) Nuclear magneton (SI). μ<0xE2><0x82><0x99> = eħ / (2mₚ)."""

NUCLEAR_MAGNETON_CGS = NUCLEAR_MAGNETON_SI * _ERGG_PER_JT # Unit: erg⋅G⁻¹
"""(μ<0xE2><0x82><0x99>) Nuclear magneton (CGS)."""

ELECTRON_G_FACTOR = -2.00231930436256 # Unit: Dimensionless
"""(gₑ) Electron g-factor (free electron) (Dimensionless). CODATA 2018. Negative sign indicates moment opposes spin angular momentum."""
# Often |gₑ| ≈ 2.0023 is used in calculations.

PROTON_G_FACTOR = 5.5856946893 # Unit: Dimensionless
"""(gₚ) Proton g-factor (Dimensionless). CODATA 2018."""


# --- Gyromagnetic Ratios (Derived) ---
# Gamma (γ) relates magnetic moment (μ) to angular momentum (J): μ = γJ
# For spin S, μ = γħS' where S' is the spin quantum number vector.
# Often refers to γ = μ / (ħ * I) where I is the spin quantum number (e.g., 1/2 for electron)
# More commonly used in NMR/ESR is γ = g * (charge*mass unit) / (2 * mass)
# For electron: γₑ = gₑ * μ<0xE2><0x82><0x8B> / ħ = gₑ * e / (2 * mₑ)
# Sign convention: If g is negative, γ is negative.

ELECTRON_GYROMAGNETIC_RATIO_SI = ELECTRON_G_FACTOR * BOHR_MAGNETON_SI / REDUCED_PLANCK_CONSTANT_SI
# Unit: rad⋅s⁻¹⋅T⁻¹
"""(γₑ) Electron gyromagnetic ratio (SI). γₑ = gₑ μ<0xE2><0x82><0x8B> / ħ. Approx -1.7608 × 10¹¹ rad⋅s⁻¹⋅T⁻¹."""

ELECTRON_GYROMAGNETIC_RATIO_MHZ_T_SI = abs(ELECTRON_GYROMAGNETIC_RATIO_SI) / (2 * math.pi * 1e6)
# Unit: MHz⋅T⁻¹
"""(|γₑ|/2π) Electron gyromagnetic ratio magnitude in MHz/T (SI). Approx 28024.95 MHz/T."""

ELECTRON_GYROMAGNETIC_RATIO_CGS = ELECTRON_GYROMAGNETIC_RATIO_SI * _RADSG_PER_RADST
# Unit: rad⋅s⁻¹⋅G⁻¹
"""(γₑ) Electron gyromagnetic ratio (CGS). Approx -1.7608 × 10⁷ rad⋅s⁻¹⋅G⁻¹."""

PROTON_GYROMAGNETIC_RATIO_SI = PROTON_G_FACTOR * NUCLEAR_MAGNETON_SI / REDUCED_PLANCK_CONSTANT_SI
# Unit: rad⋅s⁻¹⋅T⁻¹
"""(γₚ) Proton gyromagnetic ratio (SI). γₚ = gₚ μ<0xE2><0x82><0x99> / ħ. Approx 2.6752 × 10⁸ rad⋅s⁻¹⋅T⁻¹."""

PROTON_GYROMAGNETIC_RATIO_MHZ_T_SI = abs(PROTON_GYROMAGNETIC_RATIO_SI) / (2 * math.pi * 1e6)
# Unit: MHz⋅T⁻¹
"""(|γₚ|/2π) Proton gyromagnetic ratio magnitude in MHz/T (SI). Approx 42.577 MHz/T."""

PROTON_GYROMAGNETIC_RATIO_CGS = PROTON_GYROMAGNETIC_RATIO_SI * _RADSG_PER_RADST
# Unit: rad⋅s⁻¹⋅G⁻¹
"""(γₚ) Proton gyromagnetic ratio (CGS). Approx 2.6752 × 10⁴ rad⋅s⁻¹⋅G⁻¹."""


# --- Common Aliases (SI) ---
h = PLANCK_CONSTANT_SI
hbar = REDUCED_PLANCK_CONSTANT_SI
c = SPEED_OF_LIGHT_SI
mu0 = VACUUM_PERMEABILITY_SI
eps0 = VACUUM_PERMITTIVITY_SI
e = ELEMENTARY_CHARGE_SI
m_e = ELECTRON_MASS_SI
m_p = PROTON_MASS_SI
alpha = FINE_STRUCTURE_CONSTANT
mu_B = BOHR_MAGNETON_SI
mu_N = NUCLEAR_MAGNETON_SI
g_e = ELECTRON_G_FACTOR
g_p = PROTON_G_FACTOR
gamma_e = ELECTRON_GYROMAGNETIC_RATIO_SI
gamma_e_MHz_T = ELECTRON_GYROMAGNETIC_RATIO_MHZ_T_SI
gamma_p = PROTON_GYROMAGNETIC_RATIO_SI
gamma_p_MHz_T = PROTON_GYROMAGNETIC_RATIO_MHZ_T_SI

# --- Define the public API for 'from physical_constants import *' ---
__all__ = [
    # SI Constants
    'SPEED_OF_LIGHT_SI',
    'VACUUM_PERMEABILITY_SI',
    'VACUUM_PERMITTIVITY_SI',
    'PLANCK_CONSTANT_SI',
    'REDUCED_PLANCK_CONSTANT_SI',
    'ELEMENTARY_CHARGE_SI',
    'ELECTRON_MASS_SI',
    'PROTON_MASS_SI',
    'FINE_STRUCTURE_CONSTANT', # Dimensionless
    'BOHR_MAGNETON_SI',
    'NUCLEAR_MAGNETON_SI',
    'ELECTRON_G_FACTOR', # Dimensionless
    'PROTON_G_FACTOR', # Dimensionless
    'ELECTRON_GYROMAGNETIC_RATIO_SI',
    'ELECTRON_GYROMAGNETIC_RATIO_MHZ_T_SI',
    'PROTON_GYROMAGNETIC_RATIO_SI',
    'PROTON_GYROMAGNETIC_RATIO_MHZ_T_SI',
    # CGS Constants
    'SPEED_OF_LIGHT_CGS',
    'PLANCK_CONSTANT_CGS',
    'REDUCED_PLANCK_CONSTANT_CGS',
    'ELEMENTARY_CHARGE_CGS',
    'ELECTRON_MASS_CGS',
    'PROTON_MASS_CGS',
    'BOHR_MAGNETON_CGS',
    'NUCLEAR_MAGNETON_CGS',
    'ELECTRON_GYROMAGNETIC_RATIO_CGS',
    'PROTON_GYROMAGNETIC_RATIO_CGS',
    # Common Aliases (SI)
    'h', 'hbar', 'c', 'mu0', 'eps0', 'e', 'm_e', 'm_p', 'alpha',
    'mu_B', 'mu_N', 'g_e', 'g_p', 'gamma_e', 'gamma_e_MHz_T',
    'gamma_p', 'gamma_p_MHz_T',
]

# --- Optional: Verification block ---
if __name__ == "__main__":
    print("--- Physical Constants Module (Magnetism/ESR Focus) ---")
    print("Source: CODATA 2018")
    print("-" * 60)
    print("Constant                       | SI Value                     | CGS Value")
    print("-" * 60)
    print(f"Speed of Light (c)           | {SPEED_OF_LIGHT_SI:<25.8e} m/s | {SPEED_OF_LIGHT_CGS:<.8e} cm/s")
    print(f"Planck Constant (h)          | {PLANCK_CONSTANT_SI:<25.8e} J⋅s | {PLANCK_CONSTANT_CGS:<.8e} erg⋅s")
    print(f"Reduced Planck (ħ)         | {REDUCED_PLANCK_CONSTANT_SI:<25.8e} J⋅s | {REDUCED_PLANCK_CONSTANT_CGS:<.8e} erg⋅s")
    print(f"Elementary Charge (e)      | {ELEMENTARY_CHARGE_SI:<25.8e} C   | {ELEMENTARY_CHARGE_CGS:<.8e} statC")
    print(f"Electron Mass (mₑ)         | {ELECTRON_MASS_SI:<25.8e} kg  | {ELECTRON_MASS_CGS:<.8e} g")
    print(f"Proton Mass (mₚ)           | {PROTON_MASS_SI:<25.8e} kg  | {PROTON_MASS_CGS:<.8e} g")
    print(f"Bohr Magneton (μ<0xE2><0x82><0x8B>)        | {BOHR_MAGNETON_SI:<25.8e} J/T | {BOHR_MAGNETON_CGS:<.8e} erg/G")
    print(f"Nuclear Magneton (μ<0xE2><0x82><0x99>)      | {NUCLEAR_MAGNETON_SI:<25.8e} J/T | {NUCLEAR_MAGNETON_CGS:<.8e} erg/G")
    print(f"Electron g-factor (gₑ)       | {ELECTRON_G_FACTOR:<25.15f}   | (Dimensionless)")
    print(f"Proton g-factor (gₚ)         | {PROTON_G_FACTOR:<25.10f}   | (Dimensionless)")
    print(f"Electron Gyromag. (γₑ)     | {ELECTRON_GYROMAGNETIC_RATIO_SI:<25.8e} rad/s/T | {ELECTRON_GYROMAGNETIC_RATIO_CGS:<.8e} rad/s/G")
    print(f"Electron Gyromag. (|γₑ|/2π) | {ELECTRON_GYROMAGNETIC_RATIO_MHZ_T_SI:<25.4f} MHz/T |")
    print(f"Proton Gyromag. (γₚ)       | {PROTON_GYROMAGNETIC_RATIO_SI:<25.8e} rad/s/T | {PROTON_GYROMAGNETIC_RATIO_CGS:<.8e} rad/s/G")
    print(f"Proton Gyromag. (|γₚ|/2π)   | {PROTON_GYROMAGNETIC_RATIO_MHZ_T_SI:<25.4f} MHz/T |")
    print(f"Vacuum Permeability (μ₀)   | {VACUUM_PERMEABILITY_SI:<25.8e} H/m | (Handled differently)")
    print(f"Vacuum Permittivity (ε₀)   | {VACUUM_PERMITTIVITY_SI:<25.8e} F/m | (Handled differently)")
    print(f"Fine Structure (α)         | {FINE_STRUCTURE_CONSTANT:<25.8e}   | (Dimensionless)")
    print("-" * 60)