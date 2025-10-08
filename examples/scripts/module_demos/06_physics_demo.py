#!/usr/bin/env python3
"""
EPyR Tools Demo 06: Physics Constants and Unit Conversions
==========================================================

This script demonstrates the physics module which provides physical constants
and unit conversion utilities for EPR/NMR spectroscopy calculations.

Functions demonstrated:
- Physical constants (GFREE, BMAGN, PLANCK, etc.)
- EPR-specific calculations (gamma_hz, magnetic_field_to_frequency)
- Unit conversions (mhz_to_mt, mt_to_mhz, cm_inv_to_mhz)
- Energy and temperature conversions
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_physical_constants():
    """Demonstrate access to physical constants."""
    print("=== EPyR Tools Physics Demo - Physical Constants ===")
    print()

    print("1. Fundamental physical constants (2022 CODATA):")
    print("-" * 50)

    # Display constants using the constants_summary function
    try:
        epyr.physics.constants_summary()
    except:
        # Fallback to manual display
        print("SI Units:")
        print(f"  Free electron g-factor: {epyr.physics.GFREE}")
        print(f"  Bohr magneton: {epyr.physics.BMAGN:.3e} J T⁻¹")
        print(f"  Planck constant: {epyr.physics.PLANCK:.3e} J⋅s")
        print(f"  Reduced Planck (ℏ): {epyr.physics.HBAR:.3e} J⋅s")
        print(f"  Speed of light: {epyr.physics.CLIGHT:.0f} m⋅s⁻¹")
        print(f"  Boltzmann constant: {epyr.physics.BOLTZM:.3e} J⋅K⁻¹")
        print(f"  Avogadro constant: {epyr.physics.AVOGADRO:.3e} mol⁻¹")
        print(f"  Nuclear magneton: {epyr.physics.NMAGN:.3e} J⋅T⁻¹")
        print(f"  Elementary charge: {epyr.physics.ECHARGE:.3e} C")

        print("\nCGS Units:")
        print(f"  Bohr magneton (CGS): {epyr.physics.BMAGN_CGS:.3e} erg⋅G⁻¹")
        print(f"  Planck constant (CGS): {epyr.physics.PLANCK_CGS:.3e} erg⋅s")

    print()

    # Demonstrate ratio calculations
    print("2. Important ratios and derived quantities:")
    print("-" * 42)

    print(f"  γₑ/2π (electron gyromagnetic ratio): {epyr.physics.gamma_hz():.6f} MHz/T")
    print(f"  μ_N/μ_B (nuclear/Bohr magneton ratio): {epyr.physics.NMAGN/epyr.physics.BMAGN:.6f}")
    print(f"  k_BT at 300K: {epyr.physics.thermal_energy(300):.3e} J")
    print(f"  k_BT at 4K: {epyr.physics.thermal_energy(4):.3e} J")

    # Room temperature vs. EPR quantum
    kbt_300k = epyr.physics.thermal_energy(300)
    hf_xband = epyr.physics.PLANCK * 9.5e9  # X-band photon energy
    print(f"  k_BT(300K)/hf(X-band): {kbt_300k/hf_xband:.1f}")

    print()


def demo_frequency_field_conversions():
    """Demonstrate frequency to magnetic field conversions."""
    print("3. Frequency ⟷ Magnetic field conversions:")
    print("-" * 45)

    # Common EPR band frequencies
    bands = {
        'S-band': 3.0,
        'X-band': 9.5,
        'K-band': 24.0,
        'Q-band': 35.0,
        'W-band': 95.0
    }

    print("EPR frequency bands and corresponding magnetic fields:")
    print(f"{'Band':<8} {'Frequency':<12} {'Field (mT)':<12} {'Field (T)':<10}")
    print("-" * 45)

    for band, freq_ghz in bands.items():
        freq_mhz = freq_ghz * 1000
        field_mt = epyr.physics.mhz_to_mt(freq_mhz)
        field_t = field_mt / 1000

        print(f"{band:<8} {freq_ghz:>8.1f} GHz   {field_mt:>8.1f} mT    {field_t:>6.3f} T")

    print()

    # Demonstrate g-factor effects
    print("4. g-factor effects on resonance field:")
    print("-" * 40)

    freq_mhz = 9500  # X-band
    g_factors = [1.990, 2.000, 2.002, 2.010, 2.020, 2.050]

    print(f"At {freq_mhz/1000:.1f} GHz:")
    print(f"{'g-factor':<10} {'Field (mT)':<12} {'Shift (mT)':<12}")
    print("-" * 35)

    field_free_electron = epyr.physics.mhz_to_mt(freq_mhz, g_factor=epyr.physics.GFREE)

    for g in g_factors:
        field = epyr.physics.mhz_to_mt(freq_mhz, g_factor=g)
        shift = field - field_free_electron
        print(f"{g:<10.3f} {field:>8.1f} mT    {shift:>+6.1f} mT")

    print()


def demo_energy_conversions():
    """Demonstrate energy unit conversions."""
    print("5. Energy unit conversions:")
    print("-" * 27)

    # Test energy value in different units
    energy_j = 1e-23  # 10 zJ

    print(f"Energy: {energy_j:.2e} J")
    print("Equivalent values:")

    # Convert to other energy units
    energy_ev = energy_j / epyr.physics.EVOLT
    energy_mhz = energy_j / epyr.physics.PLANCK / 1e6
    energy_cm_inv = epyr.physics.mhz_to_cm_inv(energy_mhz)
    energy_kelvin = energy_j / epyr.physics.BOLTZM

    print(f"  {energy_ev:.3e} eV")
    print(f"  {energy_mhz:.3f} MHz")
    print(f"  {energy_cm_inv:.3f} cm⁻¹")
    print(f"  {energy_kelvin:.3f} K")

    print()

    # Create energy conversion table
    print("6. Energy conversion table:")
    print("-" * 28)

    try:
        epyr.physics.energy_conversion_table()
    except:
        # Manual table
        print("Common energy equivalents:")
        print("1 eV = 241.8 THz = 8065.5 cm⁻¹ = 11604 K")
        print("1 MHz = 4.136×10⁻⁹ eV = 0.033 cm⁻¹ = 4.8×10⁻⁵ K")
        print("1 cm⁻¹ = 1.240×10⁻⁴ eV = 29.98 GHz = 1.439 K")

    print()


def demo_wavelength_frequency():
    """Demonstrate wavelength and frequency conversions."""
    print("7. Wavelength ⟷ Frequency conversions:")
    print("-" * 40)

    # Common electromagnetic spectrum regions
    frequencies = [1e6, 1e9, 1e12, 1e15, 1e18]  # Hz
    regions = ['RF', 'Microwave', 'Infrared', 'Optical', 'X-ray']

    print(f"{'Region':<12} {'Frequency':<15} {'Wavelength':<15}")
    print("-" * 45)

    for region, freq in zip(regions, frequencies):
        wavelength = epyr.physics.CLIGHT / freq
        if wavelength > 1:
            wl_str = f"{wavelength:.1f} m"
        elif wavelength > 1e-3:
            wl_str = f"{wavelength*1000:.1f} mm"
        elif wavelength > 1e-6:
            wl_str = f"{wavelength*1e6:.1f} μm"
        elif wavelength > 1e-9:
            wl_str = f"{wavelength*1e9:.1f} nm"
        else:
            wl_str = f"{wavelength*1e12:.1f} pm"

        if freq >= 1e9:
            freq_str = f"{freq/1e9:.1f} GHz"
        elif freq >= 1e6:
            freq_str = f"{freq/1e6:.1f} MHz"
        else:
            freq_str = f"{freq:.0f} Hz"

        print(f"{region:<12} {freq_str:<15} {wl_str:<15}")

    print()


def demo_temperature_calculations():
    """Demonstrate temperature-related EPR calculations."""
    print("8. Temperature effects in EPR:")
    print("-" * 30)

    temperatures = [4.2, 77, 300, 500]  # K
    freq_mhz = 9500  # X-band

    print("Thermal energy vs. EPR quantum at different temperatures:")
    print(f"{'T (K)':<8} {'k_BT (J)':<15} {'k_BT (MHz)':<12} {'k_BT/hf':<10}")
    print("-" * 50)

    hf = epyr.physics.PLANCK * freq_mhz * 1e6  # EPR quantum energy

    for T in temperatures:
        kbt = epyr.physics.thermal_energy(T)
        kbt_mhz = kbt / epyr.physics.PLANCK / 1e6
        ratio = kbt / hf

        print(f"{T:<8.1f} {kbt:<15.3e} {kbt_mhz:<12.1f} {ratio:<10.3f}")

    print()
    print("Note: When k_BT/hf >> 1, classical (high-T) limit applies")
    print("      When k_BT/hf << 1, quantum (low-T) limit applies")
    print()


def demo_unit_conversion_system():
    """Demonstrate the general unit conversion system."""
    print("9. General unit conversion system:")
    print("-" * 34)

    try:
        # Show available conversions
        print("Available unit conversions:")
        epyr.physics.list_conversions()

        print("\nDemonstration of common conversions:")
        epyr.physics.demo_conversions()

    except Exception as e:
        print(f"Unit conversion system not available: {e}")
        print("Showing manual conversions:")

        # Manual conversions
        values_and_units = [
            (9500, 'MHz', 'mT'),
            (339, 'mT', 'MHz'),
            (1000, 'MHz', 'cm⁻¹'),
            (33.36, 'cm⁻¹', 'MHz')
        ]

        for value, from_unit, to_unit in values_and_units:
            if from_unit == 'MHz' and to_unit == 'mT':
                result = epyr.physics.mhz_to_mt(value)
                print(f"  {value} {from_unit} = {result:.2f} {to_unit}")
            elif from_unit == 'mT' and to_unit == 'MHz':
                result = epyr.physics.mt_to_mhz(value)
                print(f"  {value} {from_unit} = {result:.1f} {to_unit}")
            elif from_unit == 'MHz' and to_unit == 'cm⁻¹':
                result = epyr.physics.mhz_to_cm_inv(value)
                print(f"  {value} {from_unit} = {result:.3f} {to_unit}")
            elif from_unit == 'cm⁻¹' and to_unit == 'MHz':
                result = epyr.physics.cm_inv_to_mhz(value)
                print(f"  {value} {from_unit} = {result:.1f} {to_unit}")

    print()


def demo_field_sweep_calculations():
    """Demonstrate calculations for field sweep experiments."""
    print("10. EPR field sweep calculations:")
    print("-" * 33)

    # Typical X-band EPR parameters
    center_freq = 9.5  # GHz
    center_field = epyr.physics.mhz_to_mt(center_freq * 1000)
    sweep_width = 100  # mT

    print(f"X-band EPR field sweep:")
    print(f"  Center frequency: {center_freq} GHz")
    print(f"  Center field: {center_field:.1f} mT")
    print(f"  Sweep width: {sweep_width} mT")

    # Calculate field range
    field_start = center_field - sweep_width/2
    field_end = center_field + sweep_width/2

    print(f"  Field range: {field_start:.1f} to {field_end:.1f} mT")

    # Calculate corresponding frequency range (for off-resonance effects)
    freq_start = epyr.physics.mt_to_mhz(field_start)
    freq_end = epyr.physics.mt_to_mhz(field_end)

    print(f"  Frequency range: {freq_start:.1f} to {freq_end:.1f} MHz")
    print(f"  Frequency span: {freq_end - freq_start:.1f} MHz")

    # Calculate field resolution for typical digitization
    n_points = 1024
    field_resolution = sweep_width / n_points
    freq_resolution = (freq_end - freq_start) / n_points

    print(f"  Field resolution ({n_points} points): {field_resolution:.3f} mT")
    print(f"  Frequency resolution: {freq_resolution:.3f} MHz")

    print()


def create_conversion_plots():
    """Create plots showing unit conversions."""
    print("11. Creating conversion plots:")
    print("-" * 30)

    # Plot 1: Frequency vs. Field for different g-factors
    frequencies = np.linspace(1000, 100000, 100)  # 1-100 GHz
    g_factors = [1.99, 2.00, 2.002, 2.01, 2.02]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for g in g_factors:
        fields = epyr.physics.mhz_to_mt(frequencies, g_factor=g)
        ax1.plot(frequencies/1000, fields/1000, label=f'g = {g}', linewidth=2)

    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('Magnetic Field (T)')
    ax1.set_title('EPR Frequency vs. Magnetic Field')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Energy conversion chart
    energies_mhz = np.logspace(0, 6, 100)  # 1 MHz to 1 THz
    energies_cm_inv = [epyr.physics.mhz_to_cm_inv(f) for f in energies_mhz]
    energies_ev = energies_mhz * epyr.physics.PLANCK * 1e6 / epyr.physics.EVOLT

    ax2.loglog(energies_mhz/1000, energies_ev, 'b-', label='eV', linewidth=2)
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Energy (eV)')
    ax2.set_title('Frequency to Energy Conversion')
    ax2.grid(True, alpha=0.3)

    # Add some reference lines
    ax2.axvline(9.5, color='red', linestyle='--', alpha=0.7, label='X-band')
    ax2.axvline(35, color='orange', linestyle='--', alpha=0.7, label='Q-band')
    ax2.axvline(95, color='purple', linestyle='--', alpha=0.7, label='W-band')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "06_unit_conversions.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 06_unit_conversions.png")

    # Plot 3: Temperature effects
    temperatures = np.linspace(4, 300, 100)
    thermal_energies = [epyr.physics.thermal_energy(T) for T in temperatures]
    thermal_mhz = np.array(thermal_energies) / epyr.physics.PLANCK / 1e6

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(temperatures, thermal_mhz, 'g-', linewidth=2)
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Thermal Energy (MHz)')
    ax1.set_title('Thermal Energy vs. Temperature')
    ax1.grid(True, alpha=0.3)

    # Add reference lines for common EPR frequencies
    ax1.axhline(9500, color='red', linestyle='--', alpha=0.7, label='X-band')
    ax1.axhline(35000, color='orange', linestyle='--', alpha=0.7, label='Q-band')
    ax1.legend()

    # Plot thermal energy ratio
    hf_xband = 9500  # MHz
    ratio = thermal_mhz / hf_xband

    ax2.semilogy(temperatures, ratio, 'r-', linewidth=2)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('k_BT / hf (X-band)')
    ax2.set_title('Thermal vs. EPR Quantum Energy')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(1, color='black', linestyle='--', alpha=0.7, label='Classical limit')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "06_temperature_effects.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 06_temperature_effects.png")

    print()


def main():
    """Run all physics demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_physical_constants()
    demo_frequency_field_conversions()
    demo_energy_conversions()
    demo_wavelength_frequency()
    demo_temperature_calculations()
    demo_unit_conversion_system()
    demo_field_sweep_calculations()
    create_conversion_plots()

    print("=== Physics Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- Physical constants are based on 2022 CODATA recommendations")
    print("- Available in both SI and CGS unit systems")
    print("- Direct conversion functions for common EPR units (MHz, mT, cm⁻¹)")
    print("- Temperature effects are crucial for EPR saturation and population")
    print("- g-factor variations shift resonance fields significantly")
    print("- Energy scales span many orders of magnitude in spectroscopy")
    print()
    print("Generated plot files:")
    for plot_file in sorted(output_dir.glob("06_*.png")):
        print(f"  - {plot_file.name}")


if __name__ == "__main__":
    main()