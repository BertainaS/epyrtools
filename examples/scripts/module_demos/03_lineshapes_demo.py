#!/usr/bin/env python3
"""
EPyR Tools Demo 03: Lineshape Functions and Fitting
===================================================

This script demonstrates the comprehensive lineshape capabilities of EPyR Tools.
The lineshapes module provides mathematical functions for EPR spectral analysis.

Functions demonstrated:
- gaussian() - Gaussian absorption and derivative lineshapes
- lorentzian() - Lorentzian absorption and dispersion lineshapes
- voigtian() - Voigt profile (convolution of Gaussian and Lorentzian)
- pseudo_voigt() - Pseudo-Voigt approximation
- fit_epr_signal() - EPR signal fitting
- Lineshape class - Unified interface
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr
from epyr.lineshapes import gaussian, lorentzian, voigtian, pseudo_voigt, Lineshape
from epyr.lineshapes import fit_epr_signal, fit_multiple_shapes


def demo_basic_lineshapes():
    """Demonstrate basic lineshape functions."""
    print("=== EPyR Tools Lineshapes Demo ===")
    print()

    # Create field axis
    B = np.linspace(3350, 3450, 500)
    B0 = 3400  # Center field
    width = 10  # Linewidth

    print("1. Basic lineshape functions:")
    print("-" * 30)

    # Gaussian lineshape
    y_gauss = gaussian(B, center=B0, width=width)
    print(f"Gaussian: center={B0} G, width={width} G")
    print(f"  Peak value: {np.max(y_gauss):.4f}")
    print(f"  HWHM: {width/np.sqrt(2*np.log(2)):.2f} G")

    # Lorentzian lineshape
    y_lorentz = lorentzian(B, center=B0, width=width)
    print(f"Lorentzian: center={B0} G, width={width} G")
    print(f"  Peak value: {np.max(y_lorentz):.4f}")
    print(f"  HWHM: {width/2:.2f} G")

    # Voigtian lineshape
    # voigtian uses widths tuple (gaussian_fwhm, lorentzian_fwhm)
    gaussian_fwhm = 8  # Gaussian component FWHM
    lorentzian_fwhm = 12  # Lorentzian component FWHM
    y_voigt = voigtian(B, center=B0, widths=(gaussian_fwhm, lorentzian_fwhm))
    print(f"Voigtian: center={B0} G, gaussian_fwhm={gaussian_fwhm} G, lorentzian_fwhm={lorentzian_fwhm} G")
    print(f"  Peak value: {np.max(y_voigt):.4f}")

    # Plot comparison
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.plot(B, y_gauss, 'b-', label=f'Gaussian (width={width})', linewidth=2)
    plt.plot(B, y_lorentz, 'r-', label=f'Lorentzian (width={width})', linewidth=2)
    plt.plot(B, y_voigt, 'g-', label=f'Voigt (G={gaussian_fwhm}, L={lorentzian_fwhm})', linewidth=2)
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Normalized Intensity')
    plt.title('EPR Lineshape Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)

    print()


def demo_derivatives():
    """Demonstrate derivative lineshapes."""
    print("2. Derivative lineshapes:")
    print("-" * 26)

    B = np.linspace(3350, 3450, 500)
    B0 = 3400
    width = 10

    # Gaussian derivatives
    y_gauss_0 = gaussian(B, center=B0, width=width, derivative=0)
    y_gauss_1 = gaussian(B, center=B0, width=width, derivative=1)
    y_gauss_2 = gaussian(B, center=B0, width=width, derivative=2)

    print("Gaussian derivatives:")
    print(f"  0th derivative (absorption): peak = {np.max(np.abs(y_gauss_0)):.4f}")
    print(f"  1st derivative: peak-to-peak = {np.max(y_gauss_1) - np.min(y_gauss_1):.4f}")
    print(f"  2nd derivative: peak = {np.max(np.abs(y_gauss_2)):.4f}")

    # Lorentzian derivatives
    y_lorentz_0 = lorentzian(B, center=B0, width=width, derivative=0)
    y_lorentz_1 = lorentzian(B, center=B0, width=width, derivative=1)
    y_lorentz_2 = lorentzian(B, center=B0, width=width, derivative=2)

    print("Lorentzian derivatives:")
    print(f"  0th derivative (absorption): peak = {np.max(np.abs(y_lorentz_0)):.4f}")
    print(f"  1st derivative: peak-to-peak = {np.max(y_lorentz_1) - np.min(y_lorentz_1):.4f}")
    print(f"  2nd derivative: peak = {np.max(np.abs(y_lorentz_2)):.4f}")

    # Plot derivatives
    plt.subplot(2, 1, 2)
    plt.plot(B, y_gauss_0, 'b-', label='Gaussian (0th)', linewidth=2)
    plt.plot(B, y_gauss_1, 'b--', label='Gaussian (1st)', linewidth=2)
    plt.plot(B, y_lorentz_0, 'r-', label='Lorentzian (0th)', linewidth=2)
    plt.plot(B, y_lorentz_1, 'r--', label='Lorentzian (1st)', linewidth=2)
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Normalized Intensity')
    plt.title('EPR Derivative Lineshapes')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "07_lineshapes_basic.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 07_lineshapes_basic.png")
    print()


def demo_phase_effects():
    """Demonstrate phase mixing in lineshapes."""
    print("3. Phase effects:")
    print("-" * 16)

    B = np.linspace(3350, 3450, 500)
    B0 = 3400
    width = 12

    phases = [0, 30, 60, 90]  # degrees

    plt.figure(figsize=(12, 8))

    for i, phase in enumerate(phases):
        # Lorentzian with phase (creates absorption/dispersion mix)
        y_phase = lorentzian(B, center=B0, width=width, phase=np.radians(phase))

        plt.subplot(2, 2, i+1)
        plt.plot(B, y_phase, 'b-', linewidth=2)
        plt.title(f'Phase = {phase}°')
        plt.xlabel('Magnetic Field (G)')
        plt.ylabel('Intensity')
        plt.grid(True, alpha=0.3)

        # Calculate absorption/dispersion components
        absorption_frac = np.cos(np.radians(phase))**2
        dispersion_frac = np.sin(np.radians(phase))**2
        print(f"  Phase {phase}°: {absorption_frac:.1%} absorption, {dispersion_frac:.1%} dispersion")

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "08_phase_effects.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 08_phase_effects.png")
    print()


def demo_voigt_profiles():
    """Demonstrate Voigt profile with different ratios."""
    print("4. Voigt profiles:")
    print("-" * 17)

    B = np.linspace(3350, 3450, 500)
    B0 = 3400

    # Different Gaussian/Lorentzian ratios
    configurations = [
        (10, 2, "Mostly Gaussian"),
        (5, 5, "Equal mix"),
        (2, 10, "Mostly Lorentzian"),
        (0.1, 15, "Nearly pure Lorentzian")
    ]

    plt.figure(figsize=(12, 10))

    for i, (sigma, gamma, description) in enumerate(configurations):
        y_voigt = voigtian(B, center=B0, widths=(sigma*2.35, gamma*2))

        plt.subplot(2, 2, i+1)
        plt.plot(B, y_voigt, 'purple', linewidth=2)
        plt.title(f'{description}\nσ={sigma} G, γ={gamma} G')
        plt.xlabel('Magnetic Field (G)')
        plt.ylabel('Normalized Intensity')
        plt.grid(True, alpha=0.3)

        # Calculate FWHM approximation
        fwhm_approx = 0.5346 * gamma + np.sqrt(0.2166 * gamma**2 + sigma**2 * 2 * np.log(2))
        print(f"  {description}: σ={sigma}, γ={gamma}, FWHM≈{fwhm_approx:.1f} G")

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "09_voigt_profiles.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 09_voigt_profiles.png")
    print()


def demo_lineshape_class():
    """Demonstrate the unified Lineshape class."""
    print("5. Unified Lineshape class:")
    print("-" * 28)

    B = np.linspace(3350, 3450, 500)

    # Create different lineshape objects
    gauss_shape = Lineshape('gaussian', width=10)
    lorentz_shape = Lineshape('lorentzian', width=10)
    voigt_shape = Lineshape('voigt', width=(16.5, 14))

    print("Lineshape objects created:")
    print(f"  Gaussian: {gauss_shape}")
    print(f"  Lorentzian: {lorentz_shape}")
    print(f"  Voigtian: {voigt_shape}")

    # Evaluate lineshapes
    center = 3400
    y_gauss = gauss_shape(B, center=center)
    y_lorentz = lorentz_shape(B, center=center)
    y_voigt = voigt_shape(B, center=center)

    # Plot unified interface results
    plt.figure(figsize=(10, 6))
    plt.plot(B, y_gauss, 'b-', label='Gaussian (class)', linewidth=2)
    plt.plot(B, y_lorentz, 'r-', label='Lorentzian (class)', linewidth=2)
    plt.plot(B, y_voigt, 'g-', label='Voigt (class)', linewidth=2)
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Normalized Intensity')
    plt.title('Lineshape Class Interface')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(Path(__file__).parent / "10_lineshape_class.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 10_lineshape_class.png")
    print()


def demo_signal_fitting():
    """Demonstrate EPR signal fitting."""
    print("6. EPR signal fitting:")
    print("-" * 21)

    # Create synthetic EPR data with noise
    B = np.linspace(3350, 3450, 300)
    B0_true = 3398
    width_true = 12
    amplitude_true = 0.8

    # Generate noisy Lorentzian signal
    y_true = amplitude_true * lorentzian(B, center=B0_true, width=width_true)
    noise_level = 0.05
    y_noisy = y_true + noise_level * np.random.normal(size=len(B))

    print(f"Synthetic EPR signal:")
    print(f"  True center: {B0_true} G")
    print(f"  True width: {width_true} G")
    print(f"  True amplitude: {amplitude_true}")
    print(f"  Noise level: {noise_level}")
    print()

    # Fit with single lineshape
    print("Fitting with single Lorentzian:")
    try:
        result = fit_epr_signal(B, y_noisy, 'lorentzian')

        print(f"  Fitted center: {result.params['center']:.2f} ± {result.errors['center']:.2f} G")
        print(f"  Fitted width: {result.params['width']:.2f} ± {result.errors['width']:.2f} G")
        print(f"  Fitted amplitude: {result.params['amplitude']:.3f} ± {result.errors['amplitude']:.3f}")
        print(f"  R-squared: {result.r_squared:.4f}")
        print(f"  Reduced chi-squared: {result.reduced_chi_squared:.4f}")

        # Plot fit result
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.plot(B, y_noisy, 'ko', markersize=3, alpha=0.6, label='Noisy data')
        plt.plot(B, y_true, 'g-', linewidth=2, label='True signal')
        plt.plot(B, result.fitted_data, 'r-', linewidth=2, label='Fitted signal')
        plt.xlabel('Magnetic Field (G)')
        plt.ylabel('Intensity')
        plt.title(f'EPR Signal Fitting (R² = {result.r_squared:.4f})')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Plot residuals
        plt.subplot(2, 1, 2)
        residuals = y_noisy - result.fitted_data
        plt.plot(B, residuals, 'b-', linewidth=1, alpha=0.7)
        plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        plt.xlabel('Magnetic Field (G)')
        plt.ylabel('Residuals')
        plt.title('Fit Residuals')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(Path(__file__).parent / "11_signal_fitting.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 11_signal_fitting.png")

    except Exception as e:
        print(f"  Fitting failed: {e}")

    print()


def demo_multiple_fitting():
    """Demonstrate fitting with multiple models."""
    print("7. Multiple model comparison:")
    print("-" * 30)

    # Create synthetic data that's better fit by Voigt
    B = np.linspace(3350, 3450, 400)
    B0 = 3400

    # Generate Voigt data
    y_true = voigtian(B, center=B0, widths=(18.8, 12))
    y_noisy = y_true + 0.04 * np.random.normal(size=len(B))

    print("Comparing multiple lineshape models:")

    try:
        results = fit_multiple_shapes(B, y_noisy)

        print("Fit comparison results:")
        for shape_name, result in results.items():
            print(f"  {shape_name.capitalize()}:")
            print(f"    R²: {result.r_squared:.4f}")
            print(f"    Reduced χ²: {result.reduced_chi_squared:.4f}")
            print(f"    AIC: {result.aic:.2f}")

        # Find best fit
        best_shape = min(results.keys(), key=lambda k: results[k].aic)
        print(f"\nBest fit: {best_shape.capitalize()} (lowest AIC)")

        # Plot comparison
        plt.figure(figsize=(12, 10))

        plt.subplot(2, 2, 1)
        plt.plot(B, y_noisy, 'ko', markersize=2, alpha=0.5, label='Data')
        plt.plot(B, y_true, 'k-', linewidth=2, label='True (Voigt)')
        plt.xlabel('Magnetic Field (G)')
        plt.ylabel('Intensity')
        plt.title('Original Data')
        plt.legend()
        plt.grid(True, alpha=0.3)

        colors = ['red', 'blue', 'green']
        for i, (shape_name, result) in enumerate(results.items()):
            plt.subplot(2, 2, i+2)
            plt.plot(B, y_noisy, 'ko', markersize=2, alpha=0.5, label='Data')
            plt.plot(B, result.fitted_data, color=colors[i], linewidth=2,
                    label=f'{shape_name.capitalize()} fit')
            plt.xlabel('Magnetic Field (G)')
            plt.ylabel('Intensity')
            plt.title(f'{shape_name.capitalize()} (R² = {result.r_squared:.4f})')
            plt.legend()
            plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(Path(__file__).parent / "12_multiple_fitting.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 12_multiple_fitting.png")

    except Exception as e:
        print(f"  Multiple fitting failed: {e}")

    print()


def main():
    """Run all lineshape demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_basic_lineshapes()
    demo_derivatives()
    demo_phase_effects()
    demo_voigt_profiles()
    demo_lineshape_class()
    demo_signal_fitting()
    demo_multiple_fitting()

    print("=== Lineshapes Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- Gaussian, Lorentzian, and Voigt profiles available")
    print("- Derivative forms (0th, 1st, 2nd) supported")
    print("- Phase mixing for absorption/dispersion combinations")
    print("- Unified Lineshape class interface")
    print("- Comprehensive fitting capabilities with error analysis")
    print("- Multiple model comparison with statistical criteria")
    print("- All functions are area-normalized by default")
    print()
    print("Generated plot files:")
    for plot_file in sorted(output_dir.glob("*.png")):
        if plot_file.name.startswith(('07_', '08_', '09_', '10_', '11_', '12_')):
            print(f"  - {plot_file.name}")


if __name__ == "__main__":
    main()