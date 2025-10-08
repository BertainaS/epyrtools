#!/usr/bin/env python3
"""
EPyR Tools - Lineshape Analysis Example

This example demonstrates the comprehensive lineshape analysis capabilities
of EPyR Tools, including Gaussian, Lorentzian, Voigt, and pseudo-Voigt profiles
with derivatives, phase rotation, and convolution.

Author: EPyR Tools Development Team
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Import EPyR Tools lineshape functions
import epyr
from epyr.lineshapes import (
    Lineshape, gaussian, lorentzian, voigtian, pseudo_voigt, 
    create_gaussian, create_lorentzian, create_voigt, convspec
)


def main():
    """Main lineshape analysis demonstration"""
    
    print("🧲 EPyR Tools - Lineshape Analysis Example")
    print("=" * 50)
    
    # Create analysis data
    x = np.linspace(-15, 15, 1000)
    
    # 1. Basic Lineshape Comparison
    print("\n1. Comparing Basic Lineshapes...")
    compare_basic_lineshapes(x)
    
    # 2. Derivative Analysis
    print("\n2. Analyzing Derivatives...")
    analyze_derivatives(x)
    
    # 3. Phase Rotation (Absorption vs Dispersion)
    print("\n3. Phase Rotation Analysis...")
    analyze_phase_rotation(x)
    
    # 4. Voigt vs Pseudo-Voigt Comparison
    print("\n4. Voigt Profile Analysis...")
    compare_voigt_profiles(x)
    
    # 5. Using the Lineshape Class
    print("\n5. Using the Unified Lineshape Class...")
    demonstrate_lineshape_class(x)
    
    # 6. Spectrum Convolution
    print("\n6. Spectrum Convolution...")
    demonstrate_convolution(x)
    
    # 7. EPR-specific Applications
    print("\n7. EPR-specific Applications...")
    epr_applications(x)
    
    print("\n✅ Lineshape analysis complete!")
    print("📊 Check the generated plots in the current directory.")


def compare_basic_lineshapes(x):
    """Compare Gaussian, Lorentzian, and mixed profiles"""
    
    # Generate different lineshapes with same FWHM
    width = 6.0
    
    gauss = gaussian(x, 0, width)
    lorentz = lorentzian(x, 0, width)
    pv_50 = pseudo_voigt(x, 0, width, eta=0.5)  # 50/50 mix
    
    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Linear scale
    ax1.plot(x, gauss, 'b-', linewidth=2.5, label='Gaussian')
    ax1.plot(x, lorentz, 'r-', linewidth=2.5, label='Lorentzian')
    ax1.plot(x, pv_50, 'g--', linewidth=2, label='Pseudo-Voigt (50/50)')
    
    ax1.set_xlabel('Magnetic Field (mT)')
    ax1.set_ylabel('EPR Signal (normalized)')
    ax1.set_title('Lineshape Comparison (Linear Scale)')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Log scale to show tails
    ax2.semilogy(x, gauss, 'b-', linewidth=2.5, label='Gaussian')
    ax2.semilogy(x, lorentz, 'r-', linewidth=2.5, label='Lorentzian')
    ax2.semilogy(x, pv_50, 'g--', linewidth=2, label='Pseudo-Voigt (50/50)')
    
    ax2.set_xlabel('Magnetic Field (mT)')
    ax2.set_ylabel('EPR Signal (log scale)')
    ax2.set_title('Tail Behavior (Log Scale)')
    ax2.legend()
    ax2.grid(alpha=0.3)
    ax2.set_ylim(1e-6, 1)
    
    plt.tight_layout()
    plt.savefig('lineshapes_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print FWHM verification
    for name, y in [('Gaussian', gauss), ('Lorentzian', lorentz), ('Pseudo-Voigt', pv_50)]:
        half_max = np.max(y) / 2
        indices = np.where(y >= half_max)[0]
        fwhm = x[indices[-1]] - x[indices[0]]
        print(f"  {name}: FWHM = {fwhm:.2f} mT (expected: {width:.2f})")


def analyze_derivatives(x):
    """Analyze derivative lineshapes for EPR applications"""
    
    width = 5.0
    
    # Generate function and derivatives
    gauss_0 = gaussian(x, 0, width, derivative=0)
    gauss_1 = gaussian(x, 0, width, derivative=1) 
    gauss_2 = gaussian(x, 0, width, derivative=2)
    
    lorentz_0 = lorentzian(x, 0, width, derivative=0)
    lorentz_1 = lorentzian(x, 0, width, derivative=1)
    lorentz_2 = lorentzian(x, 0, width, derivative=2)
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Gaussian derivatives
    axes[0, 0].plot(x, gauss_0, 'b-', linewidth=2.5)
    axes[0, 0].set_title('Gaussian Function')
    axes[0, 0].set_ylabel('Gaussian')
    axes[0, 0].grid(alpha=0.3)
    
    axes[0, 1].plot(x, gauss_1, 'b-', linewidth=2.5)
    axes[0, 1].set_title('First Derivative')
    axes[0, 1].grid(alpha=0.3)
    
    axes[0, 2].plot(x, gauss_2, 'b-', linewidth=2.5)
    axes[0, 2].set_title('Second Derivative')
    axes[0, 2].grid(alpha=0.3)
    
    # Lorentzian derivatives
    axes[1, 0].plot(x, lorentz_0, 'r-', linewidth=2.5)
    axes[1, 0].set_title('Lorentzian Function')
    axes[1, 0].set_ylabel('Lorentzian')
    axes[1, 0].set_xlabel('Magnetic Field (mT)')
    axes[1, 0].grid(alpha=0.3)
    
    axes[1, 1].plot(x, lorentz_1, 'r-', linewidth=2.5)
    axes[1, 1].set_title('First Derivative')
    axes[1, 1].set_xlabel('Magnetic Field (mT)')
    axes[1, 1].grid(alpha=0.3)
    
    axes[1, 2].plot(x, lorentz_2, 'r-', linewidth=2.5)
    axes[1, 2].set_title('Second Derivative')
    axes[1, 2].set_xlabel('Magnetic Field (mT)')
    axes[1, 2].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('derivatives_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("  📈 Derivatives useful for:")
    print("    - First derivative: Enhanced resolution, zero-crossing detection")
    print("    - Second derivative: Peak identification, overlapping signal separation")


def analyze_phase_rotation(x):
    """Analyze absorption vs dispersion modes"""
    
    width = 6.0
    phases = np.array([0, np.pi/6, np.pi/4, np.pi/3, np.pi/2])
    phase_labels = ['0° (Absorption)', '30°', '45°', '60°', '90° (Dispersion)']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    colors = plt.cm.plasma(np.linspace(0, 1, len(phases)))
    
    # Gaussian phase rotation
    for i, (phase, label) in enumerate(zip(phases, phase_labels)):
        y = gaussian(x, 0, width, phase=phase)
        ax1.plot(x, y, color=colors[i], linewidth=2.5, label=label)
    
    ax1.set_xlabel('Magnetic Field (mT)')
    ax1.set_ylabel('EPR Signal')
    ax1.set_title('Gaussian Phase Rotation')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Lorentzian phase rotation
    for i, (phase, label) in enumerate(zip(phases, phase_labels)):
        y = lorentzian(x, 0, width, phase=phase)
        ax2.plot(x, y, color=colors[i], linewidth=2.5, label=label)
    
    ax2.set_xlabel('Magnetic Field (mT)')
    ax2.set_ylabel('EPR Signal')
    ax2.set_title('Lorentzian Phase Rotation')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('phase_rotation.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("  🌊 Phase rotation applications:")
    print("    - 0°: Pure absorption (symmetric)")
    print("    - 90°: Pure dispersion (antisymmetric)")
    print("    - Mixed phases: Complex impedance analysis")


def compare_voigt_profiles(x):
    """Compare true Voigt vs pseudo-Voigt profiles"""
    
    # Different mixing scenarios
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    scenarios = [
        {'gauss': 4, 'lorentz': 2, 'title': 'Gaussian-dominated'},
        {'gauss': 2, 'lorentz': 4, 'title': 'Lorentzian-dominated'},
        {'gauss': 3, 'lorentz': 3, 'title': 'Equal widths'},
        {'gauss': 6, 'lorentz': 1, 'title': 'Narrow Lorentzian'}
    ]
    
    for i, scenario in enumerate(scenarios):
        ax = axes[i//2, i%2]
        
        gauss_w = scenario['gauss']
        lorentz_w = scenario['lorentz']
        
        # True Voigt (convolution)
        try:
            voigt_true = voigtian(x, 0, (gauss_w, lorentz_w))
            ax.plot(x, voigt_true, 'b-', linewidth=2.5, label='True Voigt')
        except:
            print(f"    Skipping true Voigt for {scenario['title']} (not implemented)")
        
        # Pseudo-Voigt approximation with optimal mixing
        # Use empirical formula for eta parameter
        total_fwhm = np.sqrt(gauss_w**2 + lorentz_w**2)
        eta_approx = lorentz_w / total_fwhm  # Simplified mixing
        
        pv_approx = pseudo_voigt(x, 0, total_fwhm, eta=eta_approx)
        ax.plot(x, pv_approx, 'r--', linewidth=2, label=f'Pseudo-Voigt (η={eta_approx:.2f})')
        
        # Individual components
        gauss_comp = gaussian(x, 0, gauss_w) * (1 - eta_approx)
        lorentz_comp = lorentzian(x, 0, lorentz_w) * eta_approx
        
        ax.plot(x, gauss_comp, 'g:', alpha=0.7, label='Gaussian component')
        ax.plot(x, lorentz_comp, 'm:', alpha=0.7, label='Lorentzian component')
        
        ax.set_title(f'{scenario["title"]} (G={gauss_w}, L={lorentz_w})')
        ax.legend()
        ax.grid(alpha=0.3)
        ax.set_xlabel('Magnetic Field (mT)')
        ax.set_ylabel('EPR Signal')
    
    plt.tight_layout()
    plt.savefig('voigt_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()


def demonstrate_lineshape_class(x):
    """Demonstrate the unified Lineshape class"""
    
    print("  🎯 Creating lineshapes with unified interface...")
    
    # Create different lineshape objects
    gauss = Lineshape('gaussian', width=5.0)
    lorentz = Lineshape('lorentzian', width=5.0, derivative=1)
    pv = Lineshape('pseudo_voigt', width=5.0, alpha=0.3)
    
    # Generate data
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Basic lineshapes
    ax = axes[0, 0]
    ax.plot(x, gauss(x, 0), 'b-', linewidth=2.5, label='Gaussian')
    ax.plot(x, lorentz(x, 0), 'r-', linewidth=2.5, label='Lorentzian (1st deriv)')
    ax.plot(x, pv(x, 0), 'g-', linewidth=2.5, label='Pseudo-Voigt (α=0.3)')
    ax.set_title('Unified Lineshape Interface')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Method demonstrations
    ax = axes[0, 1]
    gauss_abs = gauss.absorption(x, 0)
    gauss_disp = gauss.dispersion(x, 0)
    ax.plot(x, gauss_abs, 'b-', linewidth=2.5, label='Absorption')
    ax.plot(x, gauss_disp, 'r-', linewidth=2.5, label='Dispersion')
    ax.set_title('Absorption vs Dispersion')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Parameter modification
    ax = axes[1, 0]
    gauss_narrow = gauss.set_width(3.0)
    gauss_wide = gauss.set_width(8.0)
    ax.plot(x, gauss(x, 0), 'b-', linewidth=2.5, label='Original (w=5)')
    ax.plot(x, gauss_narrow(x, 0), 'g-', linewidth=2, label='Narrow (w=3)')
    ax.plot(x, gauss_wide(x, 0), 'r-', linewidth=2, label='Wide (w=8)')
    ax.set_title('Width Modification')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Factory functions
    ax = axes[1, 1]
    factory_gauss = create_gaussian(width=5.0)
    factory_lorentz = create_lorentzian(width=5.0)
    factory_voigt = create_voigt(3.0, 2.0)
    
    ax.plot(x, factory_gauss(x, 0), 'b-', linewidth=2.5, label='Factory Gaussian')
    ax.plot(x, factory_lorentz(x, 0), 'r-', linewidth=2.5, label='Factory Lorentzian')
    ax.plot(x, factory_voigt(x, 0), 'g-', linewidth=2, label='Factory Voigt')
    ax.set_title('Factory Functions')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('lineshape_class_demo.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print object information
    print(f"    Gaussian info: {gauss}")
    print(f"    Info dict: {gauss.info()}")


def demonstrate_convolution(x):
    """Demonstrate spectrum convolution capabilities"""
    
    # Create stick spectrum (discrete peaks)
    stick_spectrum = np.zeros_like(x)
    peak_positions = [-6, -2, 3, 8]
    peak_intensities = [1.0, 0.7, 1.2, 0.5]
    
    for pos, intensity in zip(peak_positions, peak_intensities):
        idx = np.argmin(np.abs(x - pos))
        stick_spectrum[idx] = intensity
    
    # Apply different broadening
    dx = x[1] - x[0]
    
    # Gaussian broadening
    gauss_broad = convspec(stick_spectrum, dx, width=1.5, alpha=1.0)
    
    # Lorentzian broadening  
    lorentz_broad = convspec(stick_spectrum, dx, width=1.5, alpha=0.0)
    
    # Mixed broadening
    mixed_broad = convspec(stick_spectrum, dx, width=1.5, alpha=0.5)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Original stick spectrum
    ax1.stem(x, stick_spectrum, linefmt='k-', markerfmt='ko', basefmt=' ')
    ax1.set_title('Original Stick Spectrum')
    ax1.set_ylabel('Intensity')
    ax1.grid(alpha=0.3)
    ax1.set_xlim(x[0], x[-1])
    
    # Broadened spectra
    ax2.plot(x, gauss_broad, 'b-', linewidth=2.5, label='Gaussian broadening')
    ax2.plot(x, lorentz_broad, 'r-', linewidth=2.5, label='Lorentzian broadening')
    ax2.plot(x, mixed_broad, 'g--', linewidth=2, label='Mixed broadening (50/50)')
    
    ax2.set_title('Convolution Results')
    ax2.set_xlabel('Magnetic Field (mT)')
    ax2.set_ylabel('EPR Signal')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('convolution_demo.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("  📡 Convolution applications:")
    print("    - Converting stick spectra to realistic lineshapes")
    print("    - Modeling instrumental broadening")
    print("    - Simulating EPR spectra from calculated parameters")


def epr_applications(x):
    """Demonstrate EPR-specific applications"""
    
    print("  🧲 EPR-specific lineshape applications...")
    
    # 1. Multi-component spectrum fitting
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Simulate multi-component EPR spectrum
    centers = [-4, 1, 6]
    widths = [2.5, 1.8, 3.2]
    intensities = [1.0, 0.7, 1.2]
    
    # Individual components
    components = []
    for i, (center, width, intensity) in enumerate(zip(centers, widths, intensities)):
        component = intensity * lorentzian(x, center, width)
        components.append(component)
    
    total_spectrum = np.sum(components, axis=0)
    
    ax = axes[0, 0]
    for i, comp in enumerate(components):
        ax.plot(x, comp, '--', alpha=0.7, label=f'Component {i+1}')
    ax.plot(x, total_spectrum, 'k-', linewidth=2.5, label='Total spectrum')
    ax.set_title('Multi-component EPR Spectrum')
    ax.set_ylabel('EPR Signal')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 2. g-factor analysis using different lineshapes
    ax = axes[0, 1]
    
    # Simulate different radical environments
    radical_types = {
        'Organic radical': {'shape': 'gaussian', 'width': 0.5},
        'Metal complex': {'shape': 'lorentzian', 'width': 1.2},
        'Solid state': {'shape': 'pseudo_voigt', 'width': 2.0, 'alpha': 0.3}
    }
    
    g_center = 334.5  # mT, assuming X-band (9.5 GHz)
    
    for i, (name, params) in enumerate(radical_types.items()):
        if params['shape'] == 'gaussian':
            y = gaussian(x + g_center, g_center, params['width'])
        elif params['shape'] == 'lorentzian':
            y = lorentzian(x + g_center, g_center, params['width'])
        else:
            y = pseudo_voigt(x + g_center, g_center, params['width'], eta=params['alpha'])
        
        ax.plot(x, y, linewidth=2.5, label=name)
    
    ax.set_title('Lineshapes for Different Radical Types')
    ax.set_xlabel('Magnetic Field (mT)')
    ax.set_ylabel('EPR Signal')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 3. Temperature-dependent broadening
    ax = axes[1, 0]
    
    temperatures = [77, 200, 300, 400]  # Kelvin
    base_width = 1.0
    
    for T in temperatures:
        # Simple temperature broadening model
        width_T = base_width * (T / 300) ** 0.5
        y = lorentzian(x, 0, width_T)
        ax.plot(x, y, linewidth=2.5, label=f'T = {T} K')
    
    ax.set_title('Temperature-dependent Broadening')
    ax.set_xlabel('Magnetic Field (mT)')
    ax.set_ylabel('EPR Signal')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 4. Derivative spectroscopy (common in EPR)
    ax = axes[1, 1]
    
    # Overlapping signals - easier to resolve with derivatives
    signal1 = gaussian(x, -1, 3)
    signal2 = gaussian(x, 1, 3)
    overlapped = signal1 + signal2
    
    # First derivatives
    signal1_d1 = gaussian(x, -1, 3, derivative=1)
    signal2_d1 = gaussian(x, 1, 3, derivative=1)
    overlapped_d1 = signal1_d1 + signal2_d1
    
    ax.plot(x, overlapped, 'b-', linewidth=2.5, alpha=0.7, label='Overlapped absorption')
    ax.plot(x, overlapped_d1, 'r-', linewidth=2.5, label='First derivative')
    ax.axhline(0, color='k', linestyle=':', alpha=0.5)
    
    ax.set_title('Derivative Spectroscopy for Resolution')
    ax.set_xlabel('Magnetic Field (mT)')
    ax.set_ylabel('EPR Signal')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('epr_applications.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("    📊 Key EPR applications demonstrated:")
    print("      - Multi-component spectrum analysis")
    print("      - Radical type identification via lineshape")
    print("      - Temperature-dependent broadening")
    print("      - Derivative spectroscopy for enhanced resolution")


if __name__ == "__main__":
    main()