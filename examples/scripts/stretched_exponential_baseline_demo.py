#!/usr/bin/env python3
"""
Demonstration of Stretched Exponential Baseline Correction

This script shows how to use the stretched exponential baseline correction
for EPR relaxation data (T1, T2 measurements).
"""

import epyr
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("=== STRETCHED EXPONENTIAL BASELINE CORRECTION DEMO ===")
    
    # Load DMTTFI T2 relaxation data (if available)
    dmttfi_file = "../data/20210508_DMTTFI_T2EH_5p8K_10dB_20_26ns_hperpc.DSC"
    
    try:
        print("📁 Loading DMTTFI T2 relaxation data...")
        x, y, params, filepath = epyr.eprload(dmttfi_file, plot_if_possible=False)
        print(f"   Data loaded: {len(y)} points, {'complex' if np.iscomplexobj(y) else 'real'} data")
        
        # Basic stretched exponential correction
        print("\n🧪 Performing basic stretched exponential baseline correction...")
        corrected, baseline_fit = epyr.baseline_correction.baseline_stretched_exponential_1d(
            x, y, params,
            beta_range=(0.01, 5.0),  # Allow full range for stretching parameter
            use_real_part=True       # Fit real part of complex data
        )
        
        # Create comparison plot
        plt.figure(figsize=(12, 8))
        
        # Plot original data (real part)
        plt.subplot(2, 2, 1)
        plt.plot(x, np.real(y), 'b-', alpha=0.7, label='Original Data (Real)')
        plt.plot(x, baseline_fit, 'r--', label='Fitted Baseline')
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal')
        plt.title('Original Data with Fitted Baseline')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Plot corrected data
        plt.subplot(2, 2, 2)
        plt.plot(x, np.real(corrected), 'g-', label='Corrected Data (Real)')
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal')
        plt.title('Baseline-Corrected Data')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Plot magnitude comparison
        plt.subplot(2, 2, 3)
        plt.plot(x, np.abs(y), 'b-', alpha=0.7, label='Original Magnitude')
        plt.plot(x, np.abs(corrected), 'g-', label='Corrected Magnitude')
        plt.xlabel('Time (ns)')
        plt.ylabel('|Signal|')
        plt.title('Magnitude Comparison')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Semi-log plot to show exponential decay
        plt.subplot(2, 2, 4)
        positive_orig = np.abs(y)
        positive_corr = np.abs(corrected)
        positive_orig = np.maximum(positive_orig, positive_orig.max() * 1e-6)
        positive_corr = np.maximum(positive_corr, positive_corr.max() * 1e-6)
        
        plt.semilogy(x, positive_orig, 'b-', alpha=0.7, label='Original |Signal|')
        plt.semilogy(x, positive_corr, 'g-', label='Corrected |Signal|')
        plt.xlabel('Time (ns)')
        plt.ylabel('log|Signal|')
        plt.title('Semi-log View')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('stretched_exponential_baseline_demo.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print("✅ Plot saved as 'stretched_exponential_baseline_demo.png'")
        
    except FileNotFoundError:
        print("⚠️  DMTTFI data file not found, creating synthetic data...")
        
        # Create synthetic stretched exponential decay data
        x = np.linspace(0, 2000, 500)  # Time points
        
        # True parameters
        A_true = 1000
        tau_true = 500
        beta_true = 1.2  # Stretched exponential
        offset_true = 50
        
        # Generate noisy data with stretched exponential decay
        y_clean = epyr.baseline_correction.stretched_exponential_1d(x, A_true, tau_true, beta_true, offset_true)
        noise = 20 * np.random.normal(size=len(x))
        y_synthetic = y_clean + noise
        
        print(f"   Created synthetic data with β={beta_true}")
        
        # Apply correction
        print("\n🧪 Performing stretched exponential correction on synthetic data...")
        corrected, baseline_fit = epyr.baseline_correction.baseline_stretched_exponential_1d(
            x, y_synthetic, None,
            beta_range=(0.5, 2.0)  # Constrain beta range for synthetic data
        )
        
        # Create comparison plot
        plt.figure(figsize=(10, 6))
        
        plt.subplot(1, 2, 1)
        plt.plot(x, y_synthetic, 'b-', alpha=0.7, label='Synthetic Data')
        plt.plot(x, baseline_fit, 'r--', linewidth=2, label='Fitted Baseline')
        plt.plot(x, y_clean, 'k:', alpha=0.5, label='True Baseline')
        plt.xlabel('Time')
        plt.ylabel('Signal')
        plt.title('Synthetic Data with Fitted Baseline')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        plt.plot(x, corrected, 'g-', label='Corrected Data')
        plt.plot(x, noise, 'k:', alpha=0.5, label='True Noise')
        plt.xlabel('Time')
        plt.ylabel('Signal')
        plt.title('Baseline-Corrected Data')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('stretched_exponential_synthetic_demo.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print("✅ Synthetic data plot saved as 'stretched_exponential_synthetic_demo.png'")
    
    # Advanced usage examples
    print("\n📚 ADVANCED USAGE EXAMPLES:")
    print("=" * 50)
    
    print("""
# Basic usage:
corrected, baseline = epyr.baseline_correction.baseline_stretched_exponential_1d(x, y, params)

# With custom beta range:
corrected, baseline = epyr.baseline_correction.baseline_stretched_exponential_1d(
    x, y, params, beta_range=(0.1, 3.0)
)

# Interactive region selection:
corrected, baseline = epyr.baseline_correction.baseline_stretched_exponential_1d(
    x, y, params, interactive=True
)

# For complex data, use magnitude instead of real part:
corrected, baseline = epyr.baseline_correction.baseline_stretched_exponential_1d(
    x, y, params, use_real_part=False
)

# Exclude problematic points:
corrected, baseline = epyr.baseline_correction.baseline_stretched_exponential_1d(
    x, y, params, 
    exclude_initial=10,  # Skip first 10 points
    exclude_final=50     # Skip last 50 points
)

# Custom initial parameter guess:
initial_params = {'A': 1000, 'tau': 500, 'beta': 1.5, 'offset': 0}
corrected, baseline = epyr.baseline_correction.baseline_stretched_exponential_1d(
    x, y, params, initial_guess=initial_params
)
""")
    
    print("\n🎉 Stretched Exponential Baseline Correction Demo Complete!")
    print("\n📊 KEY FEATURES:")
    print("✅ Fits stretched exponential: f(x) = offset + A*exp(-(x/τ)^β)")
    print("✅ Beta parameter constrained between 0.01 and 5.0")
    print("✅ Handles complex EPR data (T1, T2 measurements)")  
    print("✅ Interactive region selection for complex baselines")
    print("✅ Smart initial parameter estimation")
    print("✅ Parameter uncertainty estimation")


if __name__ == "__main__":
    main()