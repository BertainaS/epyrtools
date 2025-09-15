#!/usr/bin/env python3
"""
Advanced Baseline Correction Demonstration

This script demonstrates the complete suite of advanced baseline correction
functions including stretched exponential, bi-exponential, and automatic
model selection for EPR data analysis.
"""

import epyr
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("=== ADVANCED BASELINE CORRECTION DEMO ===")
    print("Demonstrating stretched exponential, bi-exponential, and automatic model selection\n")
    
    # Create synthetic test data
    np.random.seed(42)
    
    # 1. Stretched exponential decay (T2 relaxation-like)
    print("📊 1. STRETCHED EXPONENTIAL BASELINE CORRECTION")
    print("-" * 50)
    
    x1 = np.linspace(0, 2000, 400)
    y1_clean = epyr.stretched_exponential_1d(x1, 800, 400, 1.3, 50)
    y1_noisy = y1_clean + 20 * np.random.normal(size=len(x1))
    
    print("Testing stretched exponential data (β=1.3)...")
    corrected1, baseline1 = epyr.baseline_stretched_exponential_1d(
        x1, y1_noisy, None, 
        beta_range=(0.5, 2.5)
    )
    
    # 2. Bi-exponential decay (complex relaxation)
    print("\n📊 2. BI-EXPONENTIAL BASELINE CORRECTION")
    print("-" * 50)
    
    x2 = np.linspace(0, 1500, 350)
    y2_clean = epyr.bi_exponential_1d(x2, 600, 150, 300, 600, 40)
    y2_noisy = y2_clean + 15 * np.random.normal(size=len(x2))
    
    print("Testing bi-exponential data (fast + slow components)...")
    corrected2, baseline2 = epyr.baseline_bi_exponential_1d(
        x2, y2_noisy, None,
        tau_ratio_min=2.5
    )
    
    # 3. Automatic model selection
    print("\n📊 3. AUTOMATIC MODEL SELECTION")
    print("-" * 50)
    
    # Test with different data types
    test_data = [
        ("Polynomial (quadratic)", 
         np.linspace(0, 100, 200),
         lambda x: 0.2 * x**2 - 8 * x + 80 + 5 * np.random.normal(size=len(x))),
        
        ("Stretched exponential", 
         x1, 
         lambda x: y1_noisy),
         
        ("Bi-exponential",
         x2,
         lambda x: y2_noisy)
    ]
    
    auto_results = []
    
    for name, x_test, y_func in test_data:
        y_test = y_func(x_test)
        
        print(f"\nTesting {name} data:")
        corrected_auto, baseline_auto, info = epyr.baseline_auto_1d(
            x_test, y_test, None, 
            verbose=False
        )
        
        print(f"  🏆 Best model: {info['best_model']}")
        print(f"  📊 Model comparison:")
        for model, aic_val in info['criteria'].items():
            marker = "🥇" if model == info['best_model'] else "  "
            print(f"    {marker} {model}: AIC={aic_val:.1f}")
        
        auto_results.append((name, info['best_model'], info['parameters']['r2']))
    
    # 4. Create visualization
    print("\n📊 4. CREATING COMPARISON PLOTS")
    print("-" * 50)
    
    plt.figure(figsize=(15, 10))
    
    # Stretched exponential example
    plt.subplot(2, 3, 1)
    plt.plot(x1, y1_noisy, 'b-', alpha=0.7, label='Noisy Data')
    plt.plot(x1, baseline1, 'r--', label='Fitted Baseline')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.title('Stretched Exponential\nBaseline Correction')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 3, 2)
    plt.plot(x1, corrected1, 'g-', label='Corrected')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.title('Baseline-Corrected Data')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Bi-exponential example
    plt.subplot(2, 3, 3)
    plt.plot(x2, y2_noisy, 'b-', alpha=0.7, label='Noisy Data')
    plt.plot(x2, baseline2, 'r--', label='Fitted Baseline')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.title('Bi-exponential\nBaseline Correction')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 3, 4)
    plt.plot(x2, corrected2, 'g-', label='Corrected')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.title('Baseline-Corrected Data')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Automatic selection summary
    plt.subplot(2, 3, 5)
    models = [result[1] for result in auto_results]
    model_counts = {model: models.count(model) for model in set(models)}
    
    plt.bar(model_counts.keys(), model_counts.values(), 
            color=['skyblue', 'lightcoral', 'lightgreen'][:len(model_counts)])
    plt.xlabel('Selected Model')
    plt.ylabel('Count')
    plt.title('Automatic Model Selection\nResults')
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3)
    
    # R² comparison
    plt.subplot(2, 3, 6)
    names = [result[0] for result in auto_results]
    r2_values = [result[2] for result in auto_results]
    
    plt.bar(range(len(names)), r2_values, 
            color=['skyblue', 'lightcoral', 'lightgreen'])
    plt.xlabel('Data Type')
    plt.ylabel('R²')
    plt.title('Fit Quality (R²)')
    plt.xticks(range(len(names)), [name.split()[0] for name in names], rotation=45)
    plt.ylim(0.98, 1.0)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('advanced_baseline_correction_demo.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("✅ Visualization saved as 'advanced_baseline_correction_demo.png'")
    
    # 5. Test with real EPR data if available
    print("\n📊 5. REAL EPR DATA ANALYSIS")
    print("-" * 50)
    
    dmttfi_file = "../data/20210508_DMTTFI_T2EH_5p8K_10dB_20_26ns_hperpc.DSC"
    try:
        print("Loading DMTTFI T2 relaxation data...")
        x_real, y_real, params_real, _ = epyr.eprload(dmttfi_file, plot_if_possible=False)
        
        print("Performing automatic baseline model selection...")
        corrected_real, baseline_real, info_real = epyr.baseline_auto_1d(
            x_real, y_real, params_real,
            use_real_part=True,
            exclude_initial=10,
            exclude_final=50,
            verbose=False
        )
        
        print(f"🏆 Best model for DMTTFI data: {info_real['best_model']}")
        print(f"📊 Fit quality: R² = {info_real['parameters']['r2']:.4f}")
        print(f"📊 Model comparison:")
        for model, aic_val in info_real['criteria'].items():
            marker = "🥇" if model == info_real['best_model'] else "  "
            print(f"    {marker} {model}: AIC={aic_val:.1f}")
        
        # Plot real data results
        plt.figure(figsize=(12, 4))
        
        plt.subplot(1, 3, 1)
        plt.plot(x_real, np.real(y_real), 'b-', alpha=0.7, label='Original')
        plt.plot(x_real, baseline_real, 'r--', label='Baseline')
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal')
        plt.title('DMTTFI T2 Data\nwith Fitted Baseline')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 3, 2)
        plt.plot(x_real, np.real(corrected_real), 'g-', label='Corrected')
        plt.xlabel('Time (ns)')
        plt.ylabel('Signal') 
        plt.title('Baseline-Corrected\nDMTTFI Data')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 3, 3)
        plt.semilogy(x_real, np.abs(y_real), 'b-', alpha=0.7, label='Original |Signal|')
        plt.semilogy(x_real, np.abs(corrected_real), 'g-', label='Corrected |Signal|')
        plt.xlabel('Time (ns)')
        plt.ylabel('log|Signal|')
        plt.title('Semi-log View')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('dmttfi_baseline_correction.png', dpi=150, bbox_inches='tight')
        plt.show()
        
        print("✅ DMTTFI analysis plot saved as 'dmttfi_baseline_correction.png'")
        
    except FileNotFoundError:
        print("⚠️  DMTTFI data file not found, skipping real data analysis")
    except Exception as e:
        print(f"⚠️  DMTTFI analysis failed: {e}")
    
    # 6. Usage examples
    print("\n📚 USAGE EXAMPLES")
    print("-" * 50)
    print("""
# Basic automatic model selection:
corrected, baseline, info = epyr.baseline_auto_1d(x, y, params)
print(f"Best model: {info['best_model']}")

# Stretched exponential with custom β range:
corrected, baseline = epyr.baseline_stretched_exponential_1d(
    x, y, params, beta_range=(0.1, 3.0)
)

# Bi-exponential with custom τ ratio constraint:
corrected, baseline = epyr.baseline_bi_exponential_1d(
    x, y, params, tau_ratio_min=3.0
)

# Interactive region selection (works in Jupyter):
corrected, baseline = epyr.baseline_auto_1d(x, y, params, interactive=True)

# Use BIC instead of AIC for model selection:
corrected, baseline, info = epyr.baseline_auto_1d(
    x, y, params, selection_criterion='bic'
)

# Restrict to specific models only:
corrected, baseline, info = epyr.baseline_auto_1d(
    x, y, params, models=['polynomial', 'stretched_exponential']
)
""")
    
    print("\n🎉 ADVANCED BASELINE CORRECTION DEMO COMPLETED!")
    print("\n🚀 NEW FEATURES SUMMARY:")
    print("✅ Stretched exponential baseline correction (β = 0.01 to 5.0)")
    print("✅ Bi-exponential baseline correction with automatic component separation") 
    print("✅ Automatic model selection using AIC/BIC/R² criteria")
    print("✅ Support for complex EPR data (T1, T2 measurements)")
    print("✅ Interactive region selection with Jupyter compatibility")
    print("✅ Parameter uncertainty estimation and model comparison")
    print("✅ Smart initial parameter guessing for robust fitting")
    
    print(f"\n📈 DEMO RESULTS SUMMARY:")
    for name, model, r2 in auto_results:
        print(f"  {name}: {model} (R² = {r2:.3f})")


if __name__ == "__main__":
    main()