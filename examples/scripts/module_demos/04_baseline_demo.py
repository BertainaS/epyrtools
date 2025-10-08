#!/usr/bin/env python3
"""
EPyR Tools Demo 04: Baseline Correction
=======================================

This script demonstrates the comprehensive baseline correction capabilities.
The baseline module provides automated and interactive baseline correction.

Functions demonstrated:
- baseline_polynomial_1d() - 1D polynomial baseline correction
- baseline_polynomial_2d() - 2D baseline correction
- interactive() - Interactive baseline selection
- baseline_auto_1d() - Automatic baseline correction
- Models and selection utilities
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr
from epyr.baseline import baseline_polynomial_1d, baseline_polynomial_2d, baseline_auto_1d


def demo_1d_baseline_correction():
    """Demonstrate 1D baseline correction."""
    print("=== EPyR Tools Baseline Demo ===")
    print()

    print("1. 1D Baseline correction:")
    print("-" * 27)

    # Load real EPR data if available
    data_dir = Path(__file__).parent.parent.parent / "data"
    cw_file = data_dir / "130406SB_CaWO4_Er_CW_5K_20.DSC"

    if cw_file.exists():
        print("Loading real EPR data...")
        x, y_raw, params, filepath = epyr.eprload(str(cw_file), plot_if_possible=False)
        print(f"Loaded: {Path(filepath).name}")
    else:
        print("Creating synthetic EPR data with baseline...")
        x = np.linspace(3350, 3450, 512)
        # EPR signal
        signal = 0.8 * np.exp(-((x - 3400)/12)**2)
        # Baseline drift
        baseline_true = 0.2 + 0.0002 * (x - 3400) + 1e-6 * (x - 3400)**2
        # Noise
        noise = 0.03 * np.random.normal(size=len(x))
        y_raw = signal + baseline_true + noise

    print(f"Data points: {len(y_raw)}")
    print(f"Signal range: {y_raw.min():.3f} to {y_raw.max():.3f}")

    # Demonstrate different polynomial orders
    orders = [1, 2, 3, 4]

    plt.figure(figsize=(15, 10))

    for i, order in enumerate(orders):
        print(f"\nTesting polynomial order {order}:")

        try:
            # Perform baseline correction
            y_corrected, baseline = baseline_polynomial_1d(x, y_raw, order=order)

            # Calculate improvement metrics
            original_baseline = np.mean([y_raw[:20].mean(), y_raw[-20:].mean()])
            corrected_baseline = np.mean([y_corrected[:20].mean(), y_corrected[-20:].mean()])
            baseline_reduction = abs(corrected_baseline / original_baseline) if original_baseline != 0 else 0

            signal_max_original = np.max(y_raw)
            signal_max_corrected = np.max(y_corrected)

            print(f"  Original baseline level: {original_baseline:.4f}")
            print(f"  Corrected baseline level: {corrected_baseline:.4f}")
            print(f"  Baseline reduction factor: {1/baseline_reduction:.1f}x" if baseline_reduction > 0 else "N/A")
            print(f"  Signal enhancement: {signal_max_corrected/signal_max_original:.2f}x")

            # Plot results
            plt.subplot(2, 2, i+1)
            plt.plot(x, y_raw, 'gray', alpha=0.7, label='Original', linewidth=1)
            plt.plot(x, y_corrected, 'blue', label=f'Corrected (order {order})', linewidth=2)
            plt.xlabel('Magnetic Field (G)')
            plt.ylabel('Intensity')
            plt.title(f'Polynomial Order {order}')
            plt.legend()
            plt.grid(True, alpha=0.3)

        except Exception as e:
            print(f"  Error with order {order}: {e}")

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "13_1d_baseline_orders.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 13_1d_baseline_orders.png")
    print()


def demo_baseline_with_signal_exclusion():
    """Demonstrate baseline correction with signal region exclusion."""
    print("2. Baseline correction with signal exclusion:")
    print("-" * 45)

    # Create EPR data with strong signal and baseline
    x = np.linspace(3300, 3500, 600)

    # Multiple EPR signals
    signal1 = 1.2 * np.exp(-((x - 3350)/8)**2)
    signal2 = 0.8 * np.exp(-((x - 3420)/12)**2)
    signal3 = 0.6 * np.exp(-((x - 3470)/10)**2)
    signal_total = signal1 + signal2 + signal3

    # Significant baseline
    baseline_true = 0.5 + 0.001 * (x - 3400) + 2e-6 * (x - 3400)**2 - 1e-9 * (x - 3400)**3

    # Noise
    noise = 0.05 * np.random.normal(size=len(x))

    y_with_baseline = signal_total + baseline_true + noise

    print(f"Created synthetic multi-line EPR spectrum")
    print(f"  3 EPR lines at 3350, 3420, 3470 G")
    print(f"  Cubic baseline polynomial")
    print(f"  SNR ≈ {np.max(signal_total) / np.std(noise):.1f}")

    # Correct without signal exclusion
    y_corrected_inclusive, baseline_inclusive = baseline_polynomial_1d(x, y_with_baseline, order=3)

    # Try to find signal regions automatically for exclusion
    try:
        # Simple signal detection: points above mean + 2*std
        threshold = np.mean(y_with_baseline) + 2 * np.std(y_with_baseline)
        signal_mask = y_with_baseline > threshold

        # Create exclusion ranges around detected signals
        signal_indices = np.where(signal_mask)[0]
        if len(signal_indices) > 0:
            # Group consecutive indices and expand ranges
            signal_ranges = []
            current_start = signal_indices[0]
            current_end = signal_indices[0]

            for idx in signal_indices[1:]:
                if idx - current_end <= 10:  # Group nearby points
                    current_end = idx
                else:
                    # Expand range
                    start_expanded = max(0, current_start - 15)
                    end_expanded = min(len(x) - 1, current_end + 15)
                    signal_ranges.append((x[start_expanded], x[end_expanded]))
                    current_start = idx
                    current_end = idx

            # Add final range
            start_expanded = max(0, current_start - 15)
            end_expanded = min(len(x) - 1, current_end + 15)
            signal_ranges.append((x[start_expanded], x[end_expanded]))

            print(f"  Detected {len(signal_ranges)} signal regions for exclusion")
            for i, (start, end) in enumerate(signal_ranges):
                print(f"    Region {i+1}: {start:.1f} - {end:.1f} G")

            # Perform correction with exclusion
            y_corrected_exclusive, baseline_exclusive = baseline_polynomial_1d(x, y_with_baseline, order=3,
                                             manual_regions=signal_ranges, region_mode='exclude')
        else:
            print("  No signal regions detected")
            y_corrected_exclusive = y_corrected_inclusive
            signal_ranges = []

    except Exception as e:
        print(f"  Error in signal exclusion: {e}")
        y_corrected_exclusive = y_corrected_inclusive
        signal_ranges = []

    # Plot comparison
    plt.figure(figsize=(15, 10))

    plt.subplot(2, 2, 1)
    plt.plot(x, signal_total, 'g-', linewidth=2, label='True signal')
    plt.plot(x, baseline_true, 'r--', linewidth=2, label='True baseline')
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Intensity')
    plt.title('True Signal and Baseline Components')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 2, 2)
    plt.plot(x, y_with_baseline, 'k-', linewidth=1, alpha=0.7, label='Raw data')
    # Highlight exclusion regions
    for start, end in signal_ranges:
        mask = (x >= start) & (x <= end)
        plt.fill_between(x[mask], y_with_baseline[mask].min(), y_with_baseline[mask].max(),
                        alpha=0.3, color='yellow', label='Excluded' if start == signal_ranges[0][0] else "")
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Intensity')
    plt.title('Raw Data with Exclusion Regions')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 2, 3)
    plt.plot(x, y_with_baseline, 'gray', alpha=0.5, label='Original')
    plt.plot(x, y_corrected_inclusive, 'blue', linewidth=2, label='Corrected (inclusive)')
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Intensity')
    plt.title('Baseline Correction - Inclusive')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 2, 4)
    plt.plot(x, y_with_baseline, 'gray', alpha=0.5, label='Original')
    plt.plot(x, y_corrected_exclusive, 'red', linewidth=2, label='Corrected (exclusive)')
    plt.plot(x, signal_total, 'g--', alpha=0.7, label='True signal')
    plt.xlabel('Magnetic Field (G)')
    plt.ylabel('Intensity')
    plt.title('Baseline Correction - Signal Excluded')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "14_baseline_exclusion.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 14_baseline_exclusion.png")

    # Calculate metrics
    if len(signal_ranges) > 0:
        mse_inclusive = np.mean((y_corrected_inclusive - signal_total)**2)
        mse_exclusive = np.mean((y_corrected_exclusive - signal_total)**2)
        improvement = mse_inclusive / mse_exclusive if mse_exclusive > 0 else 1

        print(f"  MSE (inclusive): {mse_inclusive:.6f}")
        print(f"  MSE (exclusive): {mse_exclusive:.6f}")
        print(f"  Improvement factor: {improvement:.1f}x")

    print()


def demo_2d_baseline_correction():
    """Demonstrate 2D baseline correction."""
    print("3. 2D Baseline correction:")
    print("-" * 27)

    # Load real 2D data if available
    data_dir = Path(__file__).parent.parent.parent / "data"

    # Try to load 2D data
    angular_file = data_dir / "2014_03_19_MgO_300K_111_fullrotation33dB.par"

    if angular_file.exists():
        print("Loading real 2D EPR data...")
        try:
            x_2d, y_2d_raw, params_2d, filepath = epyr.eprload(str(angular_file), plot_if_possible=False)
            print(f"Loaded: {Path(filepath).name}")
            print(f"2D data shape: {y_2d_raw.shape}")
        except Exception as e:
            print(f"Error loading 2D data: {e}")
            x_2d, y_2d_raw, params_2d = create_synthetic_2d_with_baseline()
    else:
        print("Creating synthetic 2D EPR data with baseline...")
        x_2d, y_2d_raw, params_2d = create_synthetic_2d_with_baseline()

    print(f"2D data shape: {y_2d_raw.shape}")

    # Perform 2D baseline correction
    try:
        print("Performing 2D baseline correction...")
        y_2d_corrected, baseline_2d = baseline_polynomial_2d(y_2d_raw, order=2)

        # Calculate improvement
        original_corners_mean = np.mean([
            y_2d_raw[:5, :5].mean(), y_2d_raw[-5:, :5].mean(),
            y_2d_raw[:5, -5:].mean(), y_2d_raw[-5:, -5:].mean()
        ])
        corrected_corners_mean = np.mean([
            y_2d_corrected[:5, :5].mean(), y_2d_corrected[-5:, :5].mean(),
            y_2d_corrected[:5, -5:].mean(), y_2d_corrected[-5:, -5:].mean()
        ])

        print(f"  Original baseline (corners): {original_corners_mean:.4f}")
        print(f"  Corrected baseline (corners): {corrected_corners_mean:.4f}")
        print(f"  Baseline reduction: {abs(original_corners_mean/corrected_corners_mean):.1f}x"
              if corrected_corners_mean != 0 else "N/A")

        # Plot 2D correction results
        plt.figure(figsize=(15, 6))

        # Original data
        plt.subplot(1, 3, 1)
        if isinstance(x_2d, list) and len(x_2d) == 2:
            extent = [x_2d[1][0], x_2d[1][-1], x_2d[0][0], x_2d[0][-1]]
        else:
            extent = None

        im1 = plt.imshow(y_2d_raw, aspect='auto', origin='lower', extent=extent, cmap='RdBu_r')
        plt.colorbar(im1, label='Intensity')
        plt.title('Original 2D Data')
        plt.xlabel('Field/Angle Axis')
        plt.ylabel('Other Axis')

        # Corrected data
        plt.subplot(1, 3, 2)
        im2 = plt.imshow(y_2d_corrected, aspect='auto', origin='lower', extent=extent, cmap='RdBu_r')
        plt.colorbar(im2, label='Intensity')
        plt.title('Baseline Corrected')
        plt.xlabel('Field/Angle Axis')
        plt.ylabel('Other Axis')

        # Difference (baseline removed)
        baseline_removed = y_2d_raw - y_2d_corrected
        plt.subplot(1, 3, 3)
        im3 = plt.imshow(baseline_removed, aspect='auto', origin='lower', extent=extent, cmap='viridis')
        plt.colorbar(im3, label='Intensity')
        plt.title('Removed Baseline')
        plt.xlabel('Field/Angle Axis')
        plt.ylabel('Other Axis')

        plt.tight_layout()
        plt.savefig(Path(__file__).parent / "15_2d_baseline.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 15_2d_baseline.png")

    except Exception as e:
        print(f"  2D baseline correction failed: {e}")

    print()


def create_synthetic_2d_with_baseline():
    """Create synthetic 2D EPR data with baseline for demonstration."""
    # Create 2D synthetic data
    field_axis = np.linspace(3300, 3500, 60)
    angle_axis = np.linspace(0, 180, 40)

    Field, Angle = np.meshgrid(field_axis, angle_axis)

    # EPR signal with angular dependence
    g_parallel = 2.0
    g_perp = 2.2
    theta_rad = np.radians(Angle)
    g_eff = np.sqrt(g_parallel**2 * np.cos(theta_rad)**2 + g_perp**2 * np.sin(theta_rad)**2)

    # Resonance field calculation
    h = 6.626e-34
    mu_b = 9.274e-24
    freq = 9.5e9
    B_res = h * freq / (mu_b * g_eff) * 1e4  # Convert to Gauss

    # EPR signal
    signal = np.exp(-((Field - B_res)/12)**2)

    # Add 2D baseline
    baseline = (0.2 + 0.0001 * (Field - 3400) +
               0.0002 * (Angle - 90) +
               1e-7 * (Field - 3400) * (Angle - 90))

    # Add noise
    noise = 0.03 * np.random.normal(size=signal.shape)

    y_2d = signal + baseline + noise
    x_2d = [angle_axis, field_axis]
    params_2d = {
        'XAXIS_NAME': 'Angle', 'XAXIS_UNIT': 'degrees',
        'YAXIS_NAME': 'Magnetic Field', 'YAXIS_UNIT': 'G'
    }

    return x_2d, y_2d, params_2d


def demo_baseline_quality_assessment():
    """Demonstrate baseline correction quality assessment."""
    print("4. Baseline correction quality assessment:")
    print("-" * 43)

    # Create test data with known baseline
    x = np.linspace(3350, 3450, 400)
    signal_true = 0.9 * np.exp(-((x - 3400)/15)**2)
    baseline_true = 0.3 + 0.0005 * (x - 3400) + 2e-6 * (x - 3400)**2
    noise = 0.04 * np.random.normal(size=len(x))
    y_with_baseline = signal_true + baseline_true + noise

    print("Quality metrics for different correction orders:")

    orders = [1, 2, 3, 4, 5]
    metrics = {}

    for order in orders:
        try:
            y_corrected, baseline = baseline_polynomial_1d(x, y_with_baseline, order=order)

            # Calculate quality metrics
            # 1. Baseline flatness (std of endpoints)
            endpoints = np.concatenate([y_corrected[:20], y_corrected[-20:]])
            flatness = np.std(endpoints)

            # 2. Signal preservation (correlation with true signal)
            correlation = np.corrcoef(y_corrected, signal_true)[0, 1]

            # 3. Mean baseline level
            baseline_level = np.mean(endpoints)

            # 4. Residual baseline (if we know the true baseline)
            residual_baseline = np.mean(np.abs(y_corrected - signal_true))

            metrics[order] = {
                'flatness': flatness,
                'correlation': correlation,
                'baseline_level': baseline_level,
                'residual_baseline': residual_baseline
            }

            print(f"  Order {order}:")
            print(f"    Baseline flatness (std): {flatness:.4f}")
            print(f"    Signal correlation: {correlation:.4f}")
            print(f"    Mean baseline level: {baseline_level:.4f}")
            print(f"    Residual baseline: {residual_baseline:.4f}")

        except Exception as e:
            print(f"  Order {order}: Error - {e}")

    # Find optimal order
    if metrics:
        # Score based on multiple criteria (lower is better for most)
        scores = {}
        for order, m in metrics.items():
            # Combine metrics (normalize and weight)
            score = (m['flatness'] + abs(m['baseline_level']) +
                    m['residual_baseline'] - m['correlation'])
            scores[order] = score

        best_order = min(scores.keys(), key=lambda k: scores[k])
        print(f"\n  Recommended order: {best_order} (lowest combined score: {scores[best_order]:.4f})")

    print()


def main():
    """Run all baseline correction demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_1d_baseline_correction()
    demo_baseline_with_signal_exclusion()
    demo_2d_baseline_correction()
    demo_baseline_quality_assessment()

    print("=== Baseline Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- baseline_polynomial_1d() provides polynomial baseline correction for 1D data")
    print("- baseline_polynomial_2d() handles 2D baseline correction")
    print("- Signal exclusion improves correction accuracy")
    print("- Different polynomial orders suit different baseline types")
    print("- Quality metrics help select optimal correction parameters")
    print("- Both automated and interactive correction methods available")
    print()
    print("Generated plot files:")
    for plot_file in sorted(output_dir.glob("*.png")):
        if plot_file.name.startswith(('13_', '14_', '15_')):
            print(f"  - {plot_file.name}")


if __name__ == "__main__":
    main()