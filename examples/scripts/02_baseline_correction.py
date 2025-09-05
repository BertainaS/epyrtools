#!/usr/bin/env python3
"""
EPyR Tools - Baseline Correction Example
========================================

This script demonstrates how to apply baseline correction to EPR spectra
using various algorithms available in EPyR Tools.

Requirements:
- Sample EPR data files with baseline issues in ../data/ directory
- matplotlib for plotting

Compatible with EPyR Tools v0.1.2+
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add EPyR Tools to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import epyr
from epyr.baseline import baseline_polynomial


def baseline_correction_example():
    """Demonstrate baseline correction on EPR data."""

    examples_dir = Path(__file__).parent.parent

    print("EPyR Tools - Baseline Correction Example")
    print("=" * 42)

    # Look for sample files in consolidated data directory
    data_dir = examples_dir / "data"
    sample_files = []

    # Look for BES3T and ESP files
    for file_path in data_dir.glob("*"):
        if file_path.suffix.lower() in [".dsc", ".par"]:
            sample_files.append(file_path)

    # Try to find 1D data first (baseline correction is typically for 1D spectra)
    one_d_files = []
    for file_path in sample_files[:]:
        try:
            x_test, y_test, _, _ = epyr.eprload(str(file_path), plot_if_possible=False)
            if x_test is not None and y_test is not None:
                # Check if it's 1D data
                if not isinstance(x_test, list) or len(x_test) == 1:
                    one_d_files.append(file_path)
                elif len(y_test.shape) == 1:
                    one_d_files.append(file_path)
        except Exception:
            continue

    # Prefer 1D files for baseline correction demo
    if one_d_files:
        sample_files = one_d_files

    if not sample_files:
        print("No sample files found. Creating synthetic data for demonstration.")
        create_synthetic_example(examples_dir)
        return

    # Process first available file
    file_path = sample_files[0]
    print(f"Processing file: {file_path.name}")

    try:
        # Load EPR data
        x, y, params, filepath = epyr.eprload(str(file_path), plot_if_possible=False)

        if x is None or y is None:
            print("Failed to load data. Creating synthetic example.")
            create_synthetic_example(examples_dir)
            return

        # Handle 1D vs 2D data for baseline correction
        if isinstance(x, list) and len(x) > 1 and len(y.shape) > 1:
            print(
                f"2D data detected ({y.shape}). "
                f"Using first spectrum for baseline correction demo."
            )
            # Use the first spectrum from 2D data
            x = x[0] if hasattr(x[0], "__len__") else x[0]
            y = y[0, :] if len(y.shape) == 2 else y.flatten()
        elif isinstance(x, list):
            x = x[0] if len(x) > 0 else x

        print(f"Using {len(x)} data points for baseline correction")

        # Apply different baseline corrections
        corrections = [
            ("Original", y, None),
            ("Constant Offset", *baseline_polynomial(y, x_data=x, poly_order=0)),
            ("Linear", *baseline_polynomial(y, x_data=x, poly_order=1)),
            ("Quadratic", *baseline_polynomial(y, x_data=x, poly_order=2)),
        ]

        # Create comparison plot
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        axes = axes.flatten()

        for i, (method, y_data, baseline) in enumerate(corrections):
            ax = axes[i]

            # Plot spectrum
            ax.plot(x, y_data, "b-", linewidth=1.5, label="Spectrum")

            # Plot baseline if available
            if baseline is not None:
                ax.plot(x, baseline, "r--", linewidth=1, alpha=0.7, label="Baseline")

            ax.set_xlabel("Magnetic Field (G)")
            ax.set_ylabel("EPR Signal (a.u.)")
            ax.set_title(f"{method} Correction")
            ax.grid(True, alpha=0.3)
            ax.legend()

            # Add RMS info
            if baseline is not None:
                residuals = y - baseline
                rms = np.sqrt(np.mean(residuals**2))
                ax.text(
                    0.02,
                    0.98,
                    f"RMS: {rms:.2e}",
                    transform=ax.transAxes,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                )

        plt.suptitle(f"Baseline Correction Comparison: {file_path.stem}", fontsize=14)
        plt.tight_layout()

        # Save plot
        output_file = (
            examples_dir / "scripts" / f"{file_path.stem}_baseline_correction.png"
        )
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        print(f"Comparison plot saved: {output_file.name}")

        # Demonstrate advanced baseline correction with exclusion regions
        print("\nDemonstrating baseline correction with signal exclusion...")

        # Find approximate signal regions (peaks)
        signal_threshold = np.std(y) * 2
        peak_indices = np.where(np.abs(y - np.mean(y)) > signal_threshold)[0]

        if len(peak_indices) > 0:
            # Create exclusion regions around peaks
            field_width = (x.max() - x.min()) / len(x) * 50  # ~50 point window
            exclude_regions = []

            for peak_idx in peak_indices[
                :: len(peak_indices) // 3
            ]:  # Sample some peaks
                center_field = x[peak_idx]
                exclude_regions.append(
                    (center_field - field_width / 2, center_field + field_width / 2)
                )

            print(f"Excluding {len(exclude_regions)} signal regions from baseline fit")

            # Apply correction with exclusions
            y_excluded, baseline_excluded = baseline_polynomial(
                y, x_data=x, poly_order=2, exclude_regions=exclude_regions
            )

            # Create exclusion demonstration plot
            plt.figure(figsize=(12, 8))

            plt.subplot(2, 1, 1)
            plt.plot(x, y, "b-", linewidth=1.5, label="Original")
            plt.plot(
                x,
                baseline_excluded,
                "r--",
                linewidth=2,
                label="Baseline (with exclusions)",
            )

            # Highlight excluded regions
            for start, end in exclude_regions:
                plt.axvspan(
                    start,
                    end,
                    alpha=0.2,
                    color="red",
                    label="Excluded regions" if start == exclude_regions[0][0] else "",
                )

            plt.xlabel("Magnetic Field (G)")
            plt.ylabel("EPR Signal (a.u.)")
            plt.title("Baseline Fitting with Signal Exclusion")
            plt.legend()
            plt.grid(True, alpha=0.3)

            plt.subplot(2, 1, 2)
            plt.plot(x, y_excluded, "g-", linewidth=1.5, label="Corrected")
            plt.xlabel("Magnetic Field (G)")
            plt.ylabel("EPR Signal (a.u.)")
            plt.title("Baseline-Corrected Spectrum")
            plt.legend()
            plt.grid(True, alpha=0.3)

            plt.tight_layout()

            # Save exclusion plot
            exclusion_file = (
                examples_dir / "scripts" / f"{file_path.stem}_exclusion_correction.png"
            )
            plt.savefig(exclusion_file, dpi=150, bbox_inches="tight")
            print(f"Exclusion example saved: {exclusion_file.name}")

        plt.close("all")

    except Exception as e:
        print(f"Error processing {file_path.name}: {e}")
        create_synthetic_example(examples_dir)


def create_synthetic_example(examples_dir):
    """Create synthetic EPR data for baseline correction demonstration."""

    print("Creating synthetic EPR data for demonstration...")

    # Generate synthetic EPR spectrum with baseline issues
    x = np.linspace(3200, 3400, 1000)

    # Create EPR signal (Gaussian peak with some structure)
    signal = 100 * np.exp(-((x - 3350) ** 2) / 50) + 30 * np.exp(  # Main peak
        -((x - 3370) ** 2) / 20
    )  # Secondary peak

    # Add polynomial baseline drift
    baseline_true = 0.001 * (x - 3300) ** 2 - 0.1 * (x - 3300) + 50

    # Add noise
    noise = np.random.normal(0, 2, len(x))

    # Combine everything
    y_synthetic = signal + baseline_true + noise

    # Apply baseline corrections
    corrections = [
        ("Original (with baseline)", y_synthetic, None),
        (
            "Constant Correction",
            *baseline_polynomial(y_synthetic, x_data=x, poly_order=0),
        ),
        (
            "Linear Correction",
            *baseline_polynomial(y_synthetic, x_data=x, poly_order=1),
        ),
        (
            "Quadratic Correction",
            *baseline_polynomial(y_synthetic, x_data=x, poly_order=2),
        ),
    ]

    # Create plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()

    for i, (method, y_data, baseline_fit) in enumerate(corrections):
        ax = axes[i]

        ax.plot(x, y_data, "b-", linewidth=1.5, label="Spectrum")

        if baseline_fit is not None:
            ax.plot(
                x, baseline_fit, "r--", linewidth=1, alpha=0.7, label="Fitted baseline"
            )

        # Show true baseline for original
        if i == 0:
            ax.plot(
                x, baseline_true, "k:", linewidth=1, alpha=0.7, label="True baseline"
            )
            ax.plot(x, signal, "g:", linewidth=1, alpha=0.7, label="True signal")

        ax.set_xlabel("Magnetic Field (G)")
        ax.set_ylabel("EPR Signal (a.u.)")
        ax.set_title(method)
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.suptitle("Baseline Correction on Synthetic EPR Data", fontsize=14)
    plt.tight_layout()

    # Save synthetic example
    output_file = examples_dir / "scripts" / "synthetic_baseline_correction.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Synthetic example saved: {output_file.name}")

    plt.close()


if __name__ == "__main__":
    baseline_correction_example()
