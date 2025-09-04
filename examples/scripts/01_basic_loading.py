#!/usr/bin/env python3
"""
EPyR Tools - Basic Data Loading Example
=======================================

This script demonstrates how to load EPR data from Bruker files
and perform basic visualization.

Requirements:
- Sample EPR data files in ../data/BES3T/ or ../data/ESP/
- matplotlib for plotting
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

# Add EPyR Tools to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import epyr.eprload as eprload


def load_and_plot_example():
    """Load and plot EPR data from sample files."""
    # Define data directories
    examples_dir = Path(__file__).parent.parent
    bes3t_dir = examples_dir / "data" / "BES3T"
    esp_dir = examples_dir / "data" / "ESP"

    print("EPyR Tools - Basic Data Loading Example")
    print("=" * 40)

    # Look for sample files
    sample_files = []

    # Check BES3T directory
    for dsc_file in bes3t_dir.glob("*.dsc"):
        dta_file = dsc_file.with_suffix(".dta")
        if dta_file.exists():
            sample_files.append(("BES3T", dsc_file))

    # Check ESP directory
    for par_file in esp_dir.glob("*.par"):
        spc_file = par_file.with_suffix(".spc")
        if spc_file.exists():
            sample_files.append(("ESP", par_file))

    if not sample_files:
        print("No sample EPR files found!")
        print(f"Please add sample files to:")
        print(f"  - {bes3t_dir} (BES3T format: .dsc/.dta pairs)")
        print(f"  - {esp_dir} (ESP format: .par/.spc pairs)")
        return

    # Process each sample file
    for file_format, file_path in sample_files:
        print(f"\nLoading {file_format} file: {file_path.name}")

        try:
            # Load EPR data
            x, y, params, filepath = eprload.eprload(
                str(file_path), plot_if_possible=False
            )

            if x is None or y is None:
                print(f"  Failed to load data from {file_path.name}")
                continue

            # Display basic information
            print(f"  Data points: {len(x)}")
            print(f"  Field range: {x.min():.1f} to {x.max():.1f} G")
            print(f"  Signal range: {y.min():.2e} to {y.max():.2e}")

            # Display key parameters
            key_params = {
                "MWFQ": "Microwave Frequency (Hz)",
                "MWPW": "Microwave Power (dB)",
                "AVGS": "Number of Averages",
                "HCF": "Center Field (G)",
                "HSW": "Sweep Width (G)",
                "MF": "Frequency (GHz)",
                "MP": "Power",
            }

            print("  Key Parameters:")
            for param, description in key_params.items():
                if param in params:
                    value = params[param]
                    print(f"    {description}: {value}")

            # Create plot
            plt.figure(figsize=(10, 6))
            plt.plot(x, y, "b-", linewidth=1.5)
            plt.xlabel("Magnetic Field (G)")
            plt.ylabel("EPR Signal (a.u.)")
            plt.title(f"EPR Spectrum: {file_path.stem} ({file_format} format)")
            plt.grid(True, alpha=0.3)

            # Add frequency info to plot if available
            freq = params.get("MWFQ", params.get("MF", None))
            if freq:
                if isinstance(freq, str):
                    freq_ghz = float(freq)
                else:
                    freq_ghz = freq / 1e9 if freq > 1e6 else freq
                plt.text(
                    0.02,
                    0.98,
                    f"Frequency: {freq_ghz:.3f} GHz",
                    transform=plt.gca().transAxes,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                )

            plt.tight_layout()

            # Save plot
            output_file = examples_dir / "scripts" / f"{file_path.stem}_plot.png"
            plt.savefig(output_file, dpi=150, bbox_inches="tight")
            print(f"  Plot saved: {output_file.name}")

            # Show plot (comment out for batch processing)
            # plt.show()

            plt.close()

        except Exception as e:
            print(f"  Error loading {file_path.name}: {e}")

    print(f"\nExample complete! Check the scripts/ directory for generated plots.")


if __name__ == "__main__":
    load_and_plot_example()
