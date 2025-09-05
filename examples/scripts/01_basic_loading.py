#!/usr/bin/env python3
"""
EPyR Tools - Basic Data Loading Example
=======================================

This script demonstrates how to load EPR data from Bruker files
and perform basic visualization.

Requirements:
- Sample EPR data files in ../data/ (BES3T and ESP formats)
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


def load_and_plot_example():
    """Load and plot EPR data from sample files."""
    # Define data directory (consolidated in v0.1.2)
    examples_dir = Path(__file__).parent.parent
    data_dir = examples_dir / "data"

    print("EPyR Tools - Basic Data Loading Example")
    print("=" * 40)

    # Look for sample files in consolidated data directory
    sample_files = []

    # Check for BES3T files (.dsc/.dta pairs)
    for dsc_file in data_dir.glob("*.dsc"):
        dta_file = dsc_file.with_suffix(".dta")
        if dta_file.exists():
            sample_files.append(("BES3T", dsc_file))

    # Also check for uppercase extensions
    for dsc_file in data_dir.glob("*.DSC"):
        dta_file = dsc_file.with_suffix(".DTA")
        if dta_file.exists():
            sample_files.append(("BES3T", dsc_file))

    # Check for ESP files (.par/.spc pairs)
    for par_file in data_dir.glob("*.par"):
        spc_file = par_file.with_suffix(".spc")
        if spc_file.exists():
            sample_files.append(("ESP", par_file))

    if not sample_files:
        print("No sample EPR files found!")
        print(f"Please add sample files to: {data_dir}")
        print("  - BES3T format: .dsc/.dta pairs")
        print("  - ESP format: .par/.spc pairs")
        return

    # Process each sample file
    for file_format, file_path in sample_files:
        print(f"\nLoading {file_format} file: {file_path.name}")

        try:
            # Load EPR data
            x, y, params, filepath = epyr.eprload(
                str(file_path), plot_if_possible=False
            )

            if x is None or y is None:
                print(f"  Failed to load data from {file_path.name}")
                continue

            # Handle both 1D and 2D data
            if isinstance(x, list) and len(x) > 1:
                # 2D data: x is a list of axes
                print(f"  2D Data shape: {y.shape}")
                print(f"  Dimensions: {len(x)} axes")
                if hasattr(x[0], "__len__"):
                    print(f"  Field range: {x[0].min():.1f} to {x[0].max():.1f} G")
                # Use magnitude for complex data
                y_display = np.abs(y) if np.iscomplexobj(y) else y
                print(f"  Signal range: {y_display.min():.2e} to {y_display.max():.2e}")
            else:
                # 1D data: x is a single array
                x_array = x[0] if isinstance(x, list) else x
                print(f"  Data points: {len(x_array)}")
                print(f"  Field range: {x_array.min():.1f} to {x_array.max():.1f} G")
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

            # Create plot based on data type
            plt.figure(figsize=(12, 6))

            if isinstance(x, list) and len(x) > 1:
                # 2D data plotting
                y_plot = np.abs(y) if np.iscomplexobj(y) else y

                if len(y_plot.shape) == 2:
                    # Show 2D map
                    plt.subplot(1, 2, 1)
                    plt.imshow(y_plot, aspect="auto", origin="lower", cmap="viridis")
                    plt.colorbar(label="Signal (a.u.)")
                    plt.xlabel("Field Points")
                    plt.ylabel("Parameter Points")
                    plt.title(f"2D EPR Data: {file_path.stem}")

                    # Show a few representative slices
                    plt.subplot(1, 2, 2)
                    n_slices = min(3, y_plot.shape[0])
                    for i in range(n_slices):
                        if hasattr(x[0], "__len__"):
                            plt.plot(
                                x[0], y_plot[i, :], alpha=0.7, label=f"Slice {i+1}"
                            )
                        else:
                            plt.plot(y_plot[i, :], alpha=0.7, label=f"Slice {i+1}")
                    plt.xlabel("Field (G)" if hasattr(x[0], "__len__") else "Points")
                    plt.ylabel("EPR Signal (a.u.)")
                    plt.title("Representative Spectra")
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                else:
                    # Fallback for unexpected 2D format
                    plt.plot(y_plot.flatten(), "b-", linewidth=1.5)
                    plt.xlabel("Data Points")
                    plt.ylabel("EPR Signal (a.u.)")
                    plt.title(f"EPR Data: {file_path.stem} ({file_format} format)")
                    plt.grid(True, alpha=0.3)
            else:
                # 1D data plotting
                x_array = x[0] if isinstance(x, list) else x
                plt.plot(x_array, y, "b-", linewidth=1.5)
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

                # Add to appropriate axis
                current_ax = plt.gca()
                if isinstance(x, list) and len(x) > 1 and len(y.shape) == 2:
                    # For 2D plots, add to the second subplot
                    current_ax = plt.subplot(1, 2, 2)

                current_ax.text(
                    0.02,
                    0.98,
                    f"Frequency: {freq_ghz:.3f} GHz",
                    transform=current_ax.transAxes,
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
