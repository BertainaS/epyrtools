#!/usr/bin/env python3
"""
EPyR Tools - FAIR Data Conversion Example
=========================================

This script demonstrates how to use the FAIR module in EPyR Tools to convert
Bruker EPR files into more accessible and interoperable formats (CSV, JSON, HDF5).

FAIR principles: Findable, Accessible, Interoperable, and Reusable data.

Requirements:
- Sample EPR data files in ../data/ directory
- h5py for HDF5 export
- pandas for enhanced CSV handling

Compatible with EPyR Tools v0.1.2+
"""

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add EPyR Tools to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import epyr
from epyr.fair import convert_bruker_to_fair


def fair_conversion_example():
    """Demonstrate FAIR data conversion for EPR files."""

    examples_dir = Path(__file__).parent.parent
    data_dir = examples_dir / "data"
    output_dir = examples_dir / "data" / "fair_converted"
    output_dir.mkdir(exist_ok=True)

    print("EPyR Tools - FAIR Data Conversion Example")
    print("=" * 43)
    print("Converting Bruker EPR files to FAIR formats (CSV, JSON, HDF5)")
    print()

    # Find EPR files
    epr_files = []
    for ext in ["*.dsc", "*.DSC", "*.par"]:
        epr_files.extend(data_dir.glob(ext))

    if not epr_files:
        print("❌ No EPR files found!")
        print(f"Please add EPR files to: {data_dir}")
        return

    print(f"Found {len(epr_files)} EPR files for conversion:")
    for i, file_path in enumerate(epr_files, 1):
        file_size = file_path.stat().st_size / 1024
        print(f"  {i}. {file_path.name} ({file_size:.1f} KB)")
    print()

    # Process each file
    for file_path in epr_files:
        print(f"🔄 Processing: {file_path.name}")

        try:
            # Load original data
            x, y, params, filepath = epyr.eprload(
                str(file_path), plot_if_possible=False
            )

            if x is None or y is None:
                print(f"  ❌ Failed to load {file_path.name}")
                continue

            # Determine data type
            is_2d = isinstance(x, list) and len(x) > 1 and len(y.shape) > 1
            is_complex = np.iscomplexobj(y)

            print(f"  📊 Data type: {'2D' if is_2d else '1D'}")
            print(f"  🔢 Complex data: {'Yes' if is_complex else 'No'}")

            # Use the FAIR conversion module
            base_name = file_path.stem
            output_base = output_dir / base_name

            # Call the FAIR conversion function
            convert_bruker_to_fair(str(file_path), output_dir=str(output_dir))

            # Verify converted files exist
            converted_files = []
            for ext in [".csv", ".json", ".h5"]:
                converted_file = output_dir / f"{base_name}{ext}"
                if converted_file.exists():
                    converted_files.append(converted_file)

            print(f"  ✅ Converted to {len(converted_files)} formats:")
            for conv_file in converted_files:
                file_size = conv_file.stat().st_size / 1024
                print(f"    - {conv_file.name} ({file_size:.1f} KB)")

            # Demonstrate reading converted data
            demonstrate_converted_data(output_dir, base_name, is_2d, is_complex)

        except Exception as e:
            print(f"  ❌ Error processing {file_path.name}: {e}")

        print()

    # Create comparison visualization
    create_format_comparison(output_dir)

    print("🎉 FAIR conversion complete!")
    print(f"📁 All converted files saved to: {output_dir}")
    print("\n💡 Benefits of FAIR formats:")
    print("  - CSV: Universal compatibility, easy to read in Excel/R/Python")
    print("  - JSON: Web-friendly, human-readable metadata")
    print("  - HDF5: Efficient for large datasets, preserves data types")


def demonstrate_converted_data(output_dir, base_name, is_2d, is_complex):
    """Show how to read and use the converted FAIR formats."""

    # Check CSV
    csv_file = output_dir / f"{base_name}.csv"
    if csv_file.exists():
        print(f"    📄 CSV Preview:")
        try:
            df = pd.read_csv(csv_file, comment="#")
            print(f"      Shape: {df.shape}")
            cols = list(df.columns[:3])
            suffix = "..." if len(df.columns) > 3 else ""
            print(f"      Columns: {cols}{suffix}")
        except Exception as e:
            print(f"      Error reading CSV: {e}")

    # Check JSON metadata
    json_file = output_dir / f"{base_name}.json"
    if json_file.exists():
        print(f"    📋 JSON Metadata Preview:")
        try:
            with open(json_file, "r") as f:
                metadata = json.load(f)

            keys = list(metadata.keys())[:5]
            suffix = "..." if len(metadata) > 5 else ""
            print(f"      Keys: {keys}{suffix}")

            # Show some interesting parameters
            if "measurement_info" in metadata:
                info = metadata["measurement_info"]
                for key in ["data_type", "dimensions", "total_points"]:
                    if key in info:
                        print(f"      {key}: {info[key]}")

        except Exception as e:
            print(f"      Error reading JSON: {e}")

    # Check HDF5
    h5_file = output_dir / f"{base_name}.h5"
    if h5_file.exists():
        print(f"    🗄️  HDF5 Structure:")
        try:
            import h5py

            with h5py.File(h5_file, "r") as f:

                def print_structure(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        print(f"      Dataset: {name} {obj.shape} {obj.dtype}")
                    else:
                        print(f"      Group: {name}")

                f.visititems(print_structure)

        except ImportError:
            print(f"      h5py not available for preview")
        except Exception as e:
            print(f"      Error reading HDF5: {e}")


def create_format_comparison(output_dir):
    """Create a visualization comparing the different FAIR formats."""

    print("📊 Creating format comparison visualization...")

    # Find all converted files
    csv_files = list(output_dir.glob("*.csv"))
    json_files = list(output_dir.glob("*.json"))
    h5_files = list(output_dir.glob("*.h5"))

    if not (csv_files or json_files or h5_files):
        print("  No converted files found for comparison")
        return

    # Collect file size information
    format_data = {
        "Format": [],
        "File Count": [],
        "Total Size (KB)": [],
        "Avg Size (KB)": [],
    }

    for fmt, files in [("CSV", csv_files), ("JSON", json_files), ("HDF5", h5_files)]:
        if files:
            sizes = [f.stat().st_size / 1024 for f in files]
            format_data["Format"].append(fmt)
            format_data["File Count"].append(len(files))
            format_data["Total Size (KB)"].append(sum(sizes))
            format_data["Avg Size (KB)"].append(np.mean(sizes))

    # Create comparison plot
    if format_data["Format"]:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # File count comparison
        ax1.bar(
            format_data["Format"],
            format_data["File Count"],
            color=["#2E86AB", "#A23B72", "#F18F01"],
        )
        ax1.set_title("Files Created by Format")
        ax1.set_ylabel("Number of Files")

        # Size comparison
        ax2.bar(
            format_data["Format"],
            format_data["Avg Size (KB)"],
            color=["#2E86AB", "#A23B72", "#F18F01"],
        )
        ax2.set_title("Average File Size by Format")
        ax2.set_ylabel("Size (KB)")

        plt.tight_layout()

        # Save comparison plot
        comparison_file = output_dir / "fair_format_comparison.png"
        plt.savefig(comparison_file, dpi=150, bbox_inches="tight")
        plt.close()

        print(f"  ✅ Comparison plot saved: {comparison_file.name}")

        # Print summary table
        print("\n📈 Format Comparison Summary:")
        print("  " + "-" * 50)
        print(f"  {'Format':<8} {'Files':<6} {'Total KB':<10} {'Avg KB':<8}")
        print("  " + "-" * 50)
        for i, fmt in enumerate(format_data["Format"]):
            print(
                f"  {fmt:<8} {format_data['File Count'][i]:<6} "
                f"{format_data['Total Size (KB)'][i]:<10.1f} "
                f"{format_data['Avg Size (KB)'][i]:<8.1f}"
            )
        print("  " + "-" * 50)


def demonstrate_fair_benefits():
    """Show the benefits of FAIR data formats."""

    print("\n🌟 FAIR Data Benefits Demonstration:")
    print("=" * 35)

    benefits = {
        "Findable": [
            "Standardized metadata in JSON format",
            "Consistent file naming conventions",
            "Complete parameter documentation",
        ],
        "Accessible": [
            "CSV format readable by any spreadsheet software",
            "JSON metadata human-readable",
            "No proprietary software required",
        ],
        "Interoperable": [
            "Standard formats work across platforms",
            "Easy import into R, Python, MATLAB",
            "Web-compatible JSON for online tools",
        ],
        "Reusable": [
            "Complete metadata preservation",
            "Clear data structure documentation",
            "Version information included",
        ],
    }

    for principle, items in benefits.items():
        print(f"\n📌 {principle}:")
        for item in items:
            print(f"  ✓ {item}")

    print(f"\n🔧 Usage Examples:")
    print(f"  Python: df = pd.read_csv('spectrum.csv')")
    print(f"  R: data <- read.csv('spectrum.csv')")
    print(f"  Excel: File → Open → spectrum.csv")
    print(f"  Web: fetch('spectrum.json').then(r => r.json())")


if __name__ == "__main__":
    fair_conversion_example()
    demonstrate_fair_benefits()
