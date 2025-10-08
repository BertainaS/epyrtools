#!/usr/bin/env python3
"""
EPyR Tools Demo 07: FAIR Data Conversion and Validation
=======================================================

This script demonstrates the FAIR (Findable, Accessible, Interoperable, Reusable)
data conversion and validation capabilities of EPyR Tools.

Functions demonstrated:
- convert_bruker_to_fair() - Convert Bruker EPR data to FAIR formats
- validate_fair_dataset() - Validate data for FAIR compliance
- save_to_csv_json() - Export to CSV/JSON formats with metadata
- save_to_hdf5() - Export to HDF5 format
- batch_convert_directory() - Process multiple files
"""

import sys
import json
from pathlib import Path
import numpy as np

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_basic_fair_conversion():
    """Demonstrate basic FAIR data conversion."""
    print("=== EPyR Tools FAIR Demo - Basic Data Conversion ===")
    print()

    print("1. Single file FAIR conversion:")
    print("-" * 31)

    data_dir = Path(__file__).parent.parent.parent / "data"
    output_dir = Path(__file__).parent / "fair_output"
    output_dir.mkdir(exist_ok=True)

    # Look for real EPR data files
    test_files = [
        "130406SB_CaWO4_Er_CW_5K_20.DSC",
        "2024_08_CaWO4171Yb_rabi_6K_6724G_18dB.DSC"
    ]

    converted_files = []

    for filename in test_files:
        file_path = data_dir / filename
        if file_path.exists():
            print(f"Converting: {filename}")
            try:
                # Convert to multiple FAIR formats
                success = epyr.fair.convert_bruker_to_fair(
                    str(file_path),
                    output_dir=str(output_dir),
                    formats=['csv', 'json', 'hdf5'],
                    include_metadata=True,
                    scaling=""
                )

                if success:
                    print(f"  ✓ Conversion successful")
                    converted_files.append(filename)

                    # List generated files
                    base_name = file_path.stem
                    for ext in ['csv', 'json', 'h5']:
                        output_file = output_dir / f"{base_name}.{ext}"
                        if output_file.exists():
                            size_kb = output_file.stat().st_size / 1024
                            print(f"    Generated: {output_file.name} ({size_kb:.1f} KB)")
                else:
                    print(f"  ✗ Conversion failed")

            except Exception as e:
                print(f"  ✗ Error: {e}")

            print()
            break
    else:
        # Create synthetic data if no real files available
        print("No real EPR data found, creating synthetic data for conversion...")
        create_synthetic_fair_demo(output_dir)

    return converted_files


def create_synthetic_fair_demo(output_dir):
    """Create synthetic EPR data and demonstrate FAIR conversion."""
    print("Creating synthetic EPR data for FAIR conversion demonstration:")

    # Generate synthetic 1D EPR spectrum
    field_axis = np.linspace(3300, 3500, 256)
    center = 3400
    width = 15
    signal = np.exp(-((field_axis - center)/width)**2)
    signal += 0.05 * np.random.normal(size=len(field_axis))

    # Create EPR parameters dictionary
    params = {
        'MWFQ': '9.4 GHz',
        'B0VL': center,
        'BWVL': 200.0,
        'XPTS': len(field_axis),
        'AVGS': 10,
        'SPTP': 2.048,
        'MWPW': '20 mW',
        'Temperature': '77 K',
        'Experiment': 'CW EPR',
        'XAXIS_NAME': 'Magnetic Field',
        'XAXIS_UNIT': 'G',
        'Sample': 'Synthetic radical'
    }

    # Save as synthetic Bruker-style data temporarily
    synthetic_file = output_dir / "synthetic_epr_data.txt"
    with open(synthetic_file, 'w') as f:
        f.write("# Synthetic EPR data for FAIR conversion demo\n")
        f.write("# Field(G)\tSignal\n")
        for field, sig in zip(field_axis, signal):
            f.write(f"{field:.3f}\t{sig:.6e}\n")

    print(f"  Created synthetic data: {synthetic_file.name}")

    # Manual FAIR conversion for synthetic data
    fair_data = {
        'x_axis': field_axis.tolist(),
        'y_data': signal.tolist(),
        'parameters': params,
        'metadata': {
            'source_format': 'Synthetic',
            'conversion_software': 'EPyR Tools',
            'conversion_date': '2024-01-01T00:00:00Z',
            'data_dimensions': '1D',
            'fair_principles': {
                'findable': True,
                'accessible': True,
                'interoperable': True,
                'reusable': True
            }
        }
    }

    # Save as JSON
    json_file = output_dir / "synthetic_epr_data.json"
    with open(json_file, 'w') as f:
        json.dump(fair_data, f, indent=2)
    print(f"  Saved FAIR JSON: {json_file.name}")

    # Save as CSV with metadata header
    csv_file = output_dir / "synthetic_epr_data.csv"
    with open(csv_file, 'w') as f:
        f.write("# FAIR EPR Data - Synthetic Spectrum\n")
        f.write("# Generated by EPyR Tools FAIR demo\n")
        f.write(f"# Parameters: {json.dumps(params)}\n")
        f.write("# Field(G),Signal\n")
        for field, sig in zip(field_axis, signal):
            f.write(f"{field:.3f},{sig:.6e}\n")
    print(f"  Saved FAIR CSV: {csv_file.name}")

    print()


def demo_fair_validation():
    """Demonstrate FAIR data validation."""
    print("2. FAIR data validation:")
    print("-" * 24)

    output_dir = Path(__file__).parent / "fair_output"

    # Look for generated FAIR files to validate
    json_files = list(output_dir.glob("*.json"))

    if json_files:
        for json_file in json_files[:2]:  # Validate first 2 files
            print(f"Validating: {json_file.name}")
            try:
                # Load the JSON data
                with open(json_file, 'r') as f:
                    data = json.load(f)

                # Validate using EPyR Tools validation
                result = epyr.fair.validate_fair_dataset(data)

                print(f"  Valid: {result.is_valid}")
                if result.errors:
                    print(f"  Errors: {len(result.errors)}")
                    for error in result.errors[:3]:  # Show first 3 errors
                        print(f"    - {error}")
                if result.warnings:
                    print(f"  Warnings: {len(result.warnings)}")
                    for warning in result.warnings[:3]:  # Show first 3 warnings
                        print(f"    - {warning}")
                if result.info:
                    print(f"  Info: {len(result.info)} items")

            except Exception as e:
                print(f"  Error validating {json_file.name}: {e}")

            print()
    else:
        print("No JSON files found for validation")
        print()


def demo_metadata_processing():
    """Demonstrate metadata processing and parameter mapping."""
    print("3. Metadata processing and parameter mapping:")
    print("-" * 45)

    # Show parameter mapping from Bruker to FAIR format
    try:
        param_map = epyr.fair.BRUKER_PARAM_MAP
        print("Bruker to FAIR parameter mapping (sample):")
        print(f"{'Bruker Parameter':<20} {'FAIR Parameter':<25} {'Description'}")
        print("-" * 70)

        # Show first few mappings
        for i, (bruker_key, fair_info) in enumerate(param_map.items()):
            if i >= 5:  # Show only first 5
                break
            fair_key = fair_info.get('fair_name', bruker_key)
            description = fair_info.get('description', 'No description')
            print(f"{bruker_key:<20} {fair_key:<25} {description[:20]}...")

        print(f"... and {len(param_map)-5} more mappings")

    except Exception as e:
        print(f"Error accessing parameter mapping: {e}")

    print()

    # Demonstrate parameter processing
    print("Example parameter processing:")
    bruker_params = {
        'MWFQ': '9.4 GHz',
        'B0VL': 3400.0,
        'BWVL': 100.0,
        'AVGS': 10,
        'SPTP': 2.048,
        'MWPW': '20 mW',
        'Temperature': '77 K'
    }

    try:
        processed_params = epyr.fair.process_parameters(bruker_params)
        print("Processed parameters:")
        for key, value in processed_params.items():
            print(f"  {key}: {value}")

    except Exception as e:
        print(f"Error processing parameters: {e}")

    print()


def demo_batch_conversion():
    """Demonstrate batch conversion of multiple files."""
    print("4. Batch directory conversion:")
    print("-" * 30)

    data_dir = Path(__file__).parent.parent.parent / "data"
    output_dir = Path(__file__).parent / "fair_batch_output"
    output_dir.mkdir(exist_ok=True)

    print(f"Looking for EPR files in: {data_dir}")

    # Count available files
    epr_files = []
    for pattern in ["*.DSC", "*.dsc", "*.PAR", "*.par", "*.SPC", "*.spc", "*.DTA", "*.dta"]:
        epr_files.extend(data_dir.glob(pattern))

    if epr_files:
        print(f"Found {len(epr_files)} EPR files")

        try:
            # Batch convert (limit to first 3 files for demo)
            files_to_convert = epr_files[:3]
            print(f"Converting first {len(files_to_convert)} files...")

            for i, file_path in enumerate(files_to_convert):
                print(f"  [{i+1}/{len(files_to_convert)}] {file_path.name}")
                try:
                    success = epyr.fair.convert_bruker_to_fair(
                        str(file_path),
                        output_dir=str(output_dir),
                        formats=['csv', 'json'],
                        include_metadata=True
                    )
                    print(f"    {'✓' if success else '✗'} {'Success' if success else 'Failed'}")
                except Exception as e:
                    print(f"    ✗ Error: {e}")

            # Count generated files
            generated_files = list(output_dir.glob("*"))
            print(f"\nGenerated {len(generated_files)} files in {output_dir.name}/")

        except Exception as e:
            print(f"Batch conversion error: {e}")

    else:
        print("No EPR files found for batch conversion")

    print()


def demo_format_comparison():
    """Demonstrate different output formats and their characteristics."""
    print("5. Output format comparison:")
    print("-" * 28)

    output_dir = Path(__file__).parent / "fair_output"

    # Look for files with the same base name but different extensions
    base_names = set()
    for file_path in output_dir.glob("*"):
        if file_path.suffix in ['.csv', '.json', '.h5']:
            base_names.add(file_path.stem)

    if base_names:
        base_name = list(base_names)[0]  # Take first available
        print(f"Comparing formats for: {base_name}")
        print()

        format_info = {}
        for ext in ['.csv', '.json', '.h5']:
            file_path = output_dir / f"{base_name}{ext}"
            if file_path.exists():
                size_bytes = file_path.stat().st_size
                format_info[ext] = {
                    'size_kb': size_bytes / 1024,
                    'size_bytes': size_bytes,
                    'exists': True
                }

                # Try to read and analyze content
                if ext == '.json':
                    try:
                        with open(file_path, 'r') as f:
                            data = json.load(f)
                        format_info[ext]['metadata_fields'] = len(data.get('metadata', {}))
                        format_info[ext]['parameter_count'] = len(data.get('parameters', {}))
                    except:
                        pass

        # Display comparison
        print(f"{'Format':<8} {'Size (KB)':<12} {'Features'}")
        print("-" * 40)

        for ext, info in format_info.items():
            features = []
            if ext == '.csv':
                features.append("Human readable")
                features.append("Spreadsheet compatible")
            elif ext == '.json':
                features.append("Structured metadata")
                features.append("Web compatible")
                if 'metadata_fields' in info:
                    features.append(f"{info['metadata_fields']} metadata fields")
            elif ext == '.h5':
                features.append("Binary format")
                features.append("Large data optimized")

            print(f"{ext:<8} {info['size_kb']:<12.1f} {', '.join(features)}")

    else:
        print("No converted files found for comparison")

    print()


def demo_validation_report():
    """Demonstrate comprehensive validation reporting."""
    print("6. Comprehensive validation reporting:")
    print("-" * 38)

    output_dir = Path(__file__).parent / "fair_output"
    json_files = list(output_dir.glob("*.json"))

    if json_files:
        json_file = json_files[0]  # Use first available file
        print(f"Creating validation report for: {json_file.name}")

        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

            # Create comprehensive validation report
            report = epyr.fair.create_validation_report(data)

            print("\nValidation Report:")
            print("-" * 17)
            print(f"Overall Status: {'PASS' if report['overall_status'] else 'FAIL'}")
            print(f"Validation Date: {report['validation_date']}")
            print(f"Total Checks: {report['total_checks']}")
            print(f"Passed: {report['passed_checks']}")
            print(f"Failed: {report['failed_checks']}")

            if report['errors']:
                print(f"\nErrors ({len(report['errors'])}):")
                for error in report['errors'][:3]:
                    print(f"  - {error}")

            if report['warnings']:
                print(f"\nWarnings ({len(report['warnings'])}):")
                for warning in report['warnings'][:3]:
                    print(f"  - {warning}")

            if report['recommendations']:
                print(f"\nRecommendations ({len(report['recommendations'])}):")
                for rec in report['recommendations'][:3]:
                    print(f"  - {rec}")

        except Exception as e:
            print(f"Error creating validation report: {e}")

    else:
        print("No JSON files available for validation reporting")

    print()


def demo_export_functions():
    """Demonstrate direct export functions."""
    print("7. Direct export function usage:")
    print("-" * 32)

    # Create example data
    x_data = np.linspace(3300, 3500, 100)
    y_data = np.exp(-((x_data - 3400)/20)**2) + 0.02 * np.random.normal(size=len(x_data))

    parameters = {
        'MWFQ': '9.5 GHz',
        'B0VL': 3400.0,
        'Temperature': '300 K',
        'Sample': 'Demo sample'
    }

    output_dir = Path(__file__).parent / "export_demo"
    output_dir.mkdir(exist_ok=True)

    # Demonstrate CSV/JSON export
    try:
        csv_file, json_file = epyr.fair.save_to_csv_json(
            x_data, y_data, parameters,
            output_dir=str(output_dir),
            base_filename="export_demo"
        )
        print(f"CSV export: {Path(csv_file).name}")
        print(f"JSON export: {Path(json_file).name}")

    except Exception as e:
        print(f"CSV/JSON export error: {e}")

    # Demonstrate HDF5 export
    try:
        h5_file = epyr.fair.save_to_hdf5(
            x_data, y_data, parameters,
            output_dir=str(output_dir),
            base_filename="export_demo"
        )
        print(f"HDF5 export: {Path(h5_file).name}")

    except Exception as e:
        print(f"HDF5 export error: {e}")

    print()


def main():
    """Run all FAIR data demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_basic_fair_conversion()
    demo_fair_validation()
    demo_metadata_processing()
    demo_batch_conversion()
    demo_format_comparison()
    demo_validation_report()
    demo_export_functions()

    print("=== FAIR Data Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- convert_bruker_to_fair() provides complete Bruker to FAIR conversion")
    print("- Multiple output formats (CSV, JSON, HDF5) with structured metadata")
    print("- validate_fair_dataset() ensures data meets FAIR principles")
    print("- Comprehensive parameter mapping from Bruker to standardized formats")
    print("- Batch processing capabilities for multiple file conversion")
    print("- Detailed validation reporting with recommendations")
    print("- Direct export functions for custom workflows")
    print()
    print("Generated directories:")
    for dir_path in sorted(output_dir.glob("*/")):
        if dir_path.is_dir() and dir_path.name not in ['.git', '__pycache__']:
            file_count = len(list(dir_path.glob("*")))
            print(f"  - {dir_path.name}/ ({file_count} files)")


if __name__ == "__main__":
    main()