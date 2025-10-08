#!/usr/bin/env python3
"""
EPyR Tools Demo 01: Data Loading with eprload
=============================================

This script demonstrates the core data loading functionality of EPyR Tools.
The eprload() function is the main entry point for loading Bruker EPR data files.

Functions demonstrated:
- eprload() - Main data loading function
- Format detection and parameter extraction
- 1D and 2D data handling
"""

import sys
import numpy as np
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_eprload_basic():
    """Demonstrate basic eprload functionality."""
    print("=== EPyR Tools eprload Demo ===")
    print()

    # Data directory
    data_dir = Path(__file__).parent.parent.parent / "data"

    # Demo 1: Load 1D CW EPR data (Erbium in CaWO4)
    print("1. Loading 1D CW EPR data:")
    print("-" * 30)

    cw_file = data_dir / "130406SB_CaWO4_Er_CW_5K_20.DSC"

    if cw_file.exists():
        # Load without automatic plotting
        x_field, y_signal, params, filepath = epyr.eprload(str(cw_file), plot_if_possible=False)

        print(f"File loaded: {Path(filepath).name}")
        print(f"Data type: {type(y_signal).__name__}")
        print(f"Field axis shape: {x_field.shape}")
        print(f"Signal shape: {y_signal.shape}")
        print(f"Field range: {x_field.min():.1f} to {x_field.max():.1f} G")
        print(f"Signal range: {y_signal.min():.2e} to {y_signal.max():.2e}")
        print(f"Data is complex: {np.iscomplexobj(y_signal)}")
        print(f"Number of parameters: {len(params)}")
        print()

        # Show key parameters
        print("Key parameters extracted:")
        key_params = ['MWFQ', 'B0VL', 'BWVL', 'MWPW', 'AVGS', 'SPTP', 'Temperature']
        for param in key_params:
            if param in params:
                print(f"  {param}: {params[param]}")
        print()
    else:
        print(f"File not found: {cw_file}")
        print("Creating synthetic 1D data for demonstration...")
        x_field = np.linspace(3350, 3450, 1024)
        y_signal = np.exp(-((x_field - 3400)/15)**2) + 0.05 * np.random.normal(size=len(x_field))
        params = {
            'MWFQ': '9.4 GHz',
            'B0VL': 3400.0,
            'BWVL': 100.0,
            'XPTS': len(x_field),
            'Temperature': '5.0 K'
        }
        print(f"Synthetic data shape: {y_signal.shape}")
        print()


def demo_2d_data_loading():
    """Demonstrate 2D data loading."""
    print("2. Loading 2D EPR data:")
    print("-" * 25)

    data_dir = Path(__file__).parent.parent.parent / "data"

    # Try different 2D files
    possible_2d_files = [
        "2014_03_19_MgO_300K_111_fullrotation33dB.par",
        "Rabi2D_GdCaWO4_6dB_3770G_2.DSC"
    ]

    for filename in possible_2d_files:
        file_path = data_dir / filename
        if file_path.exists():
            print(f"Loading: {filename}")
            try:
                x_axes, y_data, params, filepath = epyr.eprload(str(file_path), plot_if_possible=False)

                print(f"  Data shape: {y_data.shape}")
                print(f"  X-axis type: {type(x_axes)}")

                if isinstance(x_axes, list):
                    print(f"  Axis 1: {len(x_axes[0])} points")
                    print(f"  Axis 2: {len(x_axes[1])} points")
                    print(f"  Axis 1 range: {x_axes[0].min():.1f} to {x_axes[0].max():.1f}")
                    print(f"  Axis 2 range: {x_axes[1].min():.1f} to {x_axes[1].max():.1f}")

                print(f"  Signal is complex: {np.iscomplexobj(y_data)}")
                print(f"  Signal range: {y_data.min():.2e} to {y_data.max():.2e}")
                print()
                break

            except Exception as e:
                print(f"  Error loading {filename}: {e}")
                continue
    else:
        print("No 2D files available, creating synthetic 2D data...")
        # Create synthetic 2D data
        field_axis = np.linspace(3400, 3500, 64)
        angle_axis = np.linspace(0, 180, 32)
        Field, Angle = np.meshgrid(field_axis, angle_axis)
        y_data = np.exp(-((Field - 3450)/20)**2) * np.cos(np.radians(Angle * 2))
        y_data += 0.1 * np.random.normal(size=y_data.shape)
        x_axes = [angle_axis, field_axis]  # Note: order depends on data format
        print(f"  Synthetic 2D data shape: {y_data.shape}")
        print()


def demo_scaling_options():
    """Demonstrate data scaling options."""
    print("3. Data scaling options:")
    print("-" * 25)

    data_dir = Path(__file__).parent.parent.parent / "data"
    cw_file = data_dir / "130406SB_CaWO4_Er_CW_5K_20.DSC"

    if cw_file.exists():
        scaling_options = [
            ("", "No scaling (raw data)"),
            ("n", "Normalized by number of scans"),
            ("P", "Normalized by microwave power"),
            ("nP", "Normalized by scans and power")
        ]

        for scaling, description in scaling_options:
            try:
                x, y, params, _ = epyr.eprload(str(cw_file), scaling=scaling, plot_if_possible=False)
                signal_range = y.max() - y.min()
                print(f"  {scaling if scaling else 'Raw':<3}: {description:<30} Range: {signal_range:.2e}")
            except Exception as e:
                print(f"  {scaling if scaling else 'Raw':<3}: {description:<30} Error: {e}")
        print()
    else:
        print("  CW file not available for scaling demonstration")
        print()


def demo_parameter_extraction():
    """Demonstrate detailed parameter extraction."""
    print("4. Parameter extraction and analysis:")
    print("-" * 35)

    data_dir = Path(__file__).parent.parent.parent / "data"
    cw_file = data_dir / "130406SB_CaWO4_Er_CW_5K_20.DSC"

    if cw_file.exists():
        x, y, params, filepath = epyr.eprload(str(cw_file), plot_if_possible=False)

        print(f"Complete parameter analysis for: {Path(filepath).name}")
        print()

        # Categorize parameters
        field_params = {}
        measurement_params = {}
        acquisition_params = {}
        other_params = {}

        for key, value in params.items():
            if any(field_key in key.upper() for field_key in ['B0', 'FIELD', 'BWVL', 'XMIN', 'XMAX']):
                field_params[key] = value
            elif any(meas_key in key.upper() for meas_key in ['MWFQ', 'MWPW', 'TEMP']):
                measurement_params[key] = value
            elif any(acq_key in key.upper() for acq_key in ['AVGS', 'SPTP', 'XPTS', 'CONV']):
                acquisition_params[key] = value
            else:
                other_params[key] = value

        print("Field parameters:")
        for key, value in field_params.items():
            print(f"  {key:<15}: {value}")
        print()

        print("Measurement parameters:")
        for key, value in measurement_params.items():
            print(f"  {key:<15}: {value}")
        print()

        print("Acquisition parameters:")
        for key, value in acquisition_params.items():
            print(f"  {key:<15}: {value}")
        print()

        print(f"Other parameters: {len(other_params)} items")
        # Show first few other parameters
        for i, (key, value) in enumerate(other_params.items()):
            if i < 5:
                print(f"  {key:<15}: {value}")
            elif i == 5:
                print(f"  ... and {len(other_params)-5} more")
                break
        print()
    else:
        print("  CW file not available for parameter demonstration")
        print()


def demo_error_handling():
    """Demonstrate error handling in eprload."""
    print("5. Error handling:")
    print("-" * 18)

    # Test with non-existent file
    try:
        x, y, params, filepath = epyr.eprload("nonexistent_file.DSC", plot_if_possible=False)
    except Exception as e:
        print(f"  Non-existent file: {type(e).__name__} - {e}")

    # Test with invalid file (if any)
    try:
        x, y, params, filepath = epyr.eprload(__file__, plot_if_possible=False)  # Try to load this Python file
    except Exception as e:
        print(f"  Invalid file format: {type(e).__name__} - {str(e)[:60]}...")

    print()


def main():
    """Run all eprload demonstrations."""
    demo_eprload_basic()
    demo_2d_data_loading()
    demo_scaling_options()
    demo_parameter_extraction()
    demo_error_handling()

    print("=== eprload Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- eprload() returns (x, y, params, filepath)")
    print("- Supports both 1D and 2D data automatically")
    print("- Handles BES3T (.DSC/.DTA) and ESP (.par/.spc) formats")
    print("- Provides comprehensive parameter extraction")
    print("- Offers flexible scaling options")
    print("- Includes robust error handling")


if __name__ == "__main__":
    main()