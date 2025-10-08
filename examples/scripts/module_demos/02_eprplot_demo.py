#!/usr/bin/env python3
"""
EPyR Tools Demo 02: EPR Plotting with eprplot
==============================================

This script demonstrates the specialized EPR plotting capabilities.
The eprplot module provides functions optimized for EPR data visualization.

Functions demonstrated:
- plot_1d() - 1D EPR spectrum plotting
- plot_2d_map() - 2D color map visualization
- plot_2d_waterfall() - Waterfall plots for 2D data
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def demo_1d_plotting():
    """Demonstrate 1D EPR plotting capabilities."""
    print("=== EPyR Tools eprplot Demo - 1D Plotting ===")
    print()

    data_dir = Path(__file__).parent.parent.parent / "data"

    # Load real 1D data if available
    cw_file = data_dir / "130406SB_CaWO4_Er_CW_5K_20.DSC"

    if cw_file.exists():
        print("1. Real 1D CW EPR data plotting:")
        print("-" * 35)

        x, y, params, filepath = epyr.eprload(str(cw_file), plot_if_possible=False)

        # Basic 1D plot
        print(f"Plotting: {Path(filepath).name}")
        fig, ax = epyr.eprplot.plot_1d(x, y, params, title="Erbium in CaWO4 - CW EPR at 5K")
        plt.savefig(Path(__file__).parent / "01_cw_epr_plot.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 01_cw_epr_plot.png")
        print()

    # Create synthetic complex data for demonstration
    print("2. Complex data plotting (real and imaginary):")
    print("-" * 45)

    # Generate synthetic complex EPR data
    x_field = np.linspace(3400, 3500, 512)
    # Absorption component (real)
    absorption = np.exp(-((x_field - 3450)/10)**2)
    # Dispersion component (imaginary)
    dispersion = 0.3 * (x_field - 3450) / 10 * np.exp(-((x_field - 3450)/10)**2)
    y_complex = absorption + 1j * dispersion

    # Add noise
    y_complex += 0.02 * (np.random.normal(size=len(x_field)) +
                        1j * np.random.normal(size=len(x_field)))

    params_complex = {
        'XMIN': x_field[0], 'XMAX': x_field[-1],
        'XPTS': len(x_field), 'QuadratureDetection': 'Yes',
        'XAXIS_NAME': 'Magnetic Field', 'XAXIS_UNIT': 'G'
    }

    # Plot complex data - eprplot automatically handles real/imaginary
    fig, ax = epyr.eprplot.plot_1d(x_field, y_complex, params_complex,
                                  title="Synthetic Complex EPR Data")
    plt.savefig(Path(__file__).parent / "02_complex_epr_plot.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Complex data automatically plotted with real and imaginary components")
    print("  Saved as: 02_complex_epr_plot.png")
    print()


def demo_2d_plotting():
    """Demonstrate 2D EPR plotting capabilities."""
    print("3. 2D EPR data plotting:")
    print("-" * 25)

    data_dir = Path(__file__).parent.parent.parent / "data"

    # Try to load real 2D data
    angular_file = data_dir / "2014_03_19_MgO_300K_111_fullrotation33dB.par"

    if angular_file.exists():
        print("Loading real 2D angular study data...")
        try:
            x_2d, y_2d, params_2d, filepath = epyr.eprload(str(angular_file), plot_if_possible=False)
            print(f"Loaded: {Path(filepath).name}")
            print(f"Data shape: {y_2d.shape}")
        except Exception as e:
            print(f"Error loading real data: {e}")
            x_2d, y_2d, params_2d = create_synthetic_2d_data()
    else:
        print("Real 2D data not available, creating synthetic data...")
        x_2d, y_2d, params_2d = create_synthetic_2d_data()

    # 2D Color Map
    print("\n3a. 2D Color Map:")
    fig, ax = epyr.eprplot.plot_2d_map(x_2d, y_2d, params_2d,
                                      title="2D EPR Data - Color Map",
                                      cmap="RdBu_r")
    plt.savefig(Path(__file__).parent / "03_2d_colormap.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 03_2d_colormap.png")

    # 2D Waterfall Plot
    print("\n3b. 2D Waterfall Plot:")
    fig, ax = epyr.eprplot.plot_2d_waterfall(x_2d, y_2d, params_2d,
                                            title="2D EPR Data - Waterfall",
                                            max_traces=15)
    plt.savefig(Path(__file__).parent / "04_2d_waterfall.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 04_2d_waterfall.png")
    print()


def create_synthetic_2d_data():
    """Create synthetic 2D EPR data for demonstration."""
    # Create synthetic angular study data
    field_axis = np.linspace(3300, 3500, 80)
    angle_axis = np.linspace(0, 180, 36)

    # Create 2D data with angular dependence
    Field, Angle = np.meshgrid(field_axis, angle_axis)

    # Simulate angular-dependent g-factor effects
    g_parallel = 2.0
    g_perp = 2.2
    theta_rad = np.radians(Angle)
    g_eff = np.sqrt(g_parallel**2 * np.cos(theta_rad)**2 + g_perp**2 * np.sin(theta_rad)**2)

    # Resonance field
    h_planck = 6.626e-34
    mu_bohr = 9.274e-24
    freq = 9.5e9  # X-band frequency
    B_res = h_planck * freq / (mu_bohr * g_eff) * 1e3  # Convert to mT, then to G
    B_res *= 10  # mT to G

    # Create EPR signal
    y_2d = np.exp(-((Field - B_res)/15)**2)

    # Add noise
    y_2d += 0.05 * np.random.normal(size=y_2d.shape)

    x_2d = [angle_axis, field_axis]
    params_2d = {
        'XAXIS_NAME': 'Angle', 'XAXIS_UNIT': 'degrees',
        'YAXIS_NAME': 'Magnetic Field', 'YAXIS_UNIT': 'G',
        'Temperature': '300 K', 'Experiment': 'Angular Study'
    }

    return x_2d, y_2d, params_2d


def demo_plot_customization():
    """Demonstrate plot customization options."""
    print("4. Plot customization:")
    print("-" * 22)

    # Create test data
    x = np.linspace(3400, 3500, 200)
    y = np.exp(-((x - 3450)/8)**2) + 0.03 * np.random.normal(size=len(x))
    params = {
        'XAXIS_NAME': 'Magnetic Field', 'XAXIS_UNIT': 'G',
        'Temperature': '77 K', 'MicrowaveFrequency': '9.4 GHz'
    }

    # Custom colormap for 2D data
    x_2d, y_2d, params_2d = create_synthetic_2d_data()

    print("Different colormap options:")
    colormaps = ['viridis', 'plasma', 'RdBu_r', 'coolwarm']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, cmap in enumerate(colormaps):
        epyr.eprplot.plot_2d_map(x_2d, y_2d, params_2d,
                                title=f"Colormap: {cmap}",
                                cmap=cmap, ax=axes[i])

    plt.tight_layout()
    plt.savefig(Path(__file__).parent / "05_colormap_comparison.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("  Saved as: 05_colormap_comparison.png")
    print()


def demo_advanced_plotting():
    """Demonstrate advanced plotting features."""
    print("5. Advanced plotting features:")
    print("-" * 31)

    # Load time-domain data if available
    data_dir = Path(__file__).parent.parent.parent / "data"

    # Look for time-domain files
    time_files = [
        "2024_08_CaWO4171Yb_rabi_6K_6724G_18dB.DSC",
        "Rabi2D_GdCaWO4_6dB_3770G_2.DSC"
    ]

    for filename in time_files:
        file_path = data_dir / filename
        if file_path.exists():
            print(f"Loading time-domain data: {filename}")
            try:
                x_time, y_time, params_time, filepath = epyr.eprload(str(file_path), plot_if_possible=False)

                # Plot time-domain data
                fig, ax = epyr.eprplot.plot_1d(x_time, y_time, params_time,
                                              title=f"Time-Domain EPR: {Path(filename).stem}")
                plt.savefig(Path(__file__).parent / f"06_time_domain_{Path(filename).stem}.png",
                           dpi=150, bbox_inches='tight')
                plt.show()
                print(f"  Saved as: 06_time_domain_{Path(filename).stem}.png")
                break
            except Exception as e:
                print(f"  Error loading {filename}: {e}")
                continue
    else:
        print("No time-domain data available, creating synthetic data...")
        # Create synthetic time-domain data (Rabi oscillations)
        t = np.linspace(0, 1000, 200)  # Time in ns
        rabi_freq = 0.008  # MHz
        T2 = 300  # ns
        y_rabi = np.exp(-t/T2) * np.cos(2 * np.pi * rabi_freq * t)
        y_rabi += 0.05 * np.random.normal(size=len(t))

        params_rabi = {
            'XAXIS_NAME': 'Time', 'XAXIS_UNIT': 'ns',
            'Experiment': 'Rabi Oscillations', 'Temperature': '6 K'
        }

        fig, ax = epyr.eprplot.plot_1d(t, y_rabi, params_rabi,
                                      title="Synthetic Rabi Oscillations")
        plt.savefig(Path(__file__).parent / "06_synthetic_rabi.png", dpi=150, bbox_inches='tight')
        plt.show()
        print("  Saved as: 06_synthetic_rabi.png")

    print()


def main():
    """Run all eprplot demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    demo_1d_plotting()
    demo_2d_plotting()
    demo_plot_customization()
    demo_advanced_plotting()

    print("=== eprplot Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- plot_1d() handles both real and complex data automatically")
    print("- plot_2d_map() creates color map visualizations with proper axes")
    print("- plot_2d_waterfall() provides alternative 2D visualization")
    print("- All functions accept parameter dictionaries for automatic labeling")
    print("- Customizable colormaps and styling options available")
    print("- Functions return matplotlib figure and axes objects for further customization")
    print()
    print("Generated plot files:")
    for plot_file in sorted(output_dir.glob("*.png")):
        if plot_file.name.startswith(('01_', '02_', '03_', '04_', '05_', '06_')):
            print(f"  - {plot_file.name}")


if __name__ == "__main__":
    main()