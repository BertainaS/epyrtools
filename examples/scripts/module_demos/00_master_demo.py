#!/usr/bin/env python3
"""
EPyR Tools Master Demo - Complete Package Demonstration
=======================================================

This master script runs all individual module demonstrations to showcase
the complete functionality of EPyR Tools.

This comprehensive demonstration covers:
01. Data Loading (eprload) - Loading Bruker EPR files
02. EPR Plotting (eprplot) - Specialized EPR visualization
03. Lineshape Analysis - Mathematical lineshape functions and fitting
04. Baseline Correction - 1D and 2D baseline correction algorithms
05. Signal Processing - FFT analysis and apodization windows
06. Physics Constants - Physical constants and unit conversions
07. FAIR Data - Data conversion and validation for FAIR compliance
08. Performance - Memory optimization and caching systems
09. Configuration - Centralized settings management
10. CLI Interface - Command-line tools and workflows

Usage:
    python 00_master_demo.py [options]

Options:
    --all           Run all demonstrations (default)
    --modules N,M   Run specific module numbers (e.g., --modules 1,2,5)
    --skip N,M      Skip specific module numbers
    --no-plots      Skip plot generation (faster execution)
    --verbose       Enable verbose output
    --summary-only  Show only summary information
"""

import sys
import time
import argparse
import importlib.util
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def load_and_run_demo(demo_number, demo_name, demo_file, verbose=False, no_plots=False):
    """Load and run a specific demonstration module."""
    print(f"\n{'='*60}")
    print(f"Running Demo {demo_number:02d}: {demo_name}")
    print(f"{'='*60}")

    if verbose:
        print(f"Loading: {demo_file}")

    try:
        # Load the demo module
        spec = importlib.util.spec_from_file_location(f"demo_{demo_number:02d}", demo_file)
        if spec is None or spec.loader is None:
            print(f"✗ Error: Could not load {demo_file}")
            return False

        demo_module = importlib.util.module_from_spec(spec)

        # Modify the module to skip plots if requested
        if no_plots:
            # This is a simple approach - for a more robust solution,
            # each demo could check for a NO_PLOTS environment variable
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend

        # Execute the module
        start_time = time.time()
        spec.loader.exec_module(demo_module)

        # Run the main function if it exists
        if hasattr(demo_module, 'main'):
            demo_module.main()

        execution_time = time.time() - start_time

        print(f"\n✓ Demo {demo_number:02d} completed successfully in {execution_time:.1f} seconds")
        return True

    except Exception as e:
        print(f"✗ Error running demo {demo_number:02d}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return False


def show_demo_summary():
    """Show summary of all available demonstrations."""
    demos = [
        (1, "Data Loading (eprload)", "Loading Bruker EPR files and parameter extraction"),
        (2, "EPR Plotting (eprplot)", "Specialized visualization for 1D and 2D EPR data"),
        (3, "Lineshape Analysis", "Mathematical lineshape functions and signal fitting"),
        (4, "Baseline Correction", "1D and 2D baseline correction algorithms"),
        (5, "Signal Processing", "FFT analysis and apodization windows"),
        (6, "Physics Constants", "Physical constants and unit conversions"),
        (7, "FAIR Data Conversion", "Data conversion and validation for FAIR compliance"),
        (8, "Performance Optimization", "Memory optimization and caching systems"),
        (9, "Configuration Management", "Centralized settings and preferences"),
        (10, "CLI Interface", "Command-line tools and workflows")
    ]

    print("EPyR Tools - Complete Package Demonstration")
    print("=" * 50)
    print()
    print("Available demonstrations:")
    print("-" * 25)

    for num, name, description in demos:
        print(f"{num:2d}. {name}")
        print(f"    {description}")
        print()

    return demos


def check_dependencies():
    """Check if required dependencies are available."""
    print("Checking dependencies...")

    required_modules = [
        ('numpy', 'NumPy'),
        ('matplotlib', 'Matplotlib'),
        ('scipy', 'SciPy'),
        ('pathlib', 'Pathlib')
    ]

    optional_modules = [
        ('pandas', 'Pandas (for advanced data handling)'),
        ('h5py', 'HDF5 support'),
        ('psutil', 'Memory monitoring'),
        ('tkinter', 'GUI applications')
    ]

    missing_required = []
    missing_optional = []

    for module_name, description in required_modules:
        try:
            __import__(module_name)
            print(f"  ✓ {description}")
        except ImportError:
            print(f"  ✗ {description} - REQUIRED")
            missing_required.append(module_name)

    for module_name, description in optional_modules:
        try:
            __import__(module_name)
            print(f"  ✓ {description}")
        except ImportError:
            print(f"  ○ {description} - Optional")
            missing_optional.append(module_name)

    if missing_required:
        print(f"\nERROR: Missing required dependencies: {', '.join(missing_required)}")
        print("Please install missing dependencies before running demos.")
        return False

    if missing_optional:
        print(f"\nNote: Optional dependencies not available: {', '.join(missing_optional)}")
        print("Some features may be limited.")

    print("\nDependency check completed.")
    return True


def main():
    """Main function to run the master demonstration."""
    parser = argparse.ArgumentParser(
        description="EPyR Tools Master Demonstration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 00_master_demo.py                    # Run all demos
  python 00_master_demo.py --modules 1,2,3    # Run specific demos
  python 00_master_demo.py --skip 5,8         # Skip specific demos
  python 00_master_demo.py --no-plots         # Skip plot generation
  python 00_master_demo.py --summary-only     # Show summary only
        """
    )

    parser.add_argument('--all', action='store_true', default=True,
                       help='Run all demonstrations (default)')
    parser.add_argument('--modules', type=str,
                       help='Run specific module numbers (comma-separated)')
    parser.add_argument('--skip', type=str,
                       help='Skip specific module numbers (comma-separated)')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip plot generation for faster execution')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--summary-only', action='store_true',
                       help='Show only summary information')

    args = parser.parse_args()

    # Show summary
    demos = show_demo_summary()

    if args.summary_only:
        return

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    print()

    # Determine which demos to run
    if args.modules:
        try:
            demo_numbers = [int(x.strip()) for x in args.modules.split(',')]
        except ValueError:
            print("Error: Invalid module numbers. Use comma-separated integers.")
            sys.exit(1)
    else:
        demo_numbers = list(range(1, 11))  # All demos 1-10

    # Remove skipped demos
    if args.skip:
        try:
            skip_numbers = [int(x.strip()) for x in args.skip.split(',')]
            demo_numbers = [n for n in demo_numbers if n not in skip_numbers]
        except ValueError:
            print("Error: Invalid skip numbers. Use comma-separated integers.")
            sys.exit(1)

    # Validate demo numbers
    valid_numbers = [d[0] for d in demos]
    invalid_numbers = [n for n in demo_numbers if n not in valid_numbers]
    if invalid_numbers:
        print(f"Error: Invalid demo numbers: {invalid_numbers}")
        print(f"Valid numbers are: {valid_numbers}")
        sys.exit(1)

    print(f"Running demonstrations: {demo_numbers}")
    if args.no_plots:
        print("Plot generation disabled for faster execution")
    print()

    # Run demonstrations
    demo_dir = Path(__file__).parent
    successful_demos = []
    failed_demos = []
    total_start_time = time.time()

    for demo_num in demo_numbers:
        # Find corresponding demo info
        demo_info = next((d for d in demos if d[0] == demo_num), None)
        if not demo_info:
            continue

        _, demo_name, _ = demo_info
        demo_file = demo_dir / f"{demo_num:02d}_{demo_name.lower().replace(' ', '').replace('(', '').replace(')', '')}_demo.py"

        # Handle special naming cases
        if not demo_file.exists():
            # Try alternative naming
            if demo_num == 1:
                demo_file = demo_dir / "01_eprload_demo.py"
            elif demo_num == 2:
                demo_file = demo_dir / "02_eprplot_demo.py"
            elif demo_num == 3:
                demo_file = demo_dir / "03_lineshapes_demo.py"
            elif demo_num == 4:
                demo_file = demo_dir / "04_baseline_demo.py"
            elif demo_num == 5:
                demo_file = demo_dir / "05_signalprocessing_demo.py"
            elif demo_num == 6:
                demo_file = demo_dir / "06_physics_demo.py"
            elif demo_num == 7:
                demo_file = demo_dir / "07_fair_demo.py"
            elif demo_num == 8:
                demo_file = demo_dir / "08_performance_demo.py"
            elif demo_num == 9:
                demo_file = demo_dir / "09_config_demo.py"
            elif demo_num == 10:
                demo_file = demo_dir / "10_cli_demo.py"

        if not demo_file.exists():
            print(f"✗ Demo file not found: {demo_file}")
            failed_demos.append((demo_num, demo_name, "File not found"))
            continue

        # Run the demonstration
        success = load_and_run_demo(demo_num, demo_name, demo_file,
                                   verbose=args.verbose, no_plots=args.no_plots)

        if success:
            successful_demos.append((demo_num, demo_name))
        else:
            failed_demos.append((demo_num, demo_name, "Execution failed"))

        # Short pause between demos
        time.sleep(1)

    # Final summary
    total_time = time.time() - total_start_time

    print(f"\n{'='*60}")
    print("MASTER DEMONSTRATION COMPLETE")
    print(f"{'='*60}")
    print(f"Total execution time: {total_time:.1f} seconds")
    print(f"Successful demos: {len(successful_demos)}")
    print(f"Failed demos: {len(failed_demos)}")
    print()

    if successful_demos:
        print("✓ Successful demonstrations:")
        for demo_num, demo_name in successful_demos:
            print(f"  {demo_num:2d}. {demo_name}")
        print()

    if failed_demos:
        print("✗ Failed demonstrations:")
        for demo_num, demo_name, reason in failed_demos:
            print(f"  {demo_num:2d}. {demo_name} - {reason}")
        print()

    # Generated files summary
    output_files = []
    for pattern in ["*.png", "*.csv", "*.json", "*.txt", "*.h5"]:
        output_files.extend(demo_dir.glob(pattern))

    if output_files:
        print(f"Generated files: {len(output_files)}")
        file_types = {}
        for file_path in output_files:
            ext = file_path.suffix
            file_types[ext] = file_types.get(ext, 0) + 1

        for ext, count in sorted(file_types.items()):
            print(f"  {ext}: {count} files")
        print()

    print("EPyR Tools demonstration completed!")
    print("For individual demonstrations, run the specific demo scripts directly.")
    print("For CLI usage, see the generated cli_quick_reference.txt file.")

    if failed_demos:
        sys.exit(1)  # Exit with error if any demos failed


if __name__ == "__main__":
    main()