#!/usr/bin/env python3
"""
EPyR Tools Demo 10: Command Line Interface (CLI)
================================================

This script demonstrates the comprehensive command-line interface (CLI) of EPyR Tools.
The CLI provides professional tools for all EPR workflows without requiring Python coding.

CLI Commands demonstrated:
- epyr-convert - Convert Bruker files to FAIR formats
- epyr-baseline - Apply baseline correction
- epyr-batch-convert - Process multiple files efficiently
- epyr-config - Manage configuration settings
- epyr-info - Display system information
- epyr-plot - Interactive data visualization with measurement tools
- epyr-validate - Validate data files and FAIR compliance
- epyr-isotopes - Launch isotope database GUI
"""

import sys
import subprocess
import shutil
from pathlib import Path

# Add the project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import epyr


def check_cli_availability():
    """Check if CLI commands are available in the system."""
    print("=== EPyR Tools CLI Demo - Command Availability ===")
    print()

    print("1. Checking CLI command availability:")
    print("-" * 37)

    cli_commands = [
        'epyr-convert',
        'epyr-baseline',
        'epyr-batch-convert',
        'epyr-config',
        'epyr-info',
        'epyr-plot',
        'epyr-validate',
        'epyr-isotopes',
        'epyr'  # Main entry point
    ]

    available_commands = []

    for command in cli_commands:
        try:
            # Check if command is available
            result = subprocess.run([command, '--help'],
                                  capture_output=True,
                                  text=True,
                                  timeout=10)
            if result.returncode == 0:
                available_commands.append(command)
                print(f"  ✓ {command} - Available")
            else:
                print(f"  ✗ {command} - Not available (exit code: {result.returncode})")

        except subprocess.TimeoutExpired:
            print(f"  ✗ {command} - Timeout")
        except FileNotFoundError:
            print(f"  ✗ {command} - Command not found")
        except Exception as e:
            print(f"  ✗ {command} - Error: {e}")

    print(f"\nAvailable commands: {len(available_commands)}/{len(cli_commands)}")

    if not available_commands:
        print("\nNote: CLI commands may not be installed. Run 'pip install -e .' in the project root.")

    return available_commands


def demo_convert_command():
    """Demonstrate the epyr-convert command."""
    print("\n2. epyr-convert - Data conversion command:")
    print("-" * 42)

    data_dir = Path(__file__).parent.parent.parent / "data"
    output_dir = Path(__file__).parent / "cli_demo_output"
    output_dir.mkdir(exist_ok=True)

    # Look for test files
    test_files = list(data_dir.glob("*.DSC")) + list(data_dir.glob("*.dsc"))

    if test_files:
        test_file = test_files[0]
        print(f"Converting: {test_file.name}")

        # Show help first
        print("\nCommand help:")
        try:
            result = subprocess.run(['epyr-convert', '--help'],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                # Show first few lines of help
                help_lines = result.stdout.split('\n')[:8]
                for line in help_lines:
                    print(f"  {line}")
            else:
                print("  Help not available")
        except:
            print("  Help command failed")

        print(f"\nExample conversion command:")
        print(f"  epyr-convert {test_file.name} -o {output_dir.name} -f csv,json")

        # Simulate command execution (showing what would happen)
        print("\nThis command would:")
        print(f"  1. Load Bruker file: {test_file}")
        print(f"  2. Extract EPR data and parameters")
        print(f"  3. Convert to FAIR-compliant CSV and JSON formats")
        print(f"  4. Save to directory: {output_dir}")
        print(f"  5. Include comprehensive metadata")

    else:
        print("No EPR data files found for conversion demo")
        print("\nExample usage:")
        print("  epyr-convert spectrum.dsc")
        print("  epyr-convert data.par -o ./results -f csv,json,hdf5")
        print("  epyr-convert *.dsc --no-metadata --verbose")

    print()


def demo_baseline_command():
    """Demonstrate the epyr-baseline command."""
    print("3. epyr-baseline - Baseline correction command:")
    print("-" * 46)

    print("Command syntax and options:")
    print("  epyr-baseline input.dsc [options]")
    print("\nKey options:")
    print("  -m, --method METHOD    : polynomial, exponential, stretched_exponential")
    print("  --order INT           : Polynomial order (default: 1)")
    print("  --exclude START END   : Exclude signal region from fit")
    print("  --plot               : Generate comparison plot")
    print("  -o, --output FILE    : Output file")

    print("\nExample commands:")
    examples = [
        "epyr-baseline spectrum.dsc",
        "epyr-baseline data.dsc -m polynomial --order 2",
        "epyr-baseline spectrum.dsc --exclude 3480 3520 --plot",
        "epyr-baseline data.dsc -o corrected.csv --plot"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nThis command would:")
    print("  1. Load EPR spectrum")
    print("  2. Apply selected baseline correction algorithm")
    print("  3. Save original, baseline, and corrected data")
    print("  4. Generate comparison plot (if --plot specified)")
    print("  5. Provide quality metrics and fitting parameters")

    print()


def demo_batch_convert_command():
    """Demonstrate the epyr-batch-convert command."""
    print("4. epyr-batch-convert - Batch processing command:")
    print("-" * 49)

    print("Command for processing multiple files:")
    print("  epyr-batch-convert input_dir [options]")

    print("\nKey features:")
    print("  - Parallel processing support (-j, --jobs)")
    print("  - Custom file patterns (--pattern)")
    print("  - Progress reporting")
    print("  - Error handling per file")
    print("  - Summary statistics")

    print("\nExample commands:")
    examples = [
        "epyr-batch-convert ./data/",
        "epyr-batch-convert ./data/ --pattern '*.spc' --jobs 4",
        "epyr-batch-convert ./data/ -o ./converted -f csv,hdf5",
        "epyr-batch-convert ./data/ --jobs 2 --verbose"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nWorkflow for batch processing:")
    print("  1. Scan directory for matching files")
    print("  2. Process files in parallel (if jobs > 1)")
    print("  3. Convert each file to specified formats")
    print("  4. Handle errors gracefully")
    print("  5. Report conversion statistics")

    print()


def demo_config_command():
    """Demonstrate the epyr-config command."""
    print("5. epyr-config - Configuration management:")
    print("-" * 42)

    print("Configuration management with subcommands:")

    subcommands = {
        'show': 'Display current configuration',
        'set': 'Set configuration values',
        'reset': 'Reset to defaults',
        'export': 'Backup configuration',
        'import': 'Restore configuration'
    }

    print("\nSubcommands:")
    for cmd, desc in subcommands.items():
        print(f"  {cmd:<8} - {desc}")

    print("\nExample commands:")
    examples = [
        "epyr-config show",
        "epyr-config show plotting",
        "epyr-config set plotting.dpi 300",
        "epyr-config set plotting.figure_size '[10, 8]'",
        "epyr-config reset plotting",
        "epyr-config export my_settings.json",
        "epyr-config import my_settings.json"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nConfiguration sections:")
    sections = ['plotting', 'performance', 'data_loading', 'baseline', 'fair_conversion']
    for section in sections:
        print(f"  - {section}")

    print()


def demo_info_command():
    """Demonstrate the epyr-info command."""
    print("6. epyr-info - System information command:")
    print("-" * 42)

    print("Display comprehensive system information:")

    print("\nInformation categories:")
    categories = {
        '--config': 'Configuration details',
        '--performance': 'Performance metrics',
        '--plugins': 'Loaded plugins',
        '--all': 'Complete system report'
    }

    for flag, desc in categories.items():
        print(f"  {flag:<15} - {desc}")

    print("\nExample commands:")
    examples = [
        "epyr-info",
        "epyr-info --config",
        "epyr-info --performance",
        "epyr-info --all"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nTypical output includes:")
    info_items = [
        "EPyR Tools version",
        "Configuration file location",
        "Memory usage and system resources",
        "Loaded plugins and capabilities",
        "Performance settings",
        "Python environment details"
    ]

    for item in info_items:
        print(f"  - {item}")

    print()


def demo_plot_command():
    """Demonstrate the epyr-plot command."""
    print("7. epyr-plot - Interactive plotting command:")
    print("-" * 44)

    print("Interactive EPR data visualization with measurement tools:")

    print("\nKey features:")
    features = [
        "File dialog for data selection",
        "Interactive matplotlib backend",
        "Mouse-based measurement tool",
        "Scaling options (n=scans, P=power, G=gain)",
        "Save plots as PNG files",
        "EPR-optimized plotting functions"
    ]

    for feature in features:
        print(f"  - {feature}")

    print("\nCommand options:")
    options = [
        "-s, --scaling STRING   : Scaling string (n=scans, P=power, etc.)",
        "--interactive         : Enable interactive backend",
        "--measure            : Enable measurement tool",
        "--save              : Save plot as PNG",
        "--no-plot           : Load without plotting"
    ]

    for option in options:
        print(f"  {option}")

    print("\nExample commands:")
    examples = [
        "epyr-plot --interactive",
        "epyr-plot spectrum.dsc --interactive --measure",
        "epyr-plot data.dta -s nG --interactive --save",
        "epyr-plot --interactive --measure --verbose"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nInteractive measurement features:")
    measurements = [
        "Left-click: Select measurement points",
        "Right-click: Clear measurements",
        "Keyboard 'c': Clear measurements",
        "Keyboard 'q': Quit/close plot",
        "Display Δx, Δy, and Euclidean distance",
        "Visual annotations with lines and boxes"
    ]

    for measurement in measurements:
        print(f"  - {measurement}")

    print()


def demo_validate_command():
    """Demonstrate the epyr-validate command."""
    print("8. epyr-validate - Data validation command:")
    print("-" * 43)

    print("Validate EPR files and FAIR compliance:")

    print("\nValidation capabilities:")
    capabilities = [
        "File format integrity checking",
        "Data consistency validation",
        "FAIR metadata compliance",
        "EPR-specific parameter validation",
        "Detailed error reporting"
    ]

    for capability in capabilities:
        print(f"  - {capability}")

    print("\nCommand options:")
    options = [
        "--format FORMAT      : Expected file format",
        "--detailed          : Show detailed validation results",
        "-v, --verbose       : Enable verbose output"
    ]

    for option in options:
        print(f"  {option}")

    print("\nExample commands:")
    examples = [
        "epyr-validate spectrum.dsc",
        "epyr-validate spectrum.dsc --detailed",
        "epyr-validate *.dsc --verbose",
        "epyr-validate data.par --format ESP"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nValidation report includes:")
    report_items = [
        "File format verification",
        "Parameter completeness check",
        "Data integrity validation",
        "FAIR compliance assessment",
        "Recommendations for improvement"
    ]

    for item in report_items:
        print(f"  - {item}")

    print()


def demo_isotopes_command():
    """Demonstrate the epyr-isotopes command."""
    print("9. epyr-isotopes - Isotope database GUI:")
    print("-" * 40)

    print("Launch interactive nuclear isotope database:")

    print("\nGUI features:")
    features = [
        "Interactive periodic table visualization",
        "Color-coded elements by category",
        "Comprehensive isotope data display",
        "Real-time NMR frequency calculation",
        "EPR band presets (X, Q, W-band)",
        "Advanced filtering options",
        "Sortable data table",
        "Professional tooltips and interface"
    ]

    for feature in features:
        print(f"  - {feature}")

    print("\nData included:")
    data_items = [
        "Natural abundance percentages",
        "Nuclear spin values",
        "Gyromagnetic ratios (g-factors)",
        "Magnetic moments",
        "Quadrupole moments",
        "Both stable and radioactive isotopes"
    ]

    for item in data_items:
        print(f"  - {item}")

    print("\nUsage:")
    print("  epyr-isotopes")
    print("\nProgrammatic access:")
    print("  from epyr import isotopes")
    print("  isotopes()  # Launch GUI")

    print("\nData source:")
    print("  - EasySpin nuclear isotope database format")
    print("  - CODATA 2018 physical constants")
    print("  - Cross-platform compatibility")

    print()


def demo_main_entry_point():
    """Demonstrate the main epyr entry point."""
    print("10. epyr - Main CLI entry point:")
    print("-" * 33)

    print("Unified access to all commands through subcommands:")

    print("\nSyntax:")
    print("  epyr <command> [options]")

    print("\nAvailable subcommands:")
    subcommands = [
        "convert", "baseline", "batch-convert", "config",
        "info", "plot", "validate", "isotopes"
    ]

    for cmd in subcommands:
        print(f"  - {cmd}")

    print("\nExample usage:")
    examples = [
        "epyr convert spectrum.dsc",
        "epyr config show plotting",
        "epyr plot spectrum.dsc --interactive --measure",
        "epyr validate *.dsc --detailed",
        "epyr info --all"
    ]

    for example in examples:
        print(f"  {example}")

    print("\nHelp system:")
    help_examples = [
        "epyr --help",
        "epyr convert --help",
        "epyr baseline --help"
    ]

    for example in help_examples:
        print(f"  {example}")

    print()


def demo_cli_workflows():
    """Demonstrate common CLI workflows."""
    print("11. Common CLI workflows:")
    print("-" * 25)

    workflows = {
        "Quick data conversion": [
            "epyr-convert spectrum.dsc",
            "# Converts to CSV and JSON with metadata"
        ],

        "High-quality analysis": [
            "epyr-config set plotting.dpi 300",
            "epyr-baseline spectrum.dsc --plot",
            "epyr-plot spectrum.dsc --interactive --measure --save",
            "# Configure, correct, and analyze with measurements"
        ],

        "Batch processing": [
            "epyr-batch-convert ./data/ --jobs 4",
            "epyr-config set fair_conversion.default_formats '[\"csv\", \"hdf5\"]'",
            "# Process multiple files in parallel"
        ],

        "Quality assurance": [
            "epyr-validate *.dsc --detailed",
            "epyr-info --performance",
            "# Validate data and check system status"
        ],

        "Research workflow": [
            "epyr-isotopes  # Launch isotope database",
            "epyr-plot data.dsc --interactive --measure",
            "epyr-convert data.dsc -f csv,json,hdf5",
            "# Analyze, measure, and export for publication"
        ]
    }

    for workflow_name, commands in workflows.items():
        print(f"\n{workflow_name}:")
        for command in commands:
            if command.startswith('#'):
                print(f"  {command}")
            else:
                print(f"  $ {command}")

    print()


def create_cli_cheat_sheet():
    """Create a CLI command cheat sheet."""
    print("12. CLI Quick Reference:")
    print("-" * 24)

    cheat_sheet = """
EPyR Tools CLI Quick Reference
=============================

Data Loading & Conversion:
  epyr-convert file.dsc                    # Convert to CSV/JSON
  epyr-convert file.dsc -f csv,json,hdf5   # Multiple formats
  epyr-batch-convert ./data/ --jobs 4      # Parallel processing

Analysis:
  epyr-baseline file.dsc --plot            # Baseline correction
  epyr-plot file.dsc --interactive --measure # Interactive analysis
  epyr-validate *.dsc --detailed           # Data validation

Configuration:
  epyr-config show                         # View settings
  epyr-config set plotting.dpi 300         # Change setting
  epyr-config export settings.json         # Backup config

System:
  epyr-info --all                          # System information
  epyr-isotopes                           # Isotope database GUI

Help:
  epyr --help                             # Main help
  epyr-<command> --help                   # Command-specific help
"""

    print(cheat_sheet)

    # Save cheat sheet to file
    cheat_sheet_file = Path(__file__).parent / "cli_quick_reference.txt"
    with open(cheat_sheet_file, 'w') as f:
        f.write(cheat_sheet)

    print(f"Cheat sheet saved to: {cheat_sheet_file.name}")
    print()


def main():
    """Run all CLI demonstrations."""
    # Create output directory
    output_dir = Path(__file__).parent
    print(f"Output directory: {output_dir}")
    print()

    available_commands = check_cli_availability()

    demo_convert_command()
    demo_baseline_command()
    demo_batch_convert_command()
    demo_config_command()
    demo_info_command()
    demo_plot_command()
    demo_validate_command()
    demo_isotopes_command()
    demo_main_entry_point()
    demo_cli_workflows()
    create_cli_cheat_sheet()

    print("=== CLI Demo Complete ===")
    print()
    print("Key takeaways:")
    print("- EPyR Tools provides 9 professional CLI commands for all EPR workflows")
    print("- Commands follow consistent argument patterns and help systems")
    print("- Interactive plotting with measurement tools for precise analysis")
    print("- Batch processing capabilities for efficient multi-file workflows")
    print("- Comprehensive configuration management with persistent settings")
    print("- Data validation ensures quality and FAIR compliance")
    print("- GUI applications (isotope database) launched from command line")
    print("- Unified entry point (epyr) provides access to all subcommands")
    print()
    print("CLI advantages:")
    print("- No Python coding required for common tasks")
    print("- Scriptable and automatable workflows")
    print("- Consistent error handling and progress reporting")
    print("- Professional interface suitable for research environments")
    print("- Integration with system shell and other tools")

    if available_commands:
        print(f"\nTo try these commands, use any of the {len(available_commands)} available CLI tools.")
    else:
        print("\nTo install CLI commands: pip install -e . (from project root)")


if __name__ == "__main__":
    main()