# CLI Reference - Command Line Interface

The `epyr.cli` module provides a comprehensive command-line interface with 9 professional commands for all EPR workflows.

## Overview

EPyR Tools CLI follows modern CLI design principles with:
- Consistent argument patterns across commands
- Comprehensive help system
- Progress reporting and verbose logging
- Error handling with detailed diagnostics
- Integration with configuration system

## Command Structure

```bash
epyr <command> [options]
# or use individual commands:
epyr-<command> [options]
```

## Available Commands

### 1. `epyr-convert` - Data Conversion

Convert Bruker EPR files to FAIR-compliant formats.

```bash
epyr-convert input.dsc [options]
```

**Options:**
- `-o, --output-dir DIR`: Output directory (default: current directory)
- `-f, --formats LIST`: Output formats: csv,json,hdf5 (default: csv,json)
- `--no-metadata`: Skip metadata export
- `-v, --verbose`: Enable verbose output

**Examples:**
```bash
# Basic conversion
epyr-convert spectrum.dsc

# Multiple formats with custom output
epyr-convert spectrum.dsc -o ./results -f csv,json,hdf5

# Batch processing without metadata
epyr-convert *.dsc --no-metadata --verbose
```

**Implementation Details:**
- Uses `epyr.fair.convert_bruker_to_fair()` internally
- Validates input file existence before processing
- Creates output directory if it doesn't exist
- Provides detailed error reporting on failure
- Returns exit code 1 on errors for script integration

### 2. `epyr-baseline` - Baseline Correction

Apply baseline correction algorithms to EPR data.

```bash
epyr-baseline input.dsc [options]
```

**Options:**
- `-o, --output FILE`: Output file (default: input_baseline.csv)
- `-m, --method METHOD`: Correction method (polynomial, exponential, stretched_exponential)
- `--order INT`: Polynomial order (default: 1)
- `--exclude START END`: Exclude region from fit (can be repeated)
- `--plot`: Generate comparison plot
- `-v, --verbose`: Enable verbose output

**Examples:**
```bash
# Basic polynomial baseline
epyr-baseline spectrum.dsc

# Advanced correction with exclusions
epyr-baseline spectrum.dsc -m polynomial --order 2 --exclude 3480 3520

# Generate plot with custom output
epyr-baseline spectrum.dsc -o corrected.csv --plot
```

**Implementation Details:**
- Integrates with `epyr.baseline` module
- Saves results as CSV with original, baseline, and corrected data
- Optional matplotlib plotting with dual-panel comparison
- Supports multiple exclusion regions for signal protection
- Automatic output filename generation

### 3. `epyr-batch-convert` - Batch Processing

Convert multiple EPR files efficiently with automatic figure generation for data quality verification.

```bash
epyr-batch-convert input_dir [options]
```

**Options:**
- `-o, --output-dir DIR`: Output directory (default: input_dir/converted)
- `-f, --formats LIST`: Output formats (default: csv,json)
- `--save-jpg`: Save JPG figures of loaded data (default: True)
- `--no-jpg`: Disable JPG figure generation
- `-j, --jobs INT`: Number of parallel jobs (default: 1)
- `-v, --verbose`: Enable verbose output

**Examples:**
```bash
# Convert all .dsc and .spc files with JPG generation
epyr-batch-convert ./data/

# Disable JPG generation for faster processing
epyr-batch-convert ./data/ --no-jpg

# Custom output directory with all formats
epyr-batch-convert ./data/ -o ./converted -f csv,json,hdf5

# Verbose mode for detailed progress tracking
epyr-batch-convert ./data/ -v
```

**JPG Figure Generation:**
By default, batch conversion now generates high-quality JPG figures (150 DPI) for visual verification of converted data:

- **1D EPR Data**: Line plots with real and imaginary components
  - Proper axis labels from EPR parameters (field unit, intensity)
  - Grid overlay for easy reading
  - Legend for complex data

- **2D EPR Data**: Color map visualizations
  - Intensity colorbars with proper units
  - Optimized layout without tight_layout conflicts
  - Clear axis labels for point and scan indices

**File Discovery:**
Automatically discovers EPR files with case-insensitive patterns:
- `*.dsc`, `*.DSC` (Bruker BES3T descriptor files)
- `*.spc`, `*.SPC` (Bruker ESP spectrum files)
- Removes duplicates on case-insensitive filesystems
- Processes files in sorted order for consistency

**Implementation Details:**
- Loads each file with `eprload()` before conversion
- Generates JPG only after successful data loading
- Independent error handling for conversion and figure generation
- Progress reporting: `[N/M] Processing filename.dsc`
- Detailed summary statistics:
  - Successfully converted: X/Y files
  - Failed: X/Y files
  - Output directory location
- Comprehensive error handling per file with optional verbose traceback
- Automatic output directory creation

### 4. `epyr-config` - Configuration Management

Manage EPyR Tools configuration with subcommands.

```bash
epyr-config <subcommand> [options]
```

**Subcommands:**

#### `show` - Display Configuration
```bash
epyr-config show [section]
```
- Show all configuration or specific section
- JSON-formatted output for easy parsing

#### `set` - Set Configuration Values
```bash
epyr-config set key value
```
- Supports dot notation (e.g., `plotting.dpi`)
- Automatic JSON parsing for complex values
- Immediate save to configuration file

#### `reset` - Reset Configuration
```bash
epyr-config reset [section|all]
```
- Reset specific section or all configuration
- Restores default values

#### `export/import` - Configuration Backup
```bash
epyr-config export config.json
epyr-config import config.json
```
- Full configuration backup and restore
- JSON format for portability

**Examples:**
```bash
# View current configuration
epyr-config show

# Set plotting DPI
epyr-config set plotting.dpi 300

# Set complex value (JSON)
epyr-config set plotting.figure_size '[10, 8]'

# Reset plotting section
epyr-config reset plotting

# Backup configuration
epyr-config export my_settings.json
```

### 5. `epyr-info` - System Information

Display comprehensive system and configuration information.

```bash
epyr-info [options]
```

**Options:**
- `--config`: Show detailed configuration
- `--performance`: Show performance metrics
- `--plugins`: Show loaded plugins
- `--all`: Show all information

**Examples:**
```bash
# Basic system info
epyr-info

# Performance diagnostics
epyr-info --performance

# Complete system report
epyr-info --all
```

**Information Displayed:**
- EPyR Tools version
- Configuration file location
- Memory usage and system resources
- Loaded plugins and capabilities
- Performance settings and optimization status

### 6. `epyr-validate` - Data Validation

Validate EPR files and check FAIR compliance.

```bash
epyr-validate files... [options]
```

**Options:**
- `--format FORMAT`: Expected file format
- `--detailed`: Show detailed validation results
- `-v, --verbose`: Enable verbose output

**Examples:**
```bash
# Basic validation
epyr-validate spectrum.dsc

# Detailed FAIR compliance check
epyr-validate spectrum.dsc --detailed

# Validate multiple files
epyr-validate *.dsc --verbose
```

**Validation Features:**
- File format integrity checking
- Data consistency validation
- FAIR metadata compliance assessment
- EPR-specific parameter validation
- Detailed error and warning reports

### 7. `epyr-plot` - Interactive Data Visualization

Load and visualize EPR data with interactive plotting and measurement tools.

```bash
epyr-plot [file] [options]
```

**Options:**
- `-s, --scaling STRING`: Scaling string (n=scans, P=power, G=gain, T=temp, c=time)
- `--no-plot`: Load data without plotting
- `--interactive`: Enable interactive matplotlib backend
- `--save`: Save plot as PNG file
- `--measure`: Enable interactive measurement tool (click two points to measure distance)
- `-v, --verbose`: Enable verbose output

**Examples:**
```bash
# Interactive plot with file dialog
epyr-plot --interactive

# Load specific file with measurement tools
epyr-plot spectrum.dsc --interactive --measure

# Load with scaling and save plot
epyr-plot data.dta -s nG --interactive --save

# Measurement mode with verbose output
epyr-plot --interactive --measure -v
```

**Interactive Measurement Features:**
- **Mouse Controls:**
  - Left-click: Select measurement points (click two points to measure)
  - Right-click: Clear all measurements
  - Standard matplotlib zoom/pan navigation
- **Keyboard Shortcuts:**
  - `c`: Clear measurements
  - `q`: Quit/close plot
- **Measurements Displayed:**
  - Δx: Horizontal distance between points
  - Δy: Vertical distance between points  
  - |Δ|: Euclidean distance between points
  - Visual annotations with yellow boxes and red dashed lines

**Console Output Example:**
```
📐 Measurement Results:
  Point 1: (3340.0000, 1.2345e-01)
  Point 2: (3360.0000, 8.9123e-01)
  Δx = 20.0000
  Δy = 7.6778e-01
  Distance = 21.4567e+00
```

**Implementation Details:**
- macOS optimized with TkAgg backend for smooth performance
- Custom InteractiveMeasurementTool class for precise measurements
- EPR-aware plotting with proper field/intensity labels
- Real-time visual feedback with point markers and connecting lines
- Multiple measurements supported - click two points repeatedly
- Smart plotting logic that disables default eprload plotting in measurement mode

### 8. `epyr-isotopes` - Isotope Database GUI

Launch the interactive nuclear isotope database interface with periodic table visualization.

```bash
epyr-isotopes
```

**Features:**
- **Interactive Periodic Table**: Color-coded elements by category (main group, transition metals, rare earth)
- **Comprehensive Isotope Data**: Natural abundance, nuclear spin, g-factors, quadrupole moments
- **Real-time NMR Frequency Calculation**: Enter magnetic field strength for frequency calculations
- **EPR Band Presets**: Quick selection for X-band (340 mT), Q-band (1200 mT), W-band (3400 mT)
- **Advanced Filtering**: Show/hide unstable isotopes (*) and non-magnetic nuclei (spin=0)
- **Element-Specific Views**: Click any element to filter isotopes for that element only
- **Sortable Data Table**: Click column headers to sort by any parameter
- **Professional GUI**: Tooltips, resizable interface, cross-platform compatibility

**Data Source:**
- EasySpin nuclear isotope database format (`sub/isotopedata.txt`)
- CODATA 2018 physical constants for accurate calculations
- Both stable and selected radioactive isotopes included

**Usage Examples:**
```python
# Programmatic access
from epyr import isotopes
isotopes()  # Launch GUI

# Direct CLI launch
epyr-isotopes
```

### 9. `epyr` - Main CLI Entry Point

Unified access to all commands through subcommands.

```bash
epyr <command> [options]
```

**Examples:**
```bash
epyr convert spectrum.dsc
epyr config show plotting
epyr plot spectrum.dsc --interactive --measure
epyr validate *.dsc --detailed
```

## CLI Architecture

### Error Handling
```python
# Consistent error patterns across commands
try:
    # Command implementation
    success = perform_operation()
    if not success:
        logger.error("Operation failed")
        sys.exit(1)
except Exception as e:
    logger.error(f"Unexpected error: {e}")
    if args.verbose:
        logger.debug("Full traceback:", exc_info=True)
    sys.exit(1)
```

### Logging Integration
```python
# Verbose logging setup
if args.verbose:
    from .logging_config import setup_logging
    setup_logging('DEBUG')
```

### Configuration Integration
```python
# All commands integrate with config system
from .config import config

# Use configured defaults
default_formats = config.get('fair_conversion.default_formats')
cache_enabled = config.get('performance.cache_enabled')
```

### Progress Reporting
```python
# Consistent progress patterns
logger.info(f"Processing {file_path}")
logger.info(f"[{current}/{total}] Converting {filename}")
logger.info(f"Operation completed: {success_count}/{total_count} successful")
```

## Development Integration

### Testing CLI Commands
```python
# Example test pattern
from unittest.mock import patch
import epyr.cli

def test_command():
    with patch('sys.argv', ['epyr-convert', 'test.dsc']):
        with patch('epyr.fair.convert_bruker_to_fair') as mock_convert:
            mock_convert.return_value = True
            epyr.cli.cmd_convert()
            mock_convert.assert_called_once()
```

### Adding New Commands

1. **Implement command function:**
```python
def cmd_newcommand():
    """New command implementation."""
    parser = argparse.ArgumentParser(
        prog='epyr-newcommand',
        description='Description of new command'
    )
    # Add arguments
    args = parser.parse_args()
    # Implementation
```

2. **Register in main CLI:**
```python
# In main() function
elif args.command == 'newcommand':
    cmd_newcommand()
```

3. **Add to pyproject.toml:**
```toml
[project.scripts]
epyr-newcommand = "epyr.cli:cmd_newcommand"
```

## Best Practices

### Command Design
- Use consistent argument names across commands
- Provide sensible defaults from configuration
- Support both verbose and quiet operation modes
- Return appropriate exit codes for script integration

### Error Messages
- Provide actionable error messages
- Include suggestions for common issues
- Use verbose mode for detailed diagnostics
- Log errors for debugging while showing user-friendly messages

### Performance
- Check file existence before processing
- Use progress reporting for long operations
- Integrate with performance monitoring
- Support interruption (Ctrl+C) gracefully

### Integration
- Leverage existing EPyR Tools modules
- Use configuration system for defaults
- Integrate with logging system
- Support plugin system extensions

This CLI system provides a professional interface to all EPyR Tools capabilities, making the package accessible to users who prefer command-line workflows while maintaining full integration with the underlying Python API.