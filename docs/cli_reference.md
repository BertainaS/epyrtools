# CLI Reference - Command Line Interface

The `epyr.cli` module provides a comprehensive command-line interface with 8 professional commands for all EPR workflows.

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

Convert multiple EPR files efficiently with parallel processing.

```bash
epyr-batch-convert input_dir [options]
```

**Options:**
- `-o, --output-dir DIR`: Output directory (default: input_dir/converted)
- `-f, --formats LIST`: Output formats (default: csv,json)
- `--pattern PATTERN`: File pattern to match (default: *.dsc)
- `-j, --jobs INT`: Number of parallel jobs (default: 1)
- `-v, --verbose`: Enable verbose output

**Examples:**
```bash
# Convert all .dsc files
epyr-batch-convert ./data/

# Parallel processing with custom pattern
epyr-batch-convert ./data/ --pattern "*.spc" --jobs 4

# Custom output and formats
epyr-batch-convert ./data/ -o ./converted -f csv,hdf5
```

**Implementation Details:**
- Discovers files using glob patterns
- Progress reporting with file counter
- Parallel processing support (future enhancement)
- Comprehensive error handling per file
- Summary statistics on completion

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

### 7. `epyr-isotopes` - Isotope Database GUI

Launch the interactive isotope database interface.

```bash
epyr-isotopes
```

**Features:**
- Tkinter-based GUI application
- Nuclear isotope database browser
- Magnetic properties and parameters
- Search and filter capabilities
- Integration with EPR calculations

### 8. `epyr` - Main CLI Entry Point

Unified access to all commands through subcommands.

```bash
epyr <command> [options]
```

**Examples:**
```bash
epyr convert spectrum.dsc
epyr config show plotting
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