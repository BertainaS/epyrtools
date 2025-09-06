# EPyR Tools User Guide

This comprehensive guide covers all aspects of using EPyR Tools for EPR spectroscopy data analysis.

## Table of Contents
1. [Getting Started](#getting-started)
2. [Command Line Interface](#command-line-interface)
3. [Data Loading and Conversion](#data-loading-and-conversion)
4. [FAIR Data Standards](#fair-data-standards)
5. [Configuration Management](#configuration-management)
6. [Plugin System](#plugin-system)
7. [Performance Optimization](#performance-optimization)
8. [Troubleshooting](#troubleshooting)

## Getting Started

### Installation

Install EPyR Tools using pip:

```bash
pip install epyr-tools
```

For development installation:

```bash
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools
pip install -e ".[dev,docs]"
```

### First Steps

1. **Test your installation**:
   ```bash
   epyr --help
   epyr-info
   ```

2. **Configure EPyR Tools**:
   ```bash
   epyr-config show
   epyr-config set plotting.dpi 300
   ```

3. **Convert your first file**:
   ```bash
   epyr-convert your_data.dsc --formats csv,json
   ```

## Command Line Interface

EPyR Tools provides 8 professional CLI commands:

### Core Commands

#### `epyr-convert` - Data Conversion
Convert Bruker EPR files to FAIR formats:

```bash
# Basic conversion
epyr-convert spectrum.dsc

# Specify output directory and formats
epyr-convert spectrum.dsc --output-dir ./results --formats csv,json,hdf5

# Skip metadata export
epyr-convert spectrum.dsc --no-metadata --verbose
```

**Supported formats:**
- **CSV**: Comma-separated values with optional metadata
- **JSON**: Structured JSON with complete metadata 
- **HDF5**: Hierarchical data format for large datasets

#### `epyr-batch-convert` - Batch Processing
Convert multiple files efficiently:

```bash
# Convert all .dsc files in directory
epyr-batch-convert ./data/ --formats csv,json

# Use custom pattern and parallel processing
epyr-batch-convert ./data/ --pattern "*.spc" --jobs 4

# Specify output directory
epyr-batch-convert ./data/ --output-dir ./converted/
```

#### `epyr-baseline` - Baseline Correction
Apply baseline correction to EPR data:

```bash
# Basic polynomial baseline correction
epyr-baseline spectrum.dsc --method polynomial --order 2

# Exclude signal regions from fitting
epyr-baseline spectrum.dsc --exclude 3480 3520 --exclude 3450 3460

# Generate comparison plot
epyr-baseline spectrum.dsc --plot --output corrected_spectrum.csv
```

**Available methods:**
- `polynomial`: Polynomial baseline (default)
- `exponential`: Exponential decay baseline
- `stretched_exponential`: Stretched exponential baseline

### Utility Commands

#### `epyr-validate` - Data Validation
Validate EPR files and check FAIR compliance:

```bash
# Basic validation
epyr-validate spectrum.dsc

# Detailed validation with FAIR compliance
epyr-validate spectrum.dsc --detailed

# Validate multiple files
epyr-validate *.dsc --verbose
```

#### `epyr-config` - Configuration Management
Manage EPyR Tools settings:

```bash
# Show all configuration
epyr-config show

# Show specific section
epyr-config show plotting

# Set configuration values
epyr-config set plotting.dpi 300
epyr-config set performance.cache_enabled true

# Reset configuration
epyr-config reset plotting
epyr-config reset all

# Export/import configuration
epyr-config export my_config.json
epyr-config import my_config.json
```

#### `epyr-info` - System Information
Display system and configuration information:

```bash
# Basic info
epyr-info

# Detailed configuration
epyr-info --config

# Performance information
epyr-info --performance

# Plugin information
epyr-info --plugins

# All information
epyr-info --all
```

#### `epyr-isotopes` - Interactive GUI
Launch the isotope database GUI:

```bash
epyr-isotopes
```

### Main CLI Entry Point

The main `epyr` command provides access to all subcommands:

```bash
epyr convert spectrum.dsc
epyr config show
epyr validate *.dsc
```

## Data Loading and Conversion

### Supported File Formats

EPyR Tools supports various Bruker EPR file formats:

- **BES3T**: `.dta/.dsc` file pairs (modern Bruker format)
- **WinEPR**: `.par/.spc` file pairs (legacy format)
- **ESP**: Single `.par` files

### Loading Data Programmatically

```python
import epyr

# Load EPR data
x, y, params, file_path = epyr.eprload('spectrum.dsc')

# Load with specific scaling
x, y, params, file_path = epyr.eprload('spectrum.dsc', scaling='nPGT')

# Disable automatic plotting
x, y, params, file_path = epyr.eprload('spectrum.dsc', plot_if_possible=False)
```

### FAIR Data Conversion

```python
from epyr.fair import convert_bruker_to_fair

# Convert single file
success = convert_bruker_to_fair(
    'spectrum.dsc',
    output_dir='./results/',
    formats=['csv', 'json', 'hdf5'],
    include_metadata=True
)

# Batch conversion
from epyr.fair import batch_convert_directory

batch_convert_directory(
    './data/',
    output_directory='./converted/',
    file_extensions=['.dsc', '.spc'],
    formats=['csv', 'json'],
    recursive=True
)
```

## FAIR Data Standards

### What are FAIR Principles?

FAIR data principles ensure data is:
- **Findable**: Rich metadata and unique identifiers
- **Accessible**: Clear access protocols and metadata preservation
- **Interoperable**: Standard formats and vocabularies
- **Reusable**: Clear licenses and provenance information

### FAIR Compliance Validation

```python
from epyr.fair.validation import validate_fair_dataset, create_validation_report

# Load and validate data
x, y, params, _ = epyr.eprload('spectrum.dsc')
data_dict = {
    'x_data': x,
    'y_data': y,
    'metadata': params
}

# Comprehensive validation
result = validate_fair_dataset(data_dict)

# Generate report
report = create_validation_report(result)
print(report)
```

### Metadata Standards

EPyR Tools enforces comprehensive metadata standards:

**Required FAIR metadata:**
- `title`: Dataset title
- `description`: Dataset description  
- `creator`: Data creator/author
- `date_created`: Creation date (ISO 8601)
- `identifier`: Unique identifier
- `format`: Data format
- `license`: Usage license

**Recommended EPR metadata:**
- `instrument`: EPR spectrometer information
- `measurement_parameters`: Measurement conditions
- `sample_information`: Sample description
- `processing_history`: Data processing steps
- `units`: Data units information

## Configuration Management

### Configuration Structure

EPyR Tools uses a hierarchical configuration system:

```python
from epyr.config import config

# Plotting settings
config.get('plotting.dpi')                 # 300
config.get('plotting.default_style')       # 'publication'
config.get('plotting.figure_size')         # [8, 6]

# Performance settings
config.get('performance.cache_enabled')    # True
config.get('performance.cache_size_mb')    # 100
config.get('performance.parallel_processing') # True

# FAIR conversion settings
config.get('fair_conversion.default_formats')  # ['csv', 'json']
config.get('fair_conversion.include_metadata') # True
```

### Configuration Files

Configuration is stored in platform-appropriate locations:
- **Linux/macOS**: `~/.config/epyrtools/config.json`
- **Windows**: `%APPDATA%/EPyRTools/config.json`

### Environment Variables

Override configuration with environment variables:

```bash
export EPYR_PLOTTING_DPI=150
export EPYR_PERFORMANCE_CACHE_ENABLED=false
export EPYR_LOGGING_LEVEL=DEBUG
```

## Plugin System

### Using Plugins

EPyR Tools supports extensible plugins for:
- **File formats**: Add support for new data formats
- **Processing methods**: Custom analysis algorithms
- **Export formats**: Additional output formats

```python
from epyr.plugins import plugin_manager

# List available plugins
plugins = plugin_manager.list_plugins()

# Get file format plugin
plugin = plugin_manager.get_file_format_plugin(Path('data.xyz'))
if plugin:
    x, y, params = plugin.load(Path('data.xyz'))

# Get export plugin
exporter = plugin_manager.get_export_plugin('csv')
exporter.export(output_path, x_data, y_data, metadata)
```

### Creating Custom Plugins

```python
from epyr.plugins import FileFormatPlugin
import numpy as np

class MyFormatPlugin(FileFormatPlugin):
    plugin_name = "My Custom Format"
    format_name = "myformat"
    file_extensions = [".myf"]
    
    def initialize(self) -> bool:
        return True
    
    def can_load(self, file_path):
        return file_path.suffix.lower() == '.myf'
    
    def load(self, file_path):
        # Your loading logic
        x_data = np.array([...])
        y_data = np.array([...])
        metadata = {...}
        return x_data, y_data, metadata

# Register plugin
from epyr.plugins import plugin_manager
plugin_manager.register_plugin(MyFormatPlugin())
```

## Performance Optimization

### Memory Management

EPyR Tools includes intelligent memory management:

```python
from epyr.performance import OptimizedLoader, DataCache

# Use optimized loader for large files
loader = OptimizedLoader(chunk_size_mb=10, cache_enabled=True)
x, y, params, path = loader.load_epr_file('large_dataset.dsc')

# Manual cache management
cache = DataCache(max_size_mb=200)
cache.put(file_path, (x, y, params))
cached_data = cache.get(file_path)
```

### Performance Configuration

```python
from epyr.config import config

# Configure performance settings
config.set('performance.cache_enabled', True)
config.set('performance.cache_size_mb', 200)
config.set('performance.chunk_size_mb', 20)
config.set('performance.memory_limit_mb', 1000)
config.set('performance.parallel_processing', True)
```

### Monitoring Performance

```python
from epyr.performance import get_performance_info, MemoryMonitor

# Get system performance info
perf_info = get_performance_info()
print(f"Memory usage: {perf_info['memory']['rss_mb']:.1f} MB")

# Monitor memory during processing
if not MemoryMonitor.check_memory_limit():
    MemoryMonitor.optimize_memory()
```

## Troubleshooting

### Common Issues

#### 1. File Loading Problems

**Problem**: "Failed to load EPR file"
```bash
# Check file format and integrity
epyr-validate your_file.dsc --detailed
```

**Problem**: "Unsupported file format"
```bash
# Check supported extensions
epyr-info --plugins
```

#### 2. Memory Issues

**Problem**: "Memory limit exceeded"
```bash
# Increase memory limit
epyr-config set performance.memory_limit_mb 1000

# Enable chunked processing
epyr-config set performance.chunk_size_mb 5
```

#### 3. Conversion Failures

**Problem**: "FAIR conversion failed"
```bash
# Run with verbose logging
epyr-convert file.dsc --verbose

# Check FAIR compliance
epyr-validate file.dsc --detailed
```

#### 4. Configuration Issues

**Problem**: "Configuration not persistent"
```bash
# Check config file location
epyr-config show | head -1

# Manually save configuration
epyr-config export backup.json
epyr-config reset all
epyr-config import backup.json
```

### Debugging

Enable debug logging:

```bash
export EPYR_LOGGING_LEVEL=DEBUG
epyr-convert file.dsc --verbose
```

Or programmatically:

```python
from epyr.logging_config import setup_logging
setup_logging('DEBUG')
```

### Getting Help

1. **Command help**:
   ```bash
   epyr-convert --help
   epyr-config --help
   ```

2. **System information**:
   ```bash
   epyr-info --all
   ```

3. **GitHub Issues**: [Report bugs](https://github.com/BertainaS/epyrtools/issues)

### Performance Tips

1. **Large datasets**:
   - Enable caching for repeated access
   - Use chunked processing for very large files
   - Increase memory limits appropriately

2. **Batch processing**:
   - Use parallel processing (`--jobs` parameter)
   - Process in smaller batches if memory constrained

3. **Configuration optimization**:
   - Disable unnecessary features for production
   - Tune cache sizes based on available memory
   - Use appropriate output formats for your needs

This user guide covers the essential aspects of using EPyR Tools effectively. For more detailed information, refer to the API documentation and example notebooks.