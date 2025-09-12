# Configuration System - epyr.config

The `epyr.config` module provides a comprehensive, hierarchical configuration management system for EPyR Tools with support for defaults, user preferences, environment variables, and configuration files.

## Overview

EPyR Tools uses a centralized configuration system that:
- Provides sensible defaults for all settings
- Supports user customization through configuration files
- Allows environment variable overrides
- Enables configuration export/import for reproducibility
- Integrates seamlessly with all EPyR Tools modules

## Architecture

### EPyRConfig Class

The core `EPyRConfig` class manages configuration through a hierarchical structure:

```python
from epyr.config import EPyRConfig, config

# Global instance is automatically available
print(config.get('plotting.dpi'))  # 300
```

### Configuration Structure

```python
{
    "plotting": {
        "default_style": "publication",
        "dpi": 300,
        "figure_size": [8, 6],
        "color_scheme": "viridis",
        "font_size": 12,
        "line_width": 1.5,
        "grid_alpha": 0.3,
        "save_format": "png"
    },
    "data_loading": {
        "auto_plot": True,
        "scaling_default": "",
        "file_dialog_remember_dir": True,
        "supported_extensions": [".dta", ".dsc", ".spc", ".par"]
    },
    "baseline": {
        "default_poly_order": 1,
        "default_method": "polynomial",
        "exclusion_buffer": 0.1,
        "max_iterations": 1000
    },
    "fair_conversion": {
        "default_formats": ["csv", "json"],
        "include_metadata": True,
        "preserve_precision": True,
        "compression": "gzip"
    },
    "performance": {
        "cache_enabled": True,
        "cache_size_mb": 100,
        "chunk_size_mb": 10,
        "memory_limit_mb": 500,
        "parallel_processing": True
    },
    "logging": {
        "level": "INFO",
        "file_logging": False,
        "log_file": None,
        "console_output": True
    },
    "gui": {
        "theme": "default",
        "window_size": [800, 600],
        "remember_position": True,
        "auto_refresh": True
    },
    "advanced": {
        "debug_mode": False,
        "developer_mode": False,
        "experimental_features": False,
        "error_reporting": True
    }
}
```

## Core Methods

### Getting Configuration Values

```python
from epyr.config import config

# Get single value with dot notation
dpi = config.get('plotting.dpi')  # 300
cache_enabled = config.get('performance.cache_enabled')  # True

# Get with default value
custom_setting = config.get('custom.setting', 'default_value')

# Get entire section
plotting_config = config.get_section('plotting')
```

### Setting Configuration Values

```python
# Set single values
config.set('plotting.dpi', 150)
config.set('performance.cache_enabled', False)

# Set nested values (creates structure if needed)
config.set('custom.new.setting', 'value')

# Set entire sections
new_plotting = {
    'dpi': 300,
    'figure_size': [10, 8],
    'color_scheme': 'plasma'
}
config.set_section('plotting', new_plotting)

# Save changes to file
config.save()
```

### Configuration Persistence

```python
# Manual save
config.save()

# Export configuration
config.export_config('my_settings.json')

# Import configuration
config.import_config('my_settings.json')
config.save()  # Persist imported settings

# Reset sections or all configuration
config.reset_section('plotting')  # Reset to defaults
config.reset_all()                 # Reset everything
```

## Configuration File Locations

### Platform-Specific Paths

**Linux/macOS:**
```
~/.config/epyrtools/config.json
```

**Windows:**
```
%APPDATA%/EPyRTools/config.json
```

### Custom Configuration Directory

```python
# Get current config file path
config_path = config.get_config_file_path()
print(f"Configuration stored at: {config_path}")
```

## Environment Variable Overrides

EPyR Tools supports environment variable overrides using the `EPYR_` prefix:

```bash
# Override individual settings
export EPYR_PLOTTING_DPI=150
export EPYR_PERFORMANCE_CACHE_ENABLED=false
export EPYR_LOGGING_LEVEL=DEBUG

# Complex values as JSON
export EPYR_PLOTTING_FIGURE_SIZE='[10, 8]'
export EPYR_CUSTOM_SETTING='{"nested": "value"}'
```

**Conversion Rules:**
- Environment variable names: `EPYR_SECTION_SETTING`
- Underscores become dots: `EPYR_PLOTTING_DPI` → `plotting.dpi`
- Values are parsed as JSON first, then as strings
- Environment variables override file settings

## Configuration Sections

### 1. Plotting Configuration

Controls matplotlib plotting behavior throughout EPyR Tools:

```python
plotting_config = config.get_section('plotting')
```

**Settings:**
- `default_style`: Plot style ('publication', 'seaborn', etc.)
- `dpi`: Resolution for saved plots
- `figure_size`: Default figure dimensions [width, height]
- `color_scheme`: Color map for multi-line plots
- `font_size`: Default font size
- `line_width`: Default line width
- `grid_alpha`: Grid transparency (0-1)
- `save_format`: Default save format ('png', 'pdf', 'svg')

**Usage:**
```python
import matplotlib.pyplot as plt
from epyr.config import get_plotting_config

plot_config = get_plotting_config()
plt.rcParams['figure.dpi'] = plot_config['dpi']
plt.rcParams['font.size'] = plot_config['font_size']
```

### 2. Data Loading Configuration

Controls EPR data loading behavior:

```python
data_config = config.get_section('data_loading')
```

**Settings:**
- `auto_plot`: Whether to automatically display plots after loading
- `scaling_default`: Default scaling option for eprload
- `file_dialog_remember_dir`: Remember last directory in file dialogs
- `supported_extensions`: List of recognized EPR file extensions

### 3. Baseline Configuration

Default parameters for baseline correction:

```python
baseline_config = config.get_section('baseline')
```

**Settings:**
- `default_poly_order`: Default polynomial order
- `default_method`: Default baseline method
- `exclusion_buffer`: Buffer around peaks to exclude from fitting
- `max_iterations`: Maximum iterations for iterative methods

### 4. FAIR Conversion Configuration

Settings for FAIR data export:

```python
fair_config = config.get_section('fair_conversion')
```

**Settings:**
- `default_formats`: List of default export formats
- `include_metadata`: Include metadata in exports by default
- `preserve_precision`: Maintain full numerical precision
- `compression`: Compression method for large files

### 5. Performance Configuration

Performance optimization settings:

```python
from epyr.config import get_performance_config
perf_config = get_performance_config()
```

**Settings:**
- `cache_enabled`: Enable/disable data caching
- `cache_size_mb`: Maximum cache size in MB
- `chunk_size_mb`: Chunk size for large file processing
- `memory_limit_mb`: Memory usage warning threshold
- `parallel_processing`: Enable multi-core processing

### 6. Logging Configuration

Logging system settings:

```python
logging_config = config.get_section('logging')
```

**Settings:**
- `level`: Log level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
- `file_logging`: Enable logging to file
- `log_file`: Path to log file (None for auto-generation)
- `console_output`: Enable console logging

### 7. GUI Configuration

Settings for GUI applications:

```python
gui_config = config.get_section('gui')
```

**Settings:**
- `theme`: GUI theme ('default', 'dark', 'light')
- `window_size`: Default window size [width, height]
- `remember_position`: Remember window positions
- `auto_refresh`: Auto-refresh data displays

### 8. Advanced Configuration

Developer and experimental settings:

```python
from epyr.config import is_debug_mode
debug_enabled = is_debug_mode()
```

**Settings:**
- `debug_mode`: Enable debug features
- `developer_mode`: Enable developer tools
- `experimental_features`: Enable experimental features
- `error_reporting`: Enable anonymous error reporting

## Convenience Functions

EPyR Tools provides convenience functions for common configuration access:

```python
from epyr.config import (
    get_plotting_config,
    get_performance_config,
    is_debug_mode,
    is_cache_enabled
)

# Get configuration sections
plotting = get_plotting_config()
performance = get_performance_config()

# Boolean checks
if is_debug_mode():
    print("Debug mode enabled")

if is_cache_enabled():
    from epyr.performance import get_global_cache
    cache = get_global_cache()
```

## Integration Examples

### Module Integration

```python
# In a module that needs configuration
from epyr.config import config
from epyr.logging_config import get_logger

logger = get_logger(__name__)

def process_data():
    # Use configuration
    cache_enabled = config.get('performance.cache_enabled')
    chunk_size = config.get('performance.chunk_size_mb')
    
    if cache_enabled:
        logger.debug("Using cached data processing")
    
    # Process with configured chunk size
    process_in_chunks(chunk_size_mb=chunk_size)
```

### CLI Integration

```python
# CLI commands automatically respect configuration
def cmd_convert():
    parser = argparse.ArgumentParser()
    parser.add_argument('--formats', 
                       default=','.join(config.get('fair_conversion.default_formats')))
    # Use configured defaults
```

### Performance Integration

```python
# Performance module uses configuration automatically
from epyr.performance import OptimizedLoader

# Loader automatically uses configured settings
loader = OptimizedLoader()  # Uses config values for cache and chunk size
```

## Configuration Migration

### Version Compatibility

```python
def migrate_config_if_needed():
    """Migrate configuration from older versions."""
    current_version = config.get('_version', '0.1.0')
    
    if current_version < '0.1.6':
        # Migrate old settings
        old_cache_enabled = config.get('cache.enabled')
        if old_cache_enabled is not None:
            config.set('performance.cache_enabled', old_cache_enabled)
            config.save()
```

### Backup and Restore

```python
import datetime
from pathlib import Path

# Backup current configuration
backup_name = f"config_backup_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
config.export_config(backup_name)

# Restore from backup
config.import_config('config_backup_20250906_120000.json')
config.save()
```

## Error Handling

The configuration system handles errors gracefully:

```python
# Invalid JSON in config file
try:
    config._load_from_file()
except json.JSONDecodeError as e:
    logger.warning(f"Invalid config file, using defaults: {e}")
    # Continue with defaults

# Missing config file
if not config_file.exists():
    logger.info("Creating new configuration file")
    config.save()  # Create with defaults
```

## Best Practices

### Application Setup

```python
def setup_application():
    """Set up application with configuration."""
    from epyr.config import config
    from epyr.logging_config import setup_logging
    
    # Set up logging based on configuration
    log_level = config.get('logging.level')
    setup_logging(log_level)
    
    # Configure matplotlib
    import matplotlib.pyplot as plt
    plt.rcParams['figure.dpi'] = config.get('plotting.dpi')
    
    # Set up performance
    if config.get('performance.cache_enabled'):
        from epyr.performance import get_global_cache
        cache = get_global_cache()
        logger.info(f"Cache enabled: {cache.max_size_mb} MB")
```

### Configuration Validation

```python
def validate_config():
    """Validate configuration settings."""
    # Check memory limits
    memory_limit = config.get('performance.memory_limit_mb')
    cache_size = config.get('performance.cache_size_mb')
    
    if cache_size > memory_limit:
        logger.warning("Cache size exceeds memory limit")
        config.set('performance.cache_size_mb', memory_limit // 2)
        config.save()
    
    # Validate plot DPI
    dpi = config.get('plotting.dpi')
    if dpi < 50 or dpi > 1200:
        logger.warning(f"Unusual DPI setting: {dpi}")
```

### User Customization

```python
def customize_for_user():
    """Example user customization."""
    # High-resolution plotting for publications
    config.set('plotting.dpi', 300)
    config.set('plotting.save_format', 'pdf')
    
    # Large dataset optimization
    config.set('performance.cache_size_mb', 500)
    config.set('performance.chunk_size_mb', 50)
    
    # Verbose logging for debugging
    config.set('logging.level', 'DEBUG')
    config.set('logging.file_logging', True)
    
    config.save()
```

The configuration system provides a robust foundation for EPyR Tools, enabling user customization while maintaining sensible defaults and ensuring consistent behavior across all modules.