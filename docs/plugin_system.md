# Plugin System - epyr.plugins

The `epyr.plugins` module provides a comprehensive, extensible plugin architecture that allows users to add support for new file formats, data processing methods, and export formats without modifying the core EPyR Tools codebase.

## Overview

The plugin system enables:
- **File Format Support**: Add loaders for custom or proprietary EPR file formats
- **Processing Extensions**: Implement custom data analysis algorithms
- **Export Capabilities**: Create exporters for specialized output formats
- **Auto-Discovery**: Automatic plugin detection and loading
- **Hot Reloading**: Dynamic plugin registration and management

## Architecture

### Plugin Hierarchy

```
BasePlugin (Abstract)
├── FileFormatPlugin (File I/O)
├── ProcessingPlugin (Data Analysis)  
└── ExportPlugin (Data Export)
```

### Core Components

1. **BasePlugin**: Abstract base class defining plugin interface
2. **Specialized Plugin Types**: File format, processing, and export plugins
3. **PluginManager**: Centralized plugin registration and discovery
4. **Auto-Discovery**: Automatic plugin loading from directories

## Base Plugin System

### BasePlugin Class

All plugins inherit from the abstract `BasePlugin` class:

```python
from epyr.plugins import BasePlugin

class MyPlugin(BasePlugin):
    # Plugin metadata
    plugin_name = "My Custom Plugin"
    plugin_version = "1.0.0"
    plugin_description = "Description of what this plugin does"
    plugin_author = "Your Name"
    
    def initialize(self) -> bool:
        """Initialize the plugin. Called when plugin is loaded."""
        # Setup code here
        return True  # Return False if initialization fails
    
    def cleanup(self):
        """Cleanup plugin resources. Called when plugin is unloaded."""
        # Cleanup code here
        pass
```

### Plugin Metadata

Every plugin provides standard metadata:

```python
plugin_info = my_plugin.get_info()
# Returns:
{
    'name': 'My Custom Plugin',
    'version': '1.0.0', 
    'description': 'Description of what this plugin does',
    'author': 'Your Name',
    'class': 'MyPlugin'
}
```

## File Format Plugins

### FileFormatPlugin Class

Extend EPyR Tools to support new file formats:

```python
from epyr.plugins import FileFormatPlugin
import numpy as np
from pathlib import Path

class CustomFormatPlugin(FileFormatPlugin):
    # Plugin metadata
    plugin_name = "Custom Format Loader"
    format_name = "custom_format"
    file_extensions = [".cst", ".custom"]
    supports_loading = True
    supports_saving = False
    
    def initialize(self) -> bool:
        """Initialize the plugin."""
        return True
    
    def can_load(self, file_path: Path) -> bool:
        """Check if this plugin can load the given file."""
        return file_path.suffix.lower() in self.file_extensions
    
    def load(self, file_path: Path, **kwargs):
        """Load data from custom format file."""
        with open(file_path, 'r') as f:
            # Parse your custom format
            lines = f.readlines()
            
            # Extract data (example)
            x_data = np.array([float(line.split()[0]) for line in lines[1:]])
            y_data = np.array([float(line.split()[1]) for line in lines[1:]])
            
            # Extract metadata
            header_line = lines[0].strip()
            metadata = {
                'format': 'custom_format',
                'file_path': str(file_path),
                'header_info': header_line,
                'data_points': len(x_data)
            }
            
            return x_data, y_data, metadata
```

### File Format Plugin Features

**Format Detection:**
```python
# Plugin manager automatically finds appropriate plugin
from epyr.plugins import plugin_manager

plugin = plugin_manager.get_file_format_plugin(Path("data.cst"))
if plugin:
    x, y, params = plugin.load(Path("data.cst"))
```

**Save Support:**
```python
class CustomFormatPlugin(FileFormatPlugin):
    supports_saving = True
    
    def save(self, file_path: Path, x_data, y_data, parameters, **kwargs):
        """Save data in custom format."""
        with open(file_path, 'w') as f:
            f.write("# Custom format file\n")
            for x, y in zip(x_data, y_data):
                f.write(f"{x:.6f}\t{y:.6f}\n")
        return True
```

**Integration with EPyR Tools:**
```python
# After registration, your format works with eprload
import epyr

# This will use your custom plugin automatically
x, y, params, path = epyr.eprload("data.cst")
```

## Processing Plugins

### ProcessingPlugin Class

Add custom data analysis algorithms:

```python
from epyr.plugins import ProcessingPlugin
import numpy as np

class NoiseReductionPlugin(ProcessingPlugin):
    plugin_name = "Advanced Noise Reduction"
    processing_name = "advanced_noise_reduction"
    input_requirements = ["y_data"]  # Required input data
    output_types = ["processed_data", "noise_estimate"]
    
    def initialize(self) -> bool:
        """Initialize processing plugin."""
        return True
    
    def process(self, data, **kwargs):
        """Apply noise reduction algorithm."""
        y_data = data['y_data']
        
        # Your custom algorithm here
        # Example: Moving average noise reduction
        window_size = kwargs.get('window_size', 5)
        processed = np.convolve(y_data, np.ones(window_size)/window_size, mode='same')
        
        # Estimate noise level
        noise_estimate = np.std(y_data - processed)
        
        return {
            'processed_data': processed,
            'noise_estimate': noise_estimate,
            'parameters': {
                'window_size': window_size,
                'original_noise_std': noise_estimate
            }
        }
    
    def validate_input(self, data):
        """Validate input data meets requirements."""
        return super().validate_input(data) and len(data['y_data']) > 0
```

### Processing Plugin Usage

```python
from epyr.plugins import plugin_manager

# Get processing plugin
processor = plugin_manager.get_processing_plugin("advanced_noise_reduction")

if processor:
    # Validate input
    input_data = {'y_data': epr_spectrum}
    if processor.validate_input(input_data):
        # Process data
        result = processor.process(input_data, window_size=7)
        
        processed_spectrum = result['processed_data']
        noise_level = result['noise_estimate']
```

## Export Plugins

### ExportPlugin Class

Create custom data export formats:

```python
from epyr.plugins import ExportPlugin
import json
from pathlib import Path

class JSONExportPlugin(ExportPlugin):
    plugin_name = "Enhanced JSON Exporter"
    export_format = "enhanced_json"
    file_extension = ".json"
    supports_metadata = True
    
    def initialize(self) -> bool:
        """Initialize export plugin."""
        return True
    
    def export(self, output_path: Path, x_data, y_data, parameters, **kwargs):
        """Export data to enhanced JSON format."""
        
        # Create comprehensive data structure
        export_data = {
            'metadata': {
                'format_version': '2.0',
                'export_timestamp': datetime.datetime.now().isoformat(),
                'data_points': len(y_data),
                'x_range': [float(np.min(x_data)), float(np.max(x_data))],
                'y_range': [float(np.min(y_data)), float(np.max(y_data))]
            },
            'parameters': parameters,
            'data': {
                'x_axis': x_data.tolist(),
                'y_axis': y_data.tolist()
            },
            'statistics': {
                'mean': float(np.mean(y_data)),
                'std': float(np.std(y_data)),
                'max': float(np.max(y_data)),
                'min': float(np.min(y_data))
            }
        }
        
        # Write to file
        with open(output_path, 'w') as f:
            json.dump(export_data, f, indent=2, ensure_ascii=False)
        
        return True
```

### Export Plugin Usage

```python
from epyr.plugins import plugin_manager

# Get export plugin
exporter = plugin_manager.get_export_plugin("enhanced_json")

if exporter:
    success = exporter.export(
        output_path=Path("spectrum_enhanced.json"),
        x_data=x_values,
        y_data=y_values,
        parameters=metadata_dict
    )
```

## Plugin Manager

### PluginManager Class

The `PluginManager` handles all plugin operations:

```python
from epyr.plugins import plugin_manager

# Global plugin manager instance is automatically available
```

### Plugin Registration

**Manual Registration:**
```python
# Create and register plugin
my_plugin = CustomFormatPlugin()
success = plugin_manager.register_plugin(my_plugin)

if success:
    print("Plugin registered successfully")
else:
    print("Plugin registration failed")
```

**Bulk Registration:**
```python
plugins = [
    CustomFormatPlugin(),
    NoiseReductionPlugin(),
    JSONExportPlugin()
]

for plugin in plugins:
    plugin_manager.register_plugin(plugin)
```

### Plugin Discovery

**Auto-Discovery from Directories:**
```python
from pathlib import Path

# Discover plugins in custom directories
plugin_dirs = [
    Path("~/.epyr/plugins").expanduser(),
    Path("./my_plugins/"),
    Path("/usr/local/lib/epyr/plugins/")
]

loaded_count = plugin_manager.discover_plugins(plugin_dirs)
print(f"Loaded {loaded_count} plugins")
```

**Default Discovery:**
```python
# Discover from default locations
loaded_count = plugin_manager.discover_plugins()
# Searches:
# - User: ~/.config/epyrtools/plugins/
# - System: <epyr_package>/plugins/
```

### Plugin Management

**List All Plugins:**
```python
plugins_info = plugin_manager.list_plugins()

# Returns categorized plugin information:
{
    'file_formats': [
        {'name': 'Custom Format Loader', 'version': '1.0.0', ...}
    ],
    'processing': [
        {'name': 'Advanced Noise Reduction', 'version': '1.0.0', ...}
    ],
    'export': [
        {'name': 'Enhanced JSON Exporter', 'version': '1.0.0', ...}
    ]
}
```

**Find Specific Plugins:**
```python
# Find file format plugin
file_plugin = plugin_manager.get_file_format_plugin(Path("data.cst"))

# Find processing plugin
proc_plugin = plugin_manager.get_processing_plugin("noise_reduction")

# Find export plugin
export_plugin = plugin_manager.get_export_plugin("enhanced_json")
```

**Unregister Plugins:**
```python
success = plugin_manager.unregister_plugin("Custom Format Loader")
if success:
    print("Plugin unregistered and cleaned up")
```

**Get Supported Extensions:**
```python
extensions = plugin_manager.get_supported_extensions()
print(f"Supported file extensions: {extensions}")
# Output: ['.cst', '.custom', '.dsc', '.dta', '.par', '.spc']
```

## Plugin Development

### Development Workflow

1. **Create Plugin File:**
```python
# my_plugin.py
from epyr.plugins import FileFormatPlugin
import numpy as np

class MyFormatPlugin(FileFormatPlugin):
    plugin_name = "My Format Plugin"
    format_name = "myformat"
    file_extensions = [".myf"]
    
    def initialize(self):
        return True
    
    def can_load(self, file_path):
        return file_path.suffix.lower() == '.myf'
    
    def load(self, file_path):
        # Implementation here
        return x_data, y_data, metadata
```

2. **Test Plugin:**
```python
# test_my_plugin.py
from my_plugin import MyFormatPlugin
from epyr.plugins import plugin_manager

# Test plugin registration
plugin = MyFormatPlugin()
assert plugin_manager.register_plugin(plugin)

# Test plugin functionality
test_file = Path("test_data.myf")
assert plugin.can_load(test_file)

x, y, params = plugin.load(test_file)
assert len(x) > 0 and len(y) > 0
```

3. **Install Plugin:**
```bash
# Option 1: Copy to user plugin directory
cp my_plugin.py ~/.config/epyrtools/plugins/

# Option 2: Install as package
pip install -e ./my_plugin_package/
```

### Plugin Distribution

**As Python Package:**
```python
# setup.py
from setuptools import setup

setup(
    name="epyr-my-plugin",
    version="1.0.0",
    py_modules=["my_plugin"],
    install_requires=["epyr-tools>=0.1.6"],
    entry_points={
        'epyr_plugins': [
            'my_plugin = my_plugin:MyFormatPlugin',
        ],
    },
)
```

**Plugin Directory Structure:**
```
my_epyr_plugins/
├── __init__.py
├── format_plugins/
│   ├── __init__.py
│   └── my_format.py
├── processing_plugins/
│   ├── __init__.py
│   └── my_processor.py
└── export_plugins/
    ├── __init__.py
    └── my_exporter.py
```

### Error Handling

**Plugin Initialization Errors:**
```python
class RobustPlugin(FileFormatPlugin):
    def initialize(self):
        try:
            # Setup that might fail
            self.setup_resources()
            return True
        except Exception as e:
            logger.error(f"Plugin initialization failed: {e}")
            return False  # Plugin won't be registered
    
    def cleanup(self):
        try:
            # Cleanup resources
            self.cleanup_resources()
        except Exception as e:
            logger.warning(f"Plugin cleanup failed: {e}")
```

**Runtime Error Handling:**
```python
def load(self, file_path):
    try:
        # File loading logic
        return self._load_file(file_path)
    except FileNotFoundError:
        raise  # Re-raise expected errors
    except Exception as e:
        logger.error(f"Unexpected error loading {file_path}: {e}")
        raise RuntimeError(f"Plugin load failed: {e}")
```

## Built-in Plugins

EPyR Tools includes built-in plugins:

### CSV Export Plugin
```python
from epyr.plugins import plugin_manager

csv_plugin = plugin_manager.get_export_plugin("csv")
if csv_plugin:
    csv_plugin.export(
        output_path=Path("spectrum.csv"),
        x_data=field_values,
        y_data=intensity_values,
        parameters=epr_parameters,
        include_metadata=True
    )
```

## Advanced Features

### Plugin Configuration

```python
class ConfigurablePlugin(ProcessingPlugin):
    def initialize(self):
        from epyr.config import config
        
        # Plugin can access EPyR configuration
        self.default_window = config.get('plugins.noise_reduction.window', 5)
        return True
    
    def process(self, data, **kwargs):
        window = kwargs.get('window_size', self.default_window)
        # Use configuration in processing
```

### Plugin Dependencies

```python
class DependentPlugin(FileFormatPlugin):
    def initialize(self):
        try:
            import specialized_library
            self.library = specialized_library
            return True
        except ImportError:
            logger.error("Required library not available")
            return False
    
    def load(self, file_path):
        return self.library.load_data(file_path)
```

### Plugin Hooks

```python
class HookedPlugin(ProcessingPlugin):
    def __init__(self):
        super().__init__()
        self.pre_process_hooks = []
        self.post_process_hooks = []
    
    def add_pre_hook(self, hook_func):
        self.pre_process_hooks.append(hook_func)
    
    def process(self, data, **kwargs):
        # Run pre-processing hooks
        for hook in self.pre_process_hooks:
            data = hook(data, **kwargs)
        
        # Main processing
        result = self._main_process(data, **kwargs)
        
        # Run post-processing hooks
        for hook in self.post_process_hooks:
            result = hook(result, **kwargs)
        
        return result
```

## Best Practices

### Plugin Design
- Keep plugins focused on a single responsibility
- Provide comprehensive error handling and logging
- Use configuration system for customizable parameters
- Include thorough documentation and examples

### Performance
- Initialize expensive resources once in `initialize()`
- Use lazy loading for optional dependencies
- Cache computed results when appropriate
- Monitor memory usage for data processing plugins

### Compatibility
- Specify minimum EPyR Tools version requirements
- Handle missing dependencies gracefully
- Provide fallback behavior when possible
- Test with multiple EPyR Tools versions

### Security
- Validate all input data thoroughly
- Avoid executing arbitrary code from data files
- Use safe file handling practices
- Follow principle of least privilege

The plugin system makes EPyR Tools highly extensible, allowing the community to add support for new instruments, analysis methods, and workflows while maintaining a clean, modular architecture.