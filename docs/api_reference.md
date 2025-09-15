# API Reference

Complete API reference for EPyR Tools modules.

## Core Modules

### epyr.eprload

Main data loading function for EPR files.

```python
def eprload(
    file_name: Optional[str] = None,
    scaling: str = "",
    plot_if_possible: bool = True
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any], str]
```

**Parameters:**
- `file_name`: Path to EPR file or None for file dialog
- `scaling`: Scaling options ('nPGT' for field normalization)
- `plot_if_possible`: Whether to display plot automatically

**Returns:**
- `x_data`: Magnetic field or frequency values
- `y_data`: EPR intensity data
- `params`: Dictionary of measurement parameters
- `file_path`: Full path to loaded file

**Example:**
```python
import epyr
x, y, params, path = epyr.eprload('spectrum.dsc')
```

## Baseline Correction Package

### epyr.baseline

Complete baseline correction system with modular architecture. See [Baseline Reference](baseline_reference.md) for detailed documentation.

#### Quick Reference

```python
# Automatic model selection (recommended)
corrected, baseline, info = epyr.baseline.baseline_auto_1d(x, y, params)

# Specific methods
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params, order=3)
corrected, baseline = epyr.baseline.baseline_stretched_exponential_1d(x, y, params)
corrected, baseline = epyr.baseline.baseline_bi_exponential_1d(x, y, params)

# Backend control
epyr.baseline.setup_inline_backend()        # Static plots
epyr.baseline.setup_widget_backend()        # Interactive plots
```

#### Key Functions

- **`baseline_auto_1d()`** - Automatic model selection using AIC/BIC/R² criteria
- **`baseline_polynomial_1d()`** - Polynomial correction for CW EPR spectra
- **`baseline_stretched_exponential_1d()`** - For T2 relaxation data (β = 0.01-5.0)
- **`baseline_bi_exponential_1d()`** - For complex decay with multiple components
- **`baseline_polynomial_2d()`** - 2D surface fitting
- **`RegionSelector`** - Interactive region selection
- **`setup_*_backend()`** - Matplotlib backend control

#### Package Structure

- `models.py` - Mathematical functions
- `correction.py` - Core algorithms  
- `selection.py` - Region utilities
- `interactive.py` - Matplotlib widgets
- `auto.py` - Model selection

## Configuration Module

### epyr.config

Configuration management system.

#### EPyRConfig Class

```python
class EPyRConfig:
    def get(self, key: str, default: Any = None) -> Any
    def set(self, key: str, value: Any) -> None
    def get_section(self, section: str) -> Dict[str, Any]
    def set_section(self, section: str, config_dict: Dict[str, Any]) -> None
    def reset_section(self, section: str) -> None
    def reset_all(self) -> None
    def save(self) -> None
    def export_config(self, file_path: Union[str, Path]) -> None
    def import_config(self, file_path: Union[str, Path]) -> None
```

#### Global Configuration Instance

```python
from epyr.config import config

# Get configuration values
dpi = config.get('plotting.dpi')
cache_enabled = config.get('performance.cache_enabled')

# Set configuration values
config.set('plotting.dpi', 300)
config.save()
```

#### Convenience Functions

```python
def get_plotting_config() -> Dict[str, Any]
def get_performance_config() -> Dict[str, Any]
def is_debug_mode() -> bool
def is_cache_enabled() -> bool
```

## FAIR Data Module

### epyr.fair

FAIR data conversion and validation.

#### Main Conversion Functions

```python
def convert_bruker_to_fair(
    input_file: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    formats: List[str] = ["csv", "json"],
    include_metadata: bool = True,
    scaling: str = "",
) -> bool
```

```python
def batch_convert_directory(
    input_directory: Union[str, Path],
    output_directory: Optional[Union[str, Path]] = None,
    file_extensions: List[str] = [".dsc", ".spc", ".par"],
    scaling: str = "",
    output_formats: List[str] = ["csv_json", "hdf5"],
    recursive: bool = False,
) -> None
```

#### Validation Functions

```python
def validate_fair_dataset(
    data_dict: Dict[str, Any], 
    file_path: Optional[Path] = None
) -> ValidationResult

def validate_fair_metadata(metadata: Dict[str, Any]) -> ValidationResult

def validate_data_integrity(
    x_data: Optional[np.ndarray], 
    y_data: np.ndarray, 
    metadata: Dict[str, Any]
) -> ValidationResult

def create_validation_report(
    result: ValidationResult, 
    output_path: Optional[Path] = None
) -> str
```

#### ValidationResult Class

```python
class ValidationResult:
    is_valid: bool
    errors: List[str]
    warnings: List[str]
    info: List[str]
    
    def add_error(self, message: str) -> None
    def add_warning(self, message: str) -> None
    def add_info(self, message: str) -> None
    def get_summary(self) -> Dict[str, Any]
```

## Performance Module

### epyr.performance

Performance optimization and memory management.

#### DataCache Class

```python
class DataCache:
    def __init__(self, max_size_mb: Optional[int] = None)
    def get(self, file_path: Path) -> Optional[Tuple]
    def put(self, file_path: Path, data: Tuple) -> None
    def clear(self) -> None
    def get_stats(self) -> Dict[str, Any]
```

#### OptimizedLoader Class

```python
class OptimizedLoader:
    def __init__(
        self, 
        chunk_size_mb: Optional[int] = None, 
        cache_enabled: bool = True
    )
    
    def load_epr_file(self, file_path: Union[str, Path]) -> Tuple
    def load_chunked_data(
        self, 
        file_path: Union[str, Path], 
        chunk_processor: Callable
    ) -> Any
```

#### MemoryMonitor Class

```python
class MemoryMonitor:
    @staticmethod
    def get_memory_info() -> Dict[str, float]
    
    @staticmethod
    def check_memory_limit() -> bool
    
    @staticmethod
    def optimize_memory() -> None
```

#### Utility Functions

```python
def get_performance_info() -> Dict[str, Any]
def optimize_numpy_operations() -> None
def get_global_cache() -> DataCache
```

## Plugin System

### epyr.plugins

Extensible plugin architecture.

#### Base Plugin Classes

```python
class BasePlugin(ABC):
    plugin_name: str
    plugin_version: str
    
    @abstractmethod
    def initialize(self) -> bool
    
    def cleanup(self) -> None
    def get_info(self) -> Dict[str, Any]
```

```python
class FileFormatPlugin(BasePlugin):
    format_name: str
    file_extensions: List[str]
    supports_loading: bool
    supports_saving: bool
    
    @abstractmethod
    def can_load(self, file_path: Path) -> bool
    
    @abstractmethod
    def load(self, file_path: Path, **kwargs) -> Tuple[
        Optional[np.ndarray], 
        Optional[Union[np.ndarray, List[np.ndarray]]], 
        Optional[Dict[str, Any]]
    ]
```

```python
class ProcessingPlugin(BasePlugin):
    processing_name: str
    input_requirements: List[str]
    output_types: List[str]
    
    @abstractmethod
    def process(self, data: Dict[str, Any], **kwargs) -> Dict[str, Any]
    
    def validate_input(self, data: Dict[str, Any]) -> bool
```

```python
class ExportPlugin(BasePlugin):
    export_format: str
    file_extension: str
    supports_metadata: bool
    
    @abstractmethod
    def export(
        self, 
        output_path: Path, 
        x_data: np.ndarray, 
        y_data: np.ndarray,
        parameters: Dict[str, Any], 
        **kwargs
    ) -> bool
```

#### PluginManager Class

```python
class PluginManager:
    def register_plugin(self, plugin: BasePlugin) -> bool
    def unregister_plugin(self, plugin_name: str) -> bool
    def discover_plugins(
        self, 
        plugin_directories: Optional[List[Path]] = None
    ) -> int
    
    def get_file_format_plugin(self, file_path: Path) -> Optional[FileFormatPlugin]
    def get_export_plugin(self, format_name: str) -> Optional[ExportPlugin]
    def get_processing_plugin(self, processing_name: str) -> Optional[ProcessingPlugin]
    
    def list_plugins(self) -> Dict[str, List[Dict[str, Any]]]
    def get_supported_extensions(self) -> List[str]
```

#### Global Plugin Manager

```python
from epyr.plugins import plugin_manager

# Register custom plugin
plugin_manager.register_plugin(MyCustomPlugin())

# Find plugins
file_plugin = plugin_manager.get_file_format_plugin(Path('data.xyz'))
export_plugin = plugin_manager.get_export_plugin('csv')
```

## Logging Module

### epyr.logging_config

Centralized logging configuration.

```python
def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_string: Optional[str] = None
) -> logging.Logger

def get_logger(name: str) -> logging.Logger
```

**Example:**
```python
from epyr.logging_config import setup_logging, get_logger

# Set up logging
setup_logging('DEBUG', log_file=Path('epyr.log'))

# Get module logger
logger = get_logger(__name__)
logger.info("Processing started")
```

## CLI Module

### epyr.cli

Command-line interface functions.

```python
def cmd_convert() -> None
def cmd_baseline() -> None  
def cmd_batch_convert() -> None
def cmd_config() -> None
def cmd_info() -> None
def cmd_isotopes() -> None
def cmd_validate() -> None
def main() -> None
```

These functions are primarily used internally by the CLI system but can be called programmatically if needed.

## Baseline Correction Module

### epyr.baseline

Baseline correction algorithms.

```python
def baseline_polynomial(
    y_data: np.ndarray,
    x_data: Optional[np.ndarray] = None,
    poly_order: int = 1,
    exclude_regions: Optional[List[Tuple[float, float]]] = None
) -> Tuple[np.ndarray, np.ndarray]
```

**Parameters:**
- `y_data`: EPR intensity data
- `x_data`: Magnetic field values (optional)
- `poly_order`: Polynomial order for fitting
- `exclude_regions`: List of (start, end) regions to exclude from fit

**Returns:**
- `corrected_data`: Baseline-corrected intensity data
- `baseline`: Fitted baseline data

## Constants Module

### epyr.constants

Physical constants and unit conversions for EPR.

```python
# Fundamental constants
ELECTRON_G_FACTOR: float
BOHR_MAGNETON: float
PLANCK_CONSTANT: float

# Unit conversions
def gauss_to_tesla(value: float) -> float
def tesla_to_gauss(value: float) -> float
def frequency_to_field(frequency: float, g_factor: float = 2.0023) -> float
def field_to_frequency(field: float, g_factor: float = 2.0023) -> float
```

## Isotope GUI Module

### epyr.isotope_gui

Interactive isotope database interface.

```python
def run_gui() -> None
```

Launches the Tkinter-based isotope database GUI for exploring nuclear properties and magnetic parameters.

## Error Handling

All modules use consistent error handling:

```python
from epyr.logging_config import get_logger

logger = get_logger(__name__)

try:
    # EPyR Tools operations
    result = epyr.eprload('file.dsc')
except Exception as e:
    logger.error(f"Operation failed: {e}")
    raise
```

## Type Hints

EPyR Tools uses comprehensive type hints:

```python
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path
import numpy as np

# Example function signature
def process_data(
    x_data: np.ndarray,
    y_data: np.ndarray,
    parameters: Dict[str, Any],
    output_path: Optional[Path] = None
) -> Tuple[np.ndarray, Dict[str, Any]]:
    pass
```

## Integration Examples

### Complete Workflow Example

```python
import epyr
from epyr.fair import convert_bruker_to_fair, validate_fair_dataset
from epyr.performance import OptimizedLoader
from epyr.config import config

# Configure EPyR Tools
config.set('performance.cache_enabled', True)
config.set('plotting.dpi', 300)

# Load data with optimization
loader = OptimizedLoader(cache_enabled=True)
x, y, params, path = loader.load_epr_file('spectrum.dsc')

# Validate FAIR compliance
data_dict = {'x_data': x, 'y_data': y, 'metadata': params}
validation_result = validate_fair_dataset(data_dict)

if validation_result.is_valid:
    # Convert to FAIR formats
    success = convert_bruker_to_fair(
        'spectrum.dsc',
        formats=['csv', 'json', 'hdf5'],
        include_metadata=True
    )
    print(f"Conversion successful: {success}")
else:
    print("FAIR validation failed:")
    for error in validation_result.errors:
        print(f"  Error: {error}")
```

This API reference provides comprehensive documentation for all public interfaces in EPyR Tools.