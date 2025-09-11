# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Essential Development Commands

### Testing
```bash
# Run all tests
make test

# Run tests with coverage report
make test-cov

# Run fast tests (exclude slow markers)  
make test-fast

# Run specific test categories
pytest tests/ -m "smoke"      # Basic functionality tests (< 1 minute)
pytest tests/ -m "standard"   # Comprehensive tests (< 5 minutes) 
pytest tests/ -m "deep"       # Edge cases and performance (< 15 minutes)
pytest tests/ -m "scientific" # Scientific accuracy validation (< 30 minutes)

# Run single test file
pytest tests/test_eprload.py -v

# Run integration tests
make test-integration
```

### Code Quality
```bash
# Run all quality checks
make quality

# Format code (black + isort)
make format

# Check formatting without changes
make format-check

# Run linting only
make lint

# Run type checking
make type-check

# Run security checks
make security

# Run pre-commit hooks
make pre-commit
```

### Development Setup
```bash
# Install in development mode with all dependencies
make install-dev

# Create and set up virtual environment
make dev-setup

# Update dependencies
make update-deps
```

### CLI Tools Testing
```bash
# Test main CLI commands
epyr --help
epyr-info
epyr-config show plotting
epyr-isotopes  # GUI tool
```

### Documentation
```bash
# Build documentation
make docs

# Serve documentation locally
make docs-serve  # Available at http://localhost:8000
```

## Package Architecture

### Core Data Loading (`epyr.eprload`)
- **Primary function**: `epyr.eprload()` - Universal EPR data loader with file dialog
- **Format support**: Automatic detection of BES3T (.dsc/.dta) and ESP (.spc/.par) formats
- **Returns**: `(x_data, y_data, parameters_dict, file_path)`
- **Performance optimized**: Uses `OptimizedLoader` and memory caching

### Sub-modules for Format-Specific Loading (`epyr.sub`)
- **`loadBES3T.py`**: Bruker BES3T/Xepr format loader (.dsc/.dta files)
- **`loadESP.py`**: Bruker ESP/WinEPR format loader (.par/.spc files)  
- **`utils.py`**: File handling utilities and format detection

### Analysis Modules
- **`baseline/`**: Multi-algorithm baseline correction (1D/2D)
  - `_1d.py`: Polynomial and exponential baseline correction
  - `_2d.py`: 2D spectrum baseline correction
  - `_utils.py`: Signal exclusion and fitting utilities

- **`lineshapes/`**: Complete EPR lineshape library
  - `gaussian.py`, `lorentzian.py`, `voigtian.py`: Core lineshape functions
  - `lineshape_class.py`: Unified `Lineshape` class interface
  - `convspec.py`: Convolution and spectrum simulation tools

### FAIR Data Conversion (`epyr.fair`)
- **Purpose**: Convert Bruker files to open formats (CSV, JSON, HDF5)
- **Key function**: `convert_bruker_to_fair()` with metadata preservation
- **Plugin system**: Extensible exporters via `plugin_manager`

### Advanced Features
- **`performance.py`**: Memory optimization, caching, and performance monitoring
- **`constants.py`**: EPR-specific physical constants and conversions
- **`isotope_gui/`**: Interactive nuclear isotope database (Tkinter GUI)
- **`plot.py`**: EPR-specific plotting functions (2D maps, publication quality)

## Testing Framework

### Test Categories (use pytest markers)
```bash
pytest -m "smoke"        # < 1 minute - basic imports and sanity checks
pytest -m "standard"     # < 5 minutes - comprehensive feature testing  
pytest -m "deep"         # < 15 minutes - edge cases and performance
pytest -m "scientific"   # < 30 minutes - scientific accuracy validation
```

### Key Test Files
- **`test_eprload.py`**: Data loading functionality across formats
- **`test_lineshapes.py`**: Mathematical validation of EPR lineshapes
- **`test_performance.py`**: Memory usage and benchmark tests
- **`test_integration.py`**: End-to-end workflow testing
- **`test_deep_protocol.py`**: Comprehensive edge case testing

### Performance Testing
```bash
# Run benchmarks
make benchmark

# Memory profiling  
make memory-test

# Profile code execution
make profile
```

## Configuration System

### User Configuration (`epyr.config`)
- **Config file**: `~/.epyr/config.ini`
- **CLI management**: `epyr-config show/set/reset`
- **Sections**: `[plotting]`, `[performance]`, `[logging]`

### Logging System
- **Setup**: `epyr.logging_config.setup_logging()`
- **Usage**: `logger = get_logger(__name__)`
- **Levels**: Controlled via config or environment variables

## Plugin Architecture

### Export Plugins (`epyr.plugins`)
- **Registration**: Automatic discovery via entry points
- **Built-in**: CSV exporter in `fair/exporters.py`
- **Custom**: Implement `BaseExporter` class and register

## Development Workflow

### Code Standards
- **Formatting**: Black (88 char line length) + isort
- **Linting**: flake8 with docstring checking
- **Type hints**: mypy (gradual typing, imports ignored)
- **Security**: bandit scanning
- **Pre-commit hooks**: Automatic enforcement

### Release Process
```bash
# Check version consistency
make check-version

# Full CI simulation locally  
make ci

# Build package
make build
```

### Memory and Performance Considerations
- **Large datasets**: Use `OptimizedLoader` with caching
- **2D data**: Special handling in BES3T loader for multi-dimensional arrays
- **Memory monitoring**: `MemoryMonitor` class for tracking usage

## Common Tasks

### Adding New EPR File Format
1. Create loader in `epyr/sub/load{FORMAT}.py`
2. Add format detection to `epyr/sub/utils.py`
3. Integrate into `epyr/eprload.py` dispatcher
4. Add tests in `tests/test_eprload.py`

### Adding New Lineshape Function  
1. Implement in `epyr/lineshapes/{name}.py`
2. Add to `lineshape_class.py` registry
3. Add mathematical validation in `tests/test_lineshapes.py`
4. Update documentation

### Adding CLI Command
1. Add function to `epyr/cli.py`
2. Register in `pyproject.toml` `[project.scripts]`
3. Add tests in `tests/test_cli.py`

## Data Files and Examples

### Test Data Location
- **Tests**: Small files in `tests/` directory
- **Examples**: Full datasets in `examples/data/`
- **Notebooks**: Tutorial notebooks in `examples/notebooks/`

### Example EPR Data Formats
- **BES3T**: `*.DSC` (metadata) + `*.DTA` (binary data)
- **ESP**: `*.par` (parameters) + `*.spc` (binary data)
- **Output**: 2D arrays for angular-dependent or pulse EPR measurements

## Scientific Validation

EPyR Tools implements rigorous scientific validation:
- **Mathematical accuracy**: Numerical integration and derivative tests
- **Physical constants**: NIST-validated EPR constants
- **Format compliance**: Strict adherence to Bruker file specifications
- **Performance**: Optimized for large EPR datasets (>1M points)