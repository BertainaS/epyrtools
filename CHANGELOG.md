# Changelog

All notable changes to EPyR Tools will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.8] - 2025-09-14

### Added
- **Modular Baseline Correction Package**: Complete refactoring of baseline correction system
  - New `epyr.baseline` package with 5 specialized modules (`models.py`, `correction.py`, `selection.py`, `interactive.py`, `auto.py`)
  - 29+ functions for comprehensive baseline correction capabilities
  - Mathematical models: polynomial, stretched exponential, bi-exponential, simple exponential
- **Advanced Baseline Methods**: New correction algorithms for complex EPR data
  - `baseline_stretched_exponential_1d()` - For T2 relaxation and echo decay measurements (β = 0.01-5.0)
  - `baseline_bi_exponential_1d()` - For complex decay with multiple components
  - `baseline_auto_1d()` - Intelligent automatic model selection using AIC/BIC/R² criteria
- **Backend Control System**: User-controllable matplotlib backend selection
  - `setup_inline_backend()` - Static plots in Jupyter
  - `setup_widget_backend()` - Interactive plots in Jupyter
  - `setup_notebook_backend()` - Alternative interactive backend
- **Enhanced Interactive Selection**: Improved region selection with Jupyter compatibility
  - Cross-platform matplotlib widget support (handles version differences)
  - Multiple window closure methods (keyboard, function calls)
  - Better error handling and user guidance

### Changed
- **Architecture**: Refactored 1357-line `baseline_correction.py` into modular 5-file package structure
- **Default Behavior**: Backend selection now defaults to `'manual'` mode (user choice preserved)
- **API Design**: Clean, organized imports with comprehensive `__all__` declarations
- **Documentation**: Enhanced help system with backend control guidance and usage examples

### Improved
- **Performance**: Specialized modules for better maintainability and testing
- **Extensibility**: Easy addition of new baseline models and correction algorithms
- **User Experience**: No more forced backend changes on import - respects user preferences
- **Code Quality**: Separation of concerns (mathematical models, algorithms, UI components)

### Fixed
- **Jupyter Backend Issue**: EPyR Tools no longer overrides user's preferred matplotlib backend
- **Import Conflicts**: Resolved namespace issues between old and new baseline systems
- **Interactive Selection**: Fixed matplotlib widget compatibility across versions (props vs rectprops)
- **Backward Compatibility**: All existing baseline functions continue to work without changes

## [0.1.7] - 2025-09-14

### Added
- **EPR Plot Module**: New `eprplot` module with simple plotting functions
  - `plot_1d()` - Plot 1D EPR spectra with automatic axis detection
  - `plot_2d_map()` - Plot 2D EPR data as color maps
  - `plot_2d_waterfall()` - Plot 2D EPR data as waterfall plots
- **Enhanced Import System**: Fixed module visibility issues in Jupyter notebooks
- **Comprehensive Test Suite**: Added interactive Jupyter notebook for testing plot functions

### Changed
- **Import Architecture**: Improved `__init__.py` structure to prevent import conflicts
- **Plot Module Organization**: Added `__all__` variables for cleaner namespace management
- **Jupyter Compatibility**: Resolved module caching issues for better development experience

### Fixed
- **Module Visibility**: Resolved `AttributeError` when importing `eprplot` functions in Jupyter
- **Import Conflicts**: Fixed namespace conflicts between `plot.py` and `eprplot.py` modules
- **Automatic Imports**: Ensured all plotting functions are available via `import epyr`

## [0.1.6] - 2025-09-12

### Added
- **Comprehensive Testing Protocol**: 4-level testing framework (SMOKE, STANDARD, DEEP, SCIENTIFIC)
- **Deep Testing Infrastructure**: Automated test runners with performance benchmarking
- **Scientific Validation**: Mathematical accuracy testing against NIST standards
- **Complete Lineshape Analysis**: Extensive testing of all lineshape functions
- **Performance Benchmarking**: Speed and memory usage validation for core functions

### Changed
- **BREAKING**: Removed all numbered example notebooks for cleaner project structure
- **Documentation Cleanup**: Removed all emojis from documentation files for professional appearance
- **Version Update**: Updated to v0.1.6 across all configuration files and documentation
- **Enhanced Testing**: Added comprehensive testing documentation and protocols

### Removed
- Numbered example notebooks (01_, 04_, 05_, 06_, 07_, 08_, 09_, 10_)
- Emojis from all markdown and RST documentation files

### Fixed
- Fixed pseudo_voigt parameter handling in lineshape_class
- Corrected voigtian function parameter structure
- Updated all version references throughout the project

## [0.1.3] - 2025-09-06

### Added
- **Complete CLI System**: 8 professional command-line tools
  - `epyr-convert`: Bruker to FAIR format conversion
  - `epyr-baseline`: Baseline correction with multiple algorithms
  - `epyr-batch-convert`: High-throughput batch processing
  - `epyr-config`: Configuration management with import/export
  - `epyr-info`: System information and diagnostics
  - `epyr-isotopes`: Interactive isotope database GUI
  - `epyr-validate`: Data validation with FAIR compliance checking
  - `epyr`: Main CLI entry point with subcommands

- **FAIR Data Standards Compliance**
  - Comprehensive metadata validation system
  - Data integrity checking with detailed reports
  - EPR-specific parameter validation
  - File format validation for CSV, JSON, HDF5
  - Automated compliance reporting

- **Plugin Architecture**
  - Extensible plugin system for file formats, processing, and export
  - Base classes for FileFormatPlugin, ProcessingPlugin, ExportPlugin
  - Auto-discovery of plugins from user and system directories
  - Built-in CSV export plugin with metadata support

- **Performance Optimization System**
  - Intelligent memory monitoring and optimization
  - LRU data caching with configurable size limits
  - Optimized data loader with chunked processing support
  - NumPy operations optimization with MKL support

- **Comprehensive Configuration Management**
  - Hierarchical configuration with 8 main sections
  - Environment variable overrides (EPYR_* prefix)
  - User and system configuration file support
  - Configuration export/import functionality

- **Professional Documentation**
  - Complete User Guide (400+ lines) with CLI tutorials
  - Comprehensive API Reference with examples
  - Troubleshooting guide and best practices
  - Installation verification procedures

- **Development Infrastructure**
  - Complete testing suite with 90+ tests
  - Pre-commit hooks with Black, isort, flake8, mypy, bandit
  - Professional Makefile with 40+ development commands
  - Code quality enforcement with security scanning

### Enhanced
- **Core Data Loading**: Improved error handling and logging
- **Baseline Correction**: Integration with CLI system
- **Package Structure**: Modular architecture with clear separation
- **GUI Modernization**: Reorganized isotope GUI into proper module

### Changed
- **Package Status**: Upgraded from Alpha to Beta (Development Status :: 4 - Beta)
- **Python Support**: Added Python 3.12 support
- **Dependencies**: Updated to latest versions with comprehensive dev dependencies
- **Configuration**: Centralized configuration system replacing scattered settings

### Developer Experience
- **Quality Tools**: Black, isort, flake8, mypy, bandit, pydocstyle integration
- **Testing**: Pytest with coverage reporting and benchmark support
- **Documentation**: Automatic API documentation generation
- **CI/CD**: Complete pipeline simulation with `make ci`

## [0.1.2] - 2025-09-05

### Removed
- **BREAKING:** Removed `epyr.sub.baseline2.py` - deprecated duplicate baseline functions
- **BREAKING:** Removed `epyr.sub.processing2.py` - deprecated duplicate processing functions
- Cleaned up duplicate code and imports

### Changed
- Updated package imports to remove references to deleted modules
- All baseline correction functions now available through `epyr.baseline` module
- Streamlined package structure for better maintainability

### Fixed
- Fixed import issues in Getting Started notebook
- Consolidated all data files into single `examples/data/` directory
- Fixed complex data handling in notebooks
- Updated path detection for cross-platform compatibility

### Documentation
- Updated README with version badge
- Created comprehensive Getting Started notebook with real data examples
- Added proper error handling and troubleshooting in notebook
- Updated all version references

### Migration Guide
If you were using the removed modules:
```python
# OLD (no longer works)
from epyr.sub.baseline2 import baseline_polynomial
from epyr.sub.processing2 import baseline_polynomial

# NEW (use instead)
from epyr.baseline import baseline_polynomial
```

## [0.1.1] - 2025-09-04

### Added
- Comprehensive README with professional documentation
- Setup.py for pip package installation
- Example notebooks and tutorials
- FAIR data conversion capabilities
- Advanced plotting functionality

### Fixed
- Various import and compatibility issues
- Documentation generation
- Test coverage improvements

## [0.1.0] - 2025-09-01

### Added
- Initial release
- EPR data loading (BES3T and ESP formats)
- Basic baseline correction
- Constants and physical parameters
- Isotope GUI application
- Basic plotting capabilities
