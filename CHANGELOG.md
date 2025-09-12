# Changelog

All notable changes to EPyR Tools will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
