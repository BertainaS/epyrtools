# Changelog

All notable changes to EPyR Tools will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
