# EPyR Tools Documentation

Welcome to the comprehensive documentation for EPyR Tools - a professional-grade Python package for Electron Paramagnetic Resonance (EPR) spectroscopy data analysis.

## Documentation Overview

This documentation covers all aspects of EPyR Tools, from basic usage to advanced development topics.

### Quick Start
- **[Installation Guide](../README.md#installation)** - Get EPyR Tools installed and running
- **[User Guide](user_guide.md)** - Complete guide to using all EPyR Tools features
- **[CLI Reference](cli_reference.md)** - Command-line interface documentation

### Core Systems  
- **[API Reference](api_reference.md)** - Complete API documentation with examples
- **[Configuration System](configuration_system.md)** - Comprehensive configuration management
- **[Performance System](performance_system.md)** - Memory optimization and caching
- **[Plugin System](plugin_system.md)** - Extensible plugin architecture

### Advanced Topics
- **[FAIR Data Standards](../README.md#fair-data-standards)** - Scientific data compliance
- **[Development Guide](#development-guide)** - Contributing to EPyR Tools
- **[Testing Framework](#testing-framework)** - Quality assurance and testing

## Architecture Overview

EPyR Tools is built with a modular, professional-grade architecture:

```
EPyR Tools Architecture
├── Core Data Loading (epyr.eprload)
├── Command Line Interface (epyr.cli) 
├── Configuration Management (epyr.config)
├── Performance Optimization (epyr.performance)
├── Plugin System (epyr.plugins)
├── FAIR Data Standards (epyr.fair)
├── Baseline Correction (epyr.baseline)
├── Visualization (epyr.plot)
├── GUI Applications (epyr.isotope_gui)
└── Utilities (epyr.constants, epyr.sub)
```

## Key Features by Module

### Command Line Interface (`epyr.cli`)
Professional CLI with 8 commands for all EPR workflows:
- `epyr-convert` - Bruker to FAIR format conversion
- `epyr-baseline` - Baseline correction with multiple algorithms
- `epyr-batch-convert` - High-throughput batch processing
- `epyr-config` - Configuration management
- `epyr-info` - System diagnostics
- `epyr-isotopes` - Interactive isotope database
- `epyr-validate` - Data validation with FAIR compliance

### Configuration System (`epyr.config`)
Centralized, hierarchical configuration with:
- 8 configuration sections covering all EPyR Tools features
- Environment variable overrides (`EPYR_*` prefix)
- Configuration export/import for reproducibility
- Platform-appropriate configuration file storage

### Performance System (`epyr.performance`)
Advanced performance optimization including:
- Intelligent memory monitoring and optimization
- LRU caching for frequently accessed data
- Optimized data loading for large datasets
- Multi-core processing support
- Resource usage tracking and alerts

### Plugin System (`epyr.plugins`)
Extensible architecture supporting:
- File format plugins for custom data formats
- Processing plugins for analysis algorithms
- Export plugins for specialized output formats
- Auto-discovery and hot-loading of plugins
- Professional plugin development framework

### FAIR Data Standards (`epyr.fair`)
Complete FAIR compliance implementation:
- Comprehensive metadata validation
- Data integrity checking
- Scientific data standards compliance
- Automated validation reporting
- Export to standard formats (CSV, JSON, HDF5)

## Quick Reference

### Command Line Usage
```bash
# Convert EPR data to FAIR formats
epyr-convert spectrum.dsc --formats csv,json,hdf5

# Apply baseline correction  
epyr-baseline spectrum.dsc --method polynomial --order 2

# Batch process multiple files
epyr-batch-convert ./data/ --jobs 4

# Validate data quality
epyr-validate *.dsc --detailed

# System information
epyr-info --all
```

### Python API Usage
```python
import epyr
from epyr.config import config
from epyr.performance import OptimizedLoader
from epyr.fair import validate_fair_dataset

# Configure EPyR Tools
config.set('performance.cache_enabled', True)
config.set('plotting.dpi', 300)

# Load data with optimization
loader = OptimizedLoader()
x, y, params, path = loader.load_epr_file('spectrum.dsc')

# Validate FAIR compliance
data_dict = {'x_data': x, 'y_data': y, 'metadata': params}
validation = validate_fair_dataset(data_dict)

# Convert to FAIR formats
from epyr.fair import convert_bruker_to_fair
convert_bruker_to_fair('spectrum.dsc', formats=['csv', 'json'])
```

### Configuration Management
```python
from epyr.config import config

# View current configuration
epyr-config show

# Set values
config.set('plotting.dpi', 300)
config.set('performance.cache_size_mb', 200)

# Export/import settings
config.export_config('my_settings.json')
config.import_config('my_settings.json')
```

## Documentation Structure

### For Users
1. Start with the **[User Guide](user_guide.md)** for comprehensive usage instructions
2. Reference the **[CLI Reference](cli_reference.md)** for command-line usage
3. Use the **[API Reference](api_reference.md)** for Python programming

### For Developers
1. Read the **[Configuration System](configuration_system.md)** to understand settings
2. Study the **[Performance System](performance_system.md)** for optimization
3. Explore the **[Plugin System](plugin_system.md)** for extensibility

### For Contributors
1. Review the **Development Guide** for contribution guidelines
2. Understand the **Testing Framework** for quality assurance
3. Follow the **Code Standards** for consistency

## Development Guide

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools

# Install with development dependencies
pip install -e ".[dev,docs]"

# Set up pre-commit hooks
pre-commit install

# Run tests
make test

# Check code quality
make quality

# Build documentation
make docs
```

### Development Workflow

1. **Code Quality**: Use provided pre-commit hooks (Black, isort, flake8, mypy, bandit)
2. **Testing**: Write comprehensive tests for new features
3. **Documentation**: Update documentation for API changes
4. **Performance**: Consider performance implications of changes

### Contributing Guidelines

- Follow existing code style and patterns
- Write tests for new functionality
- Update documentation for user-facing changes
- Use type hints for better code maintainability
- Consider backward compatibility

## Testing Framework

EPyR Tools includes a comprehensive testing suite:

### Test Structure
```
tests/
├── test_cli.py              # CLI functionality
├── test_config.py           # Configuration system
├── test_performance.py      # Performance optimization
├── test_plugins.py          # Plugin system
├── test_fair_validation.py  # FAIR validation
├── test_baseline.py         # Baseline correction
├── test_eprload.py          # Data loading
└── test_integration.py      # Integration tests
```

### Running Tests
```bash
# Run all tests
make test

# Run specific test module
pytest tests/test_cli.py -v

# Run with coverage
make test-cov

# Run performance benchmarks
make benchmark
```

### Test Coverage
EPyR Tools maintains high test coverage (>90%) across all core modules:
- CLI commands and argument handling
- Configuration management and persistence
- Performance optimization and caching
- Plugin system and auto-discovery
- FAIR validation and compliance
- Data loading and processing

## Performance Considerations

EPyR Tools is designed for high performance:

### Memory Management
- Intelligent caching system for frequently accessed files
- Memory monitoring and optimization
- Graceful degradation when memory limits are approached

### Large Dataset Handling  
- Chunked processing for files > 100MB
- Optimized data loading with configurable parameters
- Multi-core processing support for batch operations

### Caching Strategy
- LRU cache with configurable size limits
- File modification tracking for cache invalidation
- Per-operation cache statistics and monitoring

## Security Features

EPyR Tools implements security best practices:

### Input Validation
- Comprehensive validation of all input data
- Safe file handling practices
- Protection against common attack vectors

### Code Quality
- Automated security scanning with Bandit
- Type checking with MyPy
- Comprehensive linting with multiple tools

### Safe Defaults
- Secure configuration defaults
- Principle of least privilege
- Error handling without information disclosure

## Getting Help

### Documentation
- **User issues**: Check the [User Guide](user_guide.md) and [CLI Reference](cli_reference.md)
- **API questions**: See the [API Reference](api_reference.md)
- **Configuration**: Read the [Configuration System](configuration_system.md)

### Community Support
- **Bug reports**: [GitHub Issues](https://github.com/BertainaS/epyrtools/issues)
- **Feature requests**: [GitHub Discussions](https://github.com/BertainaS/epyrtools/discussions)
- **Questions**: Use GitHub Discussions or contact the maintainers

### System Diagnostics
```bash
# Get comprehensive system information
epyr-info --all

# Check configuration
epyr-config show

# Validate installation
epyr-validate --help
python -c "import epyr; print(f'EPyR Tools v{epyr.__version__} ready!')"
```

## Changelog and Versioning

EPyR Tools follows [Semantic Versioning](https://semver.org/):
- **Major versions** (1.0.0): Breaking changes
- **Minor versions** (0.2.0): New features, backward compatible
- **Patch versions** (0.1.1): Bug fixes, backward compatible

See [CHANGELOG.md](../CHANGELOG.md) for detailed version history.

## License and Citation

EPyR Tools is released under the MIT License. See [LICENSE](../LICENSE) for details.

**Author:** Sylvain Bertaina (sylvain.bertaina@cnrs.fr)
**Laboratory:** [Magnetism Group (MAG), IM2NP](https://www.im2np.fr/fr/equipe-magnetisme-mag)

### Citation
If you use EPyR Tools in your research, please cite:

```
EPyR Tools: Electron Paramagnetic Resonance Tools in Python
Bertaina, S. et al.
https://github.com/BertainaS/epyrtools
```

---

**EPyR Tools v0.1.6** - Professional-grade EPR data analysis for Python

This documentation provides comprehensive coverage of all EPyR Tools features. For specific questions or issues, please refer to the appropriate section above or contact the development team.