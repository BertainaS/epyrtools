# EPyR Tools v0.3.0 Release Notes

**Release Date:** November 2025
**Release Type:** Code Quality Enhancement

## Overview

Version 0.3.0 represents a major code quality improvement, transitioning EPyR Tools from print-based output to a professional, centralized logging system. This release migrates 95.7% of all print statements (398 out of 416) to structured logging across 18 core modules, providing better traceability, debugging capabilities, and production-ready monitoring.

This is a **backward-compatible** release with no breaking changes. All functionality is preserved while significantly improving code maintainability and professional standards.

## Major Improvements

### Professional Logging Infrastructure

EPyR Tools now uses Python's standard logging framework throughout the codebase, replacing direct print() calls with structured logging:

**Migration Statistics:**
- **398 print statements** migrated to logging (95.7% coverage)
- **18 core modules** refactored with logging infrastructure
- **Zero breaking changes** - all functionality preserved
- **100+ tests passing** - no regressions introduced

### Structured Logging Hierarchy

All logging now follows a consistent hierarchy for better control and filtering:

```python
from epyr.logging_config import get_logger

logger = get_logger(__name__)

# User-facing results and progress
logger.info("Loaded EPR data: 1024 points")

# Warnings for issues and errors
logger.warning("Missing parameter: XUNI, using default unit")

# Technical details for debugging
logger.debug("Initial fit parameters: center=350.0, width=2.5")
```

**Log Levels:**
- `DEBUG`: Technical details, initial parameters, internal operations
- `INFO`: User-facing messages, progress updates, results
- `WARNING`: Errors, issues, missing dependencies, failures
- `ERROR`: Critical errors requiring user attention

### Modules Updated

The following 18 modules have been migrated to professional logging:

**Core Loading & Plotting:**
1. `epyr/eprload.py` (22 prints) - Data loading module
2. `epyr/eprplot.py` (6 prints) - EPR plotting functions

**Signal Processing:**
3. `epyr/signalprocessing/frequency_analysis.py` (79 prints) - FFT analysis
4. `epyr/signalprocessing/apowin.py` (22 prints) - Apodization windows

**Physics & Units:**
5. `epyr/physics/conversions.py` (50 prints) - Unit conversions
6. `epyr/physics/units.py` (29 prints) - Unit utilities
7. `epyr/physics/constants.py` (14 prints) - Physical constants

**Baseline Correction:**
8. `epyr/baseline/__init__.py` (27 prints) - Package initialization
9. `epyr/baseline/correction.py` (24 prints) - Correction algorithms
10. `epyr/baseline/auto.py` (21 prints) - Automatic model selection
11. `epyr/baseline/interactive.py` (16 prints) - Interactive selection

**FAIR Data Conversion:**
12. `epyr/fair/conversion.py` (15 prints) - Format conversion
13. `epyr/fair/exporters.py` (3 prints) - Export functions
14. `epyr/fair/data_processing.py` (1 print) - Metadata processing

**Lineshape Analysis:**
15. `epyr/lineshapes/fitting.py` (5 prints) - EPR signal fitting
16. `epyr/lineshapes/convspec.py` (5 prints) - Spectrum convolution

**GUI & CLI:**
17. `epyr/cli.py` (45 prints) - Command-line interface
18. `epyr/isotope_gui.py` (14 prints) - Isotope database GUI

## Benefits

### For Users

**Better Control:**
```python
import logging

# Set logging level for all EPyR modules
logging.getLogger('epyr').setLevel(logging.WARNING)  # Only warnings and errors

# Or configure specific modules
logging.getLogger('epyr.frequency_analysis').setLevel(logging.DEBUG)
```

**Clean Output:**
- Timestamps on all messages for timing analysis
- Module names for tracing message sources
- Severity levels for filtering important messages
- Structured format for automated log parsing

**Production Ready:**
- Log to files for permanent records
- Integration with logging frameworks (Sentry, Logstash, etc.)
- Better debugging with detailed trace information
- Professional monitoring and alerting capabilities

### For Developers

**Maintainability:**
- Single logging configuration point (`logging_config.py`)
- Consistent message formatting across all modules
- Easy to add new log messages with appropriate levels
- Better code organization and clarity

**Debugging:**
- Full call chain tracing with module names
- Detailed technical information at DEBUG level
- Conditional logging without code changes (via config)
- Log file analysis for issue reproduction

## Backward Compatibility

### No Breaking Changes

This release is **100% backward compatible**:

✅ **API Unchanged:** All functions return the same values with same signatures
✅ **Functionality Preserved:** All operations work identically
✅ **Tests Passing:** 100+ tests confirm no regressions
✅ **Default Behavior:** Logging configured to match previous print() output

### Migration Notes

**No action required** for most users. The default logging configuration mimics previous print() behavior.

**Optional Configuration:**

If you want to customize logging:

```python
# In your script, before importing epyr
import logging

# Configure EPyR logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='epyr_analysis.log'
)

# Now import and use EPyR Tools
from epyr import eprload

x, y, params, filepath = eprload("data.DTA")
# Logging messages written to epyr_analysis.log
```

### CLI Commands

All CLI commands work identically:

```bash
# These commands now use structured logging
epyr-convert data.DTA output/ --format csv
epyr-baseline data.DTA --order 2
epyr-plot data.DTA --interactive

# Control verbosity with standard flags (planned for future release)
# epyr-info --debug    # Show DEBUG messages
# epyr-info --quiet    # Only WARNING and ERROR
```

## Code Quality Improvements

### Centralized Configuration

New `logging_config.py` module provides:
- Single source of logging configuration
- `get_logger(__name__)` factory for consistent logger creation
- Default formatters and handlers
- Environment-based configuration support

### Consistent Patterns

All modules follow the same pattern:

```python
# At top of module
from ..logging_config import get_logger

logger = get_logger(__name__)

# Throughout module
def my_function():
    logger.debug("Starting operation...")

    result = do_work()
    logger.info(f"Processed {len(result)} items")

    if issue_detected:
        logger.warning("Issue detected, using fallback method")

    return result
```

### Professional Standards

- PEP 8 compliant logging usage
- Proper log level selection based on message importance
- Structured messages with contextual information
- No emojis in log messages (professional appearance)
- Consistent timestamp and module information

## Testing & Validation

### Test Coverage

```bash
pytest tests/ -v
# Result: 100+ tests passed, 0 regressions
# Only 1 pre-existing test failure (unrelated to logging migration)
```

**Validated Functionality:**
- All data loading operations work correctly
- Signal processing produces identical results
- Baseline correction algorithms unchanged
- FAIR conversion preserves all metadata
- CLI commands execute successfully
- GUI applications launch properly

### Import Verification

All modules import cleanly with logging infrastructure:

```python
# These all work correctly
from epyr import eprload
from epyr.signalprocessing import analyze_frequencies
from epyr.baseline import baseline_polynomial
from epyr.fair import convert_to_fair
from epyr import isotopes

# Logging configured automatically on import
```

## Known Issues

**None.** This release introduces no new issues.

**Note:** One pre-existing test failure in `baseline/models.py` (numpy dtype casting) is unrelated to the logging migration and existed before v0.3.0.

## Performance

**No measurable performance impact.** Python's logging module is highly optimized:

- Logging calls are lazy-evaluated
- Disabled log levels have minimal overhead
- No impact on computational performance
- Memory usage unchanged

## Documentation Updates

- README.md updated with v0.3.0 version badge
- docs/changelog.rst comprehensive v0.3.0 entry
- Release notes (this document) created
- Version consistency across all files verified

## Installation

### Upgrade from v0.2.5 or earlier

```bash
pip install --upgrade epyr-tools
```

### Fresh Installation

```bash
pip install epyr-tools==0.3.0
```

### Development Installation

```bash
git clone https://github.com/BertainaS/epyrtools.git
cd epyrtools
git checkout v0.3.0
pip install -e ".[dev,docs]"
```

## Future Plans

**Planned for v0.3.1:**
- CLI verbosity flags (`--debug`, `--quiet`, `--verbose`)
- Log file rotation configuration options
- Enhanced logging documentation with examples

**Planned for v0.4.0:**
- Complete migration of remaining 18 print statements (4.3%)
- Advanced logging features (custom formatters, handlers)
- Performance profiling with logging integration

## Contributors

- **Sylvain Bertaina** - Lead Developer & Maintainer
  - Email: sylvain.bertaina@cnrs.fr
  - Affiliation: Magnetism Group (MAG), IM2NP Laboratory, CNRS
  - Contribution: Complete logging migration, testing, documentation

## Feedback & Support

- **Issues**: [GitHub Issues](https://github.com/BertainaS/epyrtools/issues)
- **Discussions**: [GitHub Discussions](https://github.com/BertainaS/epyrtools/discussions)
- **Email**: sylvain.bertaina@cnrs.fr

## Summary

EPyR Tools v0.3.0 brings the codebase to professional standards with comprehensive logging infrastructure. This release demonstrates commitment to code quality and maintainability while preserving complete backward compatibility. All 398 migrated print statements now provide better debugging, monitoring, and production deployment capabilities.

**Key Takeaway:** Same functionality, better infrastructure, professional standards.

---

**EPyR Tools v0.3.0** - Professional logging for modern EPR analysis
