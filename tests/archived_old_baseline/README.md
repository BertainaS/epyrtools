# Archived Baseline Tests

This directory contains the old baseline correction tests that were used before the baseline package refactoring in v0.1.8.

## Files Archived

- `test_baseline.py` - Original baseline tests for the old `baseline_correction.py` module
- `test_baseline_advanced.py` - Advanced baseline tests for the old system
- `test_baseline_comparison.py` - Baseline comparison tests for the old system

## Replacement

These tests have been replaced by:
- `../test_baseline_refactored.py` - Comprehensive tests for the new modular `epyr.baseline` package

## Archive Date

September 14, 2025 - EPyR Tools v0.1.8 release

## Why Archived

The baseline correction system was completely refactored from a single 1357-line file (`baseline_correction.py`) into a modular package structure (`epyr.baseline/`) with 5 specialized modules. The old tests were designed for the old API and architecture, so new comprehensive tests were created for the refactored system.

The new tests in `test_baseline_refactored.py` provide:
- Complete coverage of the new modular architecture
- Tests for new advanced features (stretched exponential, bi-exponential, automatic selection)
- Tests for backend control functionality
- Tests for the new region selection system
- Backward compatibility testing

## Migration Note

If you need to reference the old test patterns, they are preserved here, but the new `test_baseline_refactored.py` should be used for all baseline correction testing going forward.