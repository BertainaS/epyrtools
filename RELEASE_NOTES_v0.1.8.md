# EPyR Tools v0.1.8 Release Notes

**Release Date:** September 14, 2025

## 🎉 Major Baseline Correction Refactoring

Version 0.1.8 introduces a comprehensive refactoring of the baseline correction system, transforming it from a monolithic 1357-line file into a clean, modular package with advanced capabilities.

## 🚀 Key Features

### 📦 Modular Architecture
- **New Package Structure**: `epyr.baseline/` with 5 specialized modules
  - `models.py` - Mathematical functions
  - `correction.py` - Core algorithms
  - `selection.py` - Region utilities  
  - `interactive.py` - Matplotlib widgets
  - `auto.py` - Model selection
- **34+ Functions**: Comprehensive baseline correction capabilities
- **Clean API**: Well-organized imports with `__all__` declarations

### 🧪 Advanced Baseline Methods
- **Stretched Exponential**: `baseline_stretched_exponential_1d()`
  - For T2 relaxation and echo decay measurements
  - Configurable β parameter (0.01-5.0)
  - Smart parameter initialization and uncertainty estimation
- **Bi-exponential**: `baseline_bi_exponential_1d()`
  - For complex decay with multiple components
  - Automatic component separation (τ₂/τ₁ ratio constraint)
- **Automatic Selection**: `baseline_auto_1d()`
  - Intelligent model choice using AIC/BIC/R² criteria
  - Tests polynomial, stretched exponential, and bi-exponential models
  - Detailed model comparison and fit quality metrics

### 🎨 Backend Control System
- **User Choice Preserved**: No more forced matplotlib backend changes
- **Convenience Functions**:
  - `epyr.setup_inline_backend()` - Static plots
  - `epyr.setup_widget_backend()` - Interactive plots
  - `epyr.setup_notebook_backend()` - Alternative interactive
- **Smart Detection**: Automatic Jupyter environment detection
- **Default Mode**: `'manual'` - respects user's existing backend

### 🖱️ Enhanced Interactive Selection
- **Cross-platform Compatibility**: Handles matplotlib version differences
- **Multiple Closure Methods**: Keyboard events, function calls, force close
- **Better Error Handling**: Clear user guidance and troubleshooting
- **Jupyter Optimized**: Works reliably in JupyterLab and Notebook

## 🔧 Technical Improvements

### Architecture
- **Separation of Concerns**: Mathematical models separate from UI components
- **Extensibility**: Easy to add new baseline models and algorithms
- **Maintainability**: Each module has focused responsibility
- **Testing**: Comprehensive test suite with 100+ test cases

### Performance
- **Optimized Algorithms**: Enhanced fitting procedures
- **Smart Initialization**: Automatic parameter guessing for robust fitting
- **Memory Efficient**: Better handling of large datasets
- **Robust Fitting**: Improved error handling and edge case management

### Documentation
- **Comprehensive Reference**: New `baseline_reference.md` with complete API
- **Enhanced Help**: Built-in help system with usage examples
- **Migration Guide**: Backward compatibility and upgrade instructions
- **Best Practices**: Workflow recommendations and troubleshooting

## 🔄 Backward Compatibility

All existing baseline correction code continues to work without changes:

```python
# Old code still works
corrected, baseline = epyr.baseline_polynomial(x, y, params)

# New recommended usage
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params)
```

## 📋 Usage Examples

### Quick Start
```python
import epyr

# Load EPR data
x, y, params, _ = epyr.eprload("data.dsc")

# Automatic baseline correction (recommended)
corrected, baseline, info = epyr.baseline.baseline_auto_1d(x, y, params)
print(f"Best model: {info['best_model']}")
```

### Backend Control
```python
# Choose your preferred backend
epyr.setup_inline_backend()        # Static plots in Jupyter
epyr.setup_widget_backend()        # Interactive plots in Jupyter

# Then use baseline correction normally
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params)
```

### Advanced Methods
```python
# T2 relaxation data
corrected, baseline = epyr.baseline.baseline_stretched_exponential_1d(
    x, y, params, beta_range=(0.5, 2.0), use_real_part=True
)

# Complex decay systems
corrected, baseline = epyr.baseline.baseline_bi_exponential_1d(
    x, y, params, tau_ratio_min=3.0
)

# Interactive region selection
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, interactive=True
)
```

## 🧪 Testing

### New Test Suite
- **Comprehensive Coverage**: `test_baseline_refactored.py` with 100+ tests
- **All Components**: Mathematical models, correction algorithms, region selection
- **Error Handling**: Edge cases and invalid input testing
- **Backend Testing**: Matplotlib backend control validation
- **Integration Testing**: End-to-end workflow validation

### Test Organization
- **Archived Old Tests**: Previous tests moved to `tests/archived_old_baseline/`
- **Clean Structure**: Focused tests for each module component
- **Synthetic Data**: Controlled test cases with known solutions
- **Real Data Testing**: Validation with actual EPR datasets

## 📚 Documentation Updates

### New Documentation
- **`baseline_reference.md`**: Complete baseline correction reference
- **Enhanced `api_reference.md`**: Updated with baseline package info
- **Updated `CHANGELOG.md`**: Detailed release notes
- **`RELEASE_NOTES_v0.1.8.md`**: This comprehensive summary

### Updated Documentation
- **Version Numbers**: All references updated to 0.1.8
- **Function References**: New function names and capabilities
- **Usage Examples**: Updated with new API and best practices
- **Migration Guide**: Instructions for upgrading from old system

## 🔧 Version Updates

All version references updated across the project:
- `epyr/__init__.py`: `__version__ = "0.1.8"`
- `pyproject.toml`: `version = "0.1.8"`
- `docs/conf.py`: `release = "0.1.8"`
- `README.md`: Version badge updated
- `epyr/lineshapes/__init__.py`: Updated to 0.1.8

## 🎯 Migration Checklist

For users upgrading from v0.1.7:

### ✅ Automatic (No Action Required)
- All existing baseline correction code continues to work
- Backward compatibility maintained through aliases
- Import statements remain the same

### 🔧 Recommended Updates
- Replace `epyr.baseline_auto()` with `epyr.baseline.baseline_auto_1d()`
- Use `epyr.setup_inline_backend()` for backend control
- Explore new advanced methods for T2/relaxation data
- Update documentation references to new function names

### 🆕 New Capabilities to Explore
- Try automatic model selection for unknown baseline types
- Use stretched exponential correction for T2 data
- Experiment with interactive region selection
- Configure backend preferences in your Jupyter workflow

## 🐛 Known Issues

### Resolved in v0.1.8
- ✅ **Jupyter Backend Override**: EPyR no longer forces backend changes
- ✅ **Interactive Selection Compatibility**: Fixed matplotlib version issues
- ✅ **Import Conflicts**: Resolved namespace issues between old/new systems
- ✅ **Widget Window Closure**: Multiple methods to close stuck windows

### Current Limitations
- Interactive selection requires compatible matplotlib backend in Jupyter
- Bi-exponential fitting can be challenging for noisy data (automatic fallback)
- 2D baseline correction limited to polynomial models (exponential methods planned)

## 🛠️ Development

### For Contributors
- **Modular Structure**: Each module can be developed/tested independently
- **Clear API**: Well-defined interfaces between modules
- **Comprehensive Tests**: Test suite covers all functionality
- **Documentation**: Inline docs and external reference materials

### Future Enhancements (Planned)
- 2D stretched exponential baseline correction
- Additional mathematical models (power law, logarithmic)
- Machine learning-based automatic selection
- GPU acceleration for large datasets
- Integration with other EPR analysis packages

## 📞 Support

### Getting Help
```python
# Built-in help system
epyr.baseline.get_help()              # Package overview
epyr.baseline.jupyter_help()          # Jupyter-specific help
help(epyr.baseline.baseline_auto_1d)  # Function-specific help
```

### Resources
- **Documentation**: `/docs/baseline_reference.md`
- **Examples**: `/examples/notebooks/16_Complete_Baseline_Correction_Demonstration.ipynb`
- **Tests**: `/tests/test_baseline_refactored.py`
- **Issues**: GitHub issue tracker for bugs and feature requests

---

## 🎉 Conclusion

EPyR Tools v0.1.8 represents a major step forward in baseline correction capabilities, providing a clean, modular architecture with advanced features while maintaining full backward compatibility. The new system is designed for both ease of use and extensibility, making it suitable for everything from quick baseline corrections to advanced research applications.

**Upgrade today and enjoy the enhanced baseline correction experience!** 🚀

---

*EPyR Tools Development Team*  
*September 14, 2025*