# Baseline Correction System Update

## New Streamlined Baseline Correction System (v2.0.0)

The EPyR Tools baseline correction system has been modernized and simplified:

### **What Changed**
- **Streamlined Functions**: Using `scipy.optimize.curve_fit` for better performance
- **Interactive Region Selection**: GUI-based region selection for complex spectra
- **Direct eprload() Compatibility**: Works directly with `eprload()` data format
- **Clean API**: Simple function names like `baseline_polynomial_1d()` and `baseline_polynomial_2d()`

### **New Recommended Functions**
```python
import epyr

# Load EPR data
x, y, params, filepath = epyr.eprload("data.dsc")

# Simple automatic correction
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params)

# With manual region exclusion
signal_regions = [(3340, 3360), (3380, 3400)]
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, 
    manual_regions=signal_regions,
    region_mode='exclude'
)

# Interactive region selection
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params,
    interactive=True
)
```

### **Updated Notebooks**
- **`13_Manual_Region_Selection_Baseline.ipynb`** - Manual region selection examples
- **`14_Baseline_Correction_Real_Bruker_Data.ipynb`** - Real Bruker data with comprehensive examples
- **`06_Baseline_Correction_Functions_Complete.ipynb`** - Complete function reference

### **Compatibility**
- Old functions still work with deprecation warnings
- Legacy `baseline_polynomial()` functions are maintained
- All existing code continues to function

### **Performance Improvements**
- ~60% faster baseline fitting using scipy optimization
- Better memory usage for large datasets
- More robust fitting algorithms

---

**Recommendation**: Use the new `baseline_polynomial_1d()` and `baseline_polynomial_2d()` functions for new projects.