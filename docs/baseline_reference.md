# Baseline Correction Reference

Complete reference for the EPyR Tools baseline correction system (v0.1.8).

## Overview

The `epyr.baseline` package provides comprehensive baseline correction capabilities for EPR spectroscopy data. The package is organized into specialized modules for maximum flexibility and maintainability.

### Package Structure

```
epyr.baseline/
├── __init__.py           # Main API with 29+ functions
├── models.py             # Mathematical functions
├── correction.py         # Core correction algorithms
├── selection.py          # Region selection utilities
├── interactive.py        # Interactive matplotlib widgets
└── auto.py              # Automatic model selection
```

## Quick Start

### Basic Usage

```python
import epyr

# Load EPR data
x, y, params, _ = epyr.eprload("data.dsc")

# Automatic baseline correction (recommended)
corrected, baseline, info = epyr.baseline.baseline_auto_1d(x, y, params)
print(f"Best model: {info['best_model']}")

# Or use specific methods
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params, order=3)
```

### Backend Control

```python
# Control matplotlib backend in Jupyter
epyr.setup_inline_backend()        # Static plots
epyr.setup_widget_backend()        # Interactive plots
epyr.setup_notebook_backend()      # Alternative interactive
```

## Core Functions

### Automatic Model Selection

#### `baseline_auto_1d(x, y, params, **kwargs)`

Intelligent automatic selection of the best baseline model using statistical criteria.

**Parameters:**
- `x`: X-axis data (field/time) or None
- `y`: Y-axis EPR data (1D array)
- `params`: Parameter dictionary from eprload()
- `models`: List of models to test `['polynomial', 'stretched_exponential', 'bi_exponential']`
- `selection_criterion`: `'aic'`, `'bic'`, or `'r2'`
- `use_real_part`: For complex data, use real part if True
- `verbose`: Print detailed comparison if True

**Returns:**
- `corrected_data`: Baseline-corrected EPR data
- `baseline`: Fitted baseline values
- `model_info`: Dictionary with model comparison results

**Example:**
```python
# Automatic selection with model comparison
corrected, baseline, info = epyr.baseline.baseline_auto_1d(
    x, y, params, 
    models=['polynomial', 'stretched_exponential'],
    selection_criterion='aic',
    verbose=True
)

print(f"Best model: {info['best_model']}")
print(f"R² = {info['parameters']['r2']:.4f}")
for model, aic in info['criteria'].items():
    print(f"{model}: AIC = {aic:.2f}")
```

### Polynomial Correction

#### `baseline_polynomial_1d(x, y, params, **kwargs)`

Polynomial baseline correction for smooth baseline drifts, ideal for CW EPR spectra.

**Parameters:**
- `order`: Polynomial order (1-4, default: 2)
- `exclude_center`: Exclude center region containing signal
- `center_fraction`: Fraction of data to exclude from center (default: 0.3)
- `manual_regions`: List of regions `[(x1, x2), ...]`
- `region_mode`: `'exclude'` or `'include'` manual regions
- `interactive`: Enable interactive region selection

**Example:**
```python
# Basic polynomial correction
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, order=3
)

# With manual region specification
regions = [(3340, 3360), (3380, 3400)]  # Signal regions to exclude
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, 
    manual_regions=regions,
    region_mode='exclude'
)

# Interactive selection
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, interactive=True
)
```

#### `baseline_polynomial_2d(x, y, params, **kwargs)`

2D polynomial surface fitting for 2D EPR datasets.

**Parameters:**
- `order`: Polynomial order `(order_x, order_y)` or int for both
- `use_real_part`: For complex 2D data
- `exclude_center`: Exclude center region
- `interactive`: Enable 2D region selection

### Exponential Methods

#### `baseline_stretched_exponential_1d(x, y, params, **kwargs)`

Stretched exponential correction for T2 relaxation and echo decay measurements.

**Mathematical Model:** `baseline = offset + A × exp(-(x/τ)^β)`

**Parameters:**
- `beta_range`: Range for stretching exponent (default: `(0.01, 5.0)`)
  - β = 1: Pure exponential decay
  - β < 1: Sub-exponential (slower than exponential)  
  - β > 1: Super-exponential (faster than exponential)
- `use_real_part`: For complex data, fit real part if True
- `exclude_initial`: Skip first N points
- `exclude_final`: Skip last N points
- `initial_guess`: Dictionary with parameter guesses

**Example:**
```python
# T2 relaxation data correction
corrected, baseline = epyr.baseline.baseline_stretched_exponential_1d(
    x, y, params,
    beta_range=(0.5, 2.0),
    use_real_part=True,
    exclude_initial=10
)

# With custom initial guess
initial_params = {'A': 1000, 'tau': 500, 'beta': 1.2, 'offset': 50}
corrected, baseline = epyr.baseline.baseline_stretched_exponential_1d(
    x, y, params, initial_guess=initial_params
)
```

#### `baseline_bi_exponential_1d(x, y, params, **kwargs)`

Bi-exponential correction for complex decay with multiple components.

**Mathematical Model:** `baseline = offset + A₁×exp(-x/τ₁) + A₂×exp(-x/τ₂)`

**Parameters:**
- `tau_ratio_min`: Minimum ratio τ₂/τ₁ for component separation (default: 2.5)
- `use_real_part`: For complex data
- `initial_guess`: Parameter guess dictionary

**Example:**
```python
# Complex decay with fast and slow components
corrected, baseline = epyr.baseline.baseline_bi_exponential_1d(
    x, y, params,
    tau_ratio_min=3.0,
    use_real_part=True
)
```

## Interactive Selection

### RegionSelector Class

Interactive region selection using matplotlib widgets.

```python
from epyr.baseline import RegionSelector

# 1D region selection
selector = RegionSelector()
regions = selector.select_regions_1d(x, y, "Select baseline regions")

# Use selected regions
corrected, baseline = epyr.baseline.baseline_polynomial_1d(
    x, y, params, 
    manual_regions=regions,
    region_mode='include'
)
```

### Utility Functions

#### `close_selector_window()`
Close stuck interactive selection windows.

```python
# If interactive selection gets stuck
epyr.baseline.close_selector_window()
```

#### `jupyter_help()`
Display comprehensive Jupyter notebook help.

```python
epyr.baseline.jupyter_help()
```

## Mathematical Models

### Available Models

Access mathematical functions directly:

```python
from epyr.baseline import models

# Stretched exponential function
y_stretched = models.stretched_exponential_1d(x, A=1000, tau=500, beta=1.2, offset=0)

# Bi-exponential function  
y_bi = models.bi_exponential_1d(x, A1=500, tau1=100, A2=300, tau2=800, offset=0)

# List available models
available_models = models.list_available_models()
print(available_models)  # ['polynomial', 'exponential', 'stretched_exponential', 'bi_exponential']
```

### Model Information

```python
# Get model information
model_info = models.MODEL_INFO['stretched_exponential']
print(model_info['description'])  # "Stretched exponential decay (KWW function)"
print(model_info['typical_use'])  # "T2 relaxation, echo decay measurements"
```

## Backend Control

### Matplotlib Backend Management

Control matplotlib backend behavior in Jupyter notebooks:

```python
# Static plots (inline backend)
epyr.baseline.setup_inline_backend()
# Equivalent to: %matplotlib inline

# Interactive plots (widget backend) 
epyr.baseline.setup_widget_backend()
# Equivalent to: %matplotlib widget

# Alternative interactive (notebook backend)
epyr.baseline.setup_notebook_backend()
# Equivalent to: %matplotlib notebook

# Check if interactive selection will work
if epyr.baseline.is_interactive_available():
    # Use interactive selection
    corrected, baseline = epyr.baseline.baseline_polynomial_1d(
        x, y, params, interactive=True
    )
```

## Configuration

### Package Settings

```python
# Configure default settings
epyr.baseline.configure(
    polynomial_order=3,
    beta_range=(0.01, 3.0),
    selection_criterion='bic',
    interactive_backend='widget'
)

# Get current settings
settings = epyr.baseline.get_configuration()
print(settings)
```

## Advanced Usage

### Model Comparison

```python
# Detailed model comparison
results = epyr.baseline.compare_models_detailed(
    x, y, params, 
    models=['polynomial', 'stretched_exponential', 'bi_exponential']
)

for model_name, result in results.items():
    metrics = result['metrics']
    print(f"{model_name}:")
    print(f"  AIC: {metrics['aic']:.2f}")
    print(f"  BIC: {metrics['bic']:.2f}")
    print(f"  R²:  {metrics['r2']:.4f}")
```

### Recommendations

```python
# Get model recommendations based on data type
recommendations = epyr.baseline.get_model_recommendations(
    data_type='pulsed', 
    experiment_type='t2'
)
print(recommendations)  # ['stretched_exponential', 'bi_exponential', 'polynomial']

# Use recommendations
corrected, baseline, info = epyr.baseline.auto_baseline_with_recommendations(
    x, y, params, 
    data_type='pulsed',
    experiment_type='t2'
)
```

## Error Handling

### Common Issues and Solutions

1. **Interactive selection not working:**
   ```python
   # Check backend compatibility
   if not epyr.baseline.is_interactive_available():
       epyr.baseline.setup_widget_backend()
   ```

2. **Selection window stuck:**
   ```python
   # Force close windows
   epyr.baseline.close_selector_window()
   ```

3. **Complex data issues:**
   ```python
   # Handle complex EPR data properly
   corrected, baseline = epyr.baseline.baseline_stretched_exponential_1d(
       x, y, params, 
       use_real_part=True  # or False for magnitude
   )
   ```

## Best Practices

### Workflow Recommendations

1. **Start with automatic selection:**
   ```python
   corrected, baseline, info = epyr.baseline.baseline_auto_1d(x, y, params)
   print(f"Suggested model: {info['best_model']}")
   ```

2. **Fine-tune if needed:**
   ```python
   if info['best_model'] == 'polynomial':
       # Try different orders
       for order in [2, 3, 4]:
           corrected, baseline = epyr.baseline.baseline_polynomial_1d(
               x, y, params, order=order
           )
   ```

3. **Always visualize results:**
   ```python
   import matplotlib.pyplot as plt
   
   plt.figure(figsize=(12, 4))
   plt.subplot(1, 2, 1)
   plt.plot(x, y, 'b-', alpha=0.7, label='Original')
   plt.plot(x, baseline, 'r--', label='Baseline')
   plt.legend()
   
   plt.subplot(1, 2, 2)
   plt.plot(x, corrected, 'g-', label='Corrected')
   plt.legend()
   plt.show()
   ```

### Data Type Guidelines

| **EPR Data Type** | **Recommended Method** | **Key Parameters** |
|-------------------|----------------------|-------------------|
| **CW EPR spectra** | `baseline_polynomial_1d()` | `order=2-3, exclude_center=True` |
| **T2 relaxation** | `baseline_stretched_exponential_1d()` | `beta_range=(0.5, 2.0), use_real_part=True` |
| **T1 relaxation** | `baseline_stretched_exponential_1d()` | `beta_range=(0.8, 1.2)` |
| **Complex decay** | `baseline_bi_exponential_1d()` | `tau_ratio_min=2.5` |
| **2D datasets** | `baseline_polynomial_2d()` | `order=(2, 2)` |
| **Unknown** | `baseline_auto_1d()` | `verbose=True` |

## Help and Support

### Getting Help

```python
# Package-level help
epyr.baseline.get_help()

# Function-specific help
help(epyr.baseline.baseline_auto_1d)

# Jupyter-specific help
epyr.baseline.jupyter_help()

# Interactive demo
epyr.baseline.demo()
```

### Version Information

```python
import epyr
print(f"EPyR Tools version: {epyr.__version__}")
print(f"Baseline package version: {epyr.baseline.__version__}")
```

---

## Migration from Old System

If migrating from the old `baseline_correction.py` system:

### Old vs New Function Names

| **Old Function** | **New Function** | **Notes** |
|------------------|------------------|-----------|
| `baseline_polynomial()` | `baseline_polynomial_1d()` | Same functionality |
| `baseline_stretched_exponential()` | `baseline_stretched_exponential_1d()` | Enhanced parameters |
| Manual fitting | `baseline_auto_1d()` | New automatic selection |

### Backward Compatibility

All old functions are still available through compatibility aliases:

```python
# These still work (deprecated but functional)
corrected, baseline = epyr.baseline_polynomial(x, y, params)
corrected, baseline = epyr.baseline_stretched_exponential(x, y, params)

# Recommended new usage
corrected, baseline = epyr.baseline.baseline_polynomial_1d(x, y, params)
corrected, baseline = epyr.baseline.baseline_stretched_exponential_1d(x, y, params)
```