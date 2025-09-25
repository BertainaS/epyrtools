# Interactive 2D EPR Viewer Guide - plot_2d_slicer

## Overview

The `plot_2d_slicer` function in EPyR Tools allows interactive visualization of 2D EPR data slice by slice with a navigation slider. This function is particularly useful for exploring complex 2D datasets like Rabi measurements, T2 experiments, or HYSCORE experiments.

## Key Features

- 🎛️ **Interactive navigation** with matplotlib slider
- 🔄 **Two directions**: horizontal or vertical slices  
- 👁️ **Overview display** with real-time position indicator
- 📏 **Automatic Y-scale adjustment** for each slice
- 🔗 **Full integration** with `eprload()` and Bruker parameters
- 📊 **Complex data support** (uses real part)

## Installation and Prerequisites

```python
# Function is included in EPyR Tools v0.1.8+
import epyr

# Prerequisites for interactivity in Jupyter
%matplotlib widget  # or %matplotlib notebook
```

## Basic Usage

```python
# Load 2D EPR data
x, y, params, _ = epyr.eprload("your_2d_file.DTA")

# Launch interactive viewer
epyr.plot_2d_slicer(x, y, params)
```

## Function Parameters

```python
plot_2d_slicer(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    slice_direction: str = 'horizontal'
) -> Dict[str, Any]
```

### Arguments

- **`x`**: Axis data from `eprload()` (can be None, array, or list)
- **`y`**: 2D spectral data array (shape `(ny, nx)`)
- **`params`**: Parameter dictionary from `eprload()`
- **`title`**: Custom title for the plot (optional)
- **`slice_direction`**: Slice direction (`'horizontal'` or `'vertical'`)

### Return Value

Dictionary containing matplotlib objects for advanced manipulation:
- `'figure'`: Matplotlib figure
- `'ax_main'`: Main plot axes
- `'ax_overview'`: Overview axes
- `'slider'`: Slider widget
- `'line'`: Main plot line
- `'slice_line'`: Position indicator line

## Visualization Modes

### 1. Horizontal Slices (default)

```python
# Navigate in Y-axis - display horizontal slices
epyr.plot_2d_slicer(x, y, params, slice_direction='horizontal')
```

- **Navigation**: Through Y direction (parameter/time/index)
- **Display**: Spectra as function of X (field/frequency)
- **Usage**: Ideal for viewing temporal or parametric evolution

### 2. Vertical Slices

```python
# Navigate in X-axis - display vertical slices  
epyr.plot_2d_slicer(x, y, params, slice_direction='vertical')
```

- **Navigation**: Through X direction (field/frequency)
- **Display**: Evolution as function of Y (parameter/time)
- **Usage**: Useful for analyzing specific resonances

## Usage Examples

### Example 1: 2D Rabi Data

```python
import epyr
import matplotlib.pyplot as plt

# Enable interactive widgets
%matplotlib widget

# Load 2D Rabi data
x, y, params, _ = epyr.eprload("Rabi2D_experiment.DTA")

# Visualize temporal slices (horizontal)
slicer = epyr.plot_2d_slicer(
    x, y, params,
    title="2D Rabi Experiment - Time Evolution",
    slice_direction='horizontal'
)
```

### Example 2: Comparative Analysis

```python
# Compare both directions
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Horizontal direction
slicer_h = epyr.plot_2d_slicer(x, y, params, slice_direction='horizontal')

# Vertical direction  
slicer_v = epyr.plot_2d_slicer(x, y, params, slice_direction='vertical')
```

### Example 3: Advanced Customization

```python
# Launch slicer and get objects
slicer = epyr.plot_2d_slicer(x, y, params)

# Customize appearance
slicer['line'].set_color('red')
slicer['line'].set_linewidth(2)
slicer['slice_line'].set_color('yellow')

# Modify title dynamically
slicer['ax_main'].set_title("Custom Title")

# Apply changes
slicer['figure'].canvas.draw_idle()
```

## Configuration for Different Environments

### Jupyter Lab / Jupyter Notebook

```python
# Option 1: Interactive widgets (recommended)
%matplotlib widget

# Option 2: Notebook mode (alternative)
%matplotlib notebook

# Then use normally
epyr.plot_2d_slicer(x, y, params)
```

### Standalone Python Scripts

```python
import matplotlib.pyplot as plt

# Enable interactive mode
plt.ion()

# Launch slicer
epyr.plot_2d_slicer(x, y, params)

# Keep window open
plt.show(block=True)
```

## Supported Data Types

### Accepted X Data Formats

```python
# Format 1: List with X and Y axes
x = [x_axis_array, y_axis_array]  # Recommended for 2D

# Format 2: Single array (X-axis only)
x = x_axis_array

# Format 3: None (uses indices)
x = None
```

### Complex Data Handling

```python
# Complex data is automatically converted
if np.iscomplexobj(y):
    plot_data = np.real(y)  # Uses real part
```

## Common Error Handling

### Error: "matplotlib.widgets required"

```python
# Solution: Enable interactive backend
%matplotlib widget  # In Jupyter

# Or install matplotlib with widgets
pip install matplotlib[widgets]
```

### Error: "Expected 2D data"

```python
# Check dimensions
print(f"Data shape: {y.shape}")  # Must be (ny, nx)

# For 1D data, use plot_1d instead
if y.ndim == 1:
    epyr.plot_1d(x, y, params)
```

### Problem: Unresponsive Slider

```python
# Ensure correct interactive mode
%matplotlib widget  # For Jupyter Lab
%matplotlib notebook  # For classic Jupyter

# Restart kernel if necessary
```

## Usage Tips

### Performance with Large Datasets

```python
# For very large datasets, consider subsampling
if y.shape[0] > 1000:  # Too many slices
    step = y.shape[0] // 500  # Keep ~500 slices
    y_sub = y[::step, :]
    x_sub = [x[0], x[1][::step]] if isinstance(x, list) else x
    
    epyr.plot_2d_slicer(x_sub, y_sub, params)
```

### Image Saving

```python
# Get figure for saving
slicer = epyr.plot_2d_slicer(x, y, params)
slicer['figure'].savefig('my_2d_visualization.png', dpi=300)
```

### Keyboard Navigation (future)

The function could be extended to support keyboard navigation:

```python
# Future functionality
# Left/right arrows for navigation
# Number keys to jump to specific positions
```

## Integration with EPyR Tools Ecosystem

### With Baseline Correction

```python
# Correct baseline before visualization
x, y, params, _ = epyr.eprload("data.DTA")

# Apply baseline correction to each slice if needed
for i in range(y.shape[0]):
    y_corrected, _ = epyr.baseline.baseline_polynomial_1d(
        x[0], y[i, :], params, order=1
    )
    y[i, :] = y_corrected

# Visualize corrected data
epyr.plot_2d_slicer(x, y, params)
```

### With Lineshape Analysis

```python
# Analyze specific lines in the slicer
slicer = epyr.plot_2d_slicer(x, y, params)

# Use current position for analysis
current_slice_idx = int(slicer['slider'].val)
current_spectrum = np.real(y[current_slice_idx, :])

# Analyze with lineshapes module
peaks = epyr.lineshapes.find_peaks(x[0], current_spectrum)
```

## Current Limitations

- **Widgets only**: Requires interactive matplotlib backend
- **2D data only**: Only works with 2D arrays
- **Sequential navigation**: No direct navigation by parameter value

## Future Developments

- Support for keyboard navigation
- Zoom and pan in individual slices
- Automatic export of selected slices
- Support for 3D datasets (plane selection)
- Integration with quantitative analysis tools

## Support and Troubleshooting

For help:

```python
# Integrated documentation
help(epyr.plot_2d_slicer)

# EPyR Tools support
epyr.baseline.get_help()

# Check installation
print(epyr.__version__)  # Must be >= 0.1.8
```

---

This function significantly enhances EPyR Tools' 2D visualization capabilities and facilitates interactive exploration of complex EPR datasets.