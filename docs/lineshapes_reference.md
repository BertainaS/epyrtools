# EPyR Tools - Lineshape Functions Reference

EPyR Tools provides a comprehensive set of lineshape functions specifically designed for EPR spectroscopy analysis. This reference covers all available functions, their parameters, and practical applications.

## Overview

The `epyr.lineshapes` module includes:

- **Individual functions**: `gaussian`, `lorentzian`, `voigtian`, `pseudo_voigt`, `lshape`
- **Unified interface**: `Lineshape` class for consistent API
- **Convolution tools**: `convspec` for spectrum broadening
- **Factory functions**: Convenient constructors for common cases

## Quick Start

```python
import epyr
from epyr.lineshapes import Lineshape, gaussian, lorentzian

# Method 1: Direct function calls
x = np.linspace(-10, 10, 1000)
y_gauss = gaussian(x, center=0, width=4.0)
y_lorentz = lorentzian(x, center=0, width=4.0)

# Method 2: Unified Lineshape class
gauss_shape = Lineshape('gaussian', width=4.0)
y_gauss = gauss_shape(x, center=0)

# Method 3: Factory functions
from epyr.lineshapes import create_gaussian
gauss_func = create_gaussian(width=4.0)
y_gauss = gauss_func(x, center=0)
```

## Core Functions

### gaussian(x, center, width, derivative=0, phase=0.0, return_both=False)

Area-normalized Gaussian lineshape representing inhomogeneous broadening.

**Parameters:**
- `x` (array): Abscissa points (magnetic field)
- `center` (float): Peak center position
- `width` (float): Full width at half maximum (FWHM)
- `derivative` (int, optional): Derivative order (0, 1, 2, -1)
- `phase` (float, optional): Phase rotation in radians
- `return_both` (bool, optional): Return (absorption, dispersion) tuple

**Examples:**
```python
# Basic Gaussian
y = gaussian(x, 0, 4.0)

# First derivative
dy = gaussian(x, 0, 4.0, derivative=1)

# Dispersion mode
disp = gaussian(x, 0, 4.0, phase=np.pi/2)

# Both components
abs_part, disp_part = gaussian(x, 0, 4.0, return_both=True)
```

**Applications:**
- Modeling inhomogeneous broadening
- Powder pattern simulations
- Temperature-dependent narrowing

### lorentzian(x, center, width, derivative=0, phase=0.0, return_both=False)

Area-normalized Lorentzian lineshape representing homogeneous broadening.

**Parameters:** Same as `gaussian`

**Examples:**
```python
# Basic Lorentzian
y = lorentzian(x, 0, 4.0)

# Second derivative
d2y = lorentzian(x, 0, 4.0, derivative=2)

# Phase rotation
mixed = lorentzian(x, 0, 4.0, phase=np.pi/4)
```

**Applications:**
- Modeling homogeneous broadening
- Lifetime broadening effects
- Collision-induced broadening

### voigtian(x, center, widths, derivative=0, phase=0.0, return_both=False)

True Voigt profile - convolution of Gaussian and Lorentzian components.

**Parameters:**
- `widths` (tuple): (gaussian_fwhm, lorentzian_fwhm)
- Other parameters same as above

**Examples:**
```python
# Equal Gaussian and Lorentzian widths
y = voigtian(x, 0, (3.0, 3.0))

# Gaussian-dominated
y_g = voigtian(x, 0, (5.0, 1.0))

# Lorentzian-dominated
y_l = voigtian(x, 0, (1.0, 5.0))
```

**Applications:**
- Realistic EPR lineshapes
- Combined broadening mechanisms
- High-precision spectroscopy

### pseudo_voigt(x, center, width, eta=0.5)

Fast pseudo-Voigt approximation: η×Lorentzian + (1-η)×Gaussian

**Parameters:**
- `eta` (float): Mixing parameter (0=Gaussian, 1=Lorentzian)

**Examples:**
```python
# 50/50 mix
y = pseudo_voigt(x, 0, 4.0, eta=0.5)

# Mostly Gaussian
y_g = pseudo_voigt(x, 0, 4.0, eta=0.2)

# Mostly Lorentzian
y_l = pseudo_voigt(x, 0, 4.0, eta=0.8)
```

**Applications:**
- Fast approximation to Voigt profiles
- Fitting experimental data
- Parameter optimization

### lshape(x, center, width, derivative=0, alpha=1.0, phase=0.0)

General lineshape function with flexible mixing.

**Parameters:**
- `alpha` (float): Shape parameter (1=Gaussian, 0=Lorentzian)
- `width` (float or tuple): FWHM, can be different for each component

**Examples:**
```python
# Pure Gaussian (same as gaussian())
y = lshape(x, 0, 4.0, alpha=1.0)

# Pure Lorentzian (same as lorentzian())
y = lshape(x, 0, 4.0, alpha=0.0)

# Mixed with different widths
y = lshape(x, 0, (3.0, 5.0), alpha=0.3)
```

## Unified Lineshape Class

### Lineshape(shape_type, width, alpha=1.0, derivative=0, phase=0.0, normalize=True)

Unified interface for all lineshape types.

**Parameters:**
- `shape_type` (str): 'gaussian', 'lorentzian', 'voigt', 'pseudo_voigt'
- Other parameters as above

**Methods:**
- `__call__(x, center, **kwargs)`: Generate lineshape
- `absorption(x, center)`: Pure absorption mode
- `dispersion(x, center)`: Pure dispersion mode
- `derivative(x, center, order)`: Derivative lineshape
- `both_components(x, center)`: (absorption, dispersion) tuple
- `set_width(width)`: Create new instance with different width
- `set_alpha(alpha)`: Create new instance with different mixing
- `info()`: Get lineshape information

**Examples:**
```python
# Create lineshape objects
gauss = Lineshape('gaussian', width=4.0)
lorentz = Lineshape('lorentzian', width=4.0, derivative=1)
voigt = Lineshape('voigt', width=(3.0, 2.0))

# Generate data
y1 = gauss(x, center=0)
y2 = lorentz(x, center=2)
y3 = voigt.absorption(x, center=-2)

# Modify parameters
wide_gauss = gauss.set_width(8.0)
gauss_deriv = gauss.set_derivative(1)

# Get information
info = gauss.info()
print(f"Shape: {info['shape_type']}, Width: {info['width']}")
```

## Convolution Function

### convspec(spectrum, step_size, width, derivative=0, alpha=1.0, phase=0.0)

Convolve spectrum with lineshape functions for broadening simulation.

**Parameters:**
- `spectrum` (array): Input spectrum to convolve
- `step_size` (float): Abscissa step size
- `width` (float): FWHM for broadening
- Other parameters as above

**Examples:**
```python
# Create stick spectrum
stick = np.zeros(1000)
stick[300] = 1.0  # Peak at position 300
stick[700] = 0.5  # Peak at position 700

# Apply Gaussian broadening
dx = 0.1  # mT per point
broadened = convspec(stick, dx, width=2.0, alpha=1.0)

# Apply Lorentzian broadening
lorentz_broad = convspec(stick, dx, width=2.0, alpha=0.0)
```

**Applications:**
- Converting calculated stick spectra to realistic lineshapes
- Modeling instrumental broadening
- Simulating temperature effects

## Factory Functions

Convenient constructors for common lineshape types:

```python
from epyr.lineshapes import (
    create_gaussian, create_lorentzian, 
    create_voigt, create_pseudo_voigt
)

# Create pre-configured lineshape functions
gauss_func = create_gaussian(width=4.0)
lorentz_func = create_lorentzian(width=4.0, derivative=1)
voigt_func = create_voigt(gaussian_width=3.0, lorentzian_width=2.0)
pv_func = create_pseudo_voigt(width=4.0, alpha=0.3)

# Use them
y = gauss_func(x, center=0)
```

## Advanced Features

### Derivative Spectroscopy

All lineshapes support analytical derivatives:

```python
# Function and its derivatives
y0 = gaussian(x, 0, 4, derivative=0)  # Function
y1 = gaussian(x, 0, 4, derivative=1)  # 1st derivative
y2 = gaussian(x, 0, 4, derivative=2)  # 2nd derivative

# Applications
# - Enhanced resolution
# - Zero-crossing detection
# - Overlapping signal separation
```

### Phase Rotation

Control absorption/dispersion character:

```python
# Pure absorption
abs_mode = gaussian(x, 0, 4, phase=0)

# Pure dispersion  
disp_mode = gaussian(x, 0, 4, phase=np.pi/2)

# Mixed phase (45°)
mixed = gaussian(x, 0, 4, phase=np.pi/4)

# Applications
# - Phase correction
# - Complex impedance analysis
# - Signal optimization
```

### Multi-component Analysis

Combine multiple lineshapes for complex spectra:

```python
# Three-component system
component1 = 1.0 * gaussian(x, -2, 1.5)
component2 = 0.8 * lorentzian(x, 1, 2.0)  
component3 = 1.2 * pseudo_voigt(x, 4, 1.8, eta=0.3)

total_spectrum = component1 + component2 + component3
```

## Performance Considerations

### Optimization Tips

1. **Use appropriate lineshape type:**
   - Gaussian: Fast, good for inhomogeneous broadening
   - Lorentzian: Fast, good for homogeneous broadening
   - Pseudo-Voigt: Fast approximation to Voigt
   - True Voigt: Most accurate but slower

2. **Vectorized operations:**
   ```python
   # Good: Vectorized
   centers = [-2, 0, 2]
   total = sum(gaussian(x, c, 2.0) for c in centers)
   
   # Better: Use broadcasting when possible
   ```

3. **Caching for repeated calculations:**
   ```python
   # Create once, use many times
   shape = Lineshape('gaussian', width=4.0)
   results = [shape(x, center=c) for c in center_list]
   ```

### Memory Usage

- All functions work with NumPy arrays
- Memory usage scales linearly with array size
- Large arrays (>10⁶ points) may need chunked processing

## Common Applications in EPR

### g-Factor Analysis

```python
# Different g-factors shift resonance position
g_values = [1.99, 2.00, 2.01]
microwave_freq = 9.5e9  # Hz
mu_B = 9.274e-24  # J/T
h = 6.626e-34  # J·s

for g in g_values:
    B_res = h * microwave_freq / (mu_B * g) * 1000  # mT
    spectrum = gaussian(field_range, B_res, linewidth)
```

### Hyperfine Coupling

```python
# Triplet from nitrogen coupling (I=1)
A_nitrogen = 1.5  # mT
intensities = [1, 1, 1]  # Equal intensities
positions = [-A_nitrogen, 0, A_nitrogen]

total = sum(I * gaussian(x, pos, width) 
           for I, pos in zip(intensities, positions))
```

### Temperature Dependence

```python
temperatures = [77, 200, 300, 400]  # K
base_width = 1.0  # mT

for T in temperatures:
    # Width increases with temperature
    width_T = base_width * (T / 300) ** 0.5
    spectrum_T = lorentzian(x, 0, width_T)
```

## Error Handling

The lineshape functions include comprehensive input validation:

```python
# These will raise ValueError:
gaussian(x, "invalid", 4)     # Non-numeric center
gaussian(x, 0, -1)            # Negative width
gaussian(x, 0, 4, derivative=-2)  # Invalid derivative
Lineshape('invalid_type')     # Unknown shape type
```

## Integration with EPyR Tools

The lineshape functions integrate seamlessly with other EPyR Tools modules:

```python
import epyr

# Load EPR data
x, y, params, path = epyr.eprload()

# Apply baseline correction
y_corrected = epyr.baseline.baseline_polynomial(y, x, order=1)

# Fit with lineshape
fitted_shape = epyr.gaussian(x, center=params.get('center', 0), 
                           width=params.get('width', 2.0))

# Convert to FAIR format with fit results
epyr.fair.convert_with_analysis(path, fitted_params={'shape': 'gaussian'})
```

## Best Practices

1. **Choose the right lineshape:**
   - Gaussian for inhomogeneous broadening
   - Lorentzian for homogeneous broadening  
   - Voigt for realistic cases

2. **Use derivatives for enhanced resolution:**
   - 1st derivative: Zero crossings mark centers
   - 2nd derivative: Better peak separation

3. **Validate with standards:**
   - Test with known samples (DPPH, TEMPO)
   - Compare with literature values

4. **Consider instrumental effects:**
   - Use convolution for broadening simulation
   - Account for phase errors

5. **Optimize parameters systematically:**
   - Start with simple models
   - Add complexity as needed
   - Use statistical validation

## See Also

- [User Guide](user_guide.md) - General EPyR Tools usage
- [API Reference](api_reference.md) - Complete function documentation  
- [Examples](../examples/) - Practical application examples
- [Tutorials](../examples/notebooks/) - Interactive learning materials