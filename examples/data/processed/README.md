# Processed Data Examples

This directory contains examples of processed EPR data in various open formats.

## File Types

### CSV Format
Simple comma-separated values for spectral data:
```csv
# EPR Spectrum Data
# Original file: sample.dsc
field_gauss,intensity_au
3200.0,0.123
3201.0,0.125
...
```

### JSON Format
Structured metadata with human-readable parameter names:
```json
{
  "original_file": "sample.dsc",
  "measurement_parameters": {
    "microwave_frequency": {"value": 9.4e9, "unit": "Hz"},
    "magnetic_field_center": {"value": 3350.0, "unit": "G"}
  },
  "data": {
    "field_axis": [3200.0, 3201.0, ...],
    "intensity": [0.123, 0.125, ...]
  }
}
```

### HDF5 Format
Self-contained hierarchical format:
```
sample.h5
├── data/
│   ├── intensity          # EPR signal
│   └── field_axis        # Magnetic field
├── metadata/
│   ├── parameters_fair/   # FAIR-mapped parameters
│   └── parameters_original/ # Original Bruker parameters
└── attributes             # Global metadata
```

## Baseline Corrected Data

Examples of EPR spectra before and after baseline correction:
- Original spectra with baseline drift
- Polynomial-corrected spectra
- Exponential-corrected spectra

## Usage Examples

### Load Processed Data
```python
import pandas as pd
import json
import h5py

# Load CSV
df = pd.read_csv('sample_corrected.csv')

# Load JSON
with open('sample_metadata.json', 'r') as f:
    data = json.load(f)

# Load HDF5
with h5py.File('sample.h5', 'r') as f:
    intensity = f['data/intensity'][:]
    field = f['data/field_axis'][:]
```

## File Organization

- `*_original.*` - Raw converted data from Bruker files
- `*_corrected.*` - Baseline-corrected spectra
- `*_metadata.*` - Parameter information in FAIR format
- `*_analysis.*` - Results from EPR analysis workflows
