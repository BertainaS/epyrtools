# BES3T Format Sample Data

This directory contains Bruker BES3T format EPR data files (.dta/.dsc pairs).

## File Format Information

**BES3T** is the modern Bruker format used by:
- EMX spectrometers
- Elexsys series (E500, E580, etc.)
- Newer Bruker EPR systems

### File Pairs
Each EPR measurement consists of two files:
- **`.dsc`** - Descriptor file (text format) containing measurement parameters
- **`.dta`** - Data file (binary format) containing the actual spectral data

Both files must be present and have the same base filename.

## Usage

```python
import epyr.eprload as eprload

# Load BES3T data (specify either .dsc or .dta)
x, y, params, filepath = eprload.eprload('filename.dsc')

# Or let the file dialog choose
x, y, params, filepath = eprload.eprload()
```

## Expected File Naming

Place your BES3T files here with descriptive names:
- `sample_name.dsc` + `sample_name.dta`
- Examples: `DPPH_powder.dsc/.dta`, `nitroxide_solution.dsc/.dta`

## Sample Data Guidelines

When adding sample data:
- Use descriptive filenames
- Keep file sizes reasonable (< 5 MB per file pair)
- Include diverse EPR samples (radicals, transition metals, etc.)
- Document measurement conditions in filename or separate notes
