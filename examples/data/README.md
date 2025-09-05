# ESP/WinEPR Format Sample Data

This directory contains Bruker ESP/WinEPR format EPR data files (.spc/.par pairs).

## File Format Information

**ESP/WinEPR** is the legacy Bruker format used by:
- ESP-300 series spectrometers
- Older Bruker EPR systems
- WinEPR software package

### File Pairs
Each EPR measurement consists of two files:
- **`.par`** - Parameter file (text format) containing measurement settings
- **`.spc`** - Spectrum file (binary format) containing the spectral data

Both files must be present and have the same base filename.

## Usage

```python
import epyr.eprload as eprload

# Load ESP data (specify either .par or .spc)
x, y, params, filepath = eprload.eprload('filename.par')

# Or let the file dialog choose
x, y, params, filepath = eprload.eprload()
```

## Parameter File Format

ESP `.par` files contain space-separated parameter values:
```
DOS  FormatANZ 75776MIN -20505.589844MAX 24769.074219JSS 4096SSX 2048...
HCF 3350.000000HSW 700.000000RES 2048REY 37MF  9.399263MP  1.001e-001...
```

Common parameters:
- `HCF`: Center field (Gauss)
- `HSW`: Sweep width (Gauss)
- `RES`: Resolution (number of points)
- `MF`: Microwave frequency (GHz)
- `MP`: Microwave power

## Expected File Naming

Place your ESP files here with descriptive names:
- `sample_name.par` + `sample_name.spc`
- Examples: `DPPH_standard.par/.spc`, `Cu_complex.par/.spc`

## Sample Data Guidelines

When adding sample data:
- Use descriptive filenames indicating sample type
- Keep file sizes reasonable (< 5 MB per file pair)
- Include variety of EPR samples and measurement conditions
- Legacy ESP data is valuable for testing compatibility
