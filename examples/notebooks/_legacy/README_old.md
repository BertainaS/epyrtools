# EPyR Tools - Jupyter Notebooks

This directory contains interactive Jupyter notebooks for learning and using EPyR Tools.

## Available Notebooks

### Getting_Started.ipynb

**Perfect for beginners!** This notebook demonstrates:

- Loading 1D and 2D EPR data from real measurements
- Handling both BES3T (.dsc/.dta) and ESP (.par/.spc) formats
- Automatic detection of complex vs real data
- Professional EPR spectrum visualization
- Basic data analysis and statistics
- Exporting data to common formats (CSV, TXT, NPZ)

**Data Used:**
- `130406SB_CaWO4_Er_CW_5K_20.DSC` - 1D CW-EPR spectrum
- `Rabi2D_GdCaWO4_13dB_3057G.DSC` - 2D pulsed EPR (complex data)
- All data files are now in a single `../data/` folder

## How to Use

### 1. Start Jupyter
```bash
# From the project root directory
cd examples/notebooks
jupyter notebook
```

### 2. Open the notebook
Click on `Getting_Started.ipynb` in the Jupyter interface

### 3. Run the cells
- Press `Shift+Enter` to run each cell
- Or use "Run All" from the Cell menu

## What You'll See

The notebook will automatically:
- Find your EPR data files
- Load and display 1D spectra with proper field ranges
- Handle 2D complex data (shows magnitude for visualization)
- Create professional plots with proper axes and labels
- Export your data for use in other software

## Data Requirements

The notebook works with the sample data in `../data/`:
- All EPR files (.DSC/.DTA and .PAR/.SPC) in a single folder
- Supports both BES3T and ESP formats

## Troubleshooting

**"No files found"**: Make sure your EPR data files are in the correct directory:
- All EPR files: `examples/data/` (single folder)

**Import errors**: Make sure EPyR Tools is installed:
```bash
pip install -e .
```

**Complex plotting issues**: The notebook automatically handles complex data by showing the magnitude.

## Next Steps

After completing the Getting Started notebook:

1. **Learn baseline correction** - Remove drift and artifacts
2. **Advanced analysis** - Peak fitting, g-factor calculations
3. **Batch processing** - Handle multiple files efficiently
4. **Custom analysis** - Adapt the notebook for your specific needs

## Contributing

To add more notebooks:
1. Follow the same structure as Getting_Started.ipynb
2. Include clear explanations and real data examples
3. Test with the actual measurement files
4. Update this README

---

**EPyR Tools** - Making EPR data analysis accessible and reproducible!
