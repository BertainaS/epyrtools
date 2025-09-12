# EPyR Tools Examples and Tutorials

This directory contains comprehensive examples, tutorials, and sample data for EPyR Tools - a Python package for EPR (Electron Paramagnetic Resonance) data analysis.

## Directory Structure

```
examples/
├── README.md                   # This file
├── data/                      # Sample EPR data files
│   ├── BES3T/                # Bruker BES3T format (.dsc/.dta pairs)
│   ├── ESP/                  # Bruker ESP format (.par/.spc pairs)
│   └── processed/            # Converted and analyzed data
├── notebooks/                # Interactive Jupyter tutorials
│   ├── 01_Getting_Started.ipynb
│   ├── 02_Baseline_Correction.ipynb
│   └── 03_Advanced_Analysis.ipynb
└── scripts/                  # Python example scripts
    ├── 01_basic_loading.py
    └── 02_baseline_correction.py
```

## Getting Started

### Prerequisites

Ensure you have EPyR Tools installed:
```bash
pip install -r requirements.txt  # From project root
```

For Jupyter notebooks, also install:
```bash
pip install jupyter matplotlib pandas
```

### Quick Start Options

1. **Interactive Notebooks** (Recommended)
   ```bash
   cd examples/notebooks
   jupyter notebook
   ```

2. **Python Scripts**
   ```bash
   cd examples/scripts
   python 01_basic_loading.py
   ```

## Tutorial Overview

### 01_Getting_Started
**Perfect for beginners**
- Loading EPR data from Bruker files
- Understanding EPR measurement parameters
- Basic data visualization and plotting
- Converting to open data formats (CSV, JSON, HDF5)
- FAIR data principles for EPR

**Time:** 15-20 minutes
**Level:** Beginner
**Files:** `notebooks/01_Getting_Started.ipynb`, `scripts/01_basic_loading.py`

### 02_Baseline_Correction
**Essential EPR data processing**
- Understanding baseline issues in EPR
- Polynomial correction methods (constant, linear, quadratic)
- Advanced correction with signal exclusion regions
- Method comparison and quality assessment
- Best practices and guidelines

**Time:** 30-40 minutes
**Level:** Intermediate
**Files:** `notebooks/02_Baseline_Correction.ipynb`, `scripts/02_baseline_correction.py`

### 03_Advanced_Analysis
**Comprehensive EPR analysis**
- Peak detection and hyperfine structure analysis
- g-factor calculations and field calibration
- Spectral integration for quantitative analysis
- Hyperfine coupling constant determination
- Pattern recognition and interpretation
- Data export for publication

**Time:** 45-60 minutes
**Level:** Advanced
**Files:** `notebooks/03_Advanced_Analysis.ipynb`

## Sample Data

### Adding Your EPR Data

1. **BES3T Format** (Modern Bruker)
   - Place `.dsc/.dta` file pairs in `data/BES3T/`
   - Used by EMX, Elexsys series spectrometers

2. **ESP Format** (Legacy Bruker)
   - Place `.par/.spc` file pairs in `data/ESP/`
   - Used by ESP-300 series and WinEPR

### Data Organization

- Use descriptive filenames (e.g., `TEMPO_solution_RT.dsc`)
- Keep file pairs together with identical base names
- Organize by sample type or measurement conditions
- Document measurement parameters in file names or metadata

## Tutorial Features

### Interactive Elements
- Real-time data visualization
- Parameter adjustment examples
- Quality assessment metrics
- Method comparison tools

### Educational Content
- EPR physics background
- Method theory and limitations
- Best practices and guidelines
- Common pitfalls and solutions

### Practical Skills
- File format handling
- Data processing workflows
- Quality control procedures
- Results interpretation

## Learning Path

**Recommended progression:**

1. **Start Here**: `01_Getting_Started.ipynb`
   - Learn basic EPyR Tools usage
   - Understand file formats and data loading
   - Practice basic visualization

2. **Essential Processing**: `02_Baseline_Correction.ipynb`
   - Master baseline correction techniques
   - Learn quality assessment
   - Understand method selection

3. **Advanced Techniques**: `03_Advanced_Analysis.ipynb`
   - Perform quantitative EPR analysis
   - Calculate physical parameters
   - Export publication-ready results

4. **Your Research**: Apply learned techniques to your EPR data!

## EPR Concepts Covered

### Experimental Techniques
- Continuous wave (CW) EPR
- Field-swept measurements
- First derivative detection
- Power and modulation effects

### Data Analysis Methods
- Baseline correction algorithms
- Peak detection and fitting
- Spectral integration
- g-factor calculations
- Hyperfine analysis

### Physical Parameters
- Electronic g-factors
- Hyperfine coupling constants
- Linewidths and relaxation
- Spin concentrations
- Anisotropy parameters

## Tips for Success

### Before Starting
- Read the EPR theory sections in each notebook
- Ensure your data files are properly formatted
- Have sample EPR data ready for practice
- Install all required dependencies

### During Tutorials
- Execute code cells in order
- Experiment with parameters
- Try your own EPR data
- Read all explanatory text

### After Completion
- Save analysis results and plots
- Document your workflows
- Share knowledge with colleagues
- Contribute improvements to EPyR Tools

## Getting Help

### Documentation
- Check the API documentation for function details
- Read docstrings and code comments
- Review example outputs and plots

### Community
- Report issues: [GitHub Issues](https://github.com/user/epyr-tools/issues)
- Ask questions: Use GitHub Discussions
- Contribute: Submit pull requests for improvements

### Common Issues
- **File not found**: Check data directory paths
- **Import errors**: Ensure EPyR Tools is installed
- **Plotting issues**: Install matplotlib and check backend
- **Memory issues**: Use smaller data files for testing

## Contributing

We welcome contributions to improve these tutorials:

- Add new example notebooks
- Provide sample EPR data
- Improve documentation
- Fix bugs and typos
- Suggest new analysis methods

See the project's contribution guidelines for details.

## References

### EPR Theory
- Weil, J.A. & Bolton, J.R. "Electron Paramagnetic Resonance" (2007)
- Schweiger, A. & Jeschke, G. "Principles of Pulse Electron Paramagnetic Resonance" (2001)

### Data Analysis
- Stoll, S. & Schweiger, A. "EasySpin, a comprehensive software package for spectral simulation and analysis in EPR" J. Magn. Reson. (2006)
- EPR data processing best practices literature

### FAIR Data
- Wilkinson, M.D. et al. "The FAIR Guiding Principles for scientific data management" Sci Data (2016)

---

**Happy EPR analyzing with EPyR Tools!**

*For questions or suggestions, please open an issue on the project repository.*
