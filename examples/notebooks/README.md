# EPyR Tools Tutorial Notebooks

A hands-on tour of the EPyR Tools public API. Each notebook is standalone, runs
top to bottom, and uses the real experimental data in `../data/`.

Start with `00_Tutorial_Series_Index.ipynb`, then work through the series in
order.

## The series

| # | Notebook | Covers |
|---|----------|--------|
| 00 | `00_Tutorial_Series_Index.ipynb` | Overview, environment check, dataset list |
| 01 | `01_Loading_and_Plotting.ipynb` | `eprload`, parameter inspection, 1D/2D plotting |
| 02 | `02_Baseline_Correction.ipynb` | Polynomial, automatic, and exponential baselines |
| 03 | `03_Lineshape_Analysis_and_Fitting.ipynb` | Gaussian/Lorentzian/Voigt, derivatives, fitting |
| 04 | `04_Relaxation_Fitting.ipynb` | T1/T2 decay/recovery models (new in v0.4.0) |
| 05 | `05_Signal_Processing_and_FFT.ipynb` | Frequency analysis of Rabi data, apodization |
| 06 | `06_FAIR_Conversion_and_Export.ipynb` | CSV/JSON/HDF5 export and validation |
| 07 | `07_Physics_Units_and_Constants.ipynb` | CODATA constants, field/frequency conversions |

## Running the notebooks

From the repository root:

```bash
pip install -e .          # make `epyr` importable
jupyter lab examples/notebooks
```

Notebooks are committed without cell outputs to keep the repository light. Run a
notebook to generate its figures and printed results.

## Notes

- Paths are written relative to this directory (`../data/...`), so launch Jupyter
  from here or from the repository root.
- The interactive features (`plot_2d_slicer`, the `epyr.isotopes()` GUI) need a
  desktop or `%matplotlib widget` backend; the notebooks describe them rather
  than executing them so they run unattended.
- Older, superseded notebooks are kept under `_legacy/` for reference.
