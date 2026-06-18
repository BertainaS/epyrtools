# Relaxation Fitting Module — Design Spec

Date: 2026-06-18
Status: Approved

## Purpose

EPyR Tools currently provides spectral lineshape fitting (`epyr.lineshapes.fitting`)
for field-domain CW EPR spectra (Gaussian, Lorentzian, Voigt, pseudo-Voigt, with
derivatives and phase). There is no equivalent for time-domain relaxation data:
T1 (saturation/inversion recovery) and T2 (echo decay) measurements.

This module adds decay/recovery model functions and a fitting layer for
extracting relaxation time constants from time-domain EPR signals, following
the architecture and conventions already established by `lineshapes/fitting.py`.

## Scope

- 1D curve fitting only. A single decay/recovery curve in, one `RelaxationFitResult`
  out. Batch fitting of 2D datasets (e.g. T2(B) maps from a 2D echo-decay dataset)
  is out of scope for this iteration; the caller can loop over rows themselves
  using the 1D API.
- Six decay/recovery models (see below). No ESEEM/oscillatory modulation models.

## Package layout

New package `epyr/relaxation/`, mirroring `epyr/lineshapes/`:

```
epyr/relaxation/
    __init__.py     # public exports
    models.py       # standalone decay/recovery functions
    fitting.py       # RelaxationFitResult, fit_relaxation, fit_multiple_decays
```

Top-level `epyr/__init__.py` gains `from . import relaxation` alongside the
existing `from . import baseline`, `from . import lineshapes, signalprocessing`.
Fitting entry points are used via `from epyr.relaxation import fit_relaxation`,
matching the existing `from epyr.lineshapes.fitting import fit_epr_signal` pattern.

## Models (`models.py`)

All functions are vectorized NumPy, float64, signature `function(t, ..., offset=0.0)`.
The time axis `t` is in whatever consistent unit the caller provides (e.g. ns or
us from `eprload`); fitted time constants are returned in that same unit. No
automatic unit detection, consistent with `fit_epr_signal` not detecting field units.

| Function | Formula | Typical use |
|---|---|---|
| `mono_exponential(t, amplitude, T, offset=0.0)` | `A * exp(-t/T) + offset` | simple T2 echo decay, or any single-rate decay |
| `stretched_exponential(t, amplitude, T, beta, offset=0.0)` | `A * exp(-(t/T)**beta) + offset` | distributed relaxation rates |
| `biexponential(t, amplitude1, tau1, amplitude2, tau2, offset=0.0)` | `A1*exp(-t/tau1) + A2*exp(-t/tau2) + offset` | two relaxation pathways/sites |
| `inversion_recovery(t, amplitude, T1, offset=0.0)` | `A * (1 - 2*exp(-t/T1)) + offset` | T1 via inversion-recovery pulse sequence |
| `saturation_recovery(t, amplitude, T1, offset=0.0)` | `A * (1 - exp(-t/T1)) + offset` | T1 via saturation-recovery pulse sequence |
| `gamma_gaussian_decay(t, amplitude, Gamma0, GammaG, offset=0.0)` | `A * exp(-Gamma0*t) * exp(-(GammaG*t)**2) + offset` | echo decay with a spectral-diffusion (Gaussian) contribution on top of the homogeneous (exponential) one |

Notes on naming:
- `biexponential` uses `tau1`/`tau2` (not `T1`/`T2`) to avoid confusion with the
  physical T1/T2 relaxation times; either component could represent a T1 or a
  T2 process depending on the experiment.
- `inversion_recovery`/`saturation_recovery` use `T1` because they are specific
  to that measurement.
- The offset parameter is always present (no toggle to disable it); bounds keep
  it free to converge near zero when the decay reaches baseline.

## Fitting layer (`fitting.py`)

### `RelaxationFitResult` dataclass

Mirrors `lineshapes.fitting.FitResult` field-for-field: `model`, `parameters`,
`parameter_errors`, `fitted_curve`, `residuals`, `r_squared`, `chi_squared`,
`success`, `message`, `covariance_matrix`, `x_fit`, plus a `.summary()` method
producing the same formatted-string style.

### `fit_relaxation`

```python
def fit_relaxation(
    t_data: np.ndarray,
    y_data: np.ndarray,
    model: str = "mono_exponential",
    initial_params: Optional[Dict[str, float]] = None,
    bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
    time_unit: str = "",
    **fit_kwargs,
) -> RelaxationFitResult:
```

Behavior parallels `fit_epr_signal`:
- Input validation (length match, >= 4 points), NaN removal, optional boolean
  mask to exclude points.
- Model name dispatches to its fit function, parameter name list, and default
  bounds (`_get_fit_function`, mirroring the lineshapes implementation).
- Initial parameter estimation per model (`_estimate_initial_params`):
  - decay-rate models (`mono_exponential`, `stretched_exponential`,
    `gamma_gaussian_decay`): offset from the tail of the data, amplitude from
    `y[0] - offset`, rate from a log-linear regression of `log(|y - offset|)`
    against `t`, restricted to points where `|y - offset|` exceeds a small
    fraction of the amplitude estimate (avoids the regression being skewed by
    near-zero values close to the decay's asymptote).
  - `biexponential`: offset from the tail, amplitude split from the total
    span, `tau1`/`tau2` seeded at roughly 10% and 100% of the time range.
  - `inversion_recovery`/`saturation_recovery`: amplitude from the data span,
    `T1` seeded from the time at which the signal crosses its midpoint.
- Bounds (`_setup_bounds`): all time constants and rates constrained positive;
  `beta` defaults to the bound `(0.05, 5.0)` (covers sub-exponential through
  compressed-exponential behavior), user-overridable via `bounds=`;
  `amplitude`/`offset` bounded by data range, matching the pattern already
  used for `baseline_slope`/`baseline_offset` in `lineshapes/fitting.py`.
- Fit via `scipy.optimize.curve_fit` with the same defaults
  (`maxfev=10000, method="trf"`), overridable through `**fit_kwargs`.
- Failures are caught and returned as `RelaxationFitResult(success=False, ...)`
  rather than raised, matching `fit_epr_signal`.
- `time_unit` is cosmetic only: if non-empty, it annotates the plot's x-axis
  label ("Time (us)") and the printed summary; it does not affect fitting.
- Plot (`plot=True`): data + fitted curve, parameter textbox (value plus
  error, matching the existing 4/2-significant-figure display convention),
  residuals panel below — same two-panel layout as `_plot_fit_results`.

### `fit_multiple_decays`

```python
def fit_multiple_decays(
    t_data: np.ndarray,
    y_data: np.ndarray,
    models: Optional[List[str]] = None,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
) -> Dict[str, RelaxationFitResult]:
```

Default `models=["mono_exponential", "stretched_exponential", "biexponential"]`.
The recovery models and `gamma_gaussian_decay` are excluded from the default
comparison set because they assume a specific data shape (recovery curves
starting negative/zero and rising, or a known spectral-diffusion contribution)
rather than being interchangeable candidates for a generic decay; callers pass
`models=[...]` explicitly to compare those.

Logs a per-model R²/chi-squared comparison and the best model by R², and (if
`plot=True`) draws all successful fits overlaid on the data with a residuals
panel — mirroring `fit_multiple_shapes`/`_plot_comparison`.

## Error handling

Same posture as `lineshapes/fitting.py`: input-shape and point-count errors
raise `ValueError` immediately (caller mistakes); unsupported model names raise
`ValueError` listing the valid choices; `curve_fit` convergence failures are
caught and reported through `RelaxationFitResult.success = False` rather than
propagating an exception, so batch/comparison callers don't need their own
try/except.

## Testing plan

New `tests/test_relaxation.py` following the project's 4-level protocol:

- **smoke**: each model function evaluates correctly at known points (e.g.
  `mono_exponential(0, A, T) == A + offset`; `inversion_recovery(0, A, T1) ==
  -A + offset`); `fit_relaxation` recovers ground-truth parameters from
  noise-free synthetic data for every model.
- **standard**: fits on synthetic data with added Gaussian noise achieve
  R² > 0.99 and recover the generating parameters within ~5% at realistic SNR;
  `fit_multiple_decays` picks the correct generating model by R² when given
  data generated from each candidate model in turn.
- **deep**: edge cases — decay time much shorter/longer than the time range,
  `beta` near its bounds, zero/negative offset, masked points, fewer than 4
  points raises `ValueError`, unsupported model name raises `ValueError`.
- **scientific**: cross-check `mono_exponential` fit results against the
  closed-form linear regression of `log(y - offset)` vs. `t` (analytic
  solution for the unstretched, single-rate case).

## Documentation and examples

- NumPy-style docstrings throughout, physical units documented explicitly
  (time axis unit is caller-defined; document this clearly since there is no
  auto-detection).
- New example script `examples/clean/06_relaxation_fitting.py`:
  - Part 1: real data — mono-exponential and stretched-exponential fits on
    `examples/data/2020_10_DMTTFBr_T2EH_28dB_6K_20ns_40ns_hperpc.DTA` (T2 echo
    decay), compared via `fit_multiple_decays`.
  - Part 2: synthetic biexponential decay with noise, recovered via
    `fit_relaxation`.
  - Part 3: synthetic `gamma_gaussian_decay` data, recovered via
    `fit_relaxation`, demonstrating the Gamma0/GammaG model.
- `CLAUDE.md` architecture section gains a `relaxation/` entry alongside
  `lineshapes/`, describing the module and its two public functions.
- Sphinx API docs (`docs/api/`) gain a page for `epyr.relaxation`, following
  the existing per-module page pattern.
- Release notes updated for the version bump below.

## Versioning

Adding a new public module/feature is a minor version bump under semver:
`0.3.9` -> `0.4.0`, updated consistently in `pyproject.toml` and
`epyr/__init__.py` (`make check-version` must pass).

## Out of scope (explicitly deferred)

- Batch/2D fitting (one decay per row/column of a 2D dataset, e.g. to produce
  a T2(B) or T1(B) map). Noted as a natural follow-up once the 1D API is in
  place; not built now to avoid speculative interface design.
- ESEEM / oscillatory modulation models.
- Automatic time-unit detection.
