# FFT Preprocessing — Design Spec
Date: 2026-07-02
Status: approved

## Context

EPR time-domain experiments (ESEEM, HYSCORE, Rabi) require a preprocessing
pipeline before FFT: baseline subtraction, apodization, and zero-padding. The
existing `signalprocessing` module bundles these steps inside
`analyze_frequencies()` with only mean-level DC removal. Users need standalone,
chainable functions with explicit control at each step.

## Scope

Three new public functions in a dedicated module
`epyr/signalprocessing/preprocessing.py`. They are exported via
`signalprocessing/__init__.py`. No changes to `analyze_frequencies()` or
`analyze_frequencies_2d()`.

## Architecture

```
epyr/signalprocessing/
    preprocessing.py      ← new file (this spec)
    apowin.py             ← unchanged, used internally
    frequency_analysis.py ← unchanged
    __init__.py           ← updated to export the three new functions
```

## Functions

### `remove_baseline`

```python
def remove_baseline(
    time: np.ndarray,
    signal: np.ndarray,
    method: str = "polynomial",
    order: int = 3,
    end_fraction: float = 0.15,
    axis: int = 1,
    plot: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
```

Fits and subtracts a slowly-varying background from a time-domain EPR signal.

**Parameters**

| Parameter | Description |
|---|---|
| `time` | Time axis, 1D array, any unit |
| `signal` | Signal array, 1D or 2D |
| `method` | `'polynomial'` or `'exponential'` |
| `order` | Polynomial degree (ignored for `'exponential'`) |
| `end_fraction` | Fraction of trailing points used to anchor the fit (default 0.15 = last 15 %) |
| `axis` | For 2D data: axis along which baseline is removed (default 1, i.e. row-wise) |
| `plot` | Show before/after figure |

**Returns** `(signal_corrected, baseline)` — both arrays have the same shape as `signal`.

**Implementation notes**

- 1D: fit polynomial/exponential on the last `round(end_fraction * N)` points, evaluated on the full time axis, then subtracted.
- Exponential fit: `A * exp(-t / tau) + C` via `scipy.optimize.curve_fit` with robust initial guess from the signal range.
- 2D: apply the 1D procedure independently to each slice along `axis`.
- If `curve_fit` fails for exponential, raise `RuntimeError` with a clear message.

---

### `apodize`

```python
def apodize(
    signal: np.ndarray,
    window: str = "hann",
    alpha: Optional[float] = None,
    half_window: Optional[str] = None,
    axis: Union[int, str] = "both",
    plot: bool = False,
) -> np.ndarray:
```

Applies an apodization window to reduce spectral leakage from signal truncation.

**Parameters**

| Parameter | Description |
|---|---|
| `signal` | Signal array, 1D or 2D |
| `window` | Window type: any key accepted by `apowin()` |
| `alpha` | Shape parameter for `kaiser` and `gaussian` windows |
| `half_window` | `None` (symmetric), `'left'`, or `'right'` |
| `axis` | For 2D: `0`, `1`, or `'both'` (outer product of two 1D windows) |
| `plot` | Show original signal, window shape, and windowed signal |

**Returns** `signal_windowed`, same shape as `signal`.

**Implementation notes**

- Delegates window generation to `apowin()` from `apowin.py`.
- For 2D with `axis='both'`: `window_2d = np.outer(w1, w2)` where `w1`, `w2` are generated for the respective axis lengths.
- For 2D with `axis=0` or `axis=1`: broadcast a 1D window across the perpendicular axis.
- `half_window` is passed directly to `apowin()`.

---

### `zero_pad`

```python
def zero_pad(
    signal: np.ndarray,
    factor: Optional[int] = None,
    n_points: Optional[int] = None,
    axis: Union[int, str] = -1,
    plot: bool = False,
) -> np.ndarray:
```

Pads a signal with trailing zeros to increase frequency resolution after FFT.

**Parameters**

| Parameter | Description |
|---|---|
| `signal` | Signal array, 1D or 2D |
| `factor` | Multiplicative factor: output length = `factor * N` along `axis` |
| `n_points` | Absolute output length along `axis` |
| `axis` | Axis to pad: `0`, `1`, `-1`, or `'both'` (2D only) |
| `plot` | Show original and padded signal |

Exactly one of `factor` or `n_points` must be provided; passing both or neither raises `ValueError`.

For `'both'` (2D only): `factor` or `n_points` applies to each axis independently (same value for both).

**Returns** `signal_padded`, same dtype as input.

---

## Chaining example

```python
from epyr import eprload
from epyr.signalprocessing import remove_baseline, apodize, zero_pad, analyze_frequencies

time, signal, params, _ = eprload("eseem.DTA")

signal_bl, baseline = remove_baseline(time, signal, method="exponential", plot=True)
signal_apo = apodize(signal_bl, window="hann", half_window="right", plot=True)
signal_zp  = zero_pad(signal_apo, factor=4, plot=True)

result = analyze_frequencies(time, signal_zp, window=None, remove_dc=False)
```

## Tests

New test file `tests/test_preprocessing.py` with marks `smoke`, `standard`, `deep`.

| Test | Mark |
|---|---|
| `test_remove_baseline_polynomial_1d` | smoke |
| `test_remove_baseline_exponential_1d` | smoke |
| `test_remove_baseline_2d` | standard |
| `test_apodize_1d_full` | smoke |
| `test_apodize_1d_half` | standard |
| `test_apodize_2d_both_axes` | standard |
| `test_zero_pad_factor` | smoke |
| `test_zero_pad_n_points` | smoke |
| `test_zero_pad_both_factor_n_points_raises` | smoke |
| `test_zero_pad_2d` | standard |
| `test_chain_eseem_pipeline` | deep |

## Out of scope

- Modification of `analyze_frequencies()` or `analyze_frequencies_2d()`.
- Interactive baseline selection (GUI).
- Spline baselines.
- Automatic baseline method selection.
