# Pseudomodulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `epyr.lineshapes.pseudo_modulation()`, a function that simulates field-modulation lock-in detection on an EPR spectrum recorded without modulation, following Hyde et al. (1990).

**Architecture:** One new module `epyr/lineshapes/pseudo_modulation.py` implementing the Bessel-transfer-function FFT filter, exported through `epyr/lineshapes/__init__.py`. Tests live in a new `tests/test_pseudo_modulation.py`, mirroring the marker style already used in `tests/test_relaxation.py` (`smoke`, `standard`, `scientific`).

**Tech Stack:** NumPy (`np.fft`), SciPy (`scipy.special.jv`), pytest with markers registered in `pytest.ini`.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-17-pseudomodulation-design.md`.
- Function signature: `pseudo_modulation(x, y, mod_amplitude, harmonic=1, pad=True) -> np.ndarray`.
- `x` must be 1D and uniformly spaced; `y` must be 1D, same length as `x`, real or complex.
- `mod_amplitude` must be > 0 (same unit as `x`, peak-to-peak).
- `harmonic` must be `1` or `2`.
- `pad=True` edge-pads `y` by `len(y)//2` samples on each side before the FFT convolution, and the result is cropped back to the original length.
- Bessel transfer function is `(1j ** harmonic) * scipy.special.jv(harmonic, k * mod_amplitude / 2)` — the `1j ** harmonic` phase factor is required so the IFFT result is real for both odd and even `harmonic` (see spec's Algorithm section for the Hermitian-symmetry derivation). Omitting it produces a purely-imaginary result for odd `harmonic`, silently discarded by `np.real`.
- Output dtype family matches input `y` (real in -> real out, complex in -> complex out).
- Docstrings follow the project style in `CLAUDE.md` (NumPy format, physical units stated, no banned hedge words, no em dash, references or formula in Notes).
- Pytest markers used: `smoke`, `standard`, `scientific` (already registered in `pytest.ini`, `--strict-markers` is enabled so no new marker names).

---

### Task 1: Module scaffold and core implementation (smoke test)

**Files:**
- Create: `epyr/lineshapes/pseudo_modulation.py`
- Create: `tests/test_pseudo_modulation.py`

**Interfaces:**
- Produces: `pseudo_modulation(x: np.ndarray, y: np.ndarray, mod_amplitude: float, harmonic: int = 1, pad: bool = True) -> np.ndarray`, importable as `from epyr.lineshapes.pseudo_modulation import pseudo_modulation` (not yet exported from the package `__init__.py` — that happens in Task 4).

- [ ] **Step 1: Write the failing smoke test**

Create `tests/test_pseudo_modulation.py`:

```python
"""Tests for epyr.lineshapes.pseudo_modulation."""

import numpy as np
import pytest

from epyr.lineshapes.gaussian import gaussian
from epyr.lineshapes.pseudo_modulation import pseudo_modulation


@pytest.mark.smoke
class TestPseudoModulationBasic:
    def test_basic_shape_and_finiteness(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)

        assert y_pm.shape == y.shape
        assert np.all(np.isfinite(y_pm))
        assert np.isrealobj(y_pm)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_pseudo_modulation.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epyr.lineshapes.pseudo_modulation'`

- [ ] **Step 3: Implement the module**

Create `epyr/lineshapes/pseudo_modulation.py`:

```python
"""
Pseudomodulation of EPR spectra.

Simulates the signal a lock-in amplifier would record under sinusoidal
field modulation, applied numerically to a spectrum recorded without field
modulation.

References:
    - Hyde, J.S., Pasenkiewicz-Gierula, M., Jesmanowicz, A., Antholine, W.E.,
      Appl. Magn. Reson. 1, 483-496 (1990).
"""

import numpy as np
from scipy import special

from ._validation import validate_abscissa


def pseudo_modulation(x, y, mod_amplitude, harmonic=1, pad=True):
    """
    Compute the pseudomodulated spectrum for simulated lock-in detection.

    Reproduces the signal a lock-in amplifier would record under sinusoidal
    field modulation of peak-to-peak amplitude `mod_amplitude`, applied to a
    spectrum recorded without field modulation.

    Parameters
    ----------
    x : np.ndarray
        Magnetic field axis, in Gauss (G) or millitesla (mT), uniformly
        spaced (ascending or descending).
    y : np.ndarray
        Spectrum values at each point in `x`, same length as `x`. Real or
        complex.
    mod_amplitude : float
        Peak-to-peak field modulation amplitude, in the same unit as `x`.
        Must be positive.
    harmonic : int, default=1
        Detection harmonic: 1 for first-harmonic (first-derivative-like)
        detection, 2 for second-harmonic detection.
    pad : bool, default=True
        Edge-pad `y` before the FFT convolution to suppress circular
        wraparound artifacts at the spectrum boundaries, then crop the
        result back to the original length.

    Returns
    -------
    np.ndarray
        Pseudomodulated spectrum, same length as `y`. Real if `y` is real,
        complex if `y` is complex.

    Raises
    ------
    ValueError
        If `x` and `y` have different lengths, `x` is not uniformly spaced,
        `mod_amplitude` is not positive, or `harmonic` is not 1 or 2.

    Notes
    -----
    For a sinusoidal modulation B(t) = B0 + (mod_amplitude/2) sin(theta),
    the signal detected at the n-th harmonic of the modulation frequency is

        S_n(B) = IFFT[ FFT(A)(k) * i^n * J_n(k * mod_amplitude / 2) ]

    where A is the input spectrum, k = 2*pi*fftfreq(N, dx) is the spatial
    frequency conjugate to the field axis, and J_n is the Bessel function of
    the first kind of order n. The i^n phase factor is required: J_n(-z) =
    (-1)^n J_n(z), so J_n(k * mod_amplitude/2) alone is odd in k for odd n,
    which breaks the Hermitian symmetry a real spectrum's FFT has and makes
    the IFFT purely imaginary. Multiplying by i^n restores Hermitian
    symmetry for every n, so the IFFT is real whenever A is real. In the
    small-amplitude limit, J_1(z) ~ z/2, so S_1(B) approaches
    (mod_amplitude/4) * dA/dB.

    Examples
    --------
    >>> x = np.linspace(-50, 50, 2000)
    >>> y = np.exp(-x**2 / (2 * 5.0**2))
    >>> y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)
    """
    x = validate_abscissa(x)
    y = np.asarray(y)
    _validate_pseudo_modulation_inputs(x, y, mod_amplitude, harmonic)

    dx = x[1] - x[0]
    is_complex_input = np.iscomplexobj(y)

    if pad:
        pad_width = len(y) // 2
        y_work = np.pad(y, pad_width, mode="edge")
    else:
        pad_width = 0
        y_work = y

    n = len(y_work)
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    transfer = (1j**harmonic) * special.jv(harmonic, k * mod_amplitude / 2)

    spectrum_fft = np.fft.fft(y_work)
    result = np.fft.ifft(spectrum_fft * transfer)

    if not is_complex_input:
        result = np.real(result)

    if pad_width:
        result = result[pad_width : pad_width + len(y)]

    return result


def _validate_pseudo_modulation_inputs(x, y, mod_amplitude, harmonic):
    """Validate pseudo_modulation input parameters."""
    if y.ndim != 1:
        raise ValueError("y must be a 1D array")
    if x.shape != y.shape:
        raise ValueError(f"Shape mismatch: x {x.shape} vs y {y.shape}")
    if len(x) < 2:
        raise ValueError("x must contain at least two points")

    steps = np.diff(x)
    if not np.allclose(steps, steps[0], rtol=1e-6):
        raise ValueError("x must be uniformly spaced")

    if mod_amplitude <= 0:
        raise ValueError("mod_amplitude must be positive")

    if harmonic not in (1, 2):
        raise ValueError("harmonic must be 1 or 2")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_pseudo_modulation.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add epyr/lineshapes/pseudo_modulation.py tests/test_pseudo_modulation.py
git commit -m "Add pseudo_modulation core implementation with smoke test"
```

---

### Task 2: Input validation and standard tests

**Files:**
- Modify: `tests/test_pseudo_modulation.py`

**Interfaces:**
- Consumes: `pseudo_modulation(x, y, mod_amplitude, harmonic=1, pad=True)` from Task 1, and `_validate_pseudo_modulation_inputs` raising `ValueError` for the four invalid-input cases.

- [ ] **Step 1: Write the failing standard tests**

Append to `tests/test_pseudo_modulation.py`:

```python
@pytest.mark.standard
class TestPseudoModulationValidation:
    def test_harmonic_2_runs_and_differs_from_harmonic_1(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm_1 = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)
        y_pm_2 = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=2)

        assert y_pm_2.shape == y.shape
        assert np.all(np.isfinite(y_pm_2))
        assert not np.allclose(y_pm_1, y_pm_2)

    def test_pad_false_same_shape(self):
        x = np.linspace(-50, 50, 2000)
        y = gaussian(x, 0.0, 8.0)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1, pad=False)

        assert y_pm.shape == y.shape
        assert np.all(np.isfinite(y_pm))

    def test_complex_input_stays_complex(self):
        x = np.linspace(-50, 50, 500)
        y = gaussian(x, 0.0, 8.0).astype(complex) * (1 + 0.3j)

        y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)

        assert np.iscomplexobj(y_pm)
        assert y_pm.shape == y.shape

    def test_nonuniform_x_raises(self):
        x = np.concatenate([np.linspace(-50, -10, 500), np.linspace(-9, 50, 300)])
        y = gaussian(x, 0.0, 8.0)

        with pytest.raises(ValueError, match="uniformly spaced"):
            pseudo_modulation(x, y, mod_amplitude=2.0)

    @pytest.mark.parametrize("bad_amplitude", [0.0, -1.0])
    def test_invalid_mod_amplitude_raises(self, bad_amplitude):
        x = np.linspace(-50, 50, 500)
        y = gaussian(x, 0.0, 8.0)

        with pytest.raises(ValueError, match="mod_amplitude"):
            pseudo_modulation(x, y, mod_amplitude=bad_amplitude)

    @pytest.mark.parametrize("bad_harmonic", [0, 3, -1])
    def test_invalid_harmonic_raises(self, bad_harmonic):
        x = np.linspace(-50, 50, 500)
        y = gaussian(x, 0.0, 8.0)

        with pytest.raises(ValueError, match="harmonic"):
            pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=bad_harmonic)

    def test_shape_mismatch_raises(self):
        x = np.linspace(-50, 50, 500)
        y = np.zeros(400)

        with pytest.raises(ValueError, match="Shape mismatch"):
            pseudo_modulation(x, y, mod_amplitude=2.0)
```

- [ ] **Step 2: Run tests to verify they fail or pass as expected**

Run: `pytest tests/test_pseudo_modulation.py -v`
Expected: All pass immediately, since `_validate_pseudo_modulation_inputs` was already implemented in Task 1 — this step is a regression check that Task 1's validation logic is correct, not a red/green cycle. If any of these fail, fix `_validate_pseudo_modulation_inputs` in `epyr/lineshapes/pseudo_modulation.py` before proceeding.

- [ ] **Step 3: Commit**

```bash
git add tests/test_pseudo_modulation.py
git commit -m "Add standard validation tests for pseudo_modulation"
```

---

### Task 3: Scientific validation against the analytic small-amplitude limit

**Files:**
- Modify: `tests/test_pseudo_modulation.py`

**Interfaces:**
- Consumes: `pseudo_modulation` from Task 1, `gaussian` from `epyr.lineshapes.gaussian` (already imported).

- [ ] **Step 1: Write the failing scientific test**

Append to `tests/test_pseudo_modulation.py`:

```python
@pytest.mark.scientific
class TestPseudoModulationAnalyticLimit:
    def test_first_harmonic_small_amplitude_matches_derivative(self):
        x = np.linspace(-50, 50, 4000)
        center, width = 0.0, 8.0
        y = gaussian(x, center, width, derivative=0)
        dy_analytic = gaussian(x, center, width, derivative=1)

        eps = 1e-3
        y_pm = pseudo_modulation(x, y, mod_amplitude=eps, harmonic=1)
        scaled = y_pm / (eps / 4)

        # Stay away from the outermost few points where the FFT treats
        # the (non-periodic) signal as if it wrapped around.
        mask = np.abs(x) < 30
        max_dy = np.max(np.abs(dy_analytic[mask]))

        assert np.allclose(scaled[mask], dy_analytic[mask], atol=0.02 * max_dy)
```

- [ ] **Step 2: Run test to verify it passes**

Run: `pytest tests/test_pseudo_modulation.py -v -m scientific`
Expected: PASS. If it fails, check the sign and scale convention: `gaussian(..., derivative=1)` returns `dA/dx` directly (no extra minus sign — verified against `epyr/lineshapes/gaussian.py`'s `_compute_gaussian_components`), and the transfer-function derivation in the module docstring's Notes section gives `S_1(B) -> (mod_amplitude/4) * dA/dB` as `mod_amplitude -> 0`.

- [ ] **Step 3: Run the full new test file**

Run: `pytest tests/test_pseudo_modulation.py -v`
Expected: All tests (smoke, standard, scientific) PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/test_pseudo_modulation.py
git commit -m "Add scientific validation test for pseudo_modulation small-amplitude limit"
```

---

### Task 4: Export from the package and final verification

**Files:**
- Modify: `epyr/lineshapes/__init__.py`

**Interfaces:**
- Produces: `pseudo_modulation` importable as `from epyr.lineshapes import pseudo_modulation` (and transitively `from epyr import lineshapes; lineshapes.pseudo_modulation`).

- [ ] **Step 1: Add the import**

In `epyr/lineshapes/__init__.py`, change:

```python
from .lshape import lshape, pseudo_voigt
from .voigtian import voigtian
```

to:

```python
from .lshape import lshape, pseudo_voigt
from .pseudo_modulation import pseudo_modulation
from .voigtian import voigtian
```

- [ ] **Step 2: Add to `__all__`**

In the same file, change:

```python
    "pseudo_voigt",
    "convspec",
```

to:

```python
    "pseudo_voigt",
    "convspec",
    "pseudo_modulation",
```

- [ ] **Step 3: Update the module docstring function list**

In the same file's top-level docstring, change:

```python
    convspec: Spectrum convolution
```

to:

```python
    convspec: Spectrum convolution
    pseudo_modulation: Simulated field-modulation lock-in detection
```

- [ ] **Step 4: Write and run an import-path regression test**

Add to `tests/test_pseudo_modulation.py`:

```python
@pytest.mark.smoke
def test_pseudo_modulation_importable_from_package():
    from epyr.lineshapes import pseudo_modulation as exported

    assert exported is pseudo_modulation
```

Run: `pytest tests/test_pseudo_modulation.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Run the full lineshapes test suite**

Run: `pytest tests/test_lineshapes.py tests/test_pseudo_modulation.py -v`
Expected: All tests PASS, no import errors, no collection errors.

- [ ] **Step 6: Run code quality checks**

Run: `make format && make quality`
Expected: `black` and `isort` report no changes needed (or auto-fix cleanly); `make quality` (lint, type-check, security) passes with no new errors attributable to `epyr/lineshapes/pseudo_modulation.py`.

- [ ] **Step 7: Run the smoke suite**

Run: `pytest -m smoke`
Expected: PASS, including the two new smoke tests.

- [ ] **Step 8: Commit**

```bash
git add epyr/lineshapes/__init__.py tests/test_pseudo_modulation.py
git commit -m "Export pseudo_modulation from epyr.lineshapes"
```

---

## Post-plan check

After Task 4, confirm the example from the spec runs end-to-end:

```bash
python -c "
import numpy as np
from epyr.lineshapes import pseudo_modulation

x = np.linspace(-50, 50, 2000)
y = np.exp(-x**2 / (2 * 5.0**2))
y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)
print(y_pm.shape, y_pm.dtype, np.all(np.isfinite(y_pm)))
"
```

Expected output: `(2000,) float64 True`.
