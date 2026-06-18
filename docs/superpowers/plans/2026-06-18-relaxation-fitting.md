# Relaxation (T1/T2) Fitting Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new `epyr/relaxation/` package providing decay/recovery model
functions and a curve-fitting layer (`fit_relaxation`, `fit_multiple_decays`)
for extracting T1/T2 relaxation times from time-domain EPR signals.

**Architecture:** Mirrors the existing `epyr/lineshapes/` package: a
`models.py` module of standalone NumPy functions, and a `fitting.py` module
with a `RelaxationFitResult` dataclass and two public entry points, built
incrementally (one model family at a time) on top of a shared `curve_fit`
based engine.

**Tech Stack:** NumPy, SciPy (`scipy.optimize.curve_fit`), Matplotlib,
pytest (smoke/standard/deep/scientific markers already registered in
`pytest.ini`).

**Reference spec:** `docs/superpowers/specs/2026-06-18-relaxation-fitting-design.md`

---

## Before you start

All new test code in this plan goes in a single new file,
`tests/test_relaxation.py`. Each task appends to it; run only the tests
relevant to that task, but the full file should pass after every task.

The numeric tolerances and synthetic parameters used in the tests below were
verified against a working prototype of this exact fitting engine before
writing this plan (noise-free synthetic fits converge to machine precision;
noisy fits and model-comparison rankings were checked with fixed RNG seeds).
Use the exact values given; do not "simplify" the synthetic data parameters.

---

### Task 1: Package skeleton and decay/recovery model functions

**Files:**
- Create: `epyr/relaxation/__init__.py`
- Create: `epyr/relaxation/models.py`
- Modify: `epyr/__init__.py:6`
- Create: `tests/test_relaxation.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_relaxation.py`:

```python
"""Tests for epyr.relaxation: T1/T2 decay and recovery models and fitting."""

import numpy as np
import pytest

from epyr.relaxation import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)


@pytest.mark.smoke
class TestModelFunctions:
    def test_mono_exponential_at_zero(self):
        v = mono_exponential(0.0, amplitude=5.0, T=10.0, offset=1.0)
        assert v == pytest.approx(6.0)

    def test_mono_exponential_decays_to_offset(self):
        t = np.array([0.0, 1e6])
        y = mono_exponential(t, amplitude=5.0, T=10.0, offset=1.0)
        assert y[-1] == pytest.approx(1.0, abs=1e-6)

    def test_stretched_exponential_reduces_to_mono_at_beta_one(self):
        t = np.linspace(0, 50, 25)
        y_stretched = stretched_exponential(
            t, amplitude=3.0, T=8.0, beta=1.0, offset=0.5
        )
        y_mono = mono_exponential(t, amplitude=3.0, T=8.0, offset=0.5)
        np.testing.assert_allclose(y_stretched, y_mono)

    def test_biexponential_at_zero(self):
        v = biexponential(
            0.0, amplitude1=2.0, tau1=5.0, amplitude2=3.0, tau2=50.0, offset=1.0
        )
        assert v == pytest.approx(6.0)

    def test_inversion_recovery_at_zero_is_negative_amplitude(self):
        v = inversion_recovery(0.0, amplitude=4.0, T1=20.0, offset=0.5)
        assert v == pytest.approx(0.5 - 4.0)

    def test_inversion_recovery_long_time_recovers_amplitude(self):
        v = inversion_recovery(np.array([1e6]), amplitude=4.0, T1=20.0, offset=0.5)
        assert v[0] == pytest.approx(4.5, abs=1e-6)

    def test_saturation_recovery_at_zero_is_offset(self):
        v = saturation_recovery(0.0, amplitude=4.0, T1=20.0, offset=0.5)
        assert v == pytest.approx(0.5)

    def test_gamma_gaussian_decay_at_zero(self):
        v = gamma_gaussian_decay(0.0, amplitude=2.0, Gamma0=0.1, GammaG=0.05, offset=0.2)
        assert v == pytest.approx(2.2)

    def test_gamma_gaussian_decay_array_shape_and_dtype(self):
        t = np.linspace(0, 100, 50)
        y = gamma_gaussian_decay(t, amplitude=1.0, Gamma0=0.02, GammaG=0.01)
        assert y.shape == t.shape
        assert y.dtype == np.float64
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_relaxation.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epyr.relaxation'`

- [ ] **Step 3: Implement the model functions**

Create `epyr/relaxation/models.py`:

```python
"""
Time-domain relaxation models for EPR T1/T2 measurements.

Mathematical Background:
    Mono-exponential decay/recovery, stretched exponential, bi-exponential,
    and a combined homogeneous/spectral-diffusion decay model used to
    extract spin relaxation times (T1, T2) from pulsed EPR time traces.

References:
    Eaton, S.S. and Eaton, G.R., "Relaxation Times of Organic Radicals and
    Transition Metal Ions", Biol. Magn. Reson., 19, 2000.
"""

import numpy as np


def mono_exponential(
    t: np.ndarray,
    amplitude: float,
    T: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a single-rate exponential decay.

    V(t) = amplitude * exp(-t / T) + offset

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller (e.g. ns or us).
    amplitude : float
        Signal amplitude at t = 0, relative to offset.
    T : float
        Decay (or recovery) time constant, in the same unit as t.
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * np.exp(-t / T) + offset


def stretched_exponential(
    t: np.ndarray,
    amplitude: float,
    T: float,
    beta: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a stretched-exponential decay.

    V(t) = amplitude * exp(-(t / T) ** beta) + offset

    beta = 1 reduces to a mono-exponential decay. beta < 1 describes a
    distribution of relaxation rates, common in pulsed EPR echo decays.

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller.
    amplitude : float
        Signal amplitude at t = 0, relative to offset.
    T : float
        Characteristic decay time constant, in the same unit as t.
    beta : float
        Stretching exponent (dimensionless).
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * np.exp(-((t / T) ** beta)) + offset


def biexponential(
    t: np.ndarray,
    amplitude1: float,
    tau1: float,
    amplitude2: float,
    tau2: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a two-component exponential decay.

    V(t) = amplitude1 * exp(-t / tau1) + amplitude2 * exp(-t / tau2) + offset

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller.
    amplitude1 : float
        Amplitude of the first component at t = 0.
    tau1 : float
        Decay time constant of the first component, same unit as t.
    amplitude2 : float
        Amplitude of the second component at t = 0.
    tau2 : float
        Decay time constant of the second component, same unit as t.
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.

    Notes
    -----
    tau1 and tau2 are not necessarily T1 or T2: either component can
    represent either relaxation process, depending on the experiment.
    """
    t = np.asarray(t, dtype=float)
    return (
        amplitude1 * np.exp(-t / tau1) + amplitude2 * np.exp(-t / tau2) + offset
    )


def inversion_recovery(
    t: np.ndarray,
    amplitude: float,
    T1: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute an inversion-recovery T1 curve.

    V(t) = amplitude * (1 - 2 * exp(-t / T1)) + offset

    Parameters
    ----------
    t : np.ndarray
        Delay time after inversion, in the unit chosen by the caller.
    amplitude : float
        Fully-recovered signal amplitude relative to offset.
    T1 : float
        Spin-lattice relaxation time, same unit as t.
    offset : float, optional
        Baseline signal level (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * (1.0 - 2.0 * np.exp(-t / T1)) + offset


def saturation_recovery(
    t: np.ndarray,
    amplitude: float,
    T1: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute a saturation-recovery T1 curve.

    V(t) = amplitude * (1 - exp(-t / T1)) + offset

    Parameters
    ----------
    t : np.ndarray
        Delay time after saturation, in the unit chosen by the caller.
    amplitude : float
        Fully-recovered signal amplitude relative to offset.
    T1 : float
        Spin-lattice relaxation time, same unit as t.
    offset : float, optional
        Baseline signal level (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return amplitude * (1.0 - np.exp(-t / T1)) + offset


def gamma_gaussian_decay(
    t: np.ndarray,
    amplitude: float,
    Gamma0: float,
    GammaG: float,
    offset: float = 0.0,
) -> np.ndarray:
    """
    Compute an echo decay with homogeneous and spectral-diffusion terms.

    V(t) = amplitude * exp(-Gamma0 * t) * exp(-(GammaG * t) ** 2) + offset

    Gamma0 captures the homogeneous (exponential) relaxation rate. GammaG
    captures an additional Gaussian-shaped decay from spectral diffusion.

    Parameters
    ----------
    t : np.ndarray
        Time axis, in the unit chosen by the caller.
    amplitude : float
        Signal amplitude at t = 0, relative to offset.
    Gamma0 : float
        Homogeneous decay rate, in inverse time units (1 / t unit).
    GammaG : float
        Gaussian (spectral-diffusion) decay rate, in inverse time units.
    offset : float, optional
        Asymptotic signal level as t approaches infinity (default: 0.0).

    Returns
    -------
    np.ndarray
        Signal values, same shape as t.
    """
    t = np.asarray(t, dtype=float)
    return (
        amplitude * np.exp(-Gamma0 * t) * np.exp(-((GammaG * t) ** 2)) + offset
    )
```

- [ ] **Step 4: Create the package `__init__.py`**

Create `epyr/relaxation/__init__.py`:

```python
"""
EPyR Tools - Relaxation Fitting Module

Time-domain decay and recovery models for T1/T2 EPR relaxation measurements,
with a curve-fitting layer mirroring epyr.lineshapes.fitting.
"""

from .models import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)

__all__ = [
    "mono_exponential",
    "stretched_exponential",
    "biexponential",
    "inversion_recovery",
    "saturation_recovery",
    "gamma_gaussian_decay",
]
```

- [ ] **Step 5: Wire the package into `epyr/__init__.py`**

In `epyr/__init__.py`, line 6 currently reads:

```python
from . import lineshapes, signalprocessing
```

Change it to:

```python
from . import lineshapes, relaxation, signalprocessing
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS (9 tests)

- [ ] **Step 7: Commit**

```bash
git add epyr/relaxation/__init__.py epyr/relaxation/models.py epyr/__init__.py tests/test_relaxation.py
git commit -m "Add relaxation model functions (mono/stretched/bi-exponential, recovery, Gamma0/GammaG)"
```

---

### Task 2: Fitting engine core, `RelaxationFitResult`, and `mono_exponential` support

**Files:**
- Create: `epyr/relaxation/fitting.py`
- Modify: `epyr/relaxation/__init__.py`
- Modify: `tests/test_relaxation.py`

- [ ] **Step 1: Append the failing tests**

Append to `tests/test_relaxation.py` (add `import matplotlib` at the very
top of the file, before the `numpy`/`pytest` imports, and call
`matplotlib.use("Agg")` immediately after, so plotting tests do not try to
open an interactive window):

```python
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from epyr.relaxation import (
    RelaxationFitResult,
    biexponential,
    fit_relaxation,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)
```

(This replaces the original `import numpy as np` / `import pytest` /
`from epyr.relaxation import (...)` block at the top of the file with the
version above, which adds the matplotlib setup and the new names.)

Then append:

```python
@pytest.mark.smoke
class TestFitRelaxationMonoExponential:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=15.0, offset=1.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.parameters["amplitude"] == pytest.approx(5.0, rel=1e-3)
        assert result.parameters["T"] == pytest.approx(15.0, rel=1e-3)
        assert result.parameters["offset"] == pytest.approx(1.0, abs=1e-3)
        assert result.r_squared > 0.999

    def test_returns_relaxation_fit_result(self):
        t = np.linspace(0, 50, 100)
        y = mono_exponential(t, amplitude=2.0, T=10.0, offset=0.0)
        result = fit_relaxation(t, y, plot=False)
        assert isinstance(result, RelaxationFitResult)
        assert result.model == "mono_exponential"

    def test_plot_true_does_not_raise(self):
        t = np.linspace(0, 50, 80)
        y = mono_exponential(t, amplitude=3.0, T=12.0, offset=0.5)
        result = fit_relaxation(t, y, plot=True)
        assert result.success
        plt.close("all")


@pytest.mark.standard
class TestFitRelaxationNoisy:
    def test_noisy_mono_exponential(self):
        rng = np.random.default_rng(0)
        t = np.linspace(0, 100, 300)
        y_true = mono_exponential(t, amplitude=5.0, T=20.0, offset=1.0)
        y = y_true + 0.02 * 5.0 * rng.standard_normal(t.size)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.r_squared > 0.99
        assert result.parameters["T"] == pytest.approx(20.0, rel=0.05)


@pytest.mark.smoke
class TestFitRelaxationValidation:
    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError):
            fit_relaxation(np.array([1, 2, 3, 4]), np.array([1, 2, 3]), plot=False)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError):
            fit_relaxation(np.array([1, 2, 3]), np.array([1, 2, 3]), plot=False)

    def test_unsupported_model_raises(self):
        t = np.linspace(0, 10, 20)
        y = np.ones_like(t)
        with pytest.raises(ValueError):
            fit_relaxation(t, y, model="not_a_model", plot=False)

    def test_mask_excludes_points(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=15.0, offset=1.0)
        mask = np.ones_like(t, dtype=bool)
        mask[50:60] = False
        result = fit_relaxation(t, y, mask=mask, plot=False)
        assert result.success
        assert result.t_fit.size == 190
        assert result.parameters["T"] == pytest.approx(15.0, rel=1e-3)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_relaxation.py -v`
Expected: FAIL with `ImportError: cannot import name 'RelaxationFitResult' from 'epyr.relaxation'`

- [ ] **Step 3: Implement the fitting engine**

Create `epyr/relaxation/fitting.py`:

```python
"""
T1/T2 relaxation fitting module.

Fit time-domain EPR relaxation data (T1 recovery, T2 echo decay) with
mono-exponential, stretched-exponential, bi-exponential, inversion- or
saturation-recovery, or combined homogeneous/spectral-diffusion models.
"""

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from ..logging_config import get_logger
from .models import mono_exponential

logger = get_logger(__name__)

SUPPORTED_MODELS = [
    "mono_exponential",
]


@dataclass
class RelaxationFitResult:
    """
    Container for relaxation fit results.

    Attributes
    ----------
    model : str
        Name of the fitted relaxation model.
    parameters : dict
        Fitted parameter values keyed by parameter name.
    parameter_errors : dict
        Standard errors of fitted parameters (square root of covariance diagonal).
    fitted_curve : np.ndarray
        Model evaluated at the fitted points (t_fit).
    residuals : np.ndarray
        Data minus model at the fitted points.
    r_squared : float
        Coefficient of determination R-squared.
    chi_squared : float
        Reduced chi-squared: sum of squared residuals divided by degrees of freedom.
    success : bool
        True if curve_fit converged.
    message : str
        Convergence message or error description.
    covariance_matrix : np.ndarray or None
        Full parameter covariance matrix returned by curve_fit.
    t_fit : np.ndarray or None
        Time values used for fitting, after NaN removal and masking.
    """

    model: str
    parameters: Dict[str, float]
    parameter_errors: Dict[str, float]
    fitted_curve: np.ndarray
    residuals: np.ndarray
    r_squared: float
    chi_squared: float
    success: bool
    message: str
    covariance_matrix: Optional[np.ndarray] = None
    t_fit: Optional[np.ndarray] = None

    def summary(self) -> str:
        """Return a formatted string summarizing fit quality and parameters."""
        lines = [f"=== Relaxation Fit Results - {self.model} ==="]
        lines.append(f"Success: {self.success}")
        if not self.success:
            lines.append(f"Error: {self.message}")
            return "\n".join(lines)

        lines.append(f"R2 = {self.r_squared:.6f}")
        lines.append(f"chi2 = {self.chi_squared:.6f}")
        lines.append("\nParameters:")

        for param, value in self.parameters.items():
            if param in self.parameter_errors:
                error = self.parameter_errors[param]
                lines.append(f"  {param}: {value:.6g} +/- {error:.6g}")
            else:
                lines.append(f"  {param}: {value:.6g}")

        return "\n".join(lines)


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
    """
    Fit time-domain relaxation data with the specified decay/recovery model.

    Parameters
    ----------
    t_data : np.ndarray
        Time axis, in the unit chosen by the caller (e.g. ns or us).
    y_data : np.ndarray
        Relaxation signal (real-valued; take np.abs() of a complex echo
        signal before calling this function).
    model : str, optional
        Relaxation model name (default: 'mono_exponential'). See
        SUPPORTED_MODELS for the full list as models are added.
    initial_params : dict, optional
        Initial parameter guesses. Auto-estimated from data if None.
    bounds : dict, optional
        Parameter bounds as {name: (lower, upper)}, overriding data-derived
        defaults.
    mask : np.ndarray of bool, optional
        Boolean array of the same length as t_data. True selects a point for
        fitting; False excludes it. If None, all non-NaN points are used.
    plot : bool, optional
        Display a fit plot with residuals panel (default: True).
    time_unit : str, optional
        Cosmetic time-unit label for the plot axis and summary (e.g. 'ns').
        Does not affect fitting (default: '').
    **fit_kwargs
        Additional keyword arguments passed to scipy.optimize.curve_fit.

    Returns
    -------
    RelaxationFitResult
        Fit parameters, errors, statistics, fitted curve, and residuals.

    Examples
    --------
    >>> from epyr.relaxation import fit_relaxation
    >>> import numpy as np
    >>> t = np.linspace(0, 100, 200)
    >>> y = 5.0 * np.exp(-t / 15.0) + 1.0
    >>> result = fit_relaxation(t, y, 'mono_exponential', plot=False)
    >>> print(result.summary())
    """

    t_data = np.asarray(t_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    if len(t_data) != len(y_data):
        raise ValueError("t_data and y_data must have same length")

    if len(t_data) < 4:
        raise ValueError("Need at least 4 data points for fitting")

    if model not in SUPPORTED_MODELS:
        raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")

    valid_mask = ~(np.isnan(t_data) | np.isnan(y_data))
    if not np.any(valid_mask):
        raise ValueError("No valid data points (all NaN)")

    t_clean = t_data[valid_mask]
    y_clean = y_data[valid_mask]

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != t_data.shape:
            raise ValueError(
                f"mask shape {mask.shape} does not match t_data shape {t_data.shape}"
            )
        fit_mask = mask[valid_mask]
        if not np.any(fit_mask):
            raise ValueError("No valid data points after applying mask")
        t_fit = t_clean[fit_mask]
        y_fit = y_clean[fit_mask]
    else:
        t_fit = t_clean
        y_fit = y_clean

    fit_func, param_names, param_bounds = _get_fit_function(model)

    if initial_params is None:
        initial_params = _estimate_initial_params(t_fit, y_fit, model)

    initial_params = _validate_initial_params(initial_params, param_names, t_fit, y_fit)

    lower_bounds, upper_bounds = _setup_bounds(
        param_names, bounds, param_bounds, initial_params, t_fit, y_fit
    )

    p0 = [initial_params[name] for name in param_names]
    bounds_tuple = (lower_bounds, upper_bounds)

    default_kwargs = {"maxfev": 10000, "method": "trf"}
    default_kwargs.update(fit_kwargs)

    try:
        popt, pcov = curve_fit(
            fit_func, t_fit, y_fit, p0=p0, bounds=bounds_tuple, **default_kwargs
        )

        y_fitted = fit_func(t_fit, *popt)
        residuals = y_fit - y_fitted

        ss_res: float = float(np.sum(residuals**2))
        ss_tot: float = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        chi_squared = ss_res / (len(y_fit) - len(popt))

        param_dict = {name: value for name, value in zip(param_names, popt)}
        param_errors = {}

        if pcov is not None:
            param_std_errors = np.sqrt(np.diag(pcov))
            param_errors = {
                name: error for name, error in zip(param_names, param_std_errors)
            }

        result = RelaxationFitResult(
            model=model,
            parameters=param_dict,
            parameter_errors=param_errors,
            fitted_curve=y_fitted,
            residuals=residuals,
            r_squared=r_squared,
            chi_squared=chi_squared,
            success=True,
            message="Fit converged successfully",
            covariance_matrix=pcov,
            t_fit=t_fit,
        )

        if plot:
            _plot_fit_results(t_fit, y_fit, result, model, time_unit)

        return result

    except Exception as e:
        return RelaxationFitResult(
            model=model,
            parameters={},
            parameter_errors={},
            fitted_curve=np.array([]),
            residuals=np.array([]),
            r_squared=0.0,
            chi_squared=np.inf,
            success=False,
            message=str(e),
        )


def _get_fit_function(
    model: str,
) -> Tuple[Callable[..., np.ndarray], List[str], Dict[str, Tuple[float, float]]]:
    """
    Build the fit function, parameter names, and default bounds for a model.

    Parameters
    ----------
    model : str
        Relaxation model name.

    Returns
    -------
    fit_func : callable
        Function compatible with scipy.optimize.curve_fit, signature
        fit_func(t, *params).
    param_names : list of str
        Ordered parameter names matching fit_func's signature.
    param_bounds : dict
        Default bounds for each parameter as {name: (lower, upper)}.
    """

    if model == "mono_exponential":

        def fit_func(t, amplitude, T, offset):
            return mono_exponential(t, amplitude, T, offset)

        param_names = ["amplitude", "T", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    else:
        raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")

    return fit_func, param_names, param_bounds


def _estimate_decay_rate(
    t: np.ndarray, y: np.ndarray, offset_guess: float, amplitude_guess: float
) -> float:
    """
    Estimate a decay rate via log-linear regression on the part of the curve
    that lies well above the asymptotic offset.

    Parameters
    ----------
    t : np.ndarray
        Time axis.
    y : np.ndarray
        Relaxation signal.
    offset_guess : float
        Estimated asymptotic offset.
    amplitude_guess : float
        Estimated amplitude at t = 0, relative to offset.

    Returns
    -------
    float
        Estimated decay rate, in inverse time units. Falls back to the
        inverse of the time range if too few points qualify for the
        regression.
    """

    deviation = np.abs(y - offset_guess)
    threshold = 0.05 * abs(amplitude_guess) if amplitude_guess != 0 else 0.0
    usable = deviation > threshold

    if np.sum(usable) >= 2:
        log_dev = np.log(deviation[usable])
        slope, _ = np.polyfit(t[usable], log_dev, 1)
        if slope < 0:
            return float(-slope)

    time_span = float(t[-1] - t[0]) or 1.0
    return 1.0 / time_span


def _estimate_initial_params(
    t: np.ndarray, y: np.ndarray, model: str
) -> Dict[str, float]:
    """
    Estimate initial fit parameters from the data.

    Parameters
    ----------
    t : np.ndarray
        Time axis.
    y : np.ndarray
        Relaxation signal.
    model : str
        Relaxation model name.

    Returns
    -------
    dict
        Initial parameter estimates keyed by parameter name.
    """

    offset_guess = float(y[-1])
    span = float(y[0] - y[-1])
    amplitude_guess = span if span != 0 else float(y.max() - y.min()) or 1.0

    if model == "mono_exponential":
        rate_guess = _estimate_decay_rate(t, y, offset_guess, amplitude_guess)
        T_guess = 1.0 / rate_guess if rate_guess > 0 else (t[-1] - t[0]) / 2
        return {"amplitude": amplitude_guess, "T": T_guess, "offset": offset_guess}

    raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")


def _validate_initial_params(
    initial_params: Dict[str, float],
    param_names: List[str],
    t: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    """
    Fill in missing initial parameters with data-derived defaults.

    Parameters
    ----------
    initial_params : dict
        Partial or complete parameter dict.
    param_names : list of str
        Full list of parameters required by the fit function.
    t : np.ndarray
        Time axis, used to derive fallback values.
    y : np.ndarray
        Relaxation signal, used to derive fallback values.

    Returns
    -------
    dict
        Complete parameter dict with all entries in param_names present.
    """

    validated = initial_params.copy()
    time_span = float(t[-1] - t[0]) or 1.0

    for name in param_names:
        if name in validated:
            continue
        if name.startswith("amplitude"):
            validated[name] = float(y.max() - y.min()) or 1.0
        elif name in ("T", "T1") or name.startswith("tau"):
            validated[name] = time_span / 2
        elif name == "beta":
            validated[name] = 1.0
        elif name in ("Gamma0", "GammaG"):
            validated[name] = 1.0 / time_span
        elif name == "offset":
            validated[name] = float(y[-1])
        else:
            validated[name] = 1.0

    return validated


def _setup_bounds(
    param_names: List[str],
    user_bounds: Optional[Dict[str, Tuple[float, float]]],
    default_bounds: Dict[str, Tuple[float, float]],
    initial_params: Dict[str, float],
    t: np.ndarray,
    y: np.ndarray,
) -> Tuple[List[float], List[float]]:
    """
    Build ordered lower/upper bound lists for scipy.optimize.curve_fit.

    Parameters
    ----------
    param_names : list of str
        Ordered parameter names matching the fit function signature.
    user_bounds : dict or None
        Caller-supplied bounds as {name: (lower, upper)}.
    default_bounds : dict
        Default bounds per parameter as {name: (lower, upper)}.
    initial_params : dict
        Current initial values; used to ensure each value lies inside its bounds.
    t : np.ndarray
        Time axis, used to derive data-range bounds.
    y : np.ndarray
        Relaxation signal, used to derive data-range bounds.

    Returns
    -------
    lower_bounds : list of float
    upper_bounds : list of float
    """

    lower_bounds = []
    upper_bounds = []
    time_span = float(t[-1] - t[0]) or 1.0
    y_range = float(y.max() - y.min()) or 1.0

    for name in param_names:
        default_lower, default_upper = default_bounds[name]

        if user_bounds and name in user_bounds:
            lower, upper = user_bounds[name]
        else:
            lower, upper = default_lower, default_upper

        if name.startswith("amplitude"):
            if lower == -np.inf:
                lower = -10 * y_range
            if upper == np.inf:
                upper = 10 * y_range
        elif name in ("T", "T1") or name.startswith("tau"):
            if upper == np.inf:
                upper = 100 * time_span
        elif name in ("Gamma0", "GammaG"):
            dt = float(np.diff(t).mean()) or (time_span / len(t))
            if upper == np.inf:
                upper = 10.0 / dt
        elif name == "offset":
            if lower == -np.inf:
                lower = float(y.min()) - 10 * y_range
            if upper == np.inf:
                upper = float(y.max()) + 10 * y_range

        init_val = initial_params[name]
        if init_val <= lower:
            lower = init_val * 0.1 if init_val > 0 else (init_val * 10 if init_val < 0 else -1.0)
        if init_val >= upper:
            upper = init_val * 10 if init_val > 0 else (init_val * 0.1 if init_val < 0 else 1.0)

        lower_bounds.append(lower)
        upper_bounds.append(upper)

    return lower_bounds, upper_bounds


def _plot_fit_results(
    t: np.ndarray,
    y: np.ndarray,
    result: RelaxationFitResult,
    model: str,
    time_unit: str = "",
) -> None:
    """
    Plot data, fitted curve, and residuals for a single relaxation fit.

    Parameters
    ----------
    t : np.ndarray
        Time values used for fitting.
    y : np.ndarray
        Relaxation signal values used for fitting.
    result : RelaxationFitResult
        Fit result containing parameters, statistics, and curves.
    model : str
        Relaxation model name, shown in the plot title.
    time_unit : str, optional
        Unit label appended to the time axis (default: '').
    """

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [3, 1]}
    )

    time_label = f"Time ({time_unit})" if time_unit else "Time"

    ax1.plot(t, y, "o", markersize=4, alpha=0.7, label="Data", color="#1f77b4")
    ax1.plot(
        t, result.fitted_curve, "-", linewidth=2, label=f"{model} fit", color="#d62728"
    )
    ax1.set_xlabel(time_label)
    ax1.set_ylabel("Intensity")
    ax1.set_title(f"Relaxation Fitting - {model}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    results_lines = [f"R2 = {result.r_squared:.4f}", f"chi2 = {result.chi_squared:.2e}"]
    for param, value in result.parameters.items():
        value_str = f"{value:.4g}"
        if param in result.parameter_errors:
            results_lines.append(f"{param}: {value_str} +/- {result.parameter_errors[param]:.2g}")
        else:
            results_lines.append(f"{param}: {value_str}")

    ax1.text(
        0.98,
        0.98,
        "\n".join(results_lines),
        transform=ax1.transAxes,
        fontsize=9,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.8),
    )

    ax2.plot(t, result.residuals, "o-", markersize=3, alpha=0.7, color="#ff7f0e")
    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.set_xlabel(time_label)
    ax2.set_ylabel("Residuals")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
```

- [ ] **Step 4: Export the new names from the package**

In `epyr/relaxation/__init__.py`, replace the full contents with:

```python
"""
EPyR Tools - Relaxation Fitting Module

Time-domain decay and recovery models for T1/T2 EPR relaxation measurements,
with a curve-fitting layer mirroring epyr.lineshapes.fitting.
"""

from .fitting import RelaxationFitResult, fit_relaxation
from .models import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)

__all__ = [
    "mono_exponential",
    "stretched_exponential",
    "biexponential",
    "inversion_recovery",
    "saturation_recovery",
    "gamma_gaussian_decay",
    "fit_relaxation",
    "RelaxationFitResult",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS (all tests so far)

- [ ] **Step 6: Commit**

```bash
git add epyr/relaxation/fitting.py epyr/relaxation/__init__.py tests/test_relaxation.py
git commit -m "Add fit_relaxation engine with mono_exponential support"
```

---

### Task 3: Add `stretched_exponential` and `gamma_gaussian_decay` support

**Files:**
- Modify: `epyr/relaxation/fitting.py`
- Modify: `tests/test_relaxation.py`

- [ ] **Step 1: Append the failing tests**

Append to `tests/test_relaxation.py`:

```python
@pytest.mark.smoke
class TestFitRelaxationStretchedExponential:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0.1, 100, 200)
        y = stretched_exponential(t, amplitude=4.0, T=25.0, beta=0.7, offset=0.5)
        result = fit_relaxation(t, y, model="stretched_exponential", plot=False)
        assert result.success
        assert result.parameters["T"] == pytest.approx(25.0, rel=1e-3)
        assert result.parameters["beta"] == pytest.approx(0.7, rel=1e-3)
        assert result.r_squared > 0.999


@pytest.mark.smoke
class TestFitRelaxationGammaGaussianDecay:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = gamma_gaussian_decay(t, amplitude=3.0, Gamma0=0.05, GammaG=0.02, offset=0.2)
        result = fit_relaxation(t, y, model="gamma_gaussian_decay", plot=False)
        assert result.success
        assert result.parameters["Gamma0"] == pytest.approx(0.05, rel=1e-3)
        assert result.parameters["GammaG"] == pytest.approx(0.02, rel=1e-3)
        assert result.r_squared > 0.999
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_relaxation.py -v -k "Stretched or GammaGaussian"`
Expected: FAIL with `ValueError: Unsupported model: stretched_exponential...`

- [ ] **Step 3: Extend the model dispatch**

In `epyr/relaxation/fitting.py`, update the import line:

```python
from .models import gamma_gaussian_decay, mono_exponential, stretched_exponential
```

Update `SUPPORTED_MODELS`:

```python
SUPPORTED_MODELS = [
    "mono_exponential",
    "stretched_exponential",
    "gamma_gaussian_decay",
]
```

Replace `_get_fit_function` with:

```python
def _get_fit_function(
    model: str,
) -> Tuple[Callable[..., np.ndarray], List[str], Dict[str, Tuple[float, float]]]:
    """
    Build the fit function, parameter names, and default bounds for a model.

    Parameters
    ----------
    model : str
        Relaxation model name.

    Returns
    -------
    fit_func : callable
        Function compatible with scipy.optimize.curve_fit, signature
        fit_func(t, *params).
    param_names : list of str
        Ordered parameter names matching fit_func's signature.
    param_bounds : dict
        Default bounds for each parameter as {name: (lower, upper)}.
    """

    if model == "mono_exponential":

        def fit_func(t, amplitude, T, offset):
            return mono_exponential(t, amplitude, T, offset)

        param_names = ["amplitude", "T", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "stretched_exponential":

        def fit_func(t, amplitude, T, beta, offset):
            return stretched_exponential(t, amplitude, T, beta, offset)

        param_names = ["amplitude", "T", "beta", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "beta": (0.05, 5.0),
            "offset": (-np.inf, np.inf),
        }

    elif model == "gamma_gaussian_decay":

        def fit_func(t, amplitude, Gamma0, GammaG, offset):
            return gamma_gaussian_decay(t, amplitude, Gamma0, GammaG, offset)

        param_names = ["amplitude", "Gamma0", "GammaG", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "Gamma0": (0.0, np.inf),
            "GammaG": (0.0, np.inf),
            "offset": (-np.inf, np.inf),
        }

    else:
        raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")

    return fit_func, param_names, param_bounds
```

Replace `_estimate_initial_params` with:

```python
def _estimate_initial_params(
    t: np.ndarray, y: np.ndarray, model: str
) -> Dict[str, float]:
    """
    Estimate initial fit parameters from the data.

    Parameters
    ----------
    t : np.ndarray
        Time axis.
    y : np.ndarray
        Relaxation signal.
    model : str
        Relaxation model name.

    Returns
    -------
    dict
        Initial parameter estimates keyed by parameter name.
    """

    offset_guess = float(y[-1])
    span = float(y[0] - y[-1])
    amplitude_guess = span if span != 0 else float(y.max() - y.min()) or 1.0

    if model in ("mono_exponential", "stretched_exponential", "gamma_gaussian_decay"):
        rate_guess = _estimate_decay_rate(t, y, offset_guess, amplitude_guess)
        T_guess = 1.0 / rate_guess if rate_guess > 0 else (t[-1] - t[0]) / 2

        if model == "mono_exponential":
            return {"amplitude": amplitude_guess, "T": T_guess, "offset": offset_guess}

        if model == "stretched_exponential":
            return {
                "amplitude": amplitude_guess,
                "T": T_guess,
                "beta": 1.0,
                "offset": offset_guess,
            }

        return {
            "amplitude": amplitude_guess,
            "Gamma0": rate_guess,
            "GammaG": rate_guess,
            "offset": offset_guess,
        }

    raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")
```

`_validate_initial_params` and `_setup_bounds` already handle `beta`,
`Gamma0`, and `GammaG` by name (see Task 2); no changes needed there.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS (all tests so far)

- [ ] **Step 5: Commit**

```bash
git add epyr/relaxation/fitting.py tests/test_relaxation.py
git commit -m "Add stretched_exponential and gamma_gaussian_decay fit support"
```

---

### Task 4: Add `biexponential`, `inversion_recovery`, `saturation_recovery` support

**Files:**
- Modify: `epyr/relaxation/fitting.py`
- Modify: `tests/test_relaxation.py`

- [ ] **Step 1: Append the failing tests**

Append to `tests/test_relaxation.py`:

```python
@pytest.mark.smoke
class TestFitRelaxationBiexponential:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 200, 300)
        y = biexponential(
            t, amplitude1=3.0, tau1=8.0, amplitude2=2.0, tau2=80.0, offset=0.3
        )
        result = fit_relaxation(t, y, model="biexponential", plot=False)
        assert result.success
        assert result.parameters["tau1"] == pytest.approx(8.0, rel=1e-2)
        assert result.parameters["tau2"] == pytest.approx(80.0, rel=1e-2)
        assert result.r_squared > 0.999


@pytest.mark.smoke
class TestFitRelaxationInversionRecovery:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = inversion_recovery(t, amplitude=4.0, T1=20.0, offset=0.5)
        result = fit_relaxation(t, y, model="inversion_recovery", plot=False)
        assert result.success
        assert result.parameters["T1"] == pytest.approx(20.0, rel=1e-3)
        assert result.parameters["amplitude"] == pytest.approx(4.0, rel=1e-3)
        assert result.r_squared > 0.999


@pytest.mark.smoke
class TestFitRelaxationSaturationRecovery:
    def test_recovers_noise_free_parameters(self):
        t = np.linspace(0, 100, 200)
        y = saturation_recovery(t, amplitude=4.0, T1=20.0, offset=0.5)
        result = fit_relaxation(t, y, model="saturation_recovery", plot=False)
        assert result.success
        assert result.parameters["T1"] == pytest.approx(20.0, rel=1e-3)
        assert result.parameters["amplitude"] == pytest.approx(4.0, rel=1e-3)
        assert result.r_squared > 0.999
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_relaxation.py -v -k "Biexponential or InversionRecovery or SaturationRecovery"`
Expected: FAIL with `ValueError: Unsupported model: biexponential...`

- [ ] **Step 3: Extend the model dispatch**

In `epyr/relaxation/fitting.py`, update the import line:

```python
from .models import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)
```

Update `SUPPORTED_MODELS`:

```python
SUPPORTED_MODELS = [
    "mono_exponential",
    "stretched_exponential",
    "biexponential",
    "inversion_recovery",
    "saturation_recovery",
    "gamma_gaussian_decay",
]
```

Replace `_get_fit_function` with:

```python
def _get_fit_function(
    model: str,
) -> Tuple[Callable[..., np.ndarray], List[str], Dict[str, Tuple[float, float]]]:
    """
    Build the fit function, parameter names, and default bounds for a model.

    Parameters
    ----------
    model : str
        Relaxation model name.

    Returns
    -------
    fit_func : callable
        Function compatible with scipy.optimize.curve_fit, signature
        fit_func(t, *params).
    param_names : list of str
        Ordered parameter names matching fit_func's signature.
    param_bounds : dict
        Default bounds for each parameter as {name: (lower, upper)}.
    """

    if model == "mono_exponential":

        def fit_func(t, amplitude, T, offset):
            return mono_exponential(t, amplitude, T, offset)

        param_names = ["amplitude", "T", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "stretched_exponential":

        def fit_func(t, amplitude, T, beta, offset):
            return stretched_exponential(t, amplitude, T, beta, offset)

        param_names = ["amplitude", "T", "beta", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T": (1e-9, np.inf),
            "beta": (0.05, 5.0),
            "offset": (-np.inf, np.inf),
        }

    elif model == "gamma_gaussian_decay":

        def fit_func(t, amplitude, Gamma0, GammaG, offset):
            return gamma_gaussian_decay(t, amplitude, Gamma0, GammaG, offset)

        param_names = ["amplitude", "Gamma0", "GammaG", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "Gamma0": (0.0, np.inf),
            "GammaG": (0.0, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "biexponential":

        def fit_func(t, amplitude1, tau1, amplitude2, tau2, offset):
            return biexponential(t, amplitude1, tau1, amplitude2, tau2, offset)

        param_names = ["amplitude1", "tau1", "amplitude2", "tau2", "offset"]
        param_bounds = {
            "amplitude1": (-np.inf, np.inf),
            "tau1": (1e-9, np.inf),
            "amplitude2": (-np.inf, np.inf),
            "tau2": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "inversion_recovery":

        def fit_func(t, amplitude, T1, offset):
            return inversion_recovery(t, amplitude, T1, offset)

        param_names = ["amplitude", "T1", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T1": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    elif model == "saturation_recovery":

        def fit_func(t, amplitude, T1, offset):
            return saturation_recovery(t, amplitude, T1, offset)

        param_names = ["amplitude", "T1", "offset"]
        param_bounds = {
            "amplitude": (-np.inf, np.inf),
            "T1": (1e-9, np.inf),
            "offset": (-np.inf, np.inf),
        }

    else:
        raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")

    return fit_func, param_names, param_bounds
```

Replace `_estimate_initial_params` with:

```python
def _estimate_initial_params(
    t: np.ndarray, y: np.ndarray, model: str
) -> Dict[str, float]:
    """
    Estimate initial fit parameters from the data.

    Parameters
    ----------
    t : np.ndarray
        Time axis.
    y : np.ndarray
        Relaxation signal.
    model : str
        Relaxation model name.

    Returns
    -------
    dict
        Initial parameter estimates keyed by parameter name.
    """

    offset_guess = float(y[-1])
    span = float(y[0] - y[-1])
    amplitude_guess = span if span != 0 else float(y.max() - y.min()) or 1.0

    if model in ("mono_exponential", "stretched_exponential", "gamma_gaussian_decay"):
        rate_guess = _estimate_decay_rate(t, y, offset_guess, amplitude_guess)
        T_guess = 1.0 / rate_guess if rate_guess > 0 else (t[-1] - t[0]) / 2

        if model == "mono_exponential":
            return {"amplitude": amplitude_guess, "T": T_guess, "offset": offset_guess}

        if model == "stretched_exponential":
            return {
                "amplitude": amplitude_guess,
                "T": T_guess,
                "beta": 1.0,
                "offset": offset_guess,
            }

        return {
            "amplitude": amplitude_guess,
            "Gamma0": rate_guess,
            "GammaG": rate_guess,
            "offset": offset_guess,
        }

    if model == "biexponential":
        time_span = float(t[-1] - t[0]) or 1.0
        return {
            "amplitude1": amplitude_guess / 2,
            "tau1": max(time_span * 0.1, 1e-6),
            "amplitude2": amplitude_guess / 2,
            "tau2": max(time_span, 1e-6),
            "offset": offset_guess,
        }

    if model in ("inversion_recovery", "saturation_recovery"):
        recovery_amplitude = float(y[-1] - y[0]) or 1.0
        midpoint = y[0] + 0.5 * (y[-1] - y[0])
        crossing_idx = int(np.argmin(np.abs(y - midpoint)))
        T1_guess = max(float(t[crossing_idx] - t[0]), float(t[1] - t[0]))
        return {
            "amplitude": abs(recovery_amplitude),
            "T1": T1_guess,
            "offset": float(y[0]),
        }

    raise ValueError(f"Unsupported model: {model}. Choose from: {SUPPORTED_MODELS}")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS (all tests so far)

- [ ] **Step 5: Commit**

```bash
git add epyr/relaxation/fitting.py tests/test_relaxation.py
git commit -m "Add biexponential, inversion_recovery, saturation_recovery fit support"
```

---

### Task 5: `fit_multiple_decays` model comparison

**Files:**
- Modify: `epyr/relaxation/fitting.py`
- Modify: `epyr/relaxation/__init__.py`
- Modify: `tests/test_relaxation.py`

- [ ] **Step 1: Append the failing tests**

Append to `tests/test_relaxation.py`:

```python
@pytest.mark.standard
class TestFitMultipleDecays:
    def test_picks_mono_exponential_when_data_is_mono_exponential(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0, 100, 300)
        y_true = mono_exponential(t, amplitude=5.0, T=20.0, offset=1.0)
        y = y_true + 0.01 * (y_true.max() - y_true.min()) * rng.standard_normal(t.size)

        results = fit_multiple_decays(t, y, plot=False)
        best = min(results, key=lambda k: results[k].chi_squared)
        assert best == "mono_exponential"

    def test_picks_stretched_exponential_when_data_is_stretched(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0.1, 100, 300)
        y_true = stretched_exponential(t, amplitude=4.0, T=25.0, beta=0.6, offset=0.5)
        y = y_true + 0.01 * (y_true.max() - y_true.min()) * rng.standard_normal(t.size)

        results = fit_multiple_decays(t, y, plot=False)
        best = min(results, key=lambda k: results[k].chi_squared)
        assert best == "stretched_exponential"

    def test_picks_biexponential_when_data_is_biexponential(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0, 200, 300)
        y_true = biexponential(
            t, amplitude1=3.0, tau1=8.0, amplitude2=2.0, tau2=80.0, offset=0.3
        )
        y = y_true + 0.01 * (y_true.max() - y_true.min()) * rng.standard_normal(t.size)

        results = fit_multiple_decays(t, y, plot=False)
        best = min(results, key=lambda k: results[k].chi_squared)
        assert best == "biexponential"

    def test_plot_true_does_not_raise(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=20.0, offset=1.0)
        results = fit_multiple_decays(t, y, plot=True)
        assert all(r.success for r in results.values())
        plt.close("all")
```

Update the `from epyr.relaxation import (...)` block at the top of the file
to add `fit_multiple_decays`:

```python
from epyr.relaxation import (
    RelaxationFitResult,
    biexponential,
    fit_multiple_decays,
    fit_relaxation,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_relaxation.py -v -k TestFitMultipleDecays`
Expected: FAIL with `ImportError: cannot import name 'fit_multiple_decays'`

- [ ] **Step 3: Implement `fit_multiple_decays` and its comparison plot**

Append to `epyr/relaxation/fitting.py`:

```python
def fit_multiple_decays(
    t_data: np.ndarray,
    y_data: np.ndarray,
    models: Optional[List[str]] = None,
    mask: Optional[np.ndarray] = None,
    plot: bool = True,
) -> Dict[str, RelaxationFitResult]:
    """
    Fit relaxation data with multiple models and compare by reduced chi-squared.

    Parameters
    ----------
    t_data : np.ndarray
        Time axis, in the unit chosen by the caller.
    y_data : np.ndarray
        Relaxation signal.
    models : list of str, optional
        Models to fit. Default: ['mono_exponential', 'stretched_exponential',
        'biexponential']. The recovery models and 'gamma_gaussian_decay' are
        excluded by default since they assume a specific data shape rather
        than being interchangeable candidates for a generic decay.
    mask : np.ndarray of bool, optional
        Boolean array selecting points to include (True = include). Passed
        unchanged to each fit_relaxation call.
    plot : bool, optional
        Display a side-by-side comparison plot (default: True).

    Returns
    -------
    dict
        Mapping of model name to RelaxationFitResult for all attempted fits.

    Notes
    -----
    Models are ranked by reduced chi-squared, not R-squared: R-squared is
    biased toward models with more free parameters, since extra parameters
    can only reduce the residual sum of squares on the same data.
    """

    if models is None:
        models = ["mono_exponential", "stretched_exponential", "biexponential"]

    results = {}

    for model in models:
        try:
            result = fit_relaxation(t_data, y_data, model=model, mask=mask, plot=False)
            results[model] = result
        except Exception as e:
            logger.warning(f"Failed to fit {model}: {e}")
            results[model] = RelaxationFitResult(
                model=model,
                parameters={},
                parameter_errors={},
                fitted_curve=np.array([]),
                residuals=np.array([]),
                r_squared=0.0,
                chi_squared=np.inf,
                success=False,
                message=str(e),
            )

    successful_fits = {k: v for k, v in results.items() if v.success}

    if successful_fits and plot:
        _plot_comparison(t_data, y_data, successful_fits)

    logger.info("=== Relaxation Model Comparison ===")
    for model, result in results.items():
        if result.success:
            logger.info(
                f"{model:22s}: chi2_red = {result.chi_squared:.4g},"
                f" R2 = {result.r_squared:.6f}"
            )
        else:
            logger.info(f"{model:22s}: FAILED - {result.message}")

    if successful_fits:
        best_model = min(
            successful_fits.keys(), key=lambda k: successful_fits[k].chi_squared
        )
        logger.info(f"\nBest fit (lowest chi2_red): {best_model}")

    return results


def _plot_comparison(
    t: np.ndarray, y: np.ndarray, results: Dict[str, RelaxationFitResult]
) -> None:
    """
    Plot fitted curves and residuals for multiple relaxation models.

    Parameters
    ----------
    t : np.ndarray
        Full time data passed to fit_multiple_decays.
    y : np.ndarray
        Full relaxation signal passed to fit_multiple_decays.
    results : dict
        Mapping of model name to RelaxationFitResult; only successful fits
        are drawn.
    """

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 10), gridspec_kw={"height_ratios": [3, 1]}
    )

    colors = ["#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

    ax1.plot(t, y, "o", markersize=4, alpha=0.7, label="Data", color="#1f77b4")

    for i, (model, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            t_plot = result.t_fit if result.t_fit is not None else t
            ax1.plot(
                t_plot,
                result.fitted_curve,
                "-",
                linewidth=2,
                label=f"{model} (chi2_red={result.chi_squared:.3g})",
                color=color,
            )

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Intensity")
    ax1.set_title("Relaxation Fitting - Model Comparison")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    for i, (model, result) in enumerate(results.items()):
        if result.success:
            color = colors[i % len(colors)]
            t_plot = result.t_fit if result.t_fit is not None else t
            ax2.plot(
                t_plot,
                result.residuals,
                "o-",
                markersize=2,
                alpha=0.7,
                label=model,
                color=color,
            )

    ax2.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Residuals")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
```

- [ ] **Step 4: Export `fit_multiple_decays`**

In `epyr/relaxation/__init__.py`, replace the full contents with:

```python
"""
EPyR Tools - Relaxation Fitting Module

Time-domain decay and recovery models for T1/T2 EPR relaxation measurements,
with a curve-fitting layer mirroring epyr.lineshapes.fitting.
"""

from .fitting import RelaxationFitResult, fit_multiple_decays, fit_relaxation
from .models import (
    biexponential,
    gamma_gaussian_decay,
    inversion_recovery,
    mono_exponential,
    saturation_recovery,
    stretched_exponential,
)

__all__ = [
    "mono_exponential",
    "stretched_exponential",
    "biexponential",
    "inversion_recovery",
    "saturation_recovery",
    "gamma_gaussian_decay",
    "fit_relaxation",
    "fit_multiple_decays",
    "RelaxationFitResult",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS (all tests so far)

- [ ] **Step 6: Commit**

```bash
git add epyr/relaxation/fitting.py epyr/relaxation/__init__.py tests/test_relaxation.py
git commit -m "Add fit_multiple_decays model comparison ranked by reduced chi-squared"
```

---

### Task 6: Deep edge-case tests and scientific cross-check

**Files:**
- Modify: `tests/test_relaxation.py`

- [ ] **Step 1: Write the deep and scientific tests**

Append to `tests/test_relaxation.py`:

```python
@pytest.mark.deep
class TestFitRelaxationEdgeCases:
    def test_decay_time_much_shorter_than_time_range(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=0.5, offset=1.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.parameters["T"] == pytest.approx(0.5, rel=1e-2)

    def test_decay_time_much_longer_than_time_range(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=10000.0, offset=1.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        # The decay is nearly linear over this range, so individual
        # parameters are not well constrained; only require convergence
        # and a good fit, not exact parameter recovery.
        assert result.success
        assert result.r_squared > 0.999

    def test_beta_near_lower_bound(self):
        t = np.linspace(0.1, 100, 200)
        y = stretched_exponential(t, amplitude=4.0, T=25.0, beta=0.1, offset=0.5)
        result = fit_relaxation(t, y, model="stretched_exponential", plot=False)
        assert result.success
        assert result.parameters["beta"] == pytest.approx(0.1, rel=1e-2)

    def test_beta_above_one_compressed_exponential(self):
        t = np.linspace(0.1, 100, 200)
        y = stretched_exponential(t, amplitude=4.0, T=25.0, beta=2.0, offset=0.5)
        result = fit_relaxation(t, y, model="stretched_exponential", plot=False)
        assert result.success
        assert result.parameters["beta"] == pytest.approx(2.0, rel=1e-2)

    def test_negative_offset(self):
        t = np.linspace(0, 100, 200)
        y = mono_exponential(t, amplitude=5.0, T=15.0, offset=-3.0)
        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success
        assert result.parameters["offset"] == pytest.approx(-3.0, rel=1e-2)


@pytest.mark.scientific
class TestFitRelaxationAnalyticCrossCheck:
    def test_mono_exponential_matches_log_linear_regression(self):
        rng = np.random.default_rng(2)
        t = np.linspace(0, 100, 200)
        T_true = 18.0
        amplitude_true = 6.0
        offset_true = 0.8
        y_true = mono_exponential(t, amplitude=amplitude_true, T=T_true, offset=offset_true)
        y = y_true + 0.01 * amplitude_true * rng.standard_normal(t.size)

        result = fit_relaxation(t, y, model="mono_exponential", plot=False)
        assert result.success

        # Closed-form solution: with the true offset known, log(|y - offset|)
        # is linear in t with slope -1 / T.
        slope, _ = np.polyfit(t, np.log(np.abs(y - offset_true)), 1)
        T_analytic = -1.0 / slope

        assert result.parameters["T"] == pytest.approx(T_analytic, rel=0.08)
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `pytest tests/test_relaxation.py -v -m "deep or scientific"`
Expected: PASS (10 tests)

- [ ] **Step 3: Run the full file once more**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS (all tests)

- [ ] **Step 4: Commit**

```bash
git add tests/test_relaxation.py
git commit -m "Add deep edge-case and scientific cross-check tests for relaxation fitting"
```

---

### Task 7: Example script

**Files:**
- Create: `examples/clean/06_relaxation_fitting.py`

- [ ] **Step 1: Write the example script**

Create `examples/clean/06_relaxation_fitting.py`:

```python
"""
T1/T2 relaxation fitting: exponential, stretched, bi-exponential, and
Gamma0/GammaG decay models.

Part 1 -- real data:   mono- and stretched-exponential fits on a T2 echo
                       decay, compared with fit_multiple_decays.
Part 2 -- synthetic:   bi-exponential decay recovered with fit_relaxation.
Part 3 -- synthetic:   Gamma0/GammaG (homogeneous + spectral-diffusion)
                       decay recovered with fit_relaxation.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from epyr import eprload
from epyr.relaxation import (
    biexponential,
    fit_multiple_decays,
    fit_relaxation,
    gamma_gaussian_decay,
)

DATA_DIR = Path(__file__).parents[1] / "data"
RNG = np.random.default_rng(42)


# --- Part 1: T2 echo decay on real data ---------------------------------

t1, y1, params1, _ = eprload(
    DATA_DIR / "2020_10_DMTTFBr_T2EH_28dB_6K_20ns_40ns_hperpc.DTA",
    plot_if_possible=False,
)
y1_magnitude = np.abs(y1)

results1 = fit_multiple_decays(
    t1, y1_magnitude, models=["mono_exponential", "stretched_exponential"], plot=True
)
print("DMTTFBr T2 echo decay -- model comparison")
for model, result in results1.items():
    status = (
        f"chi2_red = {result.chi_squared:.4g}, T = {result.parameters['T']:.1f} ns"
        if result.success
        else "failed"
    )
    print(f"  {model:<22s}: {status}")


# --- Part 2: synthetic bi-exponential decay ------------------------------

t2 = np.linspace(0.0, 200.0, 300)  # ns
y2_true = biexponential(
    t2, amplitude1=3.0, tau1=8.0, amplitude2=2.0, tau2=80.0, offset=0.3
)
noise2 = 0.01 * (y2_true.max() - y2_true.min()) * RNG.standard_normal(t2.size)
y2 = y2_true + noise2

result2 = fit_relaxation(t2, y2, model="biexponential", time_unit="ns", plot=True)
print("\nSynthetic bi-exponential decay")
print(result2.summary())


# --- Part 3: synthetic Gamma0/GammaG decay -------------------------------

t3 = np.linspace(0.0, 100.0, 200)  # us
y3_true = gamma_gaussian_decay(t3, amplitude=3.0, Gamma0=0.05, GammaG=0.02, offset=0.2)
noise3 = 0.01 * (y3_true.max() - y3_true.min()) * RNG.standard_normal(t3.size)
y3 = y3_true + noise3

result3 = fit_relaxation(t3, y3, model="gamma_gaussian_decay", time_unit="us", plot=True)
print("\nSynthetic Gamma0/GammaG echo decay")
print(result3.summary())

plt.show()
```

- [ ] **Step 2: Run the script and verify it completes without raising**

Run: `python examples/clean/06_relaxation_fitting.py`
Expected: prints the model comparison and both fit summaries, then opens
plot windows (close them to let the script exit). No traceback.

- [ ] **Step 3: Commit**

```bash
git add examples/clean/06_relaxation_fitting.py
git commit -m "Add relaxation fitting example script"
```

---

### Task 8: Documentation and version bump

**Files:**
- Modify: `CLAUDE.md:7,316`
- Create: `docs/api/epyr.relaxation.rst`
- Modify: `docs/api/epyr.rst`
- Modify: `docs/release_notes.rst`
- Create: `docs/release_notes/v0.4.0.rst`
- Modify: `pyproject.toml:7`
- Modify: `epyr/__init__.py:126`

- [ ] **Step 1: Update `CLAUDE.md`**

Line 7 currently reads:

```
EPyR Tools is a comprehensive Python package for Electron Paramagnetic Resonance (EPR) spectroscopy data analysis. The package provides tools for loading Bruker EPR files (BES3T .dsc/.dta and ESP .par/.spc formats), baseline correction, lineshape analysis, FFT-based signal processing, FAIR data conversion, and visualization. Current version: 0.3.8.
```

Replace with:

```
EPyR Tools is a comprehensive Python package for Electron Paramagnetic Resonance (EPR) spectroscopy data analysis. The package provides tools for loading Bruker EPR files (BES3T .dsc/.dta and ESP .par/.spc formats), baseline correction, lineshape analysis, T1/T2 relaxation fitting, FFT-based signal processing, FAIR data conversion, and visualization. Current version: 0.4.0.
```

Immediately after line 319 (`- `units.py` - General `unitconvert()` interface`)
and before line 320 (`- `eprplot.py` - Specialized EPR plotting functions`),
insert:

```
- `relaxation/` - T1/T2 relaxation fitting
  - `models.py` - Decay/recovery functions: mono- and stretched-exponential,
    bi-exponential, inversion/saturation recovery, Gamma0/GammaG decay
  - `fitting.py` - `fit_relaxation()` and `fit_multiple_decays()`, mirroring
    `lineshapes/fitting.py`
```

- [ ] **Step 2: Create the Sphinx API page**

Create `docs/api/epyr.relaxation.rst`:

```rst
epyr.relaxation module
=======================

.. automodule:: epyr.relaxation
   :members:
   :undoc-members:
   :show-inheritance:

Submodules
----------

.. toctree::
   :maxdepth: 2

   generated/epyr.relaxation.models
   generated/epyr.relaxation.fitting
```

- [ ] **Step 3: Link the new page from `docs/api/epyr.rst`**

In `docs/api/epyr.rst`, the first toctree currently reads:

```rst
.. toctree::
   :maxdepth: 2

   epyr.eprload
   epyr.baseline
   epyr.fair
   epyr.lineshapes
   epyr.plot
   epyr.signalprocessing
   epyr.physics
   epyr.isotope_gui
```

Add `epyr.relaxation` after `epyr.lineshapes`:

```rst
.. toctree::
   :maxdepth: 2

   epyr.eprload
   epyr.baseline
   epyr.fair
   epyr.lineshapes
   epyr.relaxation
   epyr.plot
   epyr.signalprocessing
   epyr.physics
   epyr.isotope_gui
```

The second list (`Complete Module Index`) currently reads:

```rst
.. autosummary::
   :toctree: generated/
   :recursive:

   epyr.eprload
   epyr.baseline
   epyr.fair
   epyr.lineshapes
   epyr.eprplot
   epyr.signalprocessing
   epyr.physics
   epyr.isotope_gui
   epyr.cli
   epyr.config
   epyr.performance
   epyr.plugins
```

Add `epyr.relaxation` after `epyr.lineshapes`:

```rst
.. autosummary::
   :toctree: generated/
   :recursive:

   epyr.eprload
   epyr.baseline
   epyr.fair
   epyr.lineshapes
   epyr.relaxation
   epyr.eprplot
   epyr.signalprocessing
   epyr.physics
   epyr.isotope_gui
   epyr.cli
   epyr.config
   epyr.performance
   epyr.plugins
```

- [ ] **Step 4: Create the release notes page**

Create `docs/release_notes/v0.4.0.rst`:

```rst
Version 0.4.0
=============

**Release Date:** June 2026

**Status:** Stable release

New: relaxation fitting module
-------------------------------

A new package, ``epyr.relaxation``, fits time-domain T1/T2 relaxation data,
complementing the existing field-domain lineshape fitting in
``epyr.lineshapes.fitting``.

- ``epyr.relaxation.models``: six decay/recovery functions, ``mono_exponential``,
  ``stretched_exponential``, ``biexponential``, ``inversion_recovery``,
  ``saturation_recovery``, and ``gamma_gaussian_decay`` (a combined
  homogeneous/spectral-diffusion echo-decay model).
- ``epyr.relaxation.fit_relaxation()``: fits a single decay/recovery curve
  with any of the six models, returning a ``RelaxationFitResult`` with
  parameters, errors, R-squared, reduced chi-squared, and an optional
  annotated plot. Mirrors the existing ``fit_epr_signal()`` API.
- ``epyr.relaxation.fit_multiple_decays()``: fits several candidate models
  and ranks them by reduced chi-squared (not R-squared, which is biased
  toward models with more free parameters) to suggest the best-supported
  decay model for a given dataset.
- New example script ``examples/clean/06_relaxation_fitting.py`` demonstrates
  all three entry points on a real T2 echo decay and on synthetic
  bi-exponential and Gamma0/GammaG data.

No backwards-incompatible changes.
```

- [ ] **Step 5: Link the new release notes from `docs/release_notes.rst`**

In `docs/release_notes.rst`, the toctree currently reads:

```rst
.. toctree::
   :maxdepth: 2
   :caption: Version History:

   release_notes/v0.3.9
   release_notes/v0.3.8
   release_notes/v0.3.6
   release_notes/v0.3.5
   release_notes/v0.2.0
```

Add the new entry at the top:

```rst
.. toctree::
   :maxdepth: 2
   :caption: Version History:

   release_notes/v0.4.0
   release_notes/v0.3.9
   release_notes/v0.3.8
   release_notes/v0.3.6
   release_notes/v0.3.5
   release_notes/v0.2.0
```

Immediately before the `Version 0.3.9 (Latest)` heading, insert a new
section (and remove `(Latest)` from the `0.3.9` heading since it is no
longer the latest):

```rst
Version 0.4.0 (Latest)
-----------------------

**Release Date:** June 2026

New ``epyr.relaxation`` package for T1/T2 relaxation fitting: six
decay/recovery models (mono- and stretched-exponential, bi-exponential,
inversion/saturation recovery, Gamma0/GammaG), ``fit_relaxation()``, and
``fit_multiple_decays()`` (ranked by reduced chi-squared). See
:doc:`release_notes/v0.4.0`.

Version 0.3.9
--------------
```

(The old line `Version 0.3.9 (Latest)` becomes `Version 0.3.9`, with its
underline shortened to match, and the `**Release Date:** May 2026` /
descriptive paragraph below it stays as-is.)

- [ ] **Step 6: Bump the version**

In `pyproject.toml`, line 7:

```toml
version = "0.3.9"
```

becomes:

```toml
version = "0.4.0"
```

In `epyr/__init__.py`, line 126:

```python
__version__ = "0.3.9"
```

becomes:

```python
__version__ = "0.4.0"
```

- [ ] **Step 7: Verify version consistency**

Run: `make check-version`
Expected: `Version 0.4.0 is consistent`

- [ ] **Step 8: Commit**

```bash
git add CLAUDE.md docs/api/epyr.relaxation.rst docs/api/epyr.rst docs/release_notes.rst docs/release_notes/v0.4.0.rst pyproject.toml epyr/__init__.py
git commit -m "Document epyr.relaxation and bump version to 0.4.0"
```

---

### Task 9: Final quality gate

**Files:** none (verification only)

- [ ] **Step 1: Run the full new test file with all markers**

Run: `pytest tests/test_relaxation.py -v`
Expected: PASS, all tests (model functions, fitting per model, noisy fit,
validation, multiple-decay comparison, deep edge cases, scientific
cross-check).

- [ ] **Step 2: Run the project smoke suite to check for regressions**

Run: `pytest -m smoke`
Expected: PASS (no failures introduced in other modules).

- [ ] **Step 3: Run formatting and linting**

Run: `make format`
Expected: reformats files if needed (commit again if `epyr/relaxation/*.py`
changed); exits 0.

Run: `make quality`
Expected: exits 0 (flake8, mypy, isort, black all clean on `epyr/relaxation/`).

If `mypy` or `flake8` flag anything in the new files, fix it directly
(e.g. missing type hints, unused imports) rather than suppressing the check.

- [ ] **Step 4: Run the example script end-to-end once more**

Run: `python examples/clean/06_relaxation_fitting.py`
Expected: completes without a traceback (see Task 7, Step 2).

- [ ] **Step 5: Verify version consistency one more time**

Run: `make check-version`
Expected: `Version 0.4.0 is consistent`

- [ ] **Step 6: If any fixes were needed in this task, commit them**

```bash
git add -A
git commit -m "Quality fixes for relaxation fitting module"
```

(Skip this step if Steps 1-5 required no changes.)
