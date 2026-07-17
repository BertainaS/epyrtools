# Pseudomodulation — Design Spec
Date: 2026-07-17
Status: approved

## Context

EPR spectra recorded without field modulation (slow-passage absorption, or
signals reconstructed from pulsed/rapid-scan experiments) do not have the
first-derivative lineshape that CW-EPR users expect. Pseudomodulation
numerically reproduces the effect of lock-in detection under sinusoidal field
modulation, without re-recording the spectrum, following the algorithm of
Hyde, Pasenkiewicz-Gierula, Jesmanowicz, and Antholine, *Appl. Magn. Reson.*
1, 483-496 (1990).

## Scope

One new public function in a new module `epyr/lineshapes/pseudo_modulation.py`,
exported via `epyr/lineshapes/__init__.py`. No changes to existing lineshape
or fitting functions.

## Architecture

```
epyr/lineshapes/
    pseudo_modulation.py   ← new file (this spec)
    convspec.py             ← unchanged, sibling module (FFT/convolution pattern)
    __init__.py              ← updated to export pseudo_modulation
```

## Algorithm

Field modulation at amplitude `Bm` (peak-to-peak) sweeps the instantaneous
field as `B(t) = B0 + (Bm/2) sin(theta)`. Expanding the absorption spectrum
`A(B)` as an inverse Fourier transform in spatial frequency `k` (conjugate to
the field axis) and applying the Jacobi-Anger identity

```
(1/2*pi) * integral_0^2pi  exp(i*z*sin(theta)) * exp(-i*n*theta) d theta = J_n(z)
```

gives the signal detected at the n-th harmonic of the modulation frequency as
a linear filter of the original spectrum:

```
S_n(B) = IFFT[ FFT(A)(k) * J_n(k * Bm / 2) ]
```

where `J_n` is the Bessel function of the first kind, order `n` (the
harmonic), and `k = 2*pi * fftfreq(N, dx)`. The transfer function is real, so
a real input spectrum produces a real output. In the small-amplitude limit,
`J_1(z) ~ z/2`, so `S_1(B) -> (Bm/4) * dA/dB`: pseudomodulation with a small
`Bm` reduces to the ordinary scaled derivative, which is the basis for the
scientific validation test below.

## Function

```python
def pseudo_modulation(
    x: np.ndarray,
    y: np.ndarray,
    mod_amplitude: float,
    harmonic: int = 1,
    pad: bool = True,
) -> np.ndarray:
```

**Parameters**

| Parameter | Description |
|---|---|
| `x` | Field axis, 1D array, uniform step, unit as recorded (G or mT) |
| `y` | Spectrum, 1D array, same length as `x`, real or complex |
| `mod_amplitude` | Peak-to-peak modulation amplitude, same unit as `x`, must be > 0 |
| `harmonic` | Detection harmonic, `1` or `2` |
| `pad` | If `True`, edge-pad `y` before the FFT convolution to suppress circular-wraparound artifacts at the spectrum boundaries, then crop back to the original length |

**Returns** `np.ndarray`, same shape and dtype family as `y` (real in, real
out; complex in, complex out).

**Raises** `ValueError` if `x` and `y` have different lengths, `x` is not
uniformly spaced (relative tolerance check on `np.diff(x)`), `mod_amplitude
<= 0`, or `harmonic not in (1, 2)`.

**Implementation notes**

- Validation follows the pattern in `convspec._validate_convspec_inputs`.
- Uniform spacing: `np.diff(x)` must be constant within `rtol=1e-6`.
- Padding: `np.pad(y, pad_width, mode='edge')` with `pad_width = len(y) // 2`
  on each side; crop the result back with the matching slice.
- Bessel transfer function via `scipy.special.jv(harmonic, k * mod_amplitude / 2)`.
- Cast the result to real (`np.real`) when the input `y` is real, mirroring
  the real/complex handling already used in `convspec()`.

## Example

```python
from epyr import eprload
from epyr.lineshapes import pseudo_modulation

x, y, params, _ = eprload("slow_passage_absorption.DTA")
y_pm = pseudo_modulation(x, y, mod_amplitude=2.0, harmonic=1)  # Bm = 2 G
```

## Tests

New test file `tests/test_pseudo_modulation.py` with marks `smoke`,
`standard`, `scientific`.

| Test | Mark |
|---|---|
| `test_pseudo_modulation_basic_shape` | smoke |
| `test_pseudo_modulation_harmonic_2` | standard |
| `test_pseudo_modulation_nonuniform_x_raises` | standard |
| `test_pseudo_modulation_invalid_mod_amplitude_raises` | standard |
| `test_pseudo_modulation_invalid_harmonic_raises` | standard |
| `test_pseudo_modulation_small_amplitude_matches_derivative` | scientific |

The scientific test compares `pseudo_modulation(x, gaussian_absorption, mod_amplitude=eps, harmonic=1)`
against `(eps / 4) * d/dx[gaussian_absorption]` (finite-difference or
analytic derivative) for a small `eps`, within a numerical tolerance.

## Out of scope

- Harmonics beyond 2.
- 2D EPR spectra (row-by-row support can be added later if needed).
- GUI/interactive modulation-amplitude selection.
- Automatic modulation-amplitude estimation from experimental data.
