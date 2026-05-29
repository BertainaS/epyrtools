"""
Baseline correction and lineshape fitting for CW EPR spectra.

Part 1 — real data:    polynomial baseline on CaWO4:Er (BES3T).
Part 2 — synthetic:    Lorentzian peak on a sloping background.
Part 3 — real data:    first-derivative Lorentzian fit on CuSO4 (ESP).
Part 4 — synthetic:    model comparison on a two-line spectrum.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import epyr
from epyr.baseline import baseline_polynomial_1d
from epyr.lineshapes.fitting import fit_epr_signal, fit_multiple_shapes

DATA_DIR = Path(__file__).parents[1] / "data"
RNG = np.random.default_rng(42)


# --- Part 1: polynomial baseline on real CaWO4:Er spectrum -------------------

x, y, params, _ = epyr.eprload(
    DATA_DIR / "130406SB_CaWO4_Er_CW_5K_20.DSC",
    plot_if_possible=False,
)

# exclude_center=True (default) avoids fitting the signal region as baseline
y_corr, y_bl = baseline_polynomial_1d(x, y, params, order=2)

fig1, axes1 = plt.subplots(1, 2, figsize=(12, 4))
axes1[0].plot(x, y, lw=0.8, label="raw")
axes1[0].plot(x, y_bl, lw=1.2, ls="--", label="baseline (poly, order 2)")
axes1[0].set_xlabel("Magnetic Field (G)")
axes1[0].set_ylabel("Intensity (a.u.)")
axes1[0].set_title("CaWO4:Er — raw + baseline")
axes1[0].legend()
axes1[0].grid(True, linestyle=":", alpha=0.5)

axes1[1].plot(x, y_corr, lw=0.8, color="C1")
axes1[1].set_xlabel("Magnetic Field (G)")
axes1[1].set_ylabel("Intensity (a.u.)")
axes1[1].set_title("CaWO4:Er — after correction")
axes1[1].grid(True, linestyle=":", alpha=0.5)
fig1.suptitle("Part 1 — Polynomial baseline (real data)")
fig1.tight_layout()


# --- Part 2: baseline on synthetic absorption peak with sloping background ---

x_syn = np.linspace(3300.0, 3700.0, 512)   # G
peak = epyr.lorentzian(x_syn, 3500.0, 15.0)
# Linear background typical of a poorly balanced bridge
background = 0.15 - 0.6 * (x_syn - 3300.0) / 400.0
noise = 0.012 * RNG.standard_normal(512)
y_syn = peak + background + noise

y_corr_syn, y_bl_syn = baseline_polynomial_1d(
    x_syn, y_syn, {}, order=1, exclude_center=True
)

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4))
axes2[0].plot(x_syn, y_syn, lw=0.8, label="raw (signal + background)")
axes2[0].plot(x_syn, y_bl_syn, lw=1.2, ls="--", label="baseline fit (order 1)")
axes2[0].set_xlabel("Magnetic Field (G)")
axes2[0].set_ylabel("Intensity (a.u.)")
axes2[0].set_title("Synthetic — raw and baseline")
axes2[0].legend()
axes2[0].grid(True, linestyle=":", alpha=0.5)

axes2[1].plot(x_syn, peak, lw=0.8, ls="--", label="true peak")
axes2[1].plot(x_syn, y_corr_syn, lw=0.8, label="corrected")
axes2[1].set_xlabel("Magnetic Field (G)")
axes2[1].set_ylabel("Intensity (a.u.)")
axes2[1].set_title("Synthetic — corrected vs true signal")
axes2[1].legend()
axes2[1].grid(True, linestyle=":", alpha=0.5)
fig2.suptitle("Part 2 — Polynomial baseline (synthetic data)")
fig2.tight_layout()


# --- Part 3: first-derivative fit on real CuSO4 data ------------------------

x3, y3, params3, _ = epyr.eprload(
    DATA_DIR / "CuSO4_001.par",
    plot_if_possible=False,
)

# CW EPR with field modulation gives a first-derivative lineshape
result = fit_epr_signal(x3, y3, "lorentzian", derivative=1, plot=True)
print("\nCuSO4 — Lorentzian first-derivative fit")
print(f"  R²     = {result.r_squared:.4f}")
print(f"  center = {result.parameters['center']:.1f} G")
print(f"  width  = {result.parameters['width']:.1f} G  (HWHM)")


# --- Part 4: model comparison on a synthetic pseudo-Voigt line ---------------
# A pseudo-Voigt profile (mix of Lorentzian and Gaussian) is realistic for
# inhomogeneously broadened EPR lines. Fitting it with pure Gaussian or
# Lorentzian gives systematically worse R², which model comparison reveals.

x4 = np.linspace(3300.0, 3700.0, 512)  # G
# 60 % Lorentzian / 40 % Gaussian mix, first derivative.
# Lineshape functions are area-normalized, so the amplitude of the derivative
# profile depends on width (~1/width²). Add noise at 2 % of the signal
# peak-to-peak to achieve a realistic SNR without hard-coding a unit.
y4_signal = epyr.pseudo_voigt(x4, 3500.0, 22.0, eta=0.6, derivative=1)
noise_level = 0.02 * (y4_signal.max() - y4_signal.min())
y4 = y4_signal + noise_level * RNG.standard_normal(512)

results = fit_multiple_shapes(x4, y4, derivative=1, plot=True)
print("\nSynthetic pseudo-Voigt line (η = 0.6) — model comparison:")
for shape, r in results.items():
    status = f"R² = {r.r_squared:.4f}" if r.success else "failed"
    print(f"  {shape:<14s}: {status}")

plt.show()
