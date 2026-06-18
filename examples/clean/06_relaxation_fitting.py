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

result3 = fit_relaxation(
    t3, y3, model="gamma_gaussian_decay", time_unit="us", plot=True
)
print("\nSynthetic Gamma0/GammaG echo decay")
print(result3.summary())

plt.show()
