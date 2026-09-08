"""
Pseudomodulation of EPR spectra (epyr.lineshapes.pseudo_modulation).

Part 1 — small-amplitude limit: pseudomodulation converges to the analytic
           first derivative as mod_amplitude -> 0.
Part 2 — over-modulation: increasing mod_amplitude broadens and distorts the
           lineshape, the classic instrumental artifact CW EPR users watch for.
Part 3 — second-harmonic detection (harmonic=2).

All data here is synthetic: a two-line absorption spectrum (hyperfine doublet)
recorded without field modulation, as it would come from a pulse field-swept
or rapid-scan experiment.
"""

import matplotlib.pyplot as plt
import numpy as np

import epyr
from epyr.lineshapes import pseudo_modulation

x = np.linspace(3400.0, 3600.0, 4000)  # G
linewidth = 8.0  # G, FWHM
absorption = epyr.lorentzian(x, 3480.0, linewidth) + epyr.lorentzian(
    x, 3520.0, linewidth
)


# --- Part 1: small-amplitude limit -------------------------------------------
# S_1(B) -> (mod_amplitude/4) * dA/dB as mod_amplitude -> 0 (see docstring Notes)

mod_amplitude_small = 0.5  # G, much smaller than linewidth
y_pm_small = pseudo_modulation(x, absorption, mod_amplitude_small, harmonic=1)

analytic_derivative = epyr.lorentzian(
    x, 3480.0, linewidth, derivative=1
) + epyr.lorentzian(x, 3520.0, linewidth, derivative=1)
y_analytic_small = (mod_amplitude_small / 4) * analytic_derivative

relative_error = np.max(np.abs(y_pm_small - y_analytic_small)) / np.max(
    np.abs(y_analytic_small)
)
print(f"Small-amplitude limit: max relative error = {relative_error:.2e}")

fig1, axes1 = plt.subplots(1, 2, figsize=(12, 4))
axes1[0].plot(x, absorption, lw=1.2)
axes1[0].set_title("Absorption spectrum (no modulation)")
axes1[0].set_xlabel("Magnetic Field (G)")
axes1[0].set_ylabel("Intensity (a.u.)")
axes1[0].grid(True, linestyle=":", alpha=0.5)

axes1[1].plot(x, y_analytic_small, lw=2.5, label="analytic derivative (scaled)")
axes1[1].plot(x, y_pm_small, lw=1.2, ls="--", label="pseudo_modulation, harmonic=1")
axes1[1].set_title(f"mod_amplitude = {mod_amplitude_small} G (<< linewidth)")
axes1[1].set_xlabel("Magnetic Field (G)")
axes1[1].legend()
axes1[1].grid(True, linestyle=":", alpha=0.5)
fig1.suptitle("Part 1 — Small-amplitude limit")
fig1.tight_layout()


# --- Part 2: over-modulation broadening --------------------------------------
# As mod_amplitude approaches and exceeds the linewidth, the pseudomodulated
# signal broadens and its peak-to-peak amplitude drops, mirroring what happens
# on a real lock-in amplifier when the modulation amplitude is set too high.

mod_amplitudes = [1.0, 15.0, 40.0]  # G; 40 G matches the doublet's line spacing

fig2, ax2 = plt.subplots(figsize=(7, 5))
for mod_amplitude in mod_amplitudes:
    y_pm = pseudo_modulation(x, absorption, mod_amplitude, harmonic=1)
    ax2.plot(x, y_pm, lw=1.5, label=f"mod_amplitude = {mod_amplitude} G")
ax2.set_title("Part 2 — Over-modulation broadens the derivative lineshape")
ax2.set_xlabel("Magnetic Field (G)")
ax2.set_ylabel("Intensity (a.u.)")
ax2.legend()
ax2.grid(True, linestyle=":", alpha=0.5)
fig2.tight_layout()


# --- Part 3: second-harmonic detection ---------------------------------------

y_h1 = pseudo_modulation(x, absorption, 4.0, harmonic=1)
y_h2 = pseudo_modulation(x, absorption, 4.0, harmonic=2)

fig3, axes3 = plt.subplots(1, 2, figsize=(12, 4))
axes3[0].plot(x, y_h1, lw=1.2)
axes3[0].set_title("First harmonic (dA/dB-like)")
axes3[0].set_xlabel("Magnetic Field (G)")
axes3[0].grid(True, linestyle=":", alpha=0.5)

axes3[1].plot(x, y_h2, lw=1.2, color="C1")
axes3[1].set_title("Second harmonic (d²A/dB²-like)")
axes3[1].set_xlabel("Magnetic Field (G)")
axes3[1].grid(True, linestyle=":", alpha=0.5)
fig3.suptitle("Part 3 — Harmonic detection, mod_amplitude = 4 G")
fig3.tight_layout()

plt.show()
