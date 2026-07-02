"""
FFT preprocessing pipeline: remove_baseline, apodize, zero_pad.

Demonstrates the three standalone preprocessing functions for time-domain
EPR data (ESEEM, HYSCORE, Rabi) introduced in EPyR Tools 0.5.

Datasets used
-------------
1D (Rabi / echo decay):
    2024_08_CaWO4171Yb_rabi_6K_6724G_18dB — ¹⁷¹Yb in CaWO₄, X-band,
    6 K, 18 dB attenuation.  Real component of echo height vs. pulse
    length (ns).  The signal oscillates on top of a slow exponential
    decay envelope — a canonical case for preprocessing.

2D (HYSCORE):
    Hyscore — standard HYSCORE dataset.  Full 2D FFT after 2D baseline
    removal and 2D apodization reveals nuclear modulation cross-peaks.

Functions demonstrated
----------------------
remove_baseline  : polynomial and exponential baseline subtraction
apodize          : Hann (full and right-half), Hamming, Blackman
zero_pad         : multiplicative factor and absolute point count
analyze_frequencies : FFT after manual 1D preprocessing
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import epyr
from epyr.signalprocessing import (
    apodize,
    remove_baseline,
    zero_pad,
)

DATA_DIR = Path(__file__).parents[1] / "data"

# ── Load datasets ─────────────────────────────────────────────────────────────

x_rabi, y_rabi, p_rabi, _ = epyr.eprload(
    DATA_DIR / "2024_08_CaWO4171Yb_rabi_6K_6724G_18dB.DSC",
    plot_if_possible=False,
)
x_hyscore, y_hyscore, p_hyscore, _ = epyr.eprload(
    DATA_DIR / "Hyscore.DSC",
    plot_if_possible=False,
)

# 1D Rabi: average over all traces if 2D, take real part
if np.ndim(y_rabi) == 2:
    sig1d = np.real(y_rabi).mean(axis=0)
    t_ns = x_rabi[0]
else:
    sig1d = np.real(y_rabi)
    t_ns = x_rabi

# 2D HYSCORE: real part
sig2d = np.real(y_hyscore)
t1_ns = x_hyscore[0]
t2_ns = x_hyscore[1]

print("1D Rabi:")
print(f"  {len(sig1d)} points   t = {t_ns[0]:.0f}–{t_ns[-1]:.0f} ns"
      f"  dt = {t_ns[1]-t_ns[0]:.1f} ns")
print("2D HYSCORE:")
print(f"  shape {sig2d.shape}"
      f"  t1 = {t1_ns[0]:.0f}–{t1_ns[-1]:.0f} ns"
      f"  t2 = {t2_ns[0]:.0f}–{t2_ns[-1]:.0f} ns")


# ═══════════════════════════════════════════════════════════════════════════════
# Part 1 — Baseline removal: polynomial vs. exponential
# ═══════════════════════════════════════════════════════════════════════════════
# A Rabi signal decays as exp(-t/T_2).  Subtracting this envelope leaves a
# zero-mean oscillation that FFTs cleanly.  Polynomial order-2 is a simple
# alternative when the decay is approximately linear over the acquisition window.

print("\nPart 1 — Baseline removal ...")

sig_bl_poly, baseline_poly = remove_baseline(
    t_ns, sig1d, method="polynomial", order=2, end_fraction=0.15
)
sig_bl_exp, baseline_exp = remove_baseline(
    t_ns, sig1d, method="exponential", end_fraction=0.15
)

fig1, axes = plt.subplots(2, 2, figsize=(13, 6))

axes[0, 0].plot(t_ns, sig1d, lw=1, label="Raw signal")
axes[0, 0].plot(t_ns, baseline_poly, lw=2, ls="--", label="Polynomial baseline")
axes[0, 0].set_xlabel("Pulse length (ns)")
axes[0, 0].set_ylabel("Echo height (a.u.)")
axes[0, 0].set_title("Polynomial baseline (order 2)")
axes[0, 0].legend()
axes[0, 0].grid(True, ls=":", alpha=0.5)

axes[0, 1].plot(t_ns, sig1d, lw=1, label="Raw signal")
axes[0, 1].plot(t_ns, baseline_exp, lw=2, ls="--", color="C3",
                label="Exponential baseline")
axes[0, 1].set_xlabel("Pulse length (ns)")
axes[0, 1].set_ylabel("Echo height (a.u.)")
axes[0, 1].set_title("Exponential baseline")
axes[0, 1].legend()
axes[0, 1].grid(True, ls=":", alpha=0.5)

axes[1, 0].plot(t_ns, sig_bl_poly, lw=1, color="C1")
axes[1, 0].axhline(0, color="k", lw=0.5)
axes[1, 0].set_xlabel("Pulse length (ns)")
axes[1, 0].set_ylabel("Residual (a.u.)")
axes[1, 0].set_title("After polynomial removal")
axes[1, 0].grid(True, ls=":", alpha=0.5)

axes[1, 1].plot(t_ns, sig_bl_exp, lw=1, color="C3")
axes[1, 1].axhline(0, color="k", lw=0.5)
axes[1, 1].set_xlabel("Pulse length (ns)")
axes[1, 1].set_ylabel("Residual (a.u.)")
axes[1, 1].set_title("After exponential removal")
axes[1, 1].grid(True, ls=":", alpha=0.5)

for ax in axes.flat:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig1.suptitle("Part 1 — Baseline removal: polynomial vs. exponential")
fig1.tight_layout()


# ═══════════════════════════════════════════════════════════════════════════════
# Part 2 — Apodization: window comparison
# ═══════════════════════════════════════════════════════════════════════════════
# A Rabi signal starts large and decays: the right-half window (0→1 mapping)
# is applied from the start, so it leaves the beginning intact and tapers
# toward zero — matching the natural envelope direction.
# Full symmetric windows are appropriate when the signal has no preferred
# direction (e.g., a FID reconstructed from the echo).

print("Part 2 — Apodization ...")

windows = {
    "Hann (right half)": apodize(sig_bl_exp, window="hann", half_window="right"),
    "Hann (full)": apodize(sig_bl_exp, window="hann"),
    "Hamming": apodize(sig_bl_exp, window="hamming"),
    "Blackman": apodize(sig_bl_exp, window="blackman"),
}

fig2, axes = plt.subplots(2, 2, figsize=(13, 6))
colors = ["C0", "C1", "C2", "C3"]

for ax, (name, windowed), color in zip(axes.flat, windows.items(), colors):
    ax.plot(t_ns, sig_bl_exp, lw=0.8, color="grey", alpha=0.5, label="Baseline removed")
    ax.plot(t_ns, windowed, lw=1.2, color=color, label=name)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("Pulse length (ns)")
    ax.set_ylabel("Signal (a.u.)")
    ax.set_title(name)
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig2.suptitle("Part 2 — Apodization windows applied to baseline-corrected signal")
fig2.tight_layout()


# ═══════════════════════════════════════════════════════════════════════════════
# Part 3 — Zero padding and its effect on FFT resolution
# ═══════════════════════════════════════════════════════════════════════════════
# Zero padding does not add information: it interpolates the FFT grid.
# The result is a smoother spectrum that makes peak positions easier to read,
# especially when the signal is short.  Factor 4 is a common choice.

print("Part 3 — Zero padding ...")

sig_apo = apodize(sig_bl_exp, window="hann", half_window="right")
n_orig = len(sig_apo)

pad_configs = {
    "No padding (1×)": sig_apo,
    "Factor 2": zero_pad(sig_apo, factor=2),
    "Factor 4": zero_pad(sig_apo, factor=4),
    "1024 points": zero_pad(sig_apo, n_points=1024),
}

fig3, axes = plt.subplots(1, len(pad_configs), figsize=(15, 4), sharey=True)

for ax, (label, padded) in zip(axes, pad_configs.items()):
    n = len(padded)
    dt = (t_ns[1] - t_ns[0]) * 1e-9  # seconds
    freqs = np.fft.rfftfreq(n, dt) / 1e6  # MHz
    spectrum = np.abs(np.fft.rfft(padded)) ** 2
    spectrum /= spectrum.max() if spectrum.max() > 0 else 1.0

    ax.plot(freqs, spectrum, lw=1.0)
    ax.set_xlabel("Frequency (MHz)")
    ax.set_title(f"{label}\n({n} points, Δf = {freqs[1]:.2f} MHz)")
    ax.set_xlim(0, 80)
    ax.grid(True, ls=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].set_ylabel("Normalized power")
fig3.suptitle("Part 3 — Zero padding: effect on FFT frequency resolution")
fig3.tight_layout()


# ═══════════════════════════════════════════════════════════════════════════════
# Part 4 — Full 1D pipeline: manual preprocessing then analyze_frequencies
# ═══════════════════════════════════════════════════════════════════════════════
# The three preprocessing functions are chained explicitly.  The result is
# passed to analyze_frequencies with window=None and remove_dc=False because
# apodization and baseline removal have already been applied.

print("Part 4 — Full 1D pipeline ...")

# Step 1: exponential baseline removal
sig_step1, _ = remove_baseline(t_ns, sig1d, method="exponential")

# Step 2: right-half Hann apodization (signal decays from t=0)
sig_step2 = apodize(sig_step1, window="hann", half_window="right")

# Step 3: zero-pad ×4 for smooth spectrum
sig_step3 = zero_pad(sig_step2, factor=4)

# Step 4: FFT on the full zero-padded signal
dt = (t_ns[1] - t_ns[0]) * 1e-9  # seconds
freqs_full = np.fft.rfftfreq(len(sig_step3), dt) / 1e6  # MHz
spectrum_full = np.abs(np.fft.rfft(sig_step3)) ** 2
spectrum_full /= spectrum_full.max()

# Peak detection: local maxima above 10 % of maximum
from scipy.signal import find_peaks
peak_idx, _ = find_peaks(spectrum_full, height=0.1)
dominant = sorted(freqs_full[peak_idx[spectrum_full[peak_idx].argsort()[::-1]]])
print(f"  Dominant frequencies: {[f'{f:.3f}' for f in dominant[:5]]} MHz")

fig4, axes4 = plt.subplots(1, 4, figsize=(16, 4))

axes4[0].plot(t_ns, sig1d, lw=1)
axes4[0].set_title("① Raw signal")
axes4[0].set_xlabel("Pulse length (ns)")

axes4[1].plot(t_ns, sig_step1, lw=1)
axes4[1].axhline(0, color="k", lw=0.5)
axes4[1].set_title("② Baseline removed (exp)")
axes4[1].set_xlabel("Pulse length (ns)")

axes4[2].plot(t_ns, sig_step2, lw=1, color="C2")
axes4[2].set_title("③ Apodized (Hann right)")
axes4[2].set_xlabel("Pulse length (ns)")

axes4[3].plot(freqs_full, spectrum_full, lw=1, color="C3")
for f in dominant[:3]:
    axes4[3].axvline(f, ls="--", color="k", lw=0.8, alpha=0.6)
axes4[3].set_xlim(0, 80)
axes4[3].set_xlabel("Frequency (MHz)")
axes4[3].set_ylabel("Normalized power")
axes4[3].set_title("④ FFT spectrum (×4 padding)")

for ax in axes4:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, ls=":", alpha=0.5)

axes4[0].set_ylabel("Signal (a.u.)")
fig4.suptitle("Part 4 — Full 1D pipeline: remove_baseline → apodize → zero_pad → FFT")
fig4.tight_layout()


# ═══════════════════════════════════════════════════════════════════════════════
# Part 5 — 2D HYSCORE: background normalization, apodization, and 2D FFT
# ═══════════════════════════════════════════════════════════════════════════════
# HYSCORE echo height is always negative (phase convention) and decays along
# both t1 and t2 following T_2.  The nuclear modulations are oscillations ON
# TOP of this large negative background.  Polynomial subtraction on such data
# diverges; the correct procedure is:
#   (1) divide each t1-slice (row) by its mean absolute value to isolate the
#       fractional modulation depth;
#   (2) subtract the residual row mean (DC along t2);
#   (3) subtract the residual column mean (DC along t1).
# After this normalization, apodize(axis='both') and zero_pad(axis='both')
# work exactly as in the 1D case.

print("Part 5 — 2D HYSCORE pipeline ...")

# Background normalization: isolate fractional nuclear modulation
row_bg = np.abs(sig2d.mean(axis=1, keepdims=True))  # t1-dependent background
sig2d_mod = sig2d / row_bg  # dimensionless modulation (values near ±1)
sig2d_mod -= sig2d_mod.mean(axis=1, keepdims=True)  # remove DC along t2
sig2d_mod -= sig2d_mod.mean(axis=0, keepdims=True)  # remove DC along t1
sig2d_bl = sig2d_mod

# 2D Hann apodization
sig2d_apo = apodize(sig2d_bl, window="hann", axis="both")

# Zero-pad both axes ×2
sig2d_zp = zero_pad(sig2d_apo, factor=2, axis="both")

# Full 2D FFT on the zero-padded signal
fft2d = np.fft.fftshift(np.fft.fft2(sig2d_zp))
n1, n2 = sig2d_zp.shape
dt1 = (t1_ns[1] - t1_ns[0]) * 1e-9
dt2 = (t2_ns[1] - t2_ns[0]) * 1e-9
fq1_mhz = np.fft.fftshift(np.fft.fftfreq(n1, dt1)) / 1e6
fq2_mhz = np.fft.fftshift(np.fft.fftfreq(n2, dt2)) / 1e6
spectrum2d_full = np.abs(fft2d)

fig5, axes5 = plt.subplots(1, 3, figsize=(15, 5))

im0 = axes5[0].imshow(
    sig2d, aspect="auto", cmap="RdBu_r", origin="lower",
    extent=[t2_ns[0], t2_ns[-1], t1_ns[0], t1_ns[-1]],
)
axes5[0].set_xlabel("t₂ (ns)")
axes5[0].set_ylabel("t₁ (ns)")
axes5[0].set_title("Raw HYSCORE signal")
plt.colorbar(im0, ax=axes5[0], fraction=0.046)

im1 = axes5[1].imshow(
    sig2d_apo, aspect="auto", cmap="RdBu_r", origin="lower",
    extent=[t2_ns[0], t2_ns[-1], t1_ns[0], t1_ns[-1]],
)
axes5[1].set_xlabel("t₂ (ns)")
axes5[1].set_ylabel("t₁ (ns)")
axes5[1].set_title("After baseline removal + 2D apodization")
plt.colorbar(im1, ax=axes5[1], fraction=0.046)

vmax = np.percentile(spectrum2d_full, 99)
im2 = axes5[2].imshow(
    spectrum2d_full,
    aspect="auto",
    cmap="hot",
    origin="lower",
    extent=[fq2_mhz[0], fq2_mhz[-1], fq1_mhz[0], fq1_mhz[-1]],
    vmax=vmax,
)
axes5[2].set_xlim(-30, 30)
axes5[2].set_ylim(-30, 30)
axes5[2].set_xlabel("ν₂ (MHz)")
axes5[2].set_ylabel("ν₁ (MHz)")
axes5[2].set_title("2D FFT magnitude (×2 zero-padding)")
plt.colorbar(im2, ax=axes5[2], fraction=0.046)

fig5.suptitle("Part 5 — 2D HYSCORE preprocessing pipeline")
fig5.tight_layout()

plt.show()
