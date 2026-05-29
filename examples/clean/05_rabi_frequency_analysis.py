"""
Rabi oscillation frequency analysis using frequency_analysis.py.

Datasets: Gd³⁺ in CaWO4, X-band pulse EPR, two power levels.
  - 13 dB attenuation, B₀ = 3057 G  →  lower power, slower Rabi
  - 6 dB attenuation,  B₀ = 3770 G  →  higher power, faster Rabi

The 2D datasets are (n_traces × n_time_points):
  - rows   : field / condition index (500 repetitions)
  - columns: echo height vs. pulse length (ns)  ← FFT dimension

Functions demonstrated
----------------------
analyze_frequencies       : 1D FFT on a single or averaged trace.
analyze_frequencies_2d    : row-by-row FFT → Rabi frequency map.
power_spectrum            : Welch PSD for noise characterization.
spectrogram_analysis      : time-frequency view of a single trace.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import epyr
from epyr.signalprocessing.frequency_analysis import (
    analyze_frequencies,
    analyze_frequencies_2d,
    power_spectrum,
    spectrogram_analysis,
)

DATA_DIR = Path(__file__).parents[1] / "data"

# ── Load both Rabi datasets ───────────────────────────────────────────────────

x1, y1, p1, _ = epyr.eprload(
    DATA_DIR / "Rabi2D_GdCaWO4_13dB_3057G.DSC", plot_if_possible=False
)
x2, y2, p2, _ = epyr.eprload(
    DATA_DIR / "Rabi2D_GdCaWO4_6dB_3770G_2.DSC", plot_if_possible=False
)

# Work with real part of the echo signal
sig1 = np.real(y1)     # (500, 1024)  dt = 2 ns  → max freq 250 MHz
sig2 = np.real(y2)     # (500, 1024)  dt = 1 ns  → max freq 500 MHz
t1_ns = x1[0]          # pulse length axis for dataset 1, ns
t2_ns = x2[0]          # pulse length axis for dataset 2, ns

print("Rabi dataset 1  (13 dB, 3057 G):")
print(f"  shape = {sig1.shape}  t = {t1_ns[0]:.0f}–{t1_ns[-1]:.0f} ns"
      f"  dt = {t1_ns[1]-t1_ns[0]:.0f} ns")
print("Rabi dataset 2  (6 dB,  3770 G):")
print(f"  shape = {sig2.shape}  t = {t2_ns[0]:.0f}–{t2_ns[-1]:.0f} ns"
      f"  dt = {t2_ns[1]-t2_ns[0]:.0f} ns")


# ── Part 1 — Time-domain overview ─────────────────────────────────────────────
# Average over all traces to improve SNR, then show a representative single trace.

avg1 = sig1.mean(axis=0)
avg2 = sig2.mean(axis=0)

fig, axes = plt.subplots(2, 2, figsize=(13, 6))

axes[0, 0].imshow(sig1, aspect="auto", cmap="RdBu_r",
                   extent=[t1_ns[0], t1_ns[-1], 0, sig1.shape[0]],
                   origin="lower")
axes[0, 0].set_xlabel("Pulse length (ns)")
axes[0, 0].set_ylabel("Trace index")
axes[0, 0].set_title("13 dB — raw 2D signal (real part)")

axes[0, 1].plot(t1_ns, avg1, lw=0.9)
axes[0, 1].set_xlabel("Pulse length (ns)")
axes[0, 1].set_ylabel("Echo height (a.u.)")
axes[0, 1].set_title("13 dB — average over 500 traces")
axes[0, 1].grid(True, linestyle=":", alpha=0.5)

axes[1, 0].imshow(sig2, aspect="auto", cmap="RdBu_r",
                   extent=[t2_ns[0], t2_ns[-1], 0, sig2.shape[0]],
                   origin="lower")
axes[1, 0].set_xlabel("Pulse length (ns)")
axes[1, 0].set_ylabel("Trace index")
axes[1, 0].set_title("6 dB — raw 2D signal (real part)")

axes[1, 1].plot(t2_ns, avg2, lw=0.9)
axes[1, 1].set_xlabel("Pulse length (ns)")
axes[1, 1].set_ylabel("Echo height (a.u.)")
axes[1, 1].set_title("6 dB — average over 500 traces")
axes[1, 1].grid(True, linestyle=":", alpha=0.5)

fig.suptitle("Part 1 — Time-domain overview")
fig.tight_layout()


# ── Part 2 — Single trace FFT with analyze_frequencies ───────────────────────
# analyze_frequencies handles DC removal, windowing, zero-padding automatically.
# The built-in plot shows: raw signal, power spectrum (log + linear), windowed trace.

print("\nPart 2 — FFT of averaged trace (13 dB) ...")
result1 = analyze_frequencies(
    t1_ns, avg1,
    window="hann",
    zero_padding=4,
    remove_dc=True,
    plot=True,
    freq_range=(0, 50),       # MHz — zoom on relevant range
)
freqs1 = result1["dominant_frequencies"]
print(f"  Dominant Rabi frequencies: {[f'{f:.2f}' for f in freqs1[:5]]} MHz")

print("\nPart 2 — FFT of averaged trace (6 dB) ...")
result2 = analyze_frequencies(
    t2_ns, avg2,
    window="hann",
    zero_padding=4,
    remove_dc=True,
    plot=True,
    freq_range=(0, 50),
)
freqs2 = result2["dominant_frequencies"]
print(f"  Dominant Rabi frequencies: {[f'{f:.2f}' for f in freqs2[:5]]} MHz")


# ── Part 3 — Power dependence: overlay both spectra ──────────────────────────
# Compare Rabi frequencies at two microwave power levels on the same axes.

fig3, ax3 = plt.subplots(figsize=(10, 4))
ax3.plot(result1["frequencies"], result1["power_spectrum"],
         lw=0.9, label="13 dB (B₀ = 3057 G)", color="C0")
ax3.plot(result2["frequencies"], result2["power_spectrum"],
         lw=0.9, label="6 dB  (B₀ = 3770 G)", color="C1", alpha=0.8)
ax3.set_xlim(0, 50)
ax3.set_xlabel("Frequency (MHz)")
ax3.set_ylabel("Normalized power")
ax3.set_title("Part 3 — Rabi frequency comparison: 13 dB vs 6 dB")
ax3.legend()
ax3.grid(True, linestyle=":", alpha=0.5)
fig3.tight_layout()


# ── Part 4 — 2D row-by-row FFT → Rabi frequency map ─────────────────────────
# Process each of the 500 traces independently.
# Returns: fq (frequency axis, MHz), axis2 (trace index), spectrum (500 × n_freq).
# The spectrum is centered (fftshift) → symmetric around 0. Show positive only.

print("\nPart 4 — 2D row-by-row FFT (13 dB) ...")
fq, trace_idx, spectrum, info = analyze_frequencies_2d(
    t1_ns, sig1,
    mode="row_by_row",
    window="hann",
    zero_padding=2,
    remove_dc=True,
    plot_result=True,         # built-in 4-panel plot (raw, FFT lin, FFT log, phase)
    freq_range=(-50, 50),
)
print(f"  Frequency axis: {fq[0]:.1f} – {fq[-1]:.1f} {info['freq_unit']}")
print(f"  Spectrum shape: {spectrum.shape}")

# Custom view: clip to positive frequencies, zoom 0–50 MHz
pos = fq >= 0
fq_pos = fq[pos]
spec_pos = spectrum[:, pos]

fig4, axes4 = plt.subplots(1, 2, figsize=(13, 5))

axes4[0].pcolormesh(fq_pos, trace_idx, spec_pos,
                     shading="auto", cmap="hot",
                     vmax=np.percentile(spec_pos, 99))
axes4[0].set_xlim(0, 50)
axes4[0].set_xlabel("Frequency (MHz)")
axes4[0].set_ylabel("Trace index")
axes4[0].set_title("Rabi frequency map — linear scale")

axes4[1].pcolormesh(fq_pos, trace_idx, np.log10(spec_pos + 1),
                     shading="auto", cmap="hot")
axes4[1].set_xlim(0, 50)
axes4[1].set_xlabel("Frequency (MHz)")
axes4[1].set_ylabel("Trace index")
axes4[1].set_title("Rabi frequency map — log scale")

fig4.suptitle("Part 4 — 2D Rabi frequency map (13 dB, Hann window)")
fig4.tight_layout()


# ── Part 5 — Welch PSD for noise characterization ────────────────────────────
# The Welch method averages periodograms over overlapping segments, reducing
# variance at the cost of frequency resolution. Useful for identifying
# noise floor and broadband features.

print("\nPart 5 — Welch PSD (13 dB averaged trace) ...")
psd_result = power_spectrum(
    t1_ns, avg1,
    method="welch",
    window="hann",
    nperseg=256,          # segment length (< 1024 → better averaging)
    overlap=0.75,
    remove_dc=True,
    plot=True,
)
print(f"  PSD method: {psd_result['method']}  "
      f"  freq unit: {psd_result['freq_unit']}")


# ── Part 6 — Spectrogram (time–frequency) ────────────────────────────────────
# Shows how the Rabi oscillation frequency evolves along the pulse length axis.
# For a simple Rabi signal the pattern should be stationary; deviations reveal
# dephasing or off-resonance effects.

print("\nPart 6 — Spectrogram of averaged trace (13 dB) ...")
sg_result = spectrogram_analysis(
    t1_ns, avg1,
    window="hann",
    nperseg=128,
    overlap=0.90,
    remove_dc=True,
    plot=True,
)
print(f"  Spectrogram shape: {sg_result['spectrogram'].shape}  "
      f"  freq unit: {sg_result['freq_unit']}")

plt.show()
