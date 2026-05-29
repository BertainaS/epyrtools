"""
Effect of apodization windows on Fourier spectral analysis.

Dataset: Rabi oscillations of Gd³⁺ in CaWO4 (pulse EPR, X-band).
Each row of the 2D array is the echo height as a function of pulse length (ns).
The FFT of a trace gives the Rabi frequency; the choice of apodization window
controls the trade-off between frequency resolution and spectral leakage.

Windows compared
----------------
Rectangular  — no apodization; best resolution, highest leakage.
Hann         — cosine-squared taper; good general-purpose choice.
Hamming      — similar to Hann with slightly raised pedestal; flatter passband.
Blackman     — three-term cosine; lowest sidelobes (~-74 dB), widest main lobe.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import windows

import epyr
from epyr.eprplot import plot_2d_map
from epyr.signalprocessing.frequency_analysis import analyze_frequencies

DATA_DIR = Path(__file__).parents[1] / "data"


# --- Load 2D Rabi oscillation dataset ----------------------------------------

x, y, params, _ = epyr.eprload(
    DATA_DIR / "Rabi2D_GdCaWO4_13dB_3057G.DSC",
    plot_if_possible=False,
)

time_ns = x[0]                      # pulse length axis, ns
n_traces, n_pts = y.shape
dt_ns = time_ns[1] - time_ns[0]     # sampling interval, ns

print(f"Rabi dataset   {n_traces} field points × {n_pts} time points")
print(f"Pulse length   {time_ns[0]:.0f}–{time_ns[-1]:.0f} ns  "
      f"(step {dt_ns:.0f} ns)")

# Overview: 2D map of Rabi oscillations (real part)
fig0, ax0 = plot_2d_map(
    x, np.real(y), params,
    title="Rabi oscillations — 2D map, Gd:CaWO4 (real part)",
)


# --- Average over field points to improve SNR, then remove DC ----------------

trace_raw = np.real(y).mean(axis=0)
# DC offset removal is mandatory: a non-zero mean shifts all spectral weight to
# frequency zero and obscures the oscillation peak.
trace = trace_raw - trace_raw.mean()


# --- Manual window comparison ------------------------------------------------

window_defs = {
    "Rectangular (no window)": np.ones(n_pts),
    "Hann":                    windows.hann(n_pts),
    "Hamming":                 windows.hamming(n_pts),
    "Blackman":                windows.blackman(n_pts),
}

# Zero-pad by factor 4 to get a smooth frequency-domain display without
# adding spectral information (does not improve resolution).
n_fft = 4 * n_pts
# Frequency axis in MHz (dt in ns → dt*1e-3 in µs → freq in MHz)
freq_mhz = np.fft.rfftfreq(n_fft, d=dt_ns * 1e-3)

fig1, axes = plt.subplots(2, 2, figsize=(13, 8))
axes = axes.flatten()

for ax, (name, win) in zip(axes, window_defs.items()):
    windowed = trace * win
    spectrum = np.abs(np.fft.rfft(windowed, n=n_fft))
    spectrum /= spectrum.max()

    # Inset: windowed time-domain trace
    inset = ax.inset_axes([0.54, 0.44, 0.43, 0.52])
    inset.plot(time_ns, windowed / np.abs(windowed).max(), lw=0.6, color="C0")
    inset.set_xlabel("Time (ns)", fontsize=7)
    inset.set_title("windowed trace", fontsize=7)
    inset.tick_params(labelsize=6)

    ax.plot(freq_mhz, spectrum, lw=0.9)
    ax.set_xlim(0, 60)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title(name)
    ax.grid(True, linestyle=":", alpha=0.5)

fig1.suptitle(
    "Rabi frequency spectrum — effect of apodization window\n"
    "Gd:CaWO4, 13 dB, 3057 G"
)
fig1.tight_layout()


# --- EPyR analyze_frequencies (Hann window, built-in pipeline) ---------------

print("\nRunning analyze_frequencies (Hann window) ...")
# Pass time in µs so that returned frequencies are in MHz
result = analyze_frequencies(
    time_ns * 1e-3,        # µs
    trace,
    window="hann",
    zero_padding=4,
    remove_dc=False,       # DC already removed above
    plot=True,
)
freq_unit = result.get("freq_unit", "")
dominant = result.get("dominant_frequencies", [])
if len(dominant) > 0:
    print(f"  Dominant Rabi frequencies: {[f'{f:.3f}' for f in dominant]} {freq_unit}")


# --- Leakage demo on synthetic signal with two close frequencies -------------
# Shows how a weak second peak is hidden by the sidelobes of the rectangular
# window but becomes visible with Hann apodization.

print("\nSpectral leakage demo (synthetic signal) ...")
t_demo = np.arange(512) * dt_ns                     # ns
# Strong line at 20 MHz, weak line (1/20 amplitude) at 23 MHz
f1, f2 = 20.0, 23.0                                  # MHz
sig_demo = (
    np.cos(2 * np.pi * f1 * t_demo * 1e-3)
    + 0.05 * np.cos(2 * np.pi * f2 * t_demo * 1e-3)
)

n_demo = 4 * len(t_demo)
freq_demo = np.fft.rfftfreq(n_demo, d=dt_ns * 1e-3)  # MHz

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4))
for ax, (name, win_fn) in zip(
    axes2,
    [("Rectangular", np.ones(len(t_demo))),
     ("Hann",        windows.hann(len(t_demo)))],
):
    sp = np.abs(np.fft.rfft(sig_demo * win_fn, n=n_demo))
    sp /= sp.max()
    ax.semilogy(freq_demo, sp + 1e-6, lw=0.9)
    ax.axvline(f1, color="r", ls="--", lw=0.8, label=f"f₁ = {f1} MHz")
    ax.axvline(f2, color="g", ls="--", lw=0.8, label=f"f₂ = {f2} MHz (−26 dB)")
    ax.set_xlim(10, 35)
    ax.set_ylim(1e-4, 2)
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Amplitude (log scale)")
    ax.set_title(name)
    ax.legend(fontsize=8)
    ax.grid(True, linestyle=":", alpha=0.5)

fig2.suptitle(
    "Spectral leakage: weak line (f₂) hidden by sidelobes with rectangular window,\n"
    "resolved with Hann apodization"
)
fig2.tight_layout()

plt.show()
