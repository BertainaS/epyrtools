"""
Loading and visualizing EPR data from Bruker files.

Covers BES3T (.dsc/.dta) and ESP (.par/.spc) formats for both 1D CW spectra
and 2D datasets (angle rotation map). All files are in examples/data/.
"""

from pathlib import Path

import matplotlib.pyplot as plt

import epyr
from epyr.eprplot import plot_1d, plot_2d_map, plot_2d_waterfall

DATA_DIR = Path(__file__).parents[1] / "data"


# --- 1D CW spectra -----------------------------------------------------------

# BES3T format: CaWO4 doped with Er³⁺, CW EPR at 5 K, X-band
x, y, params, _ = epyr.eprload(
    DATA_DIR / "130406SB_CaWO4_Er_CW_5K_20.DSC",
    plot_if_possible=False,
)
print(f"CaWO4:Er   {y.shape[0]} pts  "
      f"field {x.min():.0f}–{x.max():.0f} G  "
      f"ν = {params.get('MWFQ', float('nan')):.4e} Hz")

fig1, ax1 = plot_1d(x, y, params, title="CaWO4:Er³⁺ — CW EPR at 5 K (BES3T)")


# ESP format: CuSO4 powder, room-temperature reference sample
x2, y2, params2, _ = epyr.eprload(
    DATA_DIR / "CuSO4_001.par",
    plot_if_possible=False,
)
print(f"CuSO4      {y2.shape[0]} pts  "
      f"field {x2.min():.0f}–{x2.max():.0f} G")

fig2, ax2 = plot_1d(x2, y2, params2, title="CuSO4 — CW EPR (ESP)")


# --- 2D rotation map ---------------------------------------------------------

# MgO crystal, 300 K, full angular rotation (0–180°) at X-band
# x[0]: magnetic field axis (G), x[1]: rotation angle (°)
x3, y3, params3, _ = epyr.eprload(
    DATA_DIR / "2014_03_19_MgO_300K_111_fullrotation33dB.par",
    plot_if_possible=False,
)
n_angles, n_field = y3.shape
print(f"MgO rot.   {n_angles} angles × {n_field} field pts  "
      f"field {x3[0].min():.0f}–{x3[0].max():.0f} G  "
      f"angle {x3[1].min():.0f}–{x3[1].max():.0f}°")

fig3, ax3 = plot_2d_map(
    x3, y3, params3,
    title="MgO — CW EPR rotation map (300 K)",
    cmap="RdBu_r",
)
fig4, ax4 = plot_2d_waterfall(
    x3, y3, params3,
    title="MgO — CW EPR waterfall (300 K)",
    offset_factor=0.4,
    max_traces=18,
)

plt.show()
