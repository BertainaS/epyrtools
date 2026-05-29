"""
Interactive 2D EPR slicer with adjustable color scale.

Dataset: HYSCORE spectrum (256 × 256, t1 × t2 in ns).

Controls
--------
Slice        : slider — navigate through slices along the selected axis.
Vmin / Vmax  : range slider — adjust the color scale of the 2D map.
Direction    : radio buttons — switch between horizontal (fix t2) and
               vertical (fix t1) slicing.

Requirements: matplotlib >= 3.3 (RangeSlider), %matplotlib widget in Jupyter
              or a Qt/Tk interactive backend in a terminal.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RadioButtons, RangeSlider, Slider

import epyr

DATA_DIR = Path(__file__).parents[1] / "data"


# ── Load HYSCORE dataset ──────────────────────────────────────────────────────

x, y, params, _ = epyr.eprload(
    DATA_DIR / "Hyscore.DSC",
    plot_if_possible=False,
)

data  = np.real(y)      # (256, 256) — real part of the 2D FT
ny, nx = data.shape
t1_ns  = x[0]           # t1 axis, ns — shape (nx,)
t2_ns  = x[1]           # t2 axis, ns — shape (ny,)

x_label = f"t₁ ({params.get('XAXIS_UNIT', 'ns')})"
y_label = f"t₂ ({params.get('YAXIS_UNIT', 'ns')})"

print(f"HYSCORE  {ny} × {nx} points   "
      f"t1 = t2 = {t1_ns[0]:.0f}–{t1_ns[-1]:.0f} ns")


# ── Figure layout ─────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(14, 7))
fig.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.30)

ax_map   = fig.add_subplot(1, 2, 1)
ax_slice = fig.add_subplot(1, 2, 2)


# ── 2D color map ──────────────────────────────────────────────────────────────

# Use 2nd–98th percentile as default color range to avoid outlier saturation
vmin0 = float(np.percentile(data, 2))
vmax0 = float(np.percentile(data, 98))
data_min = float(data.min())
data_max = float(data.max())

mesh = ax_map.pcolormesh(
    t1_ns, t2_ns, data,
    shading="auto", cmap="RdBu_r", vmin=vmin0, vmax=vmax0,
)
fig.colorbar(mesh, ax=ax_map, label="Intensity (a.u.)", fraction=0.046)
ax_map.set_xlabel(x_label)
ax_map.set_ylabel(y_label)
ax_map.set_aspect("equal")

# Red dashed line showing the current slice position in the 2D map
(h_indicator,) = ax_map.plot([], [], "r--", lw=1.1, alpha=0.9)
(v_indicator,) = ax_map.plot([], [], "r--", lw=1.1, alpha=0.9)


# ── 1D slice panel ────────────────────────────────────────────────────────────

(slice_line,) = ax_slice.plot([], [], lw=0.9, color="C0")
ax_slice.set_ylabel("Intensity (a.u.)")
ax_slice.grid(True, linestyle=":", alpha=0.5)
ax_slice.ticklabel_format(style="sci", axis="y", scilimits=(-3, 4))


# ── Widget axes ───────────────────────────────────────────────────────────────

ax_idx   = fig.add_axes([0.07, 0.20, 0.58, 0.03])
ax_range = fig.add_axes([0.07, 0.11, 0.58, 0.03])
ax_radio = fig.add_axes([0.74, 0.04, 0.22, 0.20])

s_idx = Slider(
    ax_idx, "Slice", 0, ny - 1,
    valinit=ny // 2, valstep=1, valfmt="%d",
)
s_range = RangeSlider(
    ax_range, "Vmin / Vmax",
    data_min, data_max,
    valinit=(vmin0, vmax0),
)
RADIO_LABELS = ["Horizontal (fix t₂)", "Vertical (fix t₁)"]
radio = RadioButtons(ax_radio, RADIO_LABELS, active=0)
horizontal = True


# ── Shared update function ────────────────────────────────────────────────────

def update(_=None):
    idx        = min(int(s_idx.val), (ny if horizontal else nx) - 1)
    vmin, vmax = s_range.val

    if horizontal:
        # Fix t₂, sweep along t₁
        slice_vals = data[idx, :]
        axis_vals  = t1_ns
        t2_val     = t2_ns[idx]
        h_indicator.set_data([t1_ns[0], t1_ns[-1]], [t2_val, t2_val])
        v_indicator.set_data([], [])
        ax_slice.set_xlabel(x_label)
        ax_slice.set_title(
            f"Horizontal slice — t₂ = {t2_val:.0f} ns  (index {idx})"
        )
    else:
        # Fix t₁, sweep along t₂
        slice_vals = data[:, idx]
        axis_vals  = t2_ns
        t1_val     = t1_ns[idx]
        v_indicator.set_data([t1_val, t1_val], [t2_ns[0], t2_ns[-1]])
        h_indicator.set_data([], [])
        ax_slice.set_xlabel(y_label)
        ax_slice.set_title(
            f"Vertical slice — t₁ = {t1_val:.0f} ns  (index {idx})"
        )

    slice_line.set_data(axis_vals, slice_vals)
    ax_slice.relim()
    ax_slice.autoscale_view()

    mesh.set_clim(vmin, vmax)
    fig.canvas.draw_idle()


def on_direction(_):
    nonlocal horizontal
    horizontal = radio.value_selected == RADIO_LABELS[0]

    # Adjust slice slider range to match the new axis length
    n = ny if horizontal else nx
    s_idx.valmax = n - 1
    s_idx.ax.set_xlim(0, n - 1)
    # Keep current index within new bounds without re-triggering on_changed
    s_idx.eventson = False
    s_idx.set_val(min(int(s_idx.val), n - 1))
    s_idx.eventson = True
    update()


s_idx.on_changed(update)
s_range.on_changed(update)
radio.on_clicked(on_direction)

# Initial render
update()

fig.suptitle("Interactive HYSCORE slicer — use sliders and radio buttons",
             fontsize=11)
plt.show()
