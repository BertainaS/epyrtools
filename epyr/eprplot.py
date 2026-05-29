#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple plotting module for EPR data from eprload.

This module provides simple plotting functions for data obtained with eprload():
- plot_1d: Plot 1D EPR spectra
- plot_2d_map: Plot 2D data as color map
- plot_2d_waterfall: Plot 2D data as waterfall plot

Based on the _plot_data function from eprload.py
"""

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from .logging_config import get_logger

logger = get_logger(__name__)


def _axis_label(
    params: Optional[Dict[str, Any]],
    name_key: str,
    unit_key: str,
    default_name: str,
    default_unit: str = "a.u.",
) -> str:
    """Build an axis label string from eprload params."""
    name = params.get(name_key, default_name) if params else default_name
    unit = params.get(unit_key, default_unit) if params else default_unit
    if isinstance(unit, list):
        unit = unit[0]
    return f"{name} ({unit})"


def plot_1d(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 1D EPR data.

    Parameters
    ----------
    x : np.ndarray, list, or None
        X-axis data from eprload. Falls back to point index if None or shape mismatch.
    y : np.ndarray
        1D EPR signal array.
    params : dict, optional
        Parameter dictionary from eprload, used to extract axis labels and units.
    title : str, optional
        Plot title.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created if not provided.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes

    Examples
    --------
    >>> from epyr import eprload, plot_1d
    >>> x, y, params, _ = eprload("examples/data/130406SB_CaWO4_Er_CW_5K_20.DSC")
    >>> fig, ax = plot_1d(x, y, params, title="CaWO4:Er, 5 K")

    Reuse an existing axes (subplot composition):

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(1, 2)
    >>> plot_1d(x, y, params, ax=axes[0])  # doctest: +SKIP
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")

    if y.ndim != 1:
        raise ValueError(f"Expected 1D data, got {y.ndim}D array.")

    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # Handle x-axis data
    if x is None or not isinstance(x, np.ndarray) or x.shape != y.shape:
        if x is not None:
            warnings.warn(
                "X-axis data missing or incompatible shape. Using index for plotting.",
                stacklevel=2,
            )
        absc = np.arange(y.size)
        x_label = "Index (points)"
    else:
        absc = x
        x_label = _axis_label(params, "XAXIS_NAME", "XAXIS_UNIT", "X Axis")

    # Plot data
    if np.isrealobj(y):
        ax.plot(absc, y, label="data", lw=0.75)
    else:
        ax.plot(absc, np.real(y), label="real", lw=0.75)
        ax.plot(absc, np.imag(y), label="imag", linestyle="--", lw=0.75)
        ax.legend()

    # Set labels and formatting
    ax.set_xlabel(x_label)
    ax.set_ylabel("Intensity (a.u.)")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(-3, 4))

    if title:
        ax.set_title(title)

    fig.tight_layout()
    return fig, ax


def plot_2d_map(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    cmap: str = "magma",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D EPR data as a color map.

    Parameters
    ----------
    x : np.ndarray, list, or None
        Axis data from eprload. A list of two arrays sets both x and y axes.
    y : np.ndarray
        2D EPR signal array, shape (ny, nx).
    params : dict, optional
        Parameter dictionary from eprload, used to extract axis labels and units.
    title : str, optional
        Plot title.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created if not provided.
    cmap : str, optional
        Matplotlib colormap name (default: "magma").
    vmin : float, optional
        Lower bound of the color scale. Defaults to data minimum.
    vmax : float, optional
        Upper bound of the color scale. Defaults to data maximum.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes

    Examples
    --------
    >>> fig, ax = plot_2d_map(x, y, params, vmin=-100, vmax=100)
    >>> fig, ax = plot_2d_map(x, y, params, vmin=0)  # clip negatives
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")

    if y.ndim != 2:
        raise ValueError(f"Expected 2D data, got {y.ndim}D array.")

    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots(layout="constrained")
    else:
        fig = ax.get_figure()

    ny, nx = y.shape

    # Default coordinates
    x_coords = np.arange(nx)
    y_coords = np.arange(ny)
    x_label = f"Index ({nx} points)"
    y_label = f"Index ({ny} points)"

    # Extract axis information
    if isinstance(x, list) and len(x) >= 2:
        x_axis, y_axis = x[0], x[1]
        if isinstance(x_axis, np.ndarray) and x_axis.size == nx:
            x_coords = x_axis
            x_label = _axis_label(params, "XAXIS_NAME", "XAXIS_UNIT", "X Axis")
        if isinstance(y_axis, np.ndarray) and y_axis.size == ny:
            y_coords = y_axis
            y_label = _axis_label(params, "YAXIS_NAME", "YAXIS_UNIT", "Y Axis")
    elif isinstance(x, np.ndarray) and x.size == nx:
        x_coords = x
        x_label = _axis_label(params, "XAXIS_NAME", "XAXIS_UNIT", "X Axis")

    # Plot data (real part if complex)
    plot_data = np.real(y)

    im = ax.pcolormesh(
        x_coords, y_coords, plot_data, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax
    )
    fig.colorbar(im, ax=ax, label="Intensity (a.u.)")

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_aspect("auto")

    if title:
        ax.set_title(title)

    return fig, ax


def plot_2d_waterfall(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    offset_factor: float = 0.5,
    max_traces: int = 20,
    cmap: str = "viridis",
    lw: float = 0.75,
    clip_factor: Optional[float] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D EPR data as a waterfall plot.

    Parameters
    ----------
    x : np.ndarray, list, or None
        Axis data from eprload. A list of two arrays sets both x and y axes.
    y : np.ndarray
        2D EPR signal array, shape (ny, nx).
    params : dict, optional
        Parameter dictionary from eprload, used to extract axis labels and units.
    title : str, optional
        Plot title.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created if not provided.
    offset_factor : float, optional
        Vertical spacing between traces as a fraction of the total data range
        (default: 0.5).
    max_traces : int, optional
        Maximum number of traces to display. If ny > max_traces, traces are
        subsampled uniformly (default: 20).
    cmap : str, optional
        Matplotlib colormap used to color-code traces (default: "viridis").
    lw : float, optional
        Line width for each trace (default: 0.75).
    clip_factor : float, optional
        Clip each trace at ``clip_factor * max(|trace|)`` before adding the
        vertical offset. None disables clipping (default). A value of 0.5
        clips at 50% of the trace maximum.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes

    Examples
    --------
    >>> from epyr import eprload, plot_2d_waterfall
    >>> x, y, params, _ = eprload("examples/data/Rabi2D_GdCaWO4_13dB_3057G.DSC")
    >>> fig, ax = plot_2d_waterfall(x, y, params, max_traces=10)

    Tighter spacing with strong clipping (useful for noisy backgrounds):

    >>> fig, ax = plot_2d_waterfall(x, y, params, offset_factor=0.2, clip_factor=0.3)
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")

    if y.ndim != 2:
        raise ValueError(f"Expected 2D data, got {y.ndim}D array.")

    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots(layout="constrained")
    else:
        fig = ax.get_figure()

    ny, nx = y.shape

    # Limit number of traces if too many
    if ny > max_traces:
        step = ny // max_traces
        trace_indices = np.arange(0, ny, step)
        warnings.warn(
            f"Too many traces ({ny}), showing every "
            f"{step}th trace ({len(trace_indices)} total).",
            stacklevel=2,
        )
    else:
        trace_indices = np.arange(ny)

    # Default x-axis
    if isinstance(x, list) and len(x) >= 1:
        x_axis = x[0]
    elif isinstance(x, np.ndarray):
        x_axis = x
    else:
        x_axis = None

    if x_axis is None or not isinstance(x_axis, np.ndarray) or x_axis.size != nx:
        x_coords = np.arange(nx)
        x_label = f"Index ({nx} points)"
    else:
        x_coords = x_axis
        x_label = _axis_label(params, "XAXIS_NAME", "XAXIS_UNIT", "X Axis")

    # Get y-axis parameter name for labeling
    y_param_name = params.get("YAXIS_NAME", "Parameter") if params else "Parameter"

    # Calculate offset
    plot_data = np.real(y)  # Use real part if complex
    data_range = np.ptp(plot_data)
    offset = data_range * offset_factor

    colormap = matplotlib.colormaps[cmap]
    colors = colormap(np.linspace(0, 1, len(trace_indices)))

    # Plot traces
    for i, trace_idx in enumerate(trace_indices):
        trace_raw = plot_data[trace_idx, :]

        # Apply clipping if requested
        if clip_factor is not None:
            threshold = clip_factor * np.max(np.abs(trace_raw))
            trace_raw = np.clip(trace_raw, -threshold, threshold)

        y_offset = i * offset
        trace_data = trace_raw + y_offset

        # Create label
        if isinstance(x, list) and len(x) >= 2:
            y_axis = x[1]
            if isinstance(y_axis, np.ndarray) and trace_idx < len(y_axis):
                param_value = y_axis[trace_idx]
                label = f"{y_param_name}={param_value:.2f}"
            else:
                label = f"{y_param_name}[{trace_idx}]"
        else:
            label = f"{y_param_name}[{trace_idx}]"

        ax.plot(
            x_coords, trace_data, label=label if i < 10 else "", color=colors[i], lw=lw
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel(f"Intensity + {y_param_name} offset (a.u.)")

    # Show legend only for first few traces
    if len(trace_indices) <= 10:
        ax.legend(loc="upper right")

    if title:
        ax.set_title(title)

    return fig, ax


def plot_2d_slicer(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    slice_direction: str = "horizontal",
    cmap: str = "magma",
) -> Dict[str, Any]:
    """
    Interactive 2D EPR data slicer with a slider control.

    Parameters
    ----------
    x : np.ndarray, list, or None
        Axis data from eprload. A list of two arrays sets both x and y axes.
    y : np.ndarray
        2D EPR signal array, shape (ny, nx).
    params : dict, optional
        Parameter dictionary from eprload, used to extract axis labels and units.
    title : str, optional
        Plot title shown above the active slice panel.
    slice_direction : {'horizontal', 'vertical'}, optional
        Axis along which to slice (default: 'horizontal').
    cmap : str, optional
        Matplotlib colormap name for the overview panel (default: "magma").

    Returns
    -------
    dict
        Keys: 'figure', 'ax_main', 'ax_overview', 'slider', 'line', 'slice_line'.

    Notes
    -----
    Requires an interactive matplotlib backend. In Jupyter, activate with
    ``%matplotlib widget`` or ``%matplotlib notebook`` before calling.

    Examples
    --------
    >>> from epyr import eprload, plot_2d_slicer
    >>> x, y, params, _ = eprload("examples/data/Rabi2D_GdCaWO4_13dB_3057G.DSC")
    >>> handles = plot_2d_slicer(x, y, params, slice_direction="vertical")
    >>> handles["slider"].set_val(50)  # programmatic slider control  # doctest: +SKIP
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")

    if y.ndim != 2:
        raise ValueError(f"Expected 2D data, got {y.ndim}D array.")

    # Import widgets here to avoid errors if not available
    try:
        from matplotlib.widgets import Slider
    except ImportError:
        raise ImportError(
            "matplotlib.widgets required for interactive function. "
            "Use %matplotlib widget in Jupyter."
        )

    # Use real part if data is complex
    plot_data = np.real(y)
    ny, nx = plot_data.shape

    # Configure axes according to direction
    if slice_direction == "horizontal":
        n_slices = ny
        slice_axis_name = "Y"
    else:  # vertical
        n_slices = nx
        slice_axis_name = "X"
        plot_data = plot_data.T  # Transpose for vertical slices

    # Prepare axes
    if isinstance(x, list) and len(x) >= 1:
        x_axis = (
            x[0] if isinstance(x[0], np.ndarray) and x[0].size == nx else np.arange(nx)
        )
        y_axis = (
            x[1]
            if len(x) >= 2 and isinstance(x[1], np.ndarray) and x[1].size == ny
            else np.arange(ny)
        )
    elif isinstance(x, np.ndarray) and x.size == nx:
        x_axis = x
        y_axis = np.arange(ny)
    else:
        x_axis = np.arange(nx)
        y_axis = np.arange(ny)

    # Determine axes and labels
    x_label = _axis_label(params, "XAXIS_NAME", "XAXIS_UNIT", "Field", "G")
    y_label = _axis_label(params, "YAXIS_NAME", "YAXIS_UNIT", "Parameter")
    if slice_direction == "horizontal":
        slice_values = y_axis
        plot_axis = x_axis
        plot_label = x_label
        slice_label = y_label
    else:
        slice_values = x_axis
        plot_axis = y_axis
        plot_label = y_label
        slice_label = x_label

    # Create figure and axes
    fig, (ax_main, ax_overview) = plt.subplots(
        2, 1, gridspec_kw={"height_ratios": [3, 1]}
    )

    # Adjust space for slider
    plt.subplots_adjust(bottom=0.15)

    # Overview (2D map)
    if slice_direction == "horizontal":
        overview_data = np.real(y)
        extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
        ax_overview.imshow(
            overview_data, aspect="auto", extent=extent, origin="lower", cmap=cmap
        )
    else:
        overview_data = np.real(y).T
        extent = [y_axis[0], y_axis[-1], x_axis[0], x_axis[-1]]
        ax_overview.imshow(
            overview_data, aspect="auto", extent=extent, origin="lower", cmap=cmap
        )

    ax_overview.set_xlabel(plot_label)
    ax_overview.set_ylabel(slice_label)
    ax_overview.set_title("Overview - Current slice position shown in red")

    # Indicator line for current slice
    if slice_direction == "horizontal":
        slice_line = ax_overview.axhline(y=slice_values[0], color="red")
    else:
        slice_line = ax_overview.axvline(x=slice_values[0], color="red")

    # Initial plot of first slice
    (line,) = ax_main.plot(plot_axis, plot_data[0], "b-")
    ax_main.set_xlabel(plot_label)
    ax_main.set_ylabel("EPR Intensity")
    ax_main.grid(True, alpha=0.3)
    ax_main.set_xlim(plot_axis[0], plot_axis[-1])

    # Initial title
    initial_title = title or "Interactive 2D EPR Viewer"
    slice_value = slice_values[0]
    ax_main.set_title(f"{initial_title} - {slice_label} = {slice_value:.3f}")

    # Create slider axis
    ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])
    slider = Slider(
        ax_slider,
        f"{slice_axis_name} Index",
        0,
        n_slices - 1,
        valinit=0,
        valfmt="%d",
        valstep=1,
    )

    # Slider update function
    def update_slice(val):
        idx = int(slider.val)

        # Update main plot
        line.set_ydata(plot_data[idx])

        # Update title with parameter value
        slice_value = slice_values[idx]
        ax_main.set_title(f"{initial_title} - {slice_label} = {slice_value:.3f}")

        # Update indicator line
        if slice_direction == "horizontal":
            slice_line.set_ydata([slice_value, slice_value])
        else:
            slice_line.set_xdata([slice_value, slice_value])

        # Auto-adjust Y scale
        y_min, y_max = np.min(plot_data[idx]), np.max(plot_data[idx])
        y_range = y_max - y_min
        ax_main.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        fig.canvas.draw_idle()

    # Connect slider to update function
    slider.on_changed(update_slice)

    logger.info(
        "Interactive 2D EPR Viewer: direction=%s, slices=%d", slice_direction, n_slices
    )

    # Show plot
    plt.show()

    # Return objects for advanced manipulation if needed
    return {
        "figure": fig,
        "ax_main": ax_main,
        "ax_overview": ax_overview,
        "slider": slider,
        "line": line,
        "slice_line": slice_line,
    }


# Define public API
__all__ = ["plot_1d", "plot_2d_map", "plot_2d_waterfall", "plot_2d_slicer"]
