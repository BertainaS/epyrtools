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
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np


def plot_1d(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 1D EPR data.
    
    Args:
        x: X-axis data from eprload (can be None, array, or list)
        y: 1D spectral data array
        params: Parameter dictionary from eprload
        title: Plot title (optional)
        ax: Matplotlib axes to plot on (optional)
        
    Returns:
        Tuple of (figure, axes)
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")
    
    if y.ndim != 1:
        raise ValueError(f"Expected 1D data, got {y.ndim}D array.")
    
    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()
    
    # Handle x-axis data
    if x is None or not isinstance(x, np.ndarray) or x.shape != y.shape:
        if x is not None:
            warnings.warn("X-axis data missing or incompatible shape. Using index for plotting.")
        absc = np.arange(y.size)
        x_label = "Index (points)"
    else:
        absc = x
        x_label = params.get("XAXIS_NAME", "X Axis") if params else "X Axis"
        x_unit = params.get("XAXIS_UNIT", "a.u.") if params else "a.u."
        if isinstance(x_unit, list):
            x_unit = x_unit[0]
        if x_unit and x_unit != "a.u.":
            x_label += f" ({x_unit})"
    
    # Plot data
    if np.isrealobj(y):
        ax.plot(absc, y, label="data", linewidth=1.2)
    else:
        ax.plot(absc, np.real(y), label="real", linewidth=1.2)
        ax.plot(absc, np.imag(y), label="imag", linewidth=1.2, linestyle='--')
        ax.legend()
    
    # Set labels and formatting
    ax.set_xlabel(x_label)
    ax.set_ylabel("Intensity (a.u.)")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(-3, 4))
    
    if title:
        ax.set_title(title, fontsize=12)
    
    plt.tight_layout()
    return fig, ax


def plot_2d_map(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    cmap: str = "viridis"
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D EPR data as a color map.
    
    Args:
        x: Axis data from eprload (can be None, array, or list of arrays)
        y: 2D spectral data array
        params: Parameter dictionary from eprload
        title: Plot title (optional)
        ax: Matplotlib axes to plot on (optional)
        cmap: Colormap name
        
    Returns:
        Tuple of (figure, axes)
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")
    
    if y.ndim != 2:
        raise ValueError(f"Expected 2D data, got {y.ndim}D array.")
    
    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
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
            x_name = params.get("XAXIS_NAME", "X Axis") if params else "X Axis"
            x_unit = params.get("XAXIS_UNIT", "a.u.") if params else "a.u."
            x_label = f"{x_name} ({x_unit})"
        if isinstance(y_axis, np.ndarray) and y_axis.size == ny:
            y_coords = y_axis
            y_name = params.get("YAXIS_NAME", "Y Axis") if params else "Y Axis"
            y_unit = params.get("YAXIS_UNIT", "a.u.") if params else "a.u."
            y_label = f"{y_name} ({y_unit})"
    elif isinstance(x, np.ndarray) and x.size == nx:
        x_coords = x
        x_name = params.get("XAXIS_NAME", "X Axis") if params else "X Axis"
        x_unit = params.get("XAXIS_UNIT", "a.u.") if params else "a.u."
        x_label = f"{x_name} ({x_unit})"
    
    # Plot data (real part if complex)
    plot_data = np.real(y)
    
    im = ax.pcolormesh(x_coords, y_coords, plot_data, shading="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label="Intensity (a.u.)")
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_aspect("auto")
    
    if title:
        ax.set_title(title, fontsize=12)
    
    plt.tight_layout()
    return fig, ax


def plot_2d_waterfall(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    offset_factor: float = 0.5,
    max_traces: int = 50
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D EPR data as waterfall plot.
    
    Args:
        x: Axis data from eprload (can be None, array, or list of arrays)
        y: 2D spectral data array
        params: Parameter dictionary from eprload
        title: Plot title (optional)
        ax: Matplotlib axes to plot on (optional)
        offset_factor: Vertical offset between traces
        max_traces: Maximum number of traces to plot
        
    Returns:
        Tuple of (figure, axes)
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")
    
    if y.ndim != 2:
        raise ValueError(f"Expected 2D data, got {y.ndim}D array.")
    
    # Create figure if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig = ax.get_figure()
    
    ny, nx = y.shape
    
    # Limit number of traces if too many
    if ny > max_traces:
        step = ny // max_traces
        trace_indices = np.arange(0, ny, step)
        warnings.warn(f"Too many traces ({ny}), showing every {step}th trace ({len(trace_indices)} total).")
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
        x_name = params.get("XAXIS_NAME", "X Axis") if params else "X Axis"
        x_unit = params.get("XAXIS_UNIT", "a.u.") if params else "a.u."
        x_label = f"{x_name} ({x_unit})"
    
    # Get y-axis parameter name for labeling
    y_param_name = params.get("YAXIS_NAME", "Parameter") if params else "Parameter"
    
    # Calculate offset
    plot_data = np.real(y)  # Use real part if complex
    data_range = np.ptp(plot_data)
    offset = data_range * offset_factor
    
    # Plot traces
    for i, trace_idx in enumerate(trace_indices):
        y_offset = i * offset
        trace_data = plot_data[trace_idx, :] + y_offset
        
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
        
        ax.plot(x_coords, trace_data, linewidth=0.8, label=label if i < 10 else "")
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(f"Intensity + {y_param_name} offset (a.u.)")
    
    # Show legend only for first few traces
    if len(trace_indices) <= 10:
        ax.legend(fontsize=8, loc='upper right')
    
    if title:
        ax.set_title(title, fontsize=12)
    
    plt.tight_layout()
    return fig, ax


def plot_2d_slicer(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: np.ndarray,
    params: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    slice_direction: str = 'horizontal'
) -> None:
    """
    Interactive 2D EPR data slicer with slider control.
    
    Allows visualization of 2D EPR data slice by slice with an interactive
    slicer to navigate in both directions.
    
    Args:
        x: X-axis data from eprload (can be None, array, or list)
        y: 2D spectral data array (ny, nx)
        params: Parameter dictionary from eprload
        title: Plot title (optional)
        slice_direction: 'horizontal' for horizontal slices, 'vertical' for vertical slices
        
    Note:
        Uses matplotlib widgets for interactivity. Works in Jupyter
        with %matplotlib widget or %matplotlib notebook.
    """
    if y is None or y.size == 0:
        raise ValueError("No data available to plot.")
    
    if y.ndim != 2:
        raise ValueError(f"Expected 2D data, got {y.ndim}D array.")
    
    # Import widgets here to avoid errors if not available
    try:
        from matplotlib.widgets import Slider
    except ImportError:
        raise ImportError("matplotlib.widgets required for interactive function. "
                         "Use %matplotlib widget in Jupyter.")
    
    # Use real part if data is complex
    plot_data = np.real(y)
    ny, nx = plot_data.shape
    
    # Configure axes according to direction
    if slice_direction == 'horizontal':
        n_slices = ny
        slice_axis_name = "Y"
        plot_axis_name = "X"
    else:  # vertical
        n_slices = nx  
        slice_axis_name = "X"
        plot_axis_name = "Y"
        plot_data = plot_data.T  # Transpose for vertical slices
    
    # Prepare axes
    if isinstance(x, list) and len(x) >= 1:
        x_axis = x[0] if isinstance(x[0], np.ndarray) and x[0].size == nx else np.arange(nx)
        y_axis = x[1] if len(x) >= 2 and isinstance(x[1], np.ndarray) and x[1].size == ny else np.arange(ny)
    elif isinstance(x, np.ndarray) and x.size == nx:
        x_axis = x
        y_axis = np.arange(ny)
    else:
        x_axis = np.arange(nx)
        y_axis = np.arange(ny)
    
    # Determine axes and labels
    if slice_direction == 'horizontal':
        slice_values = y_axis
        plot_axis = x_axis
        x_name = params.get("XAXIS_NAME", "Field") if params else "Field"
        x_unit = params.get("XAXIS_UNIT", "G") if params else "G"
        y_name = params.get("YAXIS_NAME", "Parameter") if params else "Parameter"  
        y_unit = params.get("YAXIS_UNIT", "a.u.") if params else "a.u."
        plot_label = f"{x_name} ({x_unit})"
        slice_label = f"{y_name} ({y_unit})"
    else:
        slice_values = x_axis
        plot_axis = y_axis
        x_name = params.get("XAXIS_NAME", "Field") if params else "Field"
        x_unit = params.get("XAXIS_UNIT", "G") if params else "G" 
        y_name = params.get("YAXIS_NAME", "Parameter") if params else "Parameter"
        y_unit = params.get("YAXIS_UNIT", "a.u.") if params else "a.u."
        plot_label = f"{y_name} ({y_unit})"
        slice_label = f"{x_name} ({x_unit})"
    
    # Create figure and axes
    fig, (ax_main, ax_overview) = plt.subplots(2, 1, figsize=(10, 10), 
                                               gridspec_kw={'height_ratios': [3, 1]})
    
    # Adjust space for slider
    plt.subplots_adjust(bottom=0.15)
    
    # Overview (2D map)
    if slice_direction == 'horizontal':
        overview_data = np.real(y)
        extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
        ax_overview.imshow(overview_data, aspect='auto', extent=extent, 
                          origin='lower', cmap='RdBu_r')
    else:
        overview_data = np.real(y).T
        extent = [y_axis[0], y_axis[-1], x_axis[0], x_axis[-1]]
        ax_overview.imshow(overview_data, aspect='auto', extent=extent, 
                          origin='lower', cmap='RdBu_r')
    
    ax_overview.set_xlabel(plot_label)
    ax_overview.set_ylabel(slice_label)
    ax_overview.set_title("Overview - Current slice position shown in red")
    
    # Indicator line for current slice
    if slice_direction == 'horizontal':
        slice_line = ax_overview.axhline(y=slice_values[0], color='red', linewidth=2)
    else:
        slice_line = ax_overview.axvline(x=slice_values[0], color='red', linewidth=2)
    
    # Initial plot of first slice
    line, = ax_main.plot(plot_axis, plot_data[0], 'b-', linewidth=2)
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
    slider = Slider(ax_slider, f'{slice_axis_name} Index', 0, n_slices-1, 
                   valinit=0, valfmt='%d', valstep=1)
    
    # Slider update function
    def update_slice(val):
        idx = int(slider.val)
        
        # Update main plot
        line.set_ydata(plot_data[idx])
        
        # Update title with parameter value
        slice_value = slice_values[idx]
        ax_main.set_title(f"{initial_title} - {slice_label} = {slice_value:.3f}")
        
        # Update indicator line
        if slice_direction == 'horizontal':
            slice_line.set_ydata([slice_value, slice_value])
        else:
            slice_line.set_xdata([slice_value, slice_value])
        
        # Auto-adjust Y scale
        y_min, y_max = np.min(plot_data[idx]), np.max(plot_data[idx])
        y_range = y_max - y_min
        ax_main.set_ylim(y_min - 0.1*y_range, y_max + 0.1*y_range)
        
        fig.canvas.draw_idle()
    
    # Connect slider to update function
    slider.on_changed(update_slice)
    
    # User instructions
    print("🎛️  Interactive 2D EPR Viewer")
    print("="*50)
    print(f"Direction: {slice_direction}")
    print(f"Number of slices: {n_slices}")
    print(f"Use the slider to navigate between slices")
    print("Red line in overview shows current position")
    
    # Show plot
    plt.show()
    
    # Return objects for advanced manipulation if needed
    return {
        'figure': fig,
        'ax_main': ax_main,
        'ax_overview': ax_overview,
        'slider': slider,
        'line': line,
        'slice_line': slice_line
    }


# Define public API
__all__ = ['plot_1d', 'plot_2d_map', 'plot_2d_waterfall', 'plot_2d_slicer']