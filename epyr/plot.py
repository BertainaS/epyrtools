#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for specialized plotting of EPR (Electron Paramagnetic Resonance) data.

This module provides functions for creating waterfall plots and 2D maps
of EPR spectral data, commonly used in angular sweep experiments.

Author: sylvainbertaina
Created: Tue May 6 14:54:54 2025
"""

import warnings
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

# Define public API - only export these functions
__all__ = [
    'EPRPlotConfig',
    'plot_angular_sweep_waterfall',
    'plot_2d_spectral_map', 
    'plot_epr_comparison',
    'generate_example_data',
    'run_examples'
]


class EPRPlotConfig:
    """Configuration class for EPR plotting parameters."""

    # Default plotting parameters
    DEFAULT_FIGSIZE_WATERFALL = (10, 8)
    DEFAULT_FIGSIZE_2D = (10, 8)
    DEFAULT_OFFSET_FACTOR = 0.2
    DEFAULT_CMAP = "viridis"
    DEFAULT_SHADING = "auto"
    MAX_LEGEND_LINES = 10


def _validate_data_shape(
    data_2d: np.ndarray,
    axis1: np.ndarray,
    axis2: np.ndarray,
    axis1_name: str = "axis1",
    axis2_name: str = "axis2",
) -> None:
    """
    Validate that 2D data shape matches the provided axes.

    Args:
        data_2d: 2D spectral data array
        axis1: First axis array
        axis2: Second axis array
        axis1_name: Name of first axis for error messages
        axis2_name: Name of second axis for error messages

    Raises:
        ValueError: If shapes don't match
    """
    expected_shape = (len(axis1), len(axis2))
    if data_2d.shape != expected_shape:
        raise ValueError(
            f"Shape mismatch: data shape {data_2d.shape} "
            f"does not match ({axis1_name}={len(axis1)}, {axis2_name}={len(axis2)})"
        )


def _calculate_offset_parameters(
    data_2d: np.ndarray, offset_factor: float
) -> Tuple[float, float, float]:
    """
    Calculate offset parameters for waterfall plots.

    Args:
        data_2d: 2D spectral data
        offset_factor: Factor for vertical offset calculation

    Returns:
        Tuple of (data_min, data_max, vertical_offset_per_spectrum)
    """
    data_min = np.nanmin(data_2d)
    data_max = np.nanmax(data_2d)
    data_range = data_max - data_min

    # Handle flat data
    if data_range == 0:
        data_range = 1.0
        warnings.warn("Data has zero range, using default offset", UserWarning)

    return data_min, data_max, offset_factor * data_range


def _setup_colorbar_for_lines(
    fig: plt.Figure,
    ax: plt.Axes,
    angle_axis: np.ndarray,
    angle_unit: str,
    cmap_name: str,
) -> None:
    """
    Set up colorbar for line plots when there are too many lines for a legend.

    Args:
        fig: Matplotlib figure
        ax: Matplotlib axes
        angle_axis: Array of angle values
        angle_unit: Unit for angles
        cmap_name: Name of colormap
    """
    norm = plt.Normalize(vmin=angle_axis.min(), vmax=angle_axis.max())
    cmap = cm.get_cmap(cmap_name)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation="vertical", label=f"Angle ({angle_unit})")


def _format_axis_label(param_name: str, unit: str) -> str:
    """Format axis label with unit if provided."""
    return f"{param_name} ({unit})" if unit else param_name


def plot_angular_sweep_waterfall(
    field_axis: np.ndarray,
    angle_axis: np.ndarray,
    spectral_data_2d: np.ndarray,
    field_unit: str = "G",
    angle_unit: str = "deg",
    intensity_label: str = "Intensity (a.u.)",
    offset_factor: float = EPRPlotConfig.DEFAULT_OFFSET_FACTOR,
    integrate: bool = False,
    ax: Optional[plt.Axes] = None,
    title: str = "Angular Sweep Waterfall Plot",
    cmap_name: str = EPRPlotConfig.DEFAULT_CMAP,
    line_label_prefix: str = "Angle:",
    figsize: Tuple[float, float] = EPRPlotConfig.DEFAULT_FIGSIZE_WATERFALL,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a waterfall plot for angular sweep EPR data.

    Each angle's spectrum is plotted as a line with vertical offset.
    Lines are colored based on their angle value using a colormap.

    Args:
        field_axis: 1D array for magnetic field axis
        angle_axis: 1D array for angular positions
        spectral_data_2d: 2D array of spectral intensities
                         Shape: (len(angle_axis), len(field_axis))
        field_unit: Unit for field axis
        angle_unit: Unit for angle axis
        intensity_label: Label for y-axis
        offset_factor: Factor for vertical offset between spectra
        integrate: If True, plot cumulative sum (integration) of spectra
        ax: Matplotlib axes to plot on (creates new if None)
        title: Plot title
        cmap_name: Matplotlib colormap name
        line_label_prefix: Prefix for line labels (empty string disables labels)
        figsize: Figure size for new plots

    Returns:
        Tuple of (figure, axes)

    Raises:
        ValueError: If data shape doesn't match axes
    """
    # Validate input data
    _validate_data_shape(
        spectral_data_2d, angle_axis, field_axis, "angle_axis", "field_axis"
    )

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    num_angles = len(angle_axis)
    cmap = cm.get_cmap(cmap_name)

    # Calculate offset parameters
    data_min, data_max, vertical_offset = _calculate_offset_parameters(
        spectral_data_2d, offset_factor
    )

    # Plot each spectrum
    for i in range(num_angles):
        spectrum = spectral_data_2d[i, :]
        offset = i * vertical_offset

        # Normalize color index
        color_index = i / max(1, num_angles - 1)
        color = cmap(color_index)

        # Create label if requested
        label = None
        if line_label_prefix:
            label = f"{line_label_prefix} {angle_axis[i]:.1f} {angle_unit}"

        # Plot spectrum (integrated or regular)
        y_data = np.cumsum(spectrum) + offset if integrate else spectrum + offset
        ax.plot(field_axis, y_data, color=color, label=label, linewidth=1.0)

    # Format axes
    ax.set_xlabel(_format_axis_label("Magnetic Field", field_unit))
    ax.set_ylabel(intensity_label)
    ax.set_title(title)
    ax.set_xlim(field_axis.min(), field_axis.max())

    # Remove y-ticks as they represent offset + intensity
    ax.set_yticks([])

    # Set up legend or colorbar
    if line_label_prefix:
        if num_angles <= EPRPlotConfig.MAX_LEGEND_LINES:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
        else:
            _setup_colorbar_for_lines(fig, ax, angle_axis, angle_unit, cmap_name)

    # Calculate and set y-limits with padding
    min_y = np.nanmin(spectral_data_2d[0, :])
    max_y = np.nanmax(spectral_data_2d[-1, :]) + (num_angles - 1) * vertical_offset
    y_range = max_y - min_y
    if y_range == 0:
        y_range = 1

    padding = 0.1 * y_range
    # ax.set_ylim(min_y - padding, max_y + padding)

    plt.tight_layout()
    return fig, ax


def plot_2d_spectral_map(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    spectral_data_2d: np.ndarray,
    x_label: str = "Parameter 1",
    y_label: str = "Parameter 2",
    x_unit: str = "",
    y_unit: str = "",
    intensity_label: str = "Intensity (a.u.)",
    ax: Optional[plt.Axes] = None,
    title: str = "2D Spectral Map",
    cmap_name: str = EPRPlotConfig.DEFAULT_CMAP,
    normalize: bool = True,
    shading: str = EPRPlotConfig.DEFAULT_SHADING,
    figsize: Tuple[float, float] = EPRPlotConfig.DEFAULT_FIGSIZE_2D,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a 2D color map of spectral data using pcolormesh.

    Args:
        x_axis: 1D array for x-axis values
        y_axis: 1D array for y-axis values
        spectral_data_2d: 2D array of spectral intensities
                         Shape: (len(y_axis), len(x_axis))
        x_label: Label for x-axis parameter
        y_label: Label for y-axis parameter
        x_unit: Unit for x-axis
        y_unit: Unit for y-axis
        intensity_label: Label for colorbar
        ax: Matplotlib axes to plot on (creates new if None)
        title: Plot title
        cmap_name: Matplotlib colormap name
        normalize: If True, normalize data to [0, 1] range
        shading: Shading option for pcolormesh
        figsize: Figure size for new plots
        vmin: Minimum value for colorscale (None for auto)
        vmax: Maximum value for colorscale (None for auto)

    Returns:
        Tuple of (figure, axes)

    Raises:
        ValueError: If data shape doesn't match axes
    """
    # Validate input data
    _validate_data_shape(spectral_data_2d, y_axis, x_axis, "y_axis", "x_axis")

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Process data
    plot_data = spectral_data_2d.copy()
    if normalize:
        data_max = np.nanmax(plot_data)
        if data_max != 0:
            plot_data = plot_data / data_max
        else:
            warnings.warn("Cannot normalize: maximum value is zero", UserWarning)

    # Create mesh plot
    im = ax.pcolormesh(
        x_axis, y_axis, plot_data, cmap=cmap_name, shading=shading, vmin=vmin, vmax=vmax
    )

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, label=intensity_label)

    # Format axes
    ax.set_xlabel(_format_axis_label(x_label, x_unit))
    ax.set_ylabel(_format_axis_label(y_label, y_unit))
    ax.set_title(title)
    ax.set_xlim(x_axis.min(), x_axis.max())
    ax.set_ylim(y_axis.min(), y_axis.max())

    plt.tight_layout()
    return fig, ax


def plot_epr_comparison(
    field_axis: np.ndarray,
    spectra_dict: dict,
    field_unit: str = "G",
    intensity_label: str = "Intensity (a.u.)",
    normalize: bool = True,
    ax: Optional[plt.Axes] = None,
    title: str = "EPR Spectra Comparison",
    figsize: Tuple[float, float] = (10, 6),
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot multiple EPR spectra for comparison.

    Args:
        field_axis: 1D array for magnetic field axis
        spectra_dict: Dictionary with labels as keys and 1D spectra as values
        field_unit: Unit for field axis
        intensity_label: Label for y-axis
        normalize: If True, normalize each spectrum to its maximum
        ax: Matplotlib axes to plot on (creates new if None)
        title: Plot title
        figsize: Figure size for new plots

    Returns:
        Tuple of (figure, axes)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    for label, spectrum in spectra_dict.items():
        if len(spectrum) != len(field_axis):
            raise ValueError(
                f"Spectrum '{label}' length {len(spectrum)} "
                f"doesn't match field_axis length {len(field_axis)}"
            )

        plot_spectrum = spectrum.copy()
        if normalize:
            spec_max = np.nanmax(np.abs(plot_spectrum))
            if spec_max != 0:
                plot_spectrum = plot_spectrum / spec_max

        ax.plot(field_axis, plot_spectrum, label=label, linewidth=1.5)

    ax.set_xlabel(_format_axis_label("Magnetic Field", field_unit))
    ax.set_ylabel(intensity_label)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig, ax


# Convenience aliases for backward compatibility
plot_2d_map = plot_2d_spectral_map


def generate_example_data():
    """
    Generate synthetic EPR data for testing and demonstration purposes.

    Returns:
        dict: Dictionary containing example datasets with keys:
            - 'field_axis': Magnetic field values
            - 'angle_axis': Angular positions
            - 'angular_sweep_data': 2D array for angular sweep
            - 'temp_axis': Temperature values
            - 'temp_sweep_data': 2D array for temperature sweep
            - 'single_spectra': Dictionary of individual spectra
    """
    # Magnetic field axis (Gauss)
    field_axis = np.linspace(3200, 3600, 400)

    # Angular positions (degrees)
    angle_axis = np.linspace(0, 180, 19)  # 10-degree steps

    # Temperature axis (Kelvin)
    temp_axis = np.linspace(4, 300, 25)

    # Generate synthetic EPR signal with angular dependence
    angular_sweep_data = np.zeros((len(angle_axis), len(field_axis)))

    for i, angle in enumerate(angle_axis):
        # Main resonance with angular dependence
        center1 = 3400 + 50 * np.cos(np.radians(angle * 2))
        width1 = 20 + 10 * np.sin(np.radians(angle))

        # Secondary resonance
        center2 = 3350 + 30 * np.sin(np.radians(angle))
        width2 = 15

        # Create Lorentzian-like peaks with noise
        signal1 = 100 * width1**2 / ((field_axis - center1) ** 2 + width1**2)
        signal2 = 60 * width2**2 / ((field_axis - center2) ** 2 + width2**2)

        # Add derivative-like component (common in EPR)
        derivative = np.gradient(signal1)

        # Combine signals with some noise
        total_signal = signal1 + signal2 + 0.3 * derivative
        noise = np.random.normal(0, 2, len(field_axis))

        angular_sweep_data[i, :] = total_signal + noise

    # Generate temperature-dependent data
    temp_sweep_data = np.zeros((len(temp_axis), len(field_axis)))

    for i, temp in enumerate(temp_axis):
        # Temperature-dependent line width and intensity
        intensity_factor = 1 / (1 + temp / 50)  # Curie-like behavior
        linewidth = 15 + temp * 0.2  # Thermal broadening

        center = 3400
        signal = (
            intensity_factor
            * 80
            * linewidth**2
            / ((field_axis - center) ** 2 + linewidth**2)
        )
        noise = np.random.normal(0, 1, len(field_axis))

        temp_sweep_data[i, :] = signal + noise

    # Generate individual spectra for comparison
    single_spectra = {}

    # Room temperature spectrum
    rt_signal = 50 * 20**2 / ((field_axis - 3400) ** 2 + 20**2)
    single_spectra["Room Temperature"] = rt_signal + np.random.normal(
        0, 1, len(field_axis)
    )

    # Low temperature spectrum
    lt_signal = 80 * 12**2 / ((field_axis - 3405) ** 2 + 12**2)
    single_spectra["Low Temperature"] = lt_signal + np.random.normal(
        0, 0.5, len(field_axis)
    )

    # Derivative spectrum
    derivative_signal = np.gradient(rt_signal)
    single_spectra["Derivative"] = derivative_signal + np.random.normal(
        0, 0.5, len(field_axis)
    )

    return {
        "field_axis": field_axis,
        "angle_axis": angle_axis,
        "angular_sweep_data": angular_sweep_data,
        "temp_axis": temp_axis,
        "temp_sweep_data": temp_sweep_data,
        "single_spectra": single_spectra,
    }


def run_examples():
    """
    Run comprehensive examples demonstrating all plotting functions.
    """
    print("Generating example EPR data...")
    data = generate_example_data()

    # Example 1: Angular sweep waterfall plot
    print("Creating angular sweep waterfall plot...")
    fig1, ax1 = plot_angular_sweep_waterfall(
        field_axis=data["field_axis"],
        angle_axis=data["angle_axis"],
        spectral_data_2d=data["angular_sweep_data"],
        field_unit="G",
        angle_unit="deg",
        title="EPR Angular Sweep - Waterfall View",
        cmap_name="plasma",
        offset_factor=0.3,
    )
    plt.show()

    # Example 2: Angular sweep waterfall with integration
    print("Creating integrated angular sweep plot...")
    fig2, ax2 = plot_angular_sweep_waterfall(
        field_axis=data["field_axis"],
        angle_axis=data["angle_axis"],
        spectral_data_2d=data["angular_sweep_data"],
        integrate=True,
        title="EPR Angular Sweep - Integrated Spectra",
        cmap_name="viridis",
        intensity_label="Cumulative Intensity (a.u.)",
    )
    plt.show()

    # Example 3: 2D spectral map of angular sweep
    print("Creating 2D map of angular sweep data...")
    fig3, ax3 = plot_2d_spectral_map(
        x_axis=data["field_axis"],
        y_axis=data["angle_axis"],
        spectral_data_2d=data["angular_sweep_data"],
        x_label="Magnetic Field",
        y_label="Rotation Angle",
        x_unit="G",
        y_unit="deg",
        title="EPR Angular Sweep - 2D Map",
        cmap_name="RdYlBu_r",
    )
    plt.show()

    # Example 4: Temperature-dependent 2D map
    print("Creating temperature-dependent 2D map...")
    fig4, ax4 = plot_2d_spectral_map(
        x_axis=data["field_axis"],
        y_axis=data["temp_axis"],
        spectral_data_2d=data["temp_sweep_data"],
        x_label="Magnetic Field",
        y_label="Temperature",
        x_unit="G",
        y_unit="K",
        title="EPR Temperature Dependence",
        cmap_name="hot",
        normalize=True,
    )
    plt.show()

    # Example 5: Comparison plot of individual spectra
    print("Creating spectra comparison plot...")
    fig5, ax5 = plot_epr_comparison(
        field_axis=data["field_axis"],
        spectra_dict=data["single_spectra"],
        title="EPR Spectra Comparison",
        normalize=True,
    )
    plt.show()

    # Example 6: Waterfall with many angles (triggers colorbar instead of legend)
    print("Creating high-resolution angular sweep...")
    # Generate more angles for colorbar demonstration
    fine_angles = np.linspace(0, 360, 37)  # Every 10 degrees
    fine_data = np.zeros((len(fine_angles), len(data["field_axis"])))

    for i, angle in enumerate(fine_angles):
        center = 3400 + 40 * np.cos(np.radians(angle))
        width = 18
        signal = 70 * width**2 / ((data["field_axis"] - center) ** 2 + width**2)
        noise = np.random.normal(0, 1.5, len(data["field_axis"]))
        fine_data[i, :] = signal + noise

    fig6, ax6 = plot_angular_sweep_waterfall(
        field_axis=data["field_axis"],
        angle_axis=fine_angles,
        spectral_data_2d=fine_data,
        title="High Resolution Angular Sweep (360°)",
        cmap_name="hsv",  # Good for full rotation
        angle_unit="deg",
        offset_factor=0.15,
    )
    plt.show()

    # Example 7: Custom styling example
    print("Creating custom styled plot...")
    fig7, ax7 = plt.subplots(figsize=(12, 8))

    # Plot angular sweep with custom styling
    plot_angular_sweep_waterfall(
        field_axis=data["field_axis"],
        angle_axis=data["angle_axis"][::2],  # Every other angle
        spectral_data_2d=data["angular_sweep_data"][::2, :],
        ax=ax7,
        title="Custom Styled EPR Angular Sweep",
        cmap_name="coolwarm",
        line_label_prefix="θ =",
        offset_factor=0.4,
    )

    # Add custom styling
    ax7.grid(True, alpha=0.3)
    ax7.set_facecolor("lightgray")
    plt.show()

    print("All examples completed!")


if __name__ == "__main__":
    print("EPR Data Plotting Module")
    print("=" * 50)
    print()
    print("This module provides functions for plotting EPR spectroscopic data:")
    print("1. plot_angular_sweep_waterfall() - Waterfall plots for angular sweeps")
    print("2. plot_2d_spectral_map() - 2D color maps of spectral data")
    print("3. plot_epr_comparison() - Comparison plots of multiple spectra")
    print()

    # Ask user if they want to run examples
    try:
        user_input = (
            input("Would you like to run the examples? (y/n): ").lower().strip()
        )
        if user_input in ["y", "yes", "1", "true"]:
            run_examples()
        else:
            print("Examples skipped. Import this module to use the plotting functions.")
            print("\nExample usage:")
            print(
                "from epr_plotting import plot_angular_sweep_waterfall, generate_example_data"
            )
            print("data = generate_example_data()")
            print(
                "fig, ax = plot_angular_sweep_waterfall(data['field_axis'], data['angle_axis'], data['angular_sweep_data'])"
            )
    except (KeyboardInterrupt, EOFError):
        print("\nExamples skipped.")

    print("\nModule loaded successfully!")
