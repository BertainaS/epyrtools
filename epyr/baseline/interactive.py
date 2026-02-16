#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interactive region selection for baseline correction.

This module provides matplotlib-based interactive widgets for selecting
baseline regions in Jupyter notebooks and desktop environments.
"""

from typing import List, Optional, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RectangleSelector, SpanSelector

from ..logging_config import get_logger

logger = get_logger(__name__)


def _get_widget_props_param():
    """Get the correct parameter name for matplotlib widget properties based on version."""
    # Matplotlib 3.5+ uses 'props' instead of 'rectprops'
    version = tuple(map(int, matplotlib.__version__.split(".")[:2]))
    if version >= (3, 5):
        return "props"
    else:
        return "rectprops"


class RegionSelector:
    """
    Interactive region selector for baseline correction.

    This class provides matplotlib-based interactive region selection
    for both 1D and 2D EPR data. It handles matplotlib version compatibility
    and provides multiple methods to close selection windows.
    """

    def __init__(self):
        self.regions = []
        self.current_selector = None
        self.fig = None
        self.ax = None
        self.selection_done = False

    def _on_select_1d(self, xmin, xmax):
        """Callback for 1D region selection."""
        self.regions.append((xmin, xmax))
        logger.info(f"Selected region: {xmin:.2f} - {xmax:.2f}")

    def _on_key_press(self, event):
        """Handle key press events."""
        if event.key == "enter" or event.key == "escape":
            self.selection_done = True
            plt.close(self.fig)
            logger.info("✅ Region selection completed!")

    def finish_selection(self):
        """Manually finish selection and close plot."""
        self.selection_done = True
        if self.fig:
            plt.close(self.fig)
        logger.info("✅ Region selection completed!")

    def _on_select_2d(self, eclick, erelease):
        """Callback for 2D region selection."""
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        region = ((min(x1, x2), max(x1, x2)), (min(y1, y2), max(y1, y2)))
        self.regions.append(region)
        logger.info(
            f"Selected region: x={region[0][0]:.2f}-{region[0][1]:.2f}, y={region[1][0]:.2f}-{region[1][1]:.2f}"
        )

    def select_regions_1d(
        self, x, y, title="Select regions to EXCLUDE from baseline fitting"
    ):
        """
        Interactive selection of 1D regions.

        Args:
            x: X-axis data
            y: Y-axis data
            title: Plot title with instructions

        Returns:
            list: Selected regions as [(x1, x2), ...]
        """
        self.regions = []

        self.fig, self.ax = plt.subplots(figsize=(12, 6))
        self.ax.plot(x, y, "b-", alpha=0.7)
        self.ax.set_title(
            f"{title}\nClick and drag to select regions.\nPress ENTER or ESC when done, or run selector.finish_selection()"
        )
        self.ax.grid(True, alpha=0.3)

        # Add keyboard event handling
        self.fig.canvas.mpl_connect("key_press_event", self._on_key_press)

        # Create span selector - handle matplotlib version compatibility
        props_param = _get_widget_props_param()
        selector_kwargs = {
            "useblit": True,
            props_param: dict(alpha=0.3, facecolor="red"),
        }
        self.current_selector = SpanSelector(
            self.ax, self._on_select_1d, "horizontal", **selector_kwargs
        )

        plt.show()
        return self.regions

    def select_regions_2d(
        self, x, y, z, title="Select regions to EXCLUDE from baseline fitting"
    ):
        """
        Interactive selection of 2D regions.

        Args:
            x: X-axis coordinates (1D array or meshgrid)
            y: Y-axis coordinates (1D array or meshgrid)
            z: 2D data array
            title: Plot title with instructions

        Returns:
            list: Selected regions as [((x1,x2), (y1,y2)), ...]
        """
        self.regions = []

        self.fig, self.ax = plt.subplots(figsize=(10, 8))

        # Handle coordinate arrays
        if isinstance(x, np.ndarray) and x.ndim == 1:
            X, Y = np.meshgrid(x, y)
        else:
            X, Y = x, y

        im = self.ax.pcolormesh(X, Y, z, shading="auto", cmap="viridis")
        self.fig.colorbar(im, ax=self.ax)

        self.ax.set_title(
            f"{title}\nClick and drag to select rectangular regions.\nPress ENTER or ESC when done, or run selector.finish_selection()"
        )

        # Add keyboard event handling
        self.fig.canvas.mpl_connect("key_press_event", self._on_key_press)

        # Create rectangle selector - handle matplotlib version compatibility
        props_param = _get_widget_props_param()
        selector_kwargs = {
            "useblit": True,
            "button": [1],
            "minspanx": 5,
            "minspany": 5,
            props_param: dict(alpha=0.3, facecolor="red", edgecolor="red", linewidth=2),
        }
        self.current_selector = RectangleSelector(
            self.ax, self._on_select_2d, **selector_kwargs
        )

        plt.show()
        return self.regions


# Global selector instance for Jupyter compatibility
_current_selector = None


def interactive_select_regions_1d(
    x, y, title="Select regions to EXCLUDE from baseline fitting"
):
    """
    Convenience function for interactive 1D region selection.

    Args:
        x: X-axis data
        y: Y-axis data
        title: Plot title

    Returns:
        list: Selected regions as [(x1, x2), ...]
    """
    global _current_selector
    _current_selector = RegionSelector()
    return _current_selector.select_regions_1d(x, y, title)


def interactive_select_regions_2d(
    x, y, z, title="Select regions to EXCLUDE from baseline fitting"
):
    """
    Convenience function for interactive 2D region selection.

    Args:
        x: X-axis coordinates
        y: Y-axis coordinates
        z: 2D data array
        title: Plot title

    Returns:
        list: Selected regions as [((x1,x2), (y1,y2)), ...]
    """
    global _current_selector
    _current_selector = RegionSelector()
    return _current_selector.select_regions_2d(x, y, z, title)


def close_selector_window():
    """
    Utility function to close RegionSelector windows in Jupyter notebooks.

    Use this if the interactive region selector gets stuck or won't close.
    """
    try:
        global _current_selector
        if _current_selector is not None:
            _current_selector.finish_selection()
        plt.close("all")
        logger.info("✅ All selector windows closed")
    except Exception as e:
        plt.close("all")
        logger.info(f"✅ Windows closed (with warning: {e})")


def jupyter_help():
    """
    Display help for using interactive region selection in Jupyter notebooks.
    """
    help_text = """
    📋 JUPYTER NOTEBOOK - INTERACTIVE REGION SELECTION HELP
    ====================================================
    
    When you run interactive baseline correction, a plot will appear.
    
    HOW TO SELECT REGIONS:
    1. Click and drag on the plot to select regions to exclude from baseline fitting
    2. You can select multiple regions
    
    HOW TO FINISH SELECTION:
    Method 1: Press ENTER or ESC key on the plot
    Method 2: In a new cell, run: from epyr.baseline.interactive import close_selector_window; close_selector_window()
    Method 3: In a new cell, run: plt.close('all')
    
    IF STUCK:
    - Run: plt.close('all') to force close all plots
    - Run: from epyr.baseline.interactive import close_selector_window; close_selector_window()
    
    EXAMPLE:
    --------
    import epyr
    x, y, params, filepath = epyr.eprload("data.dsc", plot_if_possible=False)
    
    # Option 1: Direct baseline correction with interactive selection
    corrected, baseline = epyr.baseline_polynomial_1d(x, y, params, interactive=True)
    
    # Option 2: Manual region selection first
    from epyr.baseline.interactive import interactive_select_regions_1d
    regions = interactive_select_regions_1d(x, y, "Select baseline regions")
    corrected, baseline = epyr.baseline_polynomial_1d(x, y, params, 
                                                     manual_regions=regions,
                                                     region_mode='include')
    
    # If the plot won't close, run in a new cell:
    from epyr.baseline.interactive import close_selector_window
    close_selector_window()
    """
    logger.info(help_text)


def is_interactive_available():
    """
    Check if interactive selection is available in the current environment.

    Returns:
        bool: True if interactive selection should work
    """
    try:
        # Check if we're in a notebook
        from IPython import get_ipython

        ipython = get_ipython()

        if ipython is None:
            return True  # Assume desktop matplotlib works

        # Check for notebook backends
        backend = matplotlib.get_backend().lower()
        if "inline" in backend:
            logger.warning(
                "⚠️  Warning: Inline backend detected. Interactive selection may not work."
            )
            logger.warning("   Try: %matplotlib widget or %matplotlib notebook")
            return False
        elif "widget" in backend or "nbagg" in backend:
            return True
        else:
            return True  # Assume it works
    except ImportError:
        return True  # Not in Jupyter, assume desktop matplotlib works


def setup_interactive_backend():
    """
    Set up the best available interactive backend for the current environment.
    """
    try:
        from IPython import get_ipython

        ipython = get_ipython()

        if ipython is not None:
            # We're in Jupyter
            current_backend = matplotlib.get_backend().lower()
            if "inline" in current_backend:
                logger.info("🔧 Setting up interactive backend...")
                try:
                    ipython.magic("matplotlib widget")
                    logger.info("✅ Switched to widget backend")
                except:
                    try:
                        ipython.magic("matplotlib notebook")
                        logger.info("✅ Switched to notebook backend")
                    except:
                        logger.warning("⚠️  Could not switch to interactive backend")
                        logger.warning("   Try running: %matplotlib widget")
        else:
            logger.info(
                "✅ Desktop matplotlib detected, interactive selection should work"
            )
    except ImportError:
        logger.info("✅ Desktop environment, interactive selection should work")
