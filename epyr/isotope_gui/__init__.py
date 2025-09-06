"""
Isotope GUI Package
==================

Interactive GUI for exploring nuclear isotope data.

This package provides a modular Tkinter-based GUI application for
browsing and searching nuclear isotope data from the database.

Main Components:
- main_window: Main application window and controller
- gui_components: Individual UI components and widgets
- gui_helpers: Helper functions for GUI operations
- isotope_data: Data management and isotope database handling

Usage:
    from epyr.isotope_gui import run_gui
    run_gui()
"""

from .main_window import run_gui

__all__ = ["run_gui"]