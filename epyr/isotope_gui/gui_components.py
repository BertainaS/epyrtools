"""
Reusable GUI components for the isotope application.

This module contains reusable widgets and GUI components that can be
used throughout the isotope GUI application.
"""

import tkinter as tk
from tkinter import ttk
from typing import Any, Callable, Dict, Optional

import pandas as pd

from .gui_helpers import GUIConfig, element_class, rgb_to_hex


class ToolTip:
    """Simple tooltip implementation for Tkinter widgets."""

    def __init__(self, widget, text="widget info"):
        """Initialize tooltip for a widget.

        Args:
            widget: The Tkinter widget to attach the tooltip to
            text: The text to display in the tooltip
        """
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        self.widget.bind("<Enter>", self.enter)  # Show on mouse enter
        self.widget.bind("<Leave>", self.leave)  # Hide on mouse leave
        self.widget.bind("<ButtonPress>", self.leave)  # Hide on click

    def enter(self, event=None):
        """Display the tooltip window."""
        x_rel, y_rel, _, _ = self.widget.bbox("insert")  # Get widget bounds
        if x_rel is None:  # Handle cases where bbox might not be available yet
            x_rel, y_rel = 0, 0

        # Calculate absolute screen coordinates for the tooltip popup
        x_abs = self.widget.winfo_rootx() + x_rel + 25  # Offset from mouse
        y_abs = self.widget.winfo_rooty() + y_rel + 20

        # Create a toplevel window for the tooltip
        self.tooltip_window = tk.Toplevel(self.widget)
        # Remove window decorations (border, title bar)
        self.tooltip_window.wm_overrideredirect(True)
        self.tooltip_window.wm_geometry(f"+{int(x_abs)}+{int(y_abs)}")

        # Add label with tooltip text
        label = tk.Label(
            self.tooltip_window,
            text=self.text,
            justify="left",
            background="#ffffe0",
            relief="solid",
            borderwidth=1,
            font=("tahoma", "8", "normal"),
        )
        label.pack(ipadx=1)

    def leave(self, event=None):
        """Destroy the tooltip window."""
        if self.tooltip_window:
            self.tooltip_window.destroy()
        self.tooltip_window = None


class PeriodicTableGrid:
    """Creates and manages the periodic table element buttons."""

    def __init__(
        self,
        parent,
        config: GUIConfig,
        element_data: pd.DataFrame,
        on_element_click: Callable[[str], None],
    ):
        """Initialize the periodic table grid.

        Args:
            parent: Parent Tkinter widget
            config: GUI configuration object
            element_data: DataFrame with element data
            on_element_click: Callback function when element is clicked
        """
        self.parent = parent
        self.config = config
        self.element_data = element_data
        self.on_element_click = on_element_click
        self.element_buttons = {}

    def create_table(self) -> tk.Frame:
        """Create the periodic table frame with element buttons.

        Returns:
            tk.Frame: Frame containing the periodic table
        """
        # Create main frame for periodic table
        table_frame = tk.Frame(self.parent)

        # Create element buttons
        self._create_element_buttons(table_frame)

        return table_frame

    def _create_element_buttons(self, parent_frame: tk.Frame):
        """Create element buttons arranged in periodic table layout."""
        button_positions = self._get_element_positions()

        for _, row in self.element_data.iterrows():
            z = int(row["Z"])
            element = row["element"]

            if z in button_positions:
                period, group = button_positions[z]
                self._create_element_button(parent_frame, element, z, period, group)

    def _create_element_button(
        self, parent: tk.Frame, element: str, z: int, period: int, group: int
    ):
        """Create a single element button.

        Args:
            parent: Parent frame
            element: Element symbol
            z: Atomic number
            period: Period number
            group: Group number
        """
        # Calculate position
        x = self.config.BORDER + (group - 1) * self.config.x_spacing
        y = self.config.BORDER + (period - 1) * self.config.y_spacing

        # Adjust for lanthanides and actinides
        if 57 <= z <= 71:  # Lanthanides
            x = self.config.BORDER + (z - 54) * self.config.x_spacing
            y = (
                self.config.BORDER
                + 7 * self.config.y_spacing
                + self.config.BORDER
                + self.config.CLASS_SPACING
            )
        elif 89 <= z <= 103:  # Actinides
            x = self.config.BORDER + (z - 86) * self.config.x_spacing
            y = (
                self.config.BORDER
                + 7 * self.config.y_spacing
                + self.config.BORDER
                + self.config.y_spacing
                + self.config.CLASS_SPACING
            )

        # Get element classification for coloring
        _, _, element_category = element_class(z)

        # Color scheme
        colors = {
            0: (240, 240, 240),  # Main group - light gray
            1: (200, 230, 255),  # Transition metals - light blue
            2: (255, 200, 200),  # Rare earths - light red
        }

        bg_color = rgb_to_hex(colors.get(element_category, (240, 240, 240)))

        # Create button
        button = tk.Button(
            parent,
            text=element,
            font=("Arial", self.config.BUTTON_FONT_SIZE, "bold"),
            width=2,
            height=1,
            bg=bg_color,
            command=lambda elem=element: self.on_element_click(elem),
        )

        button.place(
            x=x, y=y, width=self.config.ELEMENT_WIDTH, height=self.config.ELEMENT_HEIGHT
        )

        # Add tooltip
        tooltip_text = f"{element} (Z={z})"
        ToolTip(button, tooltip_text)

        self.element_buttons[element] = button

    def _get_element_positions(self) -> Dict[int, tuple]:
        """Get the (period, group) positions for elements in standard periodic table.

        Returns:
            Dict mapping atomic number to (period, group) tuple
        """
        positions = {}

        # Period 1: H, He
        positions.update({1: (1, 1), 2: (1, 18)})

        # Period 2: Li-Ne
        for z in range(3, 11):
            positions[z] = (2, z - 2)

        # Period 3: Na-Ar
        for z in range(11, 19):
            positions[z] = (3, z - 10)

        # Period 4: K-Kr
        for z in range(19, 37):
            positions[z] = (4, z - 18)

        # Period 5: Rb-Xe
        for z in range(37, 55):
            positions[z] = (5, z - 36)

        # Period 6: Cs-Ba, Hf-Rn (skip lanthanides for main table)
        for z in range(55, 57):  # Cs, Ba
            positions[z] = (6, z - 54)
        for z in range(72, 87):  # Hf-Rn
            positions[z] = (6, z - 54 - 14)  # Subtract 14 for lanthanides

        # Period 7: Fr-Ra, Rf-Og (skip actinides for main table)
        for z in range(87, 89):  # Fr, Ra
            positions[z] = (7, z - 86)
        for z in range(104, 119):  # Rf-Og
            positions[z] = (7, z - 86 - 14)  # Subtract 14 for actinides

        return positions


class IsotopeDataTable:
    """Manages the isotope data table widget."""

    def __init__(self, parent, config: GUIConfig):
        """Initialize the data table.

        Args:
            parent: Parent Tkinter widget
            config: GUI configuration object
        """
        self.parent = parent
        self.config = config
        self.tree = None
        self.data = None

    def create_table(self) -> ttk.Treeview:
        """Create the data table widget.

        Returns:
            ttk.Treeview: The created table widget
        """
        # Create frame for table
        table_frame = tk.Frame(self.parent)
        table_frame.pack(fill=tk.BOTH, expand=True)

        # Define columns
        columns = ("isotope", "abundance", "spin", "gn", "gamma", "qm", "NMRfreq")

        # Create treeview widget
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings")

        # Configure column headings and widths
        column_configs = {
            "isotope": ("Isotope", 80),
            "abundance": ("Abundance [%]", 100),
            "spin": ("Spin", 60),
            "gn": ("gn", 80),
            "gamma": ("γ/(2π) [MHz/T]", 100),
            "qm": ("Q [fm²]", 80),
            "NMRfreq": ("NMR [MHz]", 90),
        }

        for col, (heading, width) in column_configs.items():
            self.tree.heading(col, text=heading)
            self.tree.column(col, width=width, minwidth=width)

        # Add scrollbars
        v_scrollbar = ttk.Scrollbar(
            table_frame, orient=tk.VERTICAL, command=self.tree.yview
        )
        h_scrollbar = ttk.Scrollbar(
            table_frame, orient=tk.HORIZONTAL, command=self.tree.xview
        )
        self.tree.configure(
            yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set
        )

        # Pack widgets
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

        return self.tree

    def update_data(self, data: pd.DataFrame, element_filter: str = ""):
        """Update the table with new data.

        Args:
            data: DataFrame containing isotope data
            element_filter: Element symbol to filter by (empty string shows all)
        """
        # Clear existing items
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Filter data if element is specified
        if element_filter and element_filter != "":
            # Extract element from isotope column
            filtered_data = data[data["isotope"].str.contains(element_filter, na=False)]
        else:
            filtered_data = data

        # Insert filtered data
        for _, row in filtered_data.iterrows():
            values = []
            for col in self.tree["columns"]:
                if col in row:
                    value = row[col]
                    # Format numeric values
                    if isinstance(value, float):
                        if col == "NMRfreq":
                            values.append(f"{value:.3f}")
                        elif col in ["abundance", "gamma"]:
                            values.append(f"{value:.2f}")
                        else:
                            values.append(f"{value:.4f}")
                    else:
                        values.append(str(value))
                else:
                    values.append("")

            self.tree.insert("", tk.END, values=values)

        self.data = filtered_data
