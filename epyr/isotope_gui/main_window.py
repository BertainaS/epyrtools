"""
Main window and application coordination for the isotope GUI.

This module contains the main application window class that coordinates
all the GUI components and handles user interactions.
"""

import platform
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Optional

from .gui_components import IsotopeDataTable, PeriodicTableGrid, ToolTip
from .gui_helpers import GUIConfig
from .isotope_data import IsotopeDataManager


class IsotopesGUI:
    """Main GUI application to display nuclear isotope data."""

    def __init__(self, root):
        """Initialize the Isotopes GUI.

        Args:
            root: The root Tkinter window (tk.Tk).

        Raises:
            FileNotFoundError: If the isotope data file cannot be found.
            Exception: If there's an error reading or processing the data file.
        """
        self.root = root
        self.root.title("Nuclear Isotopes")

        # Configuration
        self.config = GUIConfig()

        # Data manager
        self.data_manager = IsotopeDataManager()

        # GUI state
        self.current_element = ""  # Selected element symbol ('': all)
        self.field_var = None
        self.unstable_var = None
        self.nonmagnetic_var = None

        # GUI components
        self.periodic_table = None
        self.data_table = None
        self.field_entry = None

        # Load data
        try:
            self.full_data = self.data_manager.load_data()
            self.table_data = self.data_manager.get_table_data()
        except Exception as e:
            messagebox.showerror("Data Loading Error", str(e))
            raise

        # Setup window
        self._setup_window()

        # Create GUI components
        self._create_periodic_table()
        self._create_data_table()
        self._create_controls()

        # Initialize display
        self._update_table()

    def _setup_window(self):
        """Setup window geometry and properties."""
        window_width, window_height = self.config.calculate_window_size()

        # Center window on screen
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        x_pos = max(0, (screen_width - window_width) // 2)
        y_pos = max(0, (screen_height - window_height) // 2)

        self.root.geometry(f"{window_width}x{window_height}+{x_pos}+{y_pos}")
        self.root.resizable(False, False)

        # Calculate layout positions
        periodic_table_height = self.config.BORDER + 7 * self.config.y_spacing
        lan_act_block_height = (
            self.config.BORDER + 2 * self.config.y_spacing + self.config.CLASS_SPACING
        )

        self.table_y_start = (
            periodic_table_height + lan_act_block_height + self.config.LABEL_HEIGHT / 2
        )
        self.controls_y_start = self.table_y_start + self.config.TABLE_HEIGHT_PIXELS

    def _create_periodic_table(self):
        """Create the periodic table element buttons."""
        element_data = self.data_manager.get_element_data()

        self.periodic_table = PeriodicTableGrid(
            self.root, self.config, element_data, self._element_button_pushed
        )

        table_frame = self.periodic_table.create_table()
        table_frame.place(x=0, y=0)

    def _create_data_table(self):
        """Create the isotope data table."""
        # Create frame for table
        table_frame = tk.Frame(self.root)
        table_frame.place(
            x=self.config.BORDER,
            y=int(self.table_y_start),
            width=800,  # Adjust as needed
            height=self.config.TABLE_HEIGHT_PIXELS,
        )

        self.data_table = IsotopeDataTable(table_frame, self.config)
        self.tree = self.data_table.create_table()

        # Enable column sorting
        for col in self.tree["columns"]:
            self.tree.heading(col, command=lambda c=col: self._sort_table(c, False))

    def _create_controls(self):
        """Create control widgets (checkboxes, field entry, band buttons)."""
        controls_frame = ttk.Frame(self.root)
        window_width, _ = self.config.calculate_window_size()
        controls_frame.place(
            x=self.config.BORDER,
            y=int(self.controls_y_start),
            width=window_width - 2 * self.config.BORDER,
            height=self.config.BOTTOM_HEIGHT + self.config.BORDER,
        )

        control_pady = self.config.BORDER // 2

        # Checkboxes
        self.unstable_var = tk.IntVar(value=0)
        unstable_check = ttk.Checkbutton(
            controls_frame,
            text="Show unstable isotopes (*)",
            variable=self.unstable_var,
            command=self._update_table,
        )
        unstable_check.pack(side=tk.LEFT, padx=5, pady=control_pady)

        self.nonmagnetic_var = tk.IntVar(value=1)
        nonmagnetic_check = ttk.Checkbutton(
            controls_frame,
            text="Show non-magnetic isotopes (Spin=0)",
            variable=self.nonmagnetic_var,
            command=self._update_table,
        )
        nonmagnetic_check.pack(side=tk.LEFT, padx=5, pady=control_pady)

        # Band buttons
        self._create_band_buttons(controls_frame, control_pady)

        # Field entry
        self._create_field_entry(controls_frame, control_pady)

    def _create_band_buttons(self, parent, pady):
        """Create EPR band preset buttons."""
        band_button_frame = ttk.Frame(parent)
        band_button_frame.pack(side=tk.RIGHT, padx=10, pady=pady)

        button_width = 3
        style = ttk.Style()
        style.configure("Band.TButton", padding=(5, 2))

        # X-band button
        x_button = ttk.Button(
            band_button_frame,
            text="X",
            width=button_width,
            command=self._set_field_X,
            style="Band.TButton",
        )
        x_button.pack(side=tk.LEFT, padx=2)
        ToolTip(x_button, "Set field to X-band (340 mT)")

        # Q-band button
        q_button = ttk.Button(
            band_button_frame,
            text="Q",
            width=button_width,
            command=self._set_field_Q,
            style="Band.TButton",
        )
        q_button.pack(side=tk.LEFT, padx=2)
        ToolTip(q_button, "Set field to Q-band (1200 mT)")

        # W-band button
        w_button = ttk.Button(
            band_button_frame,
            text="W",
            width=button_width,
            command=self._set_field_W,
            style="Band.TButton",
        )
        w_button.pack(side=tk.LEFT, padx=2)
        ToolTip(w_button, "Set field to W-band (3400 mT)")

    def _create_field_entry(self, parent, pady):
        """Create magnetic field entry widget."""
        self.field_var = tk.DoubleVar(value=self.config.DEFAULT_FIELD)
        self.field_entry = ttk.Entry(
            parent,
            textvariable=self.field_var,
            width=10,
            justify="right",
        )
        self.field_entry.bind("<Return>", lambda event: self._update_table())
        self.field_entry.bind("<FocusOut>", lambda event: self._update_table())
        self.field_entry.pack(side=tk.RIGHT, padx=(0, 5), pady=pady)
        ToolTip(self.field_entry, "Enter magnetic field strength and press Enter")

        field_label = ttk.Label(parent, text="Field (mT):")
        field_label.pack(side=tk.RIGHT, padx=(10, 2), pady=pady)

    # --- Event Handlers ---

    def _element_button_pushed(self, element_symbol: str):
        """Handle element button click."""
        if self.current_element == element_symbol:
            # Clicking same element deselects it
            self.current_element = ""
        else:
            self.current_element = element_symbol
        self._update_table()

    def _set_field_X(self):
        """Set field to X-band value."""
        self.field_var.set(340.0)
        self._update_table()

    def _set_field_Q(self):
        """Set field to Q-band value."""
        self.field_var.set(1200.0)
        self._update_table()

    def _set_field_W(self):
        """Set field to W-band value."""
        self.field_var.set(3400.0)
        self._update_table()

    def _validate_field(self):
        """Validate and return the current field value."""
        try:
            field_value = self.field_var.get()
            if field_value <= 0:
                raise ValueError("Field must be positive")
            return field_value
        except (tk.TclError, ValueError):
            # Reset to default if invalid
            self.field_var.set(self.config.DEFAULT_FIELD)
            return self.config.DEFAULT_FIELD

    def _update_table(self):
        """Update the data table based on current filters and field."""
        field_value = self._validate_field()

        # Update NMR frequencies
        updated_data = self.data_manager.update_nmr_frequencies(
            self.table_data.copy(), field_value
        )

        # Apply filters
        filtered_data = updated_data.copy()

        # Filter by element if selected
        if self.current_element:
            element_mask = filtered_data["isotope"].str.contains(
                self.current_element, na=False
            )
            filtered_data = filtered_data[element_mask]

        # Filter unstable isotopes
        if not self.unstable_var.get():
            unstable_mask = ~filtered_data["isotope"].str.contains(r"\*", na=False)
            filtered_data = filtered_data[unstable_mask]

        # Filter non-magnetic isotopes
        if not self.nonmagnetic_var.get():
            magnetic_mask = filtered_data["spin"] != 0
            filtered_data = filtered_data[magnetic_mask]

        # Update table display
        self.data_table.update_data(filtered_data)

    def _sort_table(self, col_id: str, reverse: bool):
        """Sort the table by a column."""
        if self.data_table.data is None:
            return

        try:
            # Sort the underlying data
            if col_id in ["abundance", "spin", "gn", "gamma", "qm", "NMRfreq"]:
                # Numeric sort
                sorted_data = self.data_table.data.sort_values(
                    by=col_id, ascending=not reverse
                )
            else:
                # String sort
                sorted_data = self.data_table.data.sort_values(
                    by=col_id, ascending=not reverse
                )

            # Update display
            self.data_table.update_data(sorted_data)

            # Update column header for next click
            self.tree.heading(
                col_id, command=lambda: self._sort_table(col_id, not reverse)
            )

        except Exception as e:
            print(f"Error sorting table: {e}")


def run_gui():
    """Launch the isotope GUI application."""
    try:
        root = tk.Tk()

        # Set application icon if available
        if platform.system() == "Windows":
            try:
                root.iconbitmap("icon.ico")  # Add icon file if available
            except tk.TclError:
                pass  # Icon file not found, continue without it

        app = IsotopesGUI(root)
        root.mainloop()

    except Exception as e:
        error_message = f"Failed to initialize GUI: {str(e)}"
        print(f"Error: {error_message}")

        # Show error in message box if possible
        try:
            messagebox.showerror("Application Error", error_message)
        except:
            pass  # If messagebox fails, just print to console

        raise
