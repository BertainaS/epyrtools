#!/usr/bin/env python3
"""
A Tkinter GUI module for displaying nuclear isotope data from a file.

Provides the IsotopesGUI class and a run_gui() function to launch the application.
Expects isotope data in 'sub/isotopedata.txt' relative to this script's location
or the current working directory.
"""

import os
import platform
import traceback

import pandas as pd
import tkinter as tk
from tkinter import messagebox, ttk

from .logging_config import get_logger

logger = get_logger(__name__)

# --- Constants ---
PLANCK = 6.62607015e-34  # J⋅s (CODATA 2018)
NMAGN = 5.0507837461e-27  # J⋅T⁻¹ (Nuclear magneton, CODATA 2018)


# --- Helper Functions ---
def rgb_to_hex(rgb_tuple):
    """Converts an RGB tuple (values 0-255) to a Tkinter hex color string."""
    r = max(0, min(255, int(rgb_tuple[0])))
    g = max(0, min(255, int(rgb_tuple[1])))
    b = max(0, min(255, int(rgb_tuple[2])))
    return f"#{r:02x}{g:02x}{b:02x}"


def element_class(atomic_number):
    """
    Determines the period, group, and class of an element.

    Args:
        atomic_number (int): The atomic number (Z) of the element.

    Returns:
        tuple: (period, group, element_category)
               element_category: 0=main, 1=transition, 2=rare earth
    """
    period_limits = [0, 2, 10, 18, 36, 54, 86, 118]  # Max Z for each period start

    # Determine period
    period = 0
    for i in range(1, len(period_limits)):
        if atomic_number <= period_limits[i]:
            period = i
            break
    if period == 0 and atomic_number > 118:  # Handle elements beyond current limits
        period = 8

    group = atomic_number - period_limits[period - 1]
    element_category = 0  # Default: Main Group

    # Determine Group and Category based on period
    if period == 1:
        if group != 1:
            group = 18  # Group 1 (H) or Group 18 (He)
    elif period in [2, 3]:
        if group > 2:
            group += 10  # Shift p-block (Groups 3-8 -> 13-18)
    elif period in [4, 5]:
        if 3 <= group <= 12:
            element_category = 1  # Transition Metal
    elif period in [6, 7]:
        # Lanthanides (57-71) and Actinides (89-103)
        is_lanthanide = period == 6 and 57 <= atomic_number <= 71
        is_actinide = period == 7 and 89 <= atomic_number <= 103
        if is_lanthanide or is_actinide:
            element_category = 2  # Rare Earth
            # Group for positioning Ln/Ac block (visual column index 3-17)
            group = (
                (atomic_number - 57 + 3) if is_lanthanide else (atomic_number - 89 + 3)
            )
        elif group < 3:
            element_category = 0  # Alkali / Alkaline Earth
        else:
            # Post-Ln/Ac transition metals or main group p-block
            effective_group = atomic_number - period_limits[period - 1] - 14
            if 3 <= effective_group <= 12:
                element_category = 1  # Transition Metal
                group = effective_group
            elif effective_group > 12:
                element_category = 0  # Main Group (p-block)
                group = effective_group
            else:  # Should not happen if logic is correct
                element_category = 1  # Fallback classification
                group = effective_group

    return period, group, element_category


# --- Tooltip Class ---
class ToolTip:
    """Simple tooltip implementation for Tkinter widgets."""

    def __init__(self, widget, text="widget info"):
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        self.widget.bind("<Enter>", self.enter)  # Show on mouse enter
        self.widget.bind("<Leave>", self.leave)  # Hide on mouse leave
        self.widget.bind("<ButtonPress>", self.leave)  # Hide on click

    def enter(self, event=None):
        """Display the tooltip window."""
        x_rel, y_rel, _, _ = self.widget.bbox(
            "insert"
        )  # Get widget bounds relative to its parent
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


# --- Main GUI Class ---
class IsotopesGUI:
    """GUI application to display nuclear isotope data."""

    def __init__(self, root):
        """
        Initializes the Isotopes GUI.

        Args:
            root: The root Tkinter window (tk.Tk).

        Raises:
            FileNotFoundError: If the isotope data file cannot be found.
            Exception: If there's an error reading or processing the data file.
        """
        self.root = root
        self.root.title("EPyR Tools - Nuclear Isotopes Database")

        # --- GUI Settings ---
        self.default_field = 340.0  # mT
        self.button_font_size = 14
        self.element_width = 36
        self.element_height = self.element_width
        self.border = 10
        self.spacing = 5
        self.x_spacing = self.element_width + self.spacing
        self.y_spacing = self.element_height + self.spacing
        self.class_spacing = 5
        self.label_height = 15
        self.table_height_pixels = 200
        self.bottom_height = 30

        # --- Load Data ---
        # This will raise an exception if the file is not found or parsing fails
        self.full_data = self._read_isotope_data_file()
        # No need for explicit check/return here, exception handles failure

        # --- Prepare Data Structures ---
        self.table_display_data_cols = [
            "isotope",
            "abundance",
            "spin",
            "gn",
            "gamma",
            "qm",
        ]
        # Create a working copy for display data, add NMR frequency column
        self.table_data = self.full_data[self.table_display_data_cols].copy()
        self.table_data["NMRfreq"] = 0.0
        self.current_element = ""  # Selected element symbol ('': all)

        # --- Calculate Window Size ---
        window_width = (
            self.border + 18 * self.x_spacing + 2 * self.class_spacing + self.border
        )
        window_width += 20  # Buffer for controls fit

        periodic_table_height = self.border + 7 * self.y_spacing
        lan_act_block_height = self.border + 2 * self.y_spacing + self.class_spacing
        gap_above_table = self.label_height / 2
        table_section_height = self.table_height_pixels
        controls_section_height = self.bottom_height + self.border * 2

        window_height = (
            periodic_table_height
            + lan_act_block_height
            + gap_above_table
            + table_section_height
            + controls_section_height
        )

        # --- Set Window Properties ---
        # Store calculated dimensions
        self.calculated_window_width = int(window_width)
        self.calculated_window_height = int(window_height * 1.15)

        # Set minimum window size (80% of calculated size)
        min_width = int(window_width * 0.8)
        min_height = int(window_height * 0.8)
        self.root.minsize(min_width, min_height)

        # Make window resizable by user
        self.root.resizable(True, True)

        # Center window on screen
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        x_pos = max(0, (screen_width - window_width) // 2)
        y_pos = max(0, (screen_height - window_height) // 2)

        # --- Create Widgets First ---
        self._create_periodic_table()
        self._create_table()
        self._create_controls()

        # --- Initialize Table Display ---
        self._update_table()

        # --- Set Final Geometry After Widget Creation ---
        # Force proper window sizing after all widgets are created
        self.root.update_idletasks()

        # Schedule geometry setting after the GUI is fully realized
        def set_proper_geometry():
            self.root.geometry(
                f"{self.calculated_window_width}x{self.calculated_window_height}+{int(x_pos)}+{int(y_pos)}"
            )

        # Set geometry on next event loop iteration
        self.root.after(1, set_proper_geometry)

    def _read_isotope_data_file(self):
        """
        Reads isotope data from 'sub/isotopedata.txt'.

        Searches in the current working directory and the script's directory.

        Returns:
            pandas.DataFrame: The loaded and processed isotope data.

        Raises:
            FileNotFoundError: If 'sub/isotopedata.txt' cannot be found.
            Exception: If there is an error reading or processing the file.
        """
        script_dir = os.path.dirname(os.path.abspath(__file__))
        possible_paths = [
            os.path.join(os.getcwd(), "sub", "isotopedata.txt"),
            os.path.join(script_dir, "sub", "isotopedata.txt"),
            # You could add more search paths here if needed
        ]
        data_file = None
        found_path = None
        for path in possible_paths:
            if os.path.exists(path):
                data_file = path
                found_path = path  # Keep track of the found path for error message
                break

        if data_file is None:
            error_msg = (
                "Could not find 'sub/isotopedata.txt'. "
                f"Checked paths:\n- {possible_paths[0]}\n- {possible_paths[1]}"
                # Add others if you expanded possible_paths
            )
            logger.error(f"Error: {error_msg}")
            raise FileNotFoundError(error_msg)

        logger.info(f"Loading data from: {found_path}")  # Info message

        try:
            # Column definitions for isotopedata.txt:
            # Z: atomic number (number of protons)
            # N: mass number A (total nucleons, NOT neutron count!)
            # radioactive: * for radioactive, - for stable
            # element: chemical symbol (H, C, N, etc.)
            # name: element name
            # spin: nuclear spin quantum number
            # gn: nuclear g-factor
            # abundance: natural abundance (%)
            # qm: electric quadrupole moment (barn)
            col_names = [
                "Z",
                "N",
                "radioactive",
                "element",
                "name",
                "spin",
                "gn",
                "abundance",
                "qm",
            ]
            # Use sep='\s+' for whitespace delimiter
            data = pd.read_csv(
                data_file,
                comment="%",
                sep=r"\s+",  # Regex for one or more whitespace
                names=col_names,
                na_values=["-"],
                skipinitialspace=True,
            )

            # Convert types after reading, coercing errors
            data["Z"] = pd.to_numeric(data["Z"], errors="coerce")
            data["N"] = pd.to_numeric(data["N"], errors="coerce")
            numeric_cols = ["spin", "gn", "abundance", "qm"]
            for col in numeric_cols:
                if col in data.columns:
                    data[col] = pd.to_numeric(data[col], errors="coerce")

            # Remove rows missing critical Z or N
            data.dropna(subset=["Z", "N"], inplace=True)

            # Convert Z/N to integer *after* cleaning NaNs
            data["Z"] = data["Z"].astype(int)
            data["N"] = data["N"].astype(int)

            # Ensure string columns are strings and fill NaNs
            str_cols = ["radioactive", "element", "name"]
            for col in str_cols:
                if col in data.columns:
                    data[col] = data[col].astype(str).fillna("")

            # Calculate gamma/(2pi) [MHz/T]
            if "gn" in data.columns and pd.api.types.is_numeric_dtype(data["gn"]):
                data["gamma"] = data["gn"].apply(
                    lambda gn_val: (
                        (gn_val * NMAGN / PLANCK / 1e6) if pd.notna(gn_val) else pd.NA
                    )
                )
            else:
                data["gamma"] = pd.NA  # Ensure column exists even if gn doesn't

            # Assemble isotope symbols (e.g., 1H, 14C*, 235U)
            # Note: In the data file, column N is actually the mass number A (not neutron count)
            isotopes = []
            required_cols = ["N", "Z", "element", "radioactive"]
            if all(col in data.columns for col in required_cols):
                for _, row in data.iterrows():
                    # Check for NaN in Z/N specifically for this calculation
                    if pd.notna(row["N"]) and pd.notna(row["Z"]):
                        mass_number = int(
                            row["N"]
                        )  # N column is actually mass number A
                        iso_str = f"{mass_number}{row['element']}"
                        if row["radioactive"] == "*":
                            isotopes.append(iso_str + "*")
                        else:
                            isotopes.append(iso_str)
                    else:
                        # This case should ideally not happen after dropna above, but safer
                        isotopes.append(pd.NA)
            else:
                # Handle case where required columns might be missing in the input file
                logger.warning(
                    "Warning: Missing required columns (N, Z, element, radioactive) for isotope symbol generation."
                )
                data["isotope"] = pd.NA  # Assign NA if cannot generate

            data["isotope"] = isotopes
            # Convert to string *after* potential NA assignments
            data["isotope"] = data["isotope"].astype(str)

            # Fill remaining NaNs in numeric columns used for display AFTER calculations
            data["spin"] = data["spin"].fillna(-1.0)  # Marker for missing
            data["abundance"] = data["abundance"].fillna(0.0)
            data["qm"] = data["qm"].fillna(0.0)
            data["gn"] = data["gn"].fillna(0.0)
            data["gamma"] = data["gamma"].fillna(0.0)  # Fill calculated gamma NaNs

            # Data subset for placing elements on the periodic table
            self.element_layout_data = data.drop_duplicates(
                subset=["Z"], keep="first"
            ).copy()

            return data

        except Exception as e:
            logger.error(f"Error reading or processing data file '{data_file}': {e}")
            traceback.print_exc()
            # Re-raise the exception to be caught by the caller (__init__ -> run_gui)
            raise

    def _create_periodic_table(self):
        """Creates the element buttons in a periodic table layout."""
        x_offset = self.border
        y_start_periodic = self.border
        processed_Z = set()
        button_font = ("Arial", self.button_font_size - 2)
        is_macos = platform.system() == "Darwin"

        if not hasattr(self, "element_layout_data") or self.element_layout_data.empty:
            logger.warning(
                "Warning: No element layout data available to create periodic table."
            )
            return  # Cannot proceed without layout data

        for _, element_data in self.element_layout_data.iterrows():
            # Ensure Z is valid before proceeding
            if pd.isna(element_data["Z"]):
                continue
            ord_number = int(element_data["Z"])

            if ord_number in processed_Z:
                continue
            processed_Z.add(ord_number)

            element_symbol = element_data["element"]
            element_name = element_data["name"]
            period, group, element_category = element_class(ord_number)

            # --- Calculate element button position (x, y) ---
            x = 0
            y = 0
            if element_category == 2:  # Ln/Ac row positioning
                main_table_rows = 7
                ln_row_y = (
                    y_start_periodic
                    + main_table_rows * self.y_spacing
                    + self.class_spacing
                )
                ac_row_y = ln_row_y + self.y_spacing
                y = ln_row_y if period == 6 else ac_row_y
                # Position based on calculated 'group' (visual column 3-17)
                x = x_offset + (group - 1) * self.x_spacing
                x += self.class_spacing  # Align block start

            else:  # Main Group or Transition Metal positioning
                y = y_start_periodic + (period - 1) * self.y_spacing
                x = x_offset + (group - 1) * self.x_spacing
                if group > 2:
                    x += self.class_spacing
                if group > 12:
                    x += self.class_spacing

            # --- Determine background color based on element class ---
            bg_col_hex = "#D9D9D9"  # Default grey
            if element_category == 0:  # Main Group
                bg_col_rgb = [99, 154, 255] if group < 3 else [255, 207, 0]
                bg_col_hex = rgb_to_hex(bg_col_rgb)
            elif element_category == 1:  # Transition Metal
                bg_col_rgb = [255, 154, 156]
                bg_col_hex = rgb_to_hex(bg_col_rgb)
            elif element_category == 2:  # Lanthanide / Actinide
                bg_col_rgb = [0, 207, 49]
                bg_col_hex = rgb_to_hex(bg_col_rgb)

            # --- Check if element has isotopes with N > 0 listed ---
            try:
                # Ensure full_data is available
                if self.full_data is not None and not self.full_data.empty:
                    element_isotopes = self.full_data[self.full_data["Z"] == ord_number]
                    if not element_isotopes.empty:
                        # Check if *any* isotope for this Z has a valid N > 0
                        has_isotopes = (
                            element_isotopes["N"].notna() & (element_isotopes["N"] > 0)
                        ).any()
                    else:
                        has_isotopes = False
                else:
                    has_isotopes = False  # No data to check against
            except Exception as e:  # Catch potential errors during lookup
                logger.warning(
                    f"Warning: Error checking isotopes for Z={ord_number}: {e}"
                )
                has_isotopes = False

            final_color = (
                bg_col_hex if has_isotopes else "#E0E0E0"
            )  # Grey out if no N>0 isotopes listed

            # --- Create Button with Platform-Specific Coloring ---
            button_config = {
                "text": element_symbol,
                "font": button_font,
                "relief": "raised",
                "borderwidth": 1,
                "command": lambda sym=element_symbol: self._element_button_pushed(sym),
                "fg": "black",  # Ensure text is visible
                "state": (
                    "normal" if has_isotopes else "disabled"
                ),  # Disable button if no usable isotopes
            }

            # Use appropriate coloring method based on OS
            if is_macos:  # Use highlightbackground on macOS
                button_config["highlightbackground"] = final_color
                button_config["highlightthickness"] = 1  # Make border visible
                # Set active background for visual feedback on click
                button_config["activebackground"] = (
                    rgb_to_hex([c * 0.9 for c in bg_col_rgb])
                    if has_isotopes
                    else final_color
                )
            else:  # Use 'bg' on other platforms (Windows, Linux)
                button_config["bg"] = final_color
                # Set active background for visual feedback on click
                button_config["activebackground"] = (
                    rgb_to_hex([c * 0.9 for c in bg_col_rgb])
                    if has_isotopes
                    else final_color
                )

            button = tk.Button(self.root, **button_config)
            button.place(
                x=int(x),
                y=int(y),
                width=self.element_width,
                height=self.element_height,
            )
            if has_isotopes:  # Only add tooltip to active elements
                ToolTip(button, text=f" {element_name} ")

        # --- Add "all" button ---
        all_x = x_offset + 16 * self.x_spacing + 2 * self.class_spacing
        main_table_rows = 7
        ln_ac_rows = 2
        all_y = (
            y_start_periodic
            + main_table_rows * self.y_spacing
            + self.class_spacing
            + ln_ac_rows * self.y_spacing
            + self.spacing
        )
        all_button_width = self.element_width + self.x_spacing
        all_button_height = self.element_height + self.y_spacing

        # Define the light grey color for the 'all' button
        all_button_bg_rgb = [230, 230, 230]
        all_button_bg_hex = rgb_to_hex(all_button_bg_rgb)
        all_button_active_bg_hex = rgb_to_hex(
            [c * 0.9 for c in all_button_bg_rgb]
        )  # Slightly darker when active

        all_button_config = {
            "text": "all",
            "font": button_font,
            "relief": "raised",
            "borderwidth": 1,
            "command": lambda: self._element_button_pushed("all"),
            "fg": "black",  # Ensure text is visible
        }

        if is_macos:
            all_button_config["highlightbackground"] = all_button_bg_hex
            all_button_config["highlightthickness"] = 1
            all_button_config["activebackground"] = all_button_active_bg_hex
        else:
            all_button_config["bg"] = all_button_bg_hex
            all_button_config["activebackground"] = all_button_active_bg_hex

        all_button = tk.Button(self.root, **all_button_config)

        all_button.place(
            x=int(all_x),
            y=int(all_y),
            width=int(all_button_width),
            height=int(all_button_height),
        )
        ToolTip(all_button, text="Show all elements")

        # Y position below periodic table for the isotope table
        self.table_y_start = all_y + all_button_height + self.label_height / 2

    def _create_table(self):
        """Creates the isotope data table (Treeview)."""
        table_frame = ttk.Frame(self.root)
        table_frame_width = self.calculated_window_width - 2 * self.border
        table_frame.place(
            x=self.border,
            y=int(self.table_y_start),
            width=int(table_frame_width),
            height=self.table_height_pixels,
        )

        cols = (
            "isotope",
            "abundance",
            "spin",
            "gn",
            "gamma",
            "qm",
            "nmrfreq",
        )
        col_names = {
            "isotope": "Isotope",
            "abundance": "Abundance (%)",
            "spin": "Spin",
            "gn": "gn value",
            "gamma": "γ/2π (MHz/T)",
            "qm": "Q (barn)",
            "nmrfreq": "Frequency (MHz)",
        }
        col_widths = {  # Adjusted column widths
            "isotope": 80,
            "abundance": 110,
            "spin": 60,
            "gn": 95,
            "gamma": 115,
            "qm": 85,
            "nmrfreq": 120,
        }
        col_anchors = {  # Alignment within columns
            "isotope": "w",
            "abundance": "e",
            "spin": "e",
            "gn": "e",
            "gamma": "e",
            "qm": "e",
            "nmrfreq": "e",
        }

        self.table = ttk.Treeview(table_frame, columns=cols, show="headings")
        style = ttk.Style()
        # Configure Treeview style for row height and potentially heading font
        style.configure("Treeview", rowheight=25, font=("Arial", 10))  # Example font
        style.configure(
            "Treeview.Heading", font=("Arial", 10, "bold")
        )  # Example heading font

        for col_id in cols:
            self.table.heading(
                col_id,
                text=col_names[col_id],
                anchor="center",
                command=lambda c=col_id: self._sort_table(c, False),
            )
            self.table.column(
                col_id,
                width=col_widths[col_id],
                anchor=col_anchors[col_id],
                stretch=tk.NO,  # Fixed width columns
            )

        # Scrollbars
        vertical_scrollbar = ttk.Scrollbar(
            table_frame, orient="vertical", command=self.table.yview
        )
        horizontal_scrollbar = ttk.Scrollbar(
            table_frame, orient="horizontal", command=self.table.xview
        )
        self.table.configure(
            yscrollcommand=vertical_scrollbar.set,
            xscrollcommand=horizontal_scrollbar.set,
        )

        # Layout table and scrollbars using grid
        table_frame.grid_rowconfigure(0, weight=1)
        table_frame.grid_columnconfigure(0, weight=1)
        self.table.grid(row=0, column=0, sticky="nsew")
        vertical_scrollbar.grid(row=0, column=1, sticky="ns")
        horizontal_scrollbar.grid(row=1, column=0, sticky="ew")

        # Y position below table for controls
        self.controls_y_start = (
            self.table_y_start + self.table_height_pixels + self.border
        )

    def _create_controls(self):
        """Creates checkboxes, field entry, and band buttons."""
        controls_frame = ttk.Frame(self.root)  # Use a frame for better layout control
        controls_frame.place(
            x=self.border,
            y=int(self.controls_y_start),
            width=self.calculated_window_width - 2 * self.border,
            height=self.bottom_height + self.border,
        )  # Adjusted height

        current_x = 0  # Relative x within the frame
        control_height = 25  # Slightly increased height for better spacing/visuals
        control_pady = self.border // 2  # Vertical padding

        # Checkboxes
        self.unstable_var = tk.IntVar(value=0)
        unstable_check = ttk.Checkbutton(
            controls_frame,
            text="Show unstable isotopes (*)",
            variable=self.unstable_var,
            command=self._update_table,
        )
        # Use pack or grid within the frame for flexibility
        unstable_check.pack(side=tk.LEFT, padx=5, pady=control_pady)

        self.nonmagnetic_var = tk.IntVar(value=1)
        nonmagnetic_check = ttk.Checkbutton(
            controls_frame,
            text="Show non-magnetic isotopes (Spin=0)",
            variable=self.nonmagnetic_var,
            command=self._update_table,
        )
        nonmagnetic_check.pack(side=tk.LEFT, padx=5, pady=control_pady)

        # Spacer or flexible element might be needed here if using pack
        # For simplicity, using pack with careful padding

        # Band Buttons (Placed before field entry for visual grouping)
        band_button_frame = ttk.Frame(controls_frame)  # Sub-frame for band buttons
        band_button_frame.pack(side=tk.RIGHT, padx=10, pady=control_pady)

        button_width = 3  # ttk button width is in text units, not pixels
        style = ttk.Style()
        style.configure(
            "Band.TButton", padding=(5, 2)
        )  # Add some padding inside buttons

        x_button = ttk.Button(
            band_button_frame,
            text="X",
            width=button_width,
            command=self._set_field_X,
            style="Band.TButton",
        )
        x_button.pack(side=tk.LEFT, padx=2)
        ToolTip(x_button, text="Set field to X-band (340 mT)")

        q_button = ttk.Button(
            band_button_frame,
            text="Q",
            width=button_width,
            command=self._set_field_Q,
            style="Band.TButton",
        )
        q_button.pack(side=tk.LEFT, padx=2)
        ToolTip(q_button, text="Set field to Q-band (1200 mT)")

        w_button = ttk.Button(
            band_button_frame,
            text="W",
            width=button_width,
            command=self._set_field_W,
            style="Band.TButton",
        )
        w_button.pack(side=tk.LEFT, padx=2)
        ToolTip(w_button, text="Set field to W-band (3400 mT)")

        # Magnetic Field Input (Packed to the right, before band buttons)
        self.field_var = tk.DoubleVar(value=self.default_field)
        self.field_entry = ttk.Entry(
            controls_frame,
            textvariable=self.field_var,
            width=10,  # Adjusted width
            justify="right",
        )
        self.field_entry.bind("<Return>", lambda event: self._update_table())
        self.field_entry.bind("<FocusOut>", lambda event: self._update_table())
        style.map("TEntry", fieldbackground=[("!disabled", "white")])
        self.field_entry.pack(
            side=tk.RIGHT, padx=(0, 5), pady=control_pady
        )  # Pad right only
        ToolTip(self.field_entry, text="Enter magnetic field strength and press Enter")

        field_label = ttk.Label(controls_frame, text="Field (mT):")
        field_label.pack(
            side=tk.RIGHT, padx=(10, 2), pady=control_pady
        )  # Pad left only

    # --- Callback Functions ---
    def _element_button_pushed(self, element_symbol):
        """Handles clicks on element buttons."""
        new_element = "" if element_symbol == "all" else element_symbol
        if new_element != self.current_element:  # Update only if changed
            self.current_element = new_element
            self._update_table()

    def _set_field_X(self):
        """Sets magnetic field to X-band value (340 mT)."""
        if self.field_var.get() != 340.0:
            self.field_var.set(340.0)
            self._update_table()

    def _set_field_Q(self):
        """Sets magnetic field to Q-band value (1200 mT)."""
        if self.field_var.get() != 1200.0:
            self.field_var.set(1200.0)
            self._update_table()

    def _set_field_W(self):
        """Sets magnetic field to W-band value (3400 mT)."""
        if self.field_var.get() != 3400.0:
            self.field_var.set(3400.0)
            self._update_table()

    def _validate_field(self):
        """Validates the magnetic field entry, resets if invalid."""
        try:
            val = self.field_var.get()
            if val < 0:  # Optional: Disallow negative field
                messagebox.showwarning(
                    "Input Warning",
                    "Magnetic field cannot be negative. Resetting to default.",
                )
                self.field_var.set(self.default_field)
            return True  # Valid number obtained
        except tk.TclError:
            messagebox.showerror(
                "Input Error", "Invalid magnetic field value. Please enter a number."
            )
            self.field_var.set(self.default_field)  # Reset to default
            return False  # Invalid input

    def _update_table(self):
        """Filters data based on settings and updates the table display."""
        if not self._validate_field():
            # If field validation failed and reset the value,
            # _validate_field already showed an error.
            # We still need the B0 value for calculation, so get it again.
            B0 = self.field_var.get()  # Get the potentially reset value
        else:
            B0 = self.field_var.get()  # Get the validated value

        try:
            show_unstable = self.unstable_var.get()
            show_nonmagnetic = self.nonmagnetic_var.get()
        except tk.TclError as e:
            # Should not happen with IntVars, but good practice
            logger.error(f"Error getting checkbox values: {e}")
            return

        if self.full_data is None or self.full_data.empty:
            # Clear table if no data is loaded
            for item in self.table.get_children():
                self.table.delete(item)
            return

        # Filter data based on selected element
        if not self.current_element:
            # Make sure to copy to avoid modifying the original full_data
            filtered_df = self.full_data.copy()
        else:
            filtered_df = self.full_data[
                self.full_data["element"] == self.current_element
            ].copy()

        # Filter based on checkbox states
        if not show_unstable:
            # Check if 'radioactive' column exists before filtering
            if "radioactive" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["radioactive"] != "*"]
            else:
                logger.warning("Warning: 'radioactive' column not found for filtering.")

        if not show_nonmagnetic:
            # Check if 'spin' column exists before filtering
            if "spin" in filtered_df.columns:
                # Keep only isotopes with non-zero spin (using tolerance)
                # Also ensure spin is not the placeholder -1.0
                filtered_df = filtered_df[
                    (abs(filtered_df["spin"] - 0.0) > 1e-9)
                    & (filtered_df["spin"] >= 0.0)
                ]
            else:
                logger.warning("Warning: 'spin' column not found for filtering.")
        else:
            # If showing non-magnetic, still exclude placeholders
            if "spin" in filtered_df.columns:
                filtered_df = filtered_df[filtered_df["spin"] >= 0.0]

        # Exclude placeholder entries explicitly (redundant if handled above, but safe)
        # This step is crucial if the previous filters didn't catch the -1.0 placeholder
        # if "spin" in filtered_df.columns:
        #     filtered_df = filtered_df[filtered_df["spin"] >= 0.0]

        # Calculate NMR frequency (MHz) for the filtered isotopes
        # Ensure 'gamma' exists and B0 is valid
        if (
            "gamma" in filtered_df.columns
            and pd.api.types.is_numeric_dtype(filtered_df["gamma"])
            and B0 is not None
        ):
            # Freq[MHz] = gamma[MHz/T] * B0[mT] * 1e-3 [T/mT]
            # Ensure gamma is positive for frequency calculation if needed, or handle 0 gamma
            # Apply calculation only where gamma is not NA
            valid_gamma_mask = filtered_df["gamma"].notna()
            filtered_df.loc[valid_gamma_mask, "NMRfreq"] = (
                filtered_df.loc[valid_gamma_mask, "gamma"] * B0 * 1e-3
            )
            # Fill remaining NMRfreq with NaN or 0 where gamma was NA or calculation not applicable
            filtered_df["NMRfreq"] = filtered_df["NMRfreq"].fillna(
                pd.NA
            )  # Use pd.NA for consistency
        else:
            # Ensure NMRfreq column exists even if calculation fails or isn't applicable
            filtered_df["NMRfreq"] = pd.NA

        # --- Update Treeview Table Display ---
        # Clear previous entries efficiently
        self.table.delete(*self.table.get_children())

        if not filtered_df.empty:
            # Select only the columns intended for display in the correct order
            display_cols_ordered = [
                "isotope",
                "abundance",
                "spin",
                "gn",
                "gamma",
                "qm",
                "NMRfreq",  # Match the order defined in _create_table
            ]
            # Ensure all display columns exist in the filtered_df before selecting
            cols_to_display = [
                col for col in display_cols_ordered if col in filtered_df.columns
            ]
            display_df = filtered_df[cols_to_display]

            # Iterate and format values for display
            table_values = []
            for _, row in display_df.iterrows():
                # Format numeric values, handling potential NA/None
                abundance_str = (
                    f"{row['abundance']:.4f}"
                    if pd.notna(row["abundance"]) and "abundance" in cols_to_display
                    else ""
                )
                gn_str = (
                    f"{row['gn']:.4f}"
                    if pd.notna(row["gn"]) and "gn" in cols_to_display
                    else ""
                )
                gamma_str = (
                    f"{row['gamma']:.4f}"
                    if pd.notna(row["gamma"]) and "gamma" in cols_to_display
                    else ""
                )
                qm_val = row.get("qm", pd.NA)  # Safely get qm value
                qm_str = (
                    f"{qm_val:.4f}" if pd.notna(qm_val) and abs(qm_val) > 1e-9 else ""
                )

                nmrfreq_val = row.get("NMRfreq", pd.NA)  # Safely get NMRfreq
                nmrfreq_str = (
                    # Show frequency only if it's positive (or non-zero if preferred)
                    f"{nmrfreq_val:.4f}"
                    if pd.notna(nmrfreq_val) and nmrfreq_val > 1e-9
                    else ""
                )

                # Special formatting for spin (int, .5, or float), handling placeholder
                spin_val = row.get("spin", -1.0)  # Default to placeholder if missing
                spin_str = ""
                if (
                    pd.notna(spin_val) and spin_val >= 0.0
                ):  # Check it's valid and not placeholder
                    if spin_val == int(spin_val):
                        spin_str = str(int(spin_val))
                    # Check if it's a half-integer like 1.5, 2.5 etc.
                    elif (
                        abs(spin_val * 2 - int(spin_val * 2)) < 1e-9
                    ):  # Ends in .5 (within tolerance)
                        spin_str = f"{spin_val:.1f}"  # Format as x.5
                    else:  # Otherwise, format as float
                        spin_str = f"{spin_val:.4f}"
                # else: spin_str remains "" if spin is placeholder or NA

                # Assemble the list of values in the correct column order
                # Use .get() with default for safety if a column was missing from display_df
                values = [
                    row.get("isotope", ""),
                    abundance_str,
                    spin_str,
                    gn_str,
                    gamma_str,
                    qm_str,
                    nmrfreq_str,
                ]
                table_values.append(values)

            # Insert all rows at once (potentially faster for large datasets, though tkinter might not optimize this much)
            for values in table_values:
                self.table.insert("", tk.END, values=values)

    def _sort_table(self, col_id, reverse):
        """Sorts the table view by the selected column."""
        # Extract data with item IDs
        # Use get() to handle potential errors if item_id is somehow invalid
        data_list = []
        for item_id in self.table.get_children(""):
            try:
                value = self.table.set(item_id, col_id)
                data_list.append((value, item_id))
            except tk.TclError:
                logger.warning(
                    f"Warning: Could not get value for item {item_id}, column {col_id}. Skipping."
                )
                continue  # Skip item if value cannot be retrieved

        # Define columns to attempt numeric sort
        numeric_cols_for_sort = ["abundance", "spin", "gn", "gamma", "qm", "nmrfreq"]

        # Sorting logic
        try:
            if col_id in numeric_cols_for_sort:
                # Key function for robust numeric sort (handles errors, empty strings)
                def sort_key_numeric(item_tuple):
                    value_str = item_tuple[0]
                    try:
                        # Attempt to convert non-empty string to float
                        if value_str and isinstance(value_str, str):
                            return float(value_str)
                        elif isinstance(value_str, (int, float)):  # Already numeric
                            return float(value_str)
                        else:  # Empty string or other non-convertible type treated as lowest/highest
                            return -float("inf") if not reverse else float("inf")
                    except (ValueError, TypeError):
                        # Handle cases like "1/2" if they weren't formatted numerically, treat as lowest/highest
                        return -float("inf") if not reverse else float("inf")

                data_list.sort(key=sort_key_numeric, reverse=reverse)
            else:  # Default to case-insensitive string sort for non-numeric columns
                data_list.sort(
                    key=lambda t: str(t[0]).lower() if t[0] else "", reverse=reverse
                )

        except Exception as e:
            logger.error(f"Error during sorting column '{col_id}': {e}")
            traceback.print_exc()
            # Fallback to basic string sort if complex sort fails entirely
            try:
                data_list.sort(key=lambda t: str(t[0]), reverse=reverse)
            except Exception as fallback_e:
                logger.error(f"Fallback sort also failed: {fallback_e}")
                # If even basic sort fails, just leave the order as is.

        # Rearrange items in the Treeview
        for index, (_, item_id) in enumerate(data_list):
            try:
                self.table.move(item_id, "", index)
            except tk.TclError:
                logger.warning(f"Warning: Could not move item {item_id} during sort.")
                continue  # Skip if item cannot be moved

        # Update heading command for next sort direction
        # Ensure lambda captures the current state of 'reverse' correctly
        self.table.heading(
            col_id,
            text=self.table.heading(col_id)["text"],  # Keep original text
            command=lambda c=col_id, r=reverse: self._sort_table(c, not r),
        )


# --- Module Execution / Entry Point ---


def run_gui():
    """Creates the main Tkinter window and runs the IsotopesGUI application."""
    root = tk.Tk()
    try:
        # Pass the root window to the application class
        app = IsotopesGUI(root)
        # Start the Tkinter event loop
        root.mainloop()
    except FileNotFoundError as e:
        # Handle the specific case where the data file wasn't found during init
        messagebox.showerror(
            "Fatal Error - Data File Not Found",
            f"Could not initialize application.\n{e}\n\nPlease ensure the 'sub' directory containing 'isotopedata.txt' is accessible.",
        )
        root.destroy()  # Close the (likely empty) root window
    except Exception as e:
        # Catch any other unexpected errors during initialization
        messagebox.showerror(
            "Fatal Error - Initialization Failed",
            f"An unexpected error occurred during application startup:\n\n{e}\n\n{traceback.format_exc()}",
        )
        root.destroy()  # Close the root window


if __name__ == "__main__":
    # This block only runs when the script is executed directly
    # (e.g., python isotopes_gui.py)
    logger.info("Running Isotopes GUI as main script...")
    run_gui()
