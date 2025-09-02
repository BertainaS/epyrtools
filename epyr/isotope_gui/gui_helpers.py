"""
Helper functions and utilities for the isotope GUI application.

This module contains utility functions for color conversion, element classification,
and other helper functions used throughout the isotope GUI.
"""

# --- Constants ---
PLANCK = 6.62607015e-34  # J⋅s (CODATA 2018)
NMAGN = 5.0507837461e-27  # J⋅T⁻¹ (Nuclear magneton, CODATA 2018)


def rgb_to_hex(rgb_tuple):
    """Convert an RGB tuple (values 0-255) to a Tkinter hex color string.

    Args:
        rgb_tuple: Tuple of (r, g, b) values from 0-255

    Returns:
        str: Hex color string in format #RRGGBB
    """
    r = max(0, min(255, int(rgb_tuple[0])))
    g = max(0, min(255, int(rgb_tuple[1])))
    b = max(0, min(255, int(rgb_tuple[2])))
    return f"#{r:02x}{g:02x}{b:02x}"


def element_class(atomic_number):
    """
    Determine the period, group, and class of an element.

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

    if period == 0:
        return 0, 0, 0  # Unknown element

    # Determine group and element class
    group = 0
    element_category = 0  # 0=main, 1=transition, 2=rare earth

    if period == 1:  # H, He
        if atomic_number == 1:
            group = 1
        else:  # He
            group = 18
    elif period == 2:  # Li-Ne
        group = atomic_number - 2
    elif period == 3:  # Na-Ar
        group = atomic_number - 10
    elif period == 4:  # K-Kr
        if atomic_number <= 20:  # K, Ca
            group = atomic_number - 18
        elif atomic_number <= 30:  # Sc-Zn (transition metals)
            group = atomic_number - 18
            element_category = 1
        else:  # Ga-Kr
            group = atomic_number - 18
    elif period == 5:  # Rb-Xe
        if atomic_number <= 38:  # Rb, Sr
            group = atomic_number - 36
        elif atomic_number <= 48:  # Y-Cd (transition metals)
            group = atomic_number - 36
            element_category = 1
        else:  # In-Xe
            group = atomic_number - 36
    elif period == 6:  # Cs-Rn
        if atomic_number <= 56:  # Cs, Ba
            group = atomic_number - 54
        elif 57 <= atomic_number <= 71:  # La-Lu (lanthanides)
            group = 3  # All lanthanides are considered group 3
            element_category = 2
        elif atomic_number <= 80:  # Hf-Hg (transition metals)
            group = atomic_number - 54 - 14  # Subtract 14 for lanthanides
            element_category = 1
        else:  # Tl-Rn
            group = atomic_number - 54 - 14
    elif period == 7:  # Fr-Og
        if atomic_number <= 88:  # Fr, Ra
            group = atomic_number - 86
        elif 89 <= atomic_number <= 103:  # Ac-Lr (actinides)
            group = 3  # All actinides are considered group 3
            element_category = 2
        elif atomic_number <= 112:  # Rf-Cn (transition metals)
            group = atomic_number - 86 - 14  # Subtract 14 for actinides
            element_category = 1
        else:  # Nh-Og
            group = atomic_number - 86 - 14

    return period, group, element_category


class GUIConfig:
    """Configuration settings for the isotope GUI."""

    # Default settings
    DEFAULT_FIELD = 340.0  # mT
    BUTTON_FONT_SIZE = 14
    ELEMENT_WIDTH = 36
    ELEMENT_HEIGHT = 36
    BORDER = 10
    SPACING = 5
    CLASS_SPACING = 5
    LABEL_HEIGHT = 15
    TABLE_HEIGHT_PIXELS = 200
    BOTTOM_HEIGHT = 30

    # Calculated spacing
    @property
    def x_spacing(self):
        return self.ELEMENT_WIDTH + self.SPACING

    @property
    def y_spacing(self):
        return self.ELEMENT_HEIGHT + self.SPACING

    def calculate_window_size(self):
        """Calculate the required window dimensions."""
        window_width = (
            self.BORDER + 18 * self.x_spacing + 2 * self.CLASS_SPACING + self.BORDER
        )
        window_width += 20  # Buffer for controls fit

        periodic_table_height = self.BORDER + 7 * self.y_spacing
        lan_act_block_height = self.BORDER + 2 * self.y_spacing + self.CLASS_SPACING
        gap_above_table = self.LABEL_HEIGHT / 2
        table_section_height = self.TABLE_HEIGHT_PIXELS
        controls_section_height = self.BOTTOM_HEIGHT + self.BORDER * 2

        window_height = (
            periodic_table_height
            + lan_act_block_height
            + gap_above_table
            + table_section_height
            + controls_section_height
        )

        return int(window_width), int(window_height)
