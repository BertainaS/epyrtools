# eprload.py

import os
import sys
import tkinter as tk
import warnings
from pathlib import Path
from tkinter import filedialog
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

# Make sure the 'sub' directory is importable
# This assumes 'eprload.py' is in the parent directory of 'sub'
try:
    from .sub import loadBES3T, loadESP
except ImportError:
    # This works if eprload.py is in the same directory as the 'sub' folder
    # or if 'sub' is in the Python path.
    # This is often the case when running a script directly.
    try:
        from sub import loadBES3T, loadESP
    except ImportError:
        # If running from a different structure, adjust sys.path if needed
        # Or ensure the package containing 'sub' is installed
        print(
            "Error: Could not import loading modules from 'sub' directory.",
            file=sys.stderr,
        )
        print(
            "Ensure eprload.py is in the parent directory of 'sub', or adjust Python's path.",
            file=sys.stderr,
        )
        sys.exit(1)


def eprload(file_name=None, scaling="", plot_if_possible=True, save_if_possible=False):
    """
    Load experimental EPR data from Bruker BES3T or ESP formats.

    Args:
        file_name (str or Path, optional): Path to the data file (.dta, .dsc, .spc, .par) or a directory.
            If None or a directory, a file browser is shown. Defaults to None (opens browser in cwd).
        scaling (str, optional): String of characters specifying scaling operations (only for Bruker files).
            Each character enables a scaling operation:
                'n': Divide by number of scans (AVGS/JSD).
                'P': Divide by sqrt of MW power in mW (MWPW/MP).
                'G': Divide by receiver gain (RCAG/RRG).
                'T': Multiply by temperature in Kelvin (STMP/TE).
                'c': Divide by conversion/sampling time in ms (SPTP/RCT).
            Defaults to "" (no scaling).
        plot_if_possible (bool, optional): If True and data is loaded successfully, a simple plot is generated using matplotlib.
            Defaults to True.

    Returns:
        tuple:
            - x (np.ndarray or list of np.ndarray): Abscissa data (or list for 2D).
            - y (np.ndarray): Ordinate data.
            - pars (dict): Dictionary of parameters from the descriptor/parameter file.
            - file_path (str): The full path of the loaded file.
            - On failure (e.g., user cancel, file error): (None, None, None, None)

    Raises:
        FileNotFoundError: If the specified file or directory does not exist.
        ValueError: If the file format is unsupported, scaling is invalid, or parameter inconsistencies are found.
        IOError: If there are problems reading files.
    """
    x, y, pars, loaded_file_path = None, None, None, None  # Initialize outputs

    if file_name is None:
        file_name = Path.cwd()  # Use current directory if none provided
    else:
        file_name = Path(file_name)  # Convert string to Path object

    # --- File/Directory Handling ---
    if file_name.is_dir():
        root = tk.Tk()
        root.withdraw()  # Hide the main tkinter window
        ui_file_path = filedialog.askopenfilename(
            title="Load EPR data file...",
            initialdir=str(file_name),
            filetypes=[
                ("Bruker BES3T", "*.DTA *.dta *.DSC *.dsc"),
                ("Bruker ESP/WinEPR", "*.spc *.SPC *.par *.PAR"),
                ("All files", "*.*"),
            ],
        )
        root.destroy()  # Close the hidden window

        if not ui_file_path:  # User cancelled
            warnings.warn("File selection cancelled.")
            return None, None, None, None  # Return None tuple on cancel
        else:
            file_path = Path(ui_file_path)
    elif file_name.is_file():
        file_path = file_name
    else:
        raise FileNotFoundError(f"The file or directory '{file_name}' does not exist!")

    # --- Determine File Format ---
    full_base_name = file_path.with_suffix("")
    file_extension = file_path.suffix

    # Handle case where extension might be missing but file exists (e.g., JEOL)
    # For this version, we only support Bruker, so require extension.
    if not file_extension:
        # Check for common Bruker extensions if none provided
        found_ext = None
        for ext in [".dta", ".DTA", ".dsc", ".DSC", ".spc", ".SPC", ".par", ".PAR"]:
            potential_file = full_base_name.with_suffix(ext)
            if potential_file.is_file():
                found_ext = ext
                file_path = (
                    potential_file  # Update file_path to include found extension
                )
                break
        if found_ext:
            file_extension = found_ext
            warnings.warn(
                f"No extension given, assuming '{found_ext}' based on existing file."
            )
        else:
            raise ValueError(
                f"File '{full_base_name}' lacks a recognized extension (.dta, .dsc, .spc, .par)."
            )

    # Determine format based on extension (case-insensitive)
    ext_upper = file_extension.upper()
    file_format = None
    if ext_upper in [".DTA", ".DSC"]:
        file_format = "BrukerBES3T"
    elif ext_upper in [".PAR", ".SPC"]:
        file_format = "BrukerESP"
    # Add other formats here if needed in the future
    # elif ...
    else:
        raise ValueError(
            f"Unsupported file extension '{file_extension}'. Only Bruker formats (.dta, .dsc, .spc, .par) supported."
        )

    # --- Validate Scaling ---
    if scaling:
        valid_scaling_chars = "nPGTc"
        invalid_chars = set(scaling) - set(valid_scaling_chars)
        if invalid_chars:
            raise ValueError(
                f"Scaling string contains invalid characters: {invalid_chars}. Allowed: '{valid_scaling_chars}'."
            )
        # Note: Scaling applicability is checked within the specific load functions

    # --- Load Data using appropriate module ---
    loaded_file_path = str(file_path.resolve())  # Store full path string

    try:
        if file_format == "BrukerBES3T":
            y, x, pars = loadBES3T.load(full_base_name, file_extension, scaling)
        elif file_format == "BrukerESP":
            y, x, pars = loadESP.load(full_base_name, file_extension, scaling)
        # No 'else' needed due to prior format check

    except (FileNotFoundError, ValueError, IOError, NotImplementedError) as e:
        # Catch errors from loading functions
        warnings.warn(f"Failed to load data from {file_path}: {e}")
        return None, None, None, None
    except Exception as e:
        # Catch unexpected errors
        warnings.warn(
            f"An unexpected error occurred loading {file_path}: {type(e).__name__}: {e}"
        )
        import traceback

        traceback.print_exc()
        return None, None, None, None

    # --- Plotting (Optional) ---
    if plot_if_possible and y is not None and x is not None:
        _plot_data(x, y, loaded_file_path, pars, save_if_possible=save_if_possible)

    # --- Return Results ---
    # Always return the full tuple, even if some elements are None (e.g., on error)
    return x, y, pars, loaded_file_path


def _plot_data(
    x: Union[np.ndarray, List[np.ndarray], None],
    y: Optional[np.ndarray],
    file_name: str,
    params: Optional[Dict[str, Any]] = None,
    save_if_possible: bool = False,
) -> None:
    """Helper function to plot the loaded EPR data using Matplotlib."""
    if y is None or y.size == 0:
        print("No data available to plot.")
        return

    fig, ax = plt.subplots()
    plot_title = Path(file_name).name
    ax.set_title(plot_title, fontsize=10)  # Use smaller font for title

    # --- 1D Data Plotting ---
    if y.ndim == 1:
        absc = x
        if absc is None or not isinstance(absc, np.ndarray) or absc.shape != y.shape:
            warnings.warn(
                "Abscissa data (x) missing or incompatible shape. Using index for plotting."
            )
            absc = np.arange(y.size)
            x_label = "Index"
            x_unit = "points"
        else:
            x_label = params.get("XAXIS_NAME")
            x_unit = params.get("XAXIS_UNIT") if params else "a.u."
            if isinstance(x_unit, list):  # Use first unit if list provided for 1D
                x_unit = x_unit[0]
            if x_unit:
                x_label += f" ({x_unit})"

        if np.isrealobj(y):
            ax.plot(absc, y, label="data")
        else:
            ax.plot(absc, np.real(y), label="real")
            ax.plot(absc, np.imag(y), label="imag")  # , linestyle='--')
            ax.legend()

        ax.set_xlabel(x_label)
        ax.set_ylabel("Intensity (a.u.)")
        ax.grid(True, linestyle=":", alpha=0.6)
        ax.ticklabel_format(
            style="sci", axis="y", scilimits=(-3, 4)
        )  # Use scientific notation if needed

    # --- 2D Data Plotting ---
    elif y.ndim == 2:
        # y shape is typically (ny, nx) after loading
        ny, nx = y.shape
        aspect_ratio = nx / ny if ny > 0 else 1.0

        # Determine x and y axes for the plot
        x_coords = np.arange(nx)  # Default x: index
        y_coords = np.arange(ny)  # Default y: index
        x_label = f"Index ({nx} points)"
        y_label = f"Index ({ny} points)"

        x_units_list = (
            [params.get("XAXIS_UNIT"), params.get("YAXIS_UNIT")] if params else None
        )

        if isinstance(x, list) and len(x) >= 2:
            # Assume x = [x_axis, y_axis, ...]
            x_axis, y_axis = x[0], x[1]
            if isinstance(x_axis, np.ndarray) and x_axis.size == nx:
                x_coords = x_axis
                x_unit = params.get(
                    "XAXIS_UNIT"
                )  # if isinstance(x_units_list, list) and len(x_units_list) > 0 else 'a.u.'
                x_label = params.get("XAXIS_NAME") + f"({x_unit})"
            if isinstance(y_axis, np.ndarray) and y_axis.size == ny:
                y_coords = y_axis
                y_unit = (
                    params.get("YAXIS_UNIT")
                    if isinstance(x_units_list, list) and len(x_units_list) > 1
                    else "a.u."
                )
                y_label = params.get("YAXIS_NAME") + f"({y_unit})"

        elif isinstance(x, np.ndarray) and x.size == nx:
            # Only x-axis provided for 2D plot
            x_coords = x
            x_unit = x_units_list if isinstance(x_units_list, str) else "a.u."
            x_label = f"X Axis ({x_unit})"
            # y-axis remains index

        print("Plotting 2D data (real part) using pcolormesh.")
        # Take real part if complex
        plot_data = np.real(y)

        # Use pcolormesh for potentially non-uniform grids
        # Shading='auto' tries to guess best behavior for pixel vs grid centers
        im = ax.pcolormesh(
            x_coords, y_coords, plot_data, shading="auto", cmap="viridis"
        )
        fig.colorbar(im, ax=ax, label="Intensity (real part, a.u.)")

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_aspect("auto")  # Fallback if ny is 0

        ax.autoscale(tight=True)  # Fit axes to data

    else:
        print(f"Warning: Cannot plot data with {y.ndim} dimensions.")
        return  # Don't show empty plot

    plt.tight_layout()
    if save_if_possible:
        # Save the figure if requested
        print(file_name)
        save_path = Path(file_name).with_suffix(".png")  # Save as PNG
        fig.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")
    plt.show()


# --- Example Usage ---
if __name__ == "__main__":
    print("Running eprload example...")
    print("This will open a file dialog to select a Bruker EPR file.")

    # Example 1: Load data using file dialog and plot
    print("\n--- Example 1: Load with dialog, default scaling, plotting enabled ---")
    x_data, y_data, parameters, file = eprload()  # plot_if_possible is True by default

    if y_data is not None:
        print(f"Successfully loaded: {file}")
        print(f"Data shape: {y_data.shape}")
        if isinstance(x_data, np.ndarray):
            print(f"Abscissa shape: {x_data.shape}")
        elif isinstance(x_data, (list, tuple)):
            print(f"Abscissa shapes: {[ax.shape for ax in x_data]}")
        print(f"Number of parameters loaded: {len(parameters)}")
        # print("Parameters:", parameters) # Uncomment to see all parameters
        # Example accessing a parameter:
        mw_freq = parameters.get("MWFQ", "N/A")  # Get MW frequency if available
        print(f"Microwave Frequency (MWFQ): {mw_freq}")
    else:
        print("Loading failed or was cancelled.")

    # Example 2: Specify a file, apply scaling, suppress plotting
    # Replace 'path/to/your/datafile.DTA' with an actual file path
    # test_file = Path('path/to/your/datafile.DTA')
    # if test_file.exists():
    #      print("\n--- Example 2: Load specific file, scaling='nG', no plot ---")
    #      x_scaled, y_scaled, pars_scaled, f_scaled = eprload(
    #          test_file,
    #          scaling='nG',
    #          plot_if_possible=False
    #      )
    #      if y_scaled is not None:
    #          print(f"Successfully loaded and scaled: {f_scaled}")
    #          print(f"Scaled data shape: {y_scaled.shape}")
    #          # You can plot manually here if needed:
    #          # import matplotlib.pyplot as plt
    #          # plt.figure()
    #          # plt.plot(x_scaled, y_scaled)
    #          # plt.title("Manually Plotted Scaled Data")
    #          # plt.show()
    #      else:
    #          print("Loading/scaling failed for specific file.")
    # else:
    #     print("\nSkipping Example 2: Test file path not found or not set.")

    print("\neprload example finished.")
