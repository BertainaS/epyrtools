"""
Isotope data loading and processing for the GUI application.

This module handles loading isotope data from files, processing it,
and providing data access methods for the GUI components.
"""

import os
import traceback
from typing import Optional, Tuple

import pandas as pd

from .gui_helpers import NMAGN, PLANCK


class IsotopeDataManager:
    """Manages loading and processing of isotope data."""

    def __init__(self):
        self.full_data: Optional[pd.DataFrame] = None
        self.element_layout_data: Optional[pd.DataFrame] = None
        self.table_display_data_cols = [
            "isotope",
            "abundance",
            "spin",
            "gn",
            "gamma",
            "qm",
        ]

    def load_data(self) -> pd.DataFrame:
        """Load isotope data from the data file.

        Returns:
            pandas.DataFrame: The loaded and processed isotope data.

        Raises:
            FileNotFoundError: If the isotope data file cannot be found.
            Exception: If there's an error reading or processing the file.
        """
        data_file_path = self._find_data_file()
        return self._read_and_process_data(data_file_path)

    def _find_data_file(self) -> str:
        """Find the isotope data file in possible locations.

        Returns:
            str: Path to the data file

        Raises:
            FileNotFoundError: If the file cannot be found
        """
        script_dir = os.path.dirname(os.path.abspath(__file__))
        possible_paths = [
            os.path.join(os.getcwd(), "sub", "isotopedata.txt"),
            os.path.join(
                script_dir, "..", "sub", "isotopedata.txt"
            ),  # Adjusted for new structure
            os.path.join(script_dir, "sub", "isotopedata.txt"),
        ]

        for path in possible_paths:
            if os.path.exists(path):
                print(f"Loading data from: {path}")
                return path

        error_msg = (
            "Could not find 'sub/isotopedata.txt'. "
            f"Checked paths:\n- {chr(10).join(possible_paths)}"
        )
        print(f"Error: {error_msg}")
        raise FileNotFoundError(error_msg)

    def _read_and_process_data(self, data_file_path: str) -> pd.DataFrame:
        """Read and process the isotope data file.

        Args:
            data_file_path: Path to the data file

        Returns:
            pandas.DataFrame: Processed isotope data

        Raises:
            Exception: If there's an error reading or processing the file
        """
        try:
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

            # Read the CSV file
            data = pd.read_csv(
                data_file_path,
                comment="%",
                sep=r"\s+",  # Regex for one or more whitespace
                names=col_names,
                na_values=["-"],
                skipinitialspace=True,
            )

            # Process and clean the data
            data = self._clean_and_convert_data(data)
            data = self._calculate_derived_values(data)
            data = self._generate_isotope_symbols(data)
            data = self._fill_missing_values(data)

            # Store element layout data
            self.element_layout_data = data.drop_duplicates(
                subset=["Z"], keep="first"
            ).copy()

            self.full_data = data
            return data

        except Exception as e:
            print(f"Error reading or processing data file '{data_file_path}': {e}")
            traceback.print_exc()
            raise

    def _clean_and_convert_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """Clean and convert data types."""
        # Convert numeric columns
        data["Z"] = pd.to_numeric(data["Z"], errors="coerce")
        data["N"] = pd.to_numeric(data["N"], errors="coerce")

        numeric_cols = ["spin", "gn", "abundance", "qm"]
        for col in numeric_cols:
            if col in data.columns:
                data[col] = pd.to_numeric(data[col], errors="coerce")

        # Remove rows missing critical Z or N
        data.dropna(subset=["Z", "N"], inplace=True)

        # Convert Z/N to integer after cleaning NaNs
        data["Z"] = data["Z"].astype(int)
        data["N"] = data["N"].astype(int)

        # Ensure string columns are strings and fill NaNs
        str_cols = ["radioactive", "element", "name"]
        for col in str_cols:
            if col in data.columns:
                data[col] = data[col].astype(str).fillna("")

        return data

    def _calculate_derived_values(self, data: pd.DataFrame) -> pd.DataFrame:
        """Calculate gamma/(2pi) [MHz/T] from gn values."""
        if "gn" in data.columns and pd.api.types.is_numeric_dtype(data["gn"]):
            data["gamma"] = data["gn"].apply(
                lambda gn_val: (
                    (gn_val * NMAGN / PLANCK / 1e6) if pd.notna(gn_val) else pd.NA
                )
            )
        else:
            data["gamma"] = pd.NA

        return data

    def _generate_isotope_symbols(self, data: pd.DataFrame) -> pd.DataFrame:
        """Generate isotope symbols (e.g., 1H, 14C*, 235U)."""
        isotopes = []
        required_cols = ["N", "Z", "element", "radioactive"]

        if all(col in data.columns for col in required_cols):
            for _, row in data.iterrows():
                if pd.notna(row["N"]) and pd.notna(row["Z"]):
                    mass_number = int(row["N"] + row["Z"])  # A = N + Z
                    iso_str = f"{mass_number}{row['element']}"
                    if row["radioactive"] == "*":
                        isotopes.append(iso_str + "*")
                    else:
                        isotopes.append(iso_str)
                else:
                    isotopes.append(pd.NA)
        else:
            print("Warning: Missing required columns for isotope symbol generation.")
            isotopes = [pd.NA] * len(data)

        data["isotope"] = isotopes
        data["isotope"] = data["isotope"].astype(str)
        return data

    def _fill_missing_values(self, data: pd.DataFrame) -> pd.DataFrame:
        """Fill remaining NaN values with appropriate defaults."""
        data["spin"] = data["spin"].fillna(-1.0)  # Marker for missing
        data["abundance"] = data["abundance"].fillna(0.0)
        data["qm"] = data["qm"].fillna(0.0)
        data["gn"] = data["gn"].fillna(0.0)
        data["gamma"] = data["gamma"].fillna(0.0)
        return data

    def get_table_data(self) -> pd.DataFrame:
        """Get data formatted for table display.

        Returns:
            pandas.DataFrame: Table data with NMRfreq column initialized to 0.0
        """
        if self.full_data is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        table_data = self.full_data[self.table_display_data_cols].copy()
        table_data["NMRfreq"] = 0.0
        return table_data

    def get_element_data(self) -> pd.DataFrame:
        """Get element data for periodic table layout.

        Returns:
            pandas.DataFrame: Element layout data
        """
        if self.element_layout_data is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        return self.element_layout_data

    def update_nmr_frequencies(
        self, table_data: pd.DataFrame, field_value: float
    ) -> pd.DataFrame:
        """Update NMR frequencies based on field value.

        Args:
            table_data: DataFrame with isotope data
            field_value: Magnetic field in mT

        Returns:
            pandas.DataFrame: Updated table data with calculated NMR frequencies
        """
        # Convert field from mT to T
        field_tesla = field_value / 1000.0

        # Calculate NMR frequencies: f = gamma * B / (2*pi)
        table_data["NMRfreq"] = table_data["gamma"] * field_tesla

        return table_data
