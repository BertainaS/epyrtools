"""Pytest configuration and shared fixtures for EPyR Tools tests."""

import tempfile
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pytest


@pytest.fixture
def sample_1d_data():
    """Generate sample 1D EPR data for testing."""
    x = np.linspace(-50, 50, 1000)  # Magnetic field range
    # Gaussian peak with some baseline
    y = 10 + 0.1 * x + 100 * np.exp(-(x**2) / 50)
    return x, y


@pytest.fixture
def sample_2d_data():
    """Generate sample 2D EPR data for testing."""
    nx, ny = 100, 50
    x_axis = np.linspace(-10, 10, nx)
    y_axis = np.linspace(0, 100, ny)
    # Create a 2D dataset with some structure
    X, Y = np.meshgrid(x_axis, y_axis)
    Z = 5 + 0.1 * Y + 50 * np.exp(-((X**2) / 4 + (Y - 50) ** 2 / 100))
    return x_axis, y_axis, Z


@pytest.fixture
def sample_epr_params():
    """Generate sample EPR parameter dictionary."""
    return {
        "MWFQ": 9.4e9,  # Microwave frequency in Hz
        "MWPW": 20.0,  # Microwave power in dB
        "AVGS": 1000,  # Number of averages
        "RCAG": 30.0,  # Receiver gain in dB
        "STMP": 298.15,  # Temperature in K
        "SPTP": 81.92,  # Sweep time in ms
        "XMIN": -50.0,  # Field minimum
        "XWID": 100.0,  # Field width
        "XAXIS_NAME": "Magnetic Field",
        "XAXIS_UNIT": "mT",
        "YAXIS_NAME": "Intensity",
        "YAXIS_UNIT": "a.u.",
    }


@pytest.fixture
def temp_data_files():
    """Create temporary test data files."""
    temp_dir = Path(tempfile.mkdtemp())

    # Create some dummy files with correct extensions
    test_files = {
        "test.dta": b"dummy_bruker_bes3t_data",
        "test.dsc": b"dummy_bruker_bes3t_descriptor",
        "test.spc": b"dummy_bruker_esp_data",
        "test.par": b"dummy_bruker_esp_parameters",
    }

    file_paths = {}
    for filename, content in test_files.items():
        file_path = temp_dir / filename
        file_path.write_bytes(content)
        file_paths[filename] = file_path

    yield file_paths

    # Cleanup
    import shutil

    shutil.rmtree(temp_dir)


@pytest.fixture
def mock_bruker_files(tmp_path):
    """Create mock Bruker EPR data files with proper structure."""
    # Create a simple .dsc file (text-based descriptor)
    dsc_content = """#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
#SPL	1	Sample
#OPER    Test Operator
#DATE    01/01/24 10:00:00
#CMNT    Test EPR measurement
*
#BSEQ	0
#IKKF	REAL
#XTYP	IDX
#YTYP	NODATA
#ZTYP	NODATA
*
#XPTS	1024
#XMIN	320.000000
#XWID	10.000000
#XNAM	Field
#XUNI	G
*
#MWFQ	9.400000e+009
#MWPW	2.000000e+001
#AVGS	1000
#SPTP	8.192000e+001
#RCAG	3.000000e+001
#RCHM	1
#B0MA	3.250000e+002
#B0MF	1.000000e+001
*
#RRES	1024
#BSEQ	0
#RESO	1024,1
#WDTH	1.000000e+003
#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
"""

    # Create .dta file (binary data)
    dta_data = np.random.random(1024).astype(np.float64)

    dsc_file = tmp_path / "test.dsc"
    dta_file = tmp_path / "test.dta"

    dsc_file.write_text(dsc_content)
    dta_file.write_bytes(dta_data.tobytes())

    return {"dsc_path": dsc_file, "dta_path": dta_file, "base_path": tmp_path / "test"}


@pytest.fixture
def baseline_test_data():
    """Generate test data with known baseline for correction testing."""
    x = np.linspace(0, 100, 200)

    # Known polynomial baseline
    true_baseline = 2 * x + 10

    # Add some signal (peaks)
    signal = 50 * np.exp(-((x - 30) ** 2) / 20) + 30 * np.exp(-((x - 70) ** 2) / 15)

    # Combine baseline + signal
    noisy_data = true_baseline + signal + np.random.normal(0, 1, len(x))

    return {
        "x": x,
        "y_with_baseline": noisy_data,
        "true_baseline": true_baseline,
        "pure_signal": signal,
        "signal_regions": [(25, 35), (65, 75)],  # Where the peaks are
    }
