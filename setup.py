#!/usr/bin/env python3
"""
Setup script for EPyR Tools - Electron Paramagnetic Resonance Tools in Python.

This setup.py provides backward compatibility and alternative installation method
alongside the modern pyproject.toml configuration.

Installation:
    pip install .                    # Standard installation
    pip install -e .                 # Development installation (editable)
    pip install -e .[dev]           # Development installation with dev tools
    pip install -e .[docs]          # Installation with documentation tools

Usage after installation:
    import epyr
    from epyr import eprload, baseline, constants
    from epyr.fair import convert_bruker_to_fair
    from epyr.isotope_gui import run_gui
"""

import os
from pathlib import Path

from setuptools import find_packages, setup

# Read long description from README
here = Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")


# Read version from epyr/__init__.py
def get_version():
    """Extract version from epyr/__init__.py"""
    version_file = here / "epyr" / "__init__.py"
    with open(version_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split('"')[1]
    raise RuntimeError("Unable to find version string")


# Define package data and include requirements
def get_package_data():
    """Get package data files"""
    data_files = []

    # Include example data files
    examples_dir = here / "examples"
    if examples_dir.exists():
        for root, dirs, files in os.walk(examples_dir):
            for file in files:
                if file.endswith((".dta", ".dsc", ".spc", ".par", ".ipynb", ".py")):
                    rel_path = os.path.relpath(os.path.join(root, file), here)
                    data_files.append(rel_path)

    return data_files


setup(
    # Basic package information
    name="epyr-tools",
    version=get_version(),
    author="Sylvain Bertaina",
    author_email="sylvain.bertaina@cnrs.fr",
    maintainer="Sylvain Bertaina",
    maintainer_email="sylvain.bertaina@cnrs.fr",
    # Package description
    description="Electron Paramagnetic Resonance (EPR) Tools in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # URLs
    url="https://github.com/BertainaS/epyrtools",
    project_urls={
        "Homepage": "https://github.com/BertainaS/epyrtools",
        "Repository": "https://github.com/BertainaS/epyrtools.git",
        "Issues": "https://github.com/BertainaS/epyrtools/issues",
        "Documentation": "https://github.com/BertainaS/epyrtools/tree/main/docs",
        "Laboratory": "https://www.im2np.fr/fr/equipe-magnetisme-mag",
    },
    # Package discovery and content
    packages=find_packages(exclude=["tests", "tests.*"]),
    package_data={
        "epyr": ["*.py"],
        "": get_package_data(),
    },
    include_package_data=True,
    # Requirements
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "h5py>=3.1.0",
        "pandas>=1.3.0",
    ],
    # Optional dependencies
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
            "black>=21.0.0",
            "flake8>=3.8.0",
            "mypy>=0.812",
            "pre-commit>=2.10.0",
        ],
        "docs": [
            "sphinx>=6.0.0",
            "sphinx-rtd-theme>=1.3.0",
            "myst-parser>=2.0.0",
            "scipy>=1.7.0",
        ],
        "all": [
            # Include all optional dependencies
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
            "black>=21.0.0",
            "flake8>=3.8.0",
            "mypy>=0.812",
            "pre-commit>=2.10.0",
            "sphinx>=6.0.0",
            "sphinx-rtd-theme>=1.3.0",
            "myst-parser>=2.0.0",
            "scipy>=1.7.0",
        ],
    },
    # Entry points for command-line tools
    entry_points={
        "console_scripts": [
            "epyr=epyr.cli:main",
            "epyr-convert=epyr.cli:cmd_convert",
            "epyr-baseline=epyr.cli:cmd_baseline",
            "epyr-batch-convert=epyr.cli:cmd_batch_convert",
            "epyr-config=epyr.cli:cmd_config",
            "epyr-info=epyr.cli:cmd_info",
            "epyr-isotopes=epyr.cli:cmd_isotopes",
            "epyr-plot=epyr.cli:cmd_plot",
            "epyr-validate=epyr.cli:cmd_validate",
        ],
    },
    # Classification
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    # License
    license="MIT",
    # Keywords
    keywords=[
        "EPR",
        "electron paramagnetic resonance",
        "spectroscopy",
        "Bruker",
        "data analysis",
        "FAIR",
        "scientific computing",
    ],
    # Additional metadata
    zip_safe=False,  # Required for package data access
)
