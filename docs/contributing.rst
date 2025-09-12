Contributing to EPyR Tools
==========================

We welcome contributions to EPyR Tools! This document provides guidelines for contributing code, documentation, and reporting issues.

Getting Started
---------------

Setting up Development Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Fork and Clone**

.. code-block:: bash

   git clone https://github.com/yourusername/epyrtools.git
   cd epyrtools

2. **Install Development Dependencies**

.. code-block:: bash

   pip install -e .[dev]

   # Or install all optional dependencies
   pip install -e .[all]

3. **Set up Pre-commit Hooks**

.. code-block:: bash

   pre-commit install

This ensures code formatting and quality checks run automatically.

Development Workflow
~~~~~~~~~~~~~~~~~~~~

1. **Create a Feature Branch**

.. code-block:: bash

   git checkout -b feature/your-feature-name

2. **Make Changes**

   * Write code following project conventions
   * Add tests for new functionality
   * Update documentation as needed

3. **Test Your Changes**

.. code-block:: bash

   # Run tests
   pytest

   # Check code formatting
   black --check .
   isort --check-only .
   flake8

4. **Commit and Push**

.. code-block:: bash

   git add .
   git commit -m "Add feature: description"
   git push origin feature/your-feature-name

5. **Create Pull Request**

   * Describe your changes clearly
   * Reference any related issues
   * Ensure all CI checks pass

Code Standards
--------------

Style Guidelines
~~~~~~~~~~~~~~~~

We follow PEP 8 with some specific conventions:

* **Line Length**: Maximum 88 characters (Black default)
* **Imports**: Organized with isort
* **Formatting**: Automatic with Black
* **Linting**: flake8 for code quality

**Example:**

.. code-block:: python

   """Module for EPR data processing."""

   import numpy as np
   from pathlib import Path
   from typing import Tuple, Optional, List

   from epyr.constants import PHYSICAL_CONSTANTS


   def process_epr_data(
       x_data: np.ndarray,
       y_data: np.ndarray,
       correction_order: int = 1
   ) -> Tuple[np.ndarray, np.ndarray]:
       """Process EPR data with baseline correction.

       Args:
           x_data: Field values in Gauss
           y_data: EPR signal intensity
           correction_order: Polynomial order for baseline

       Returns:
           Tuple of (corrected_data, baseline)

       Raises:
           ValueError: If data arrays have different lengths
       """
       if len(x_data) != len(y_data):
           raise ValueError("Data arrays must have same length")

       # Implementation here
       return corrected_data, baseline

Documentation
~~~~~~~~~~~~~

* **Docstrings**: Use Google/NumPy style
* **Type Hints**: Add for all public functions
* **Examples**: Include usage examples in docstrings
* **Comments**: Explain complex algorithms and physics

Testing
-------

Test Structure
~~~~~~~~~~~~~~

Tests are organized in the ``tests/`` directory:

.. code-block:: text

   tests/
   ├── conftest.py          # Pytest configuration
   ├── test_eprload.py      # Data loading tests
   ├── test_baseline.py     # Baseline correction tests
   ├── test_fair.py         # FAIR conversion tests
   └── data/                # Test data files
       ├── sample.dsc
       └── sample.dta

Writing Tests
~~~~~~~~~~~~~

Use pytest for all tests:

.. code-block:: python

   import pytest
   import numpy as np
   from epyr.baseline import baseline_polynomial


   class TestBaselineCorrection:
       """Test baseline correction functionality."""

       def test_polynomial_correction(self):
           """Test polynomial baseline correction."""
           # Create test data with known baseline
           x = np.linspace(0, 100, 1000)
           true_signal = np.exp(-((x - 50)**2) / 100)
           true_baseline = 0.1 * x + 5
           noisy_data = true_signal + true_baseline + np.random.normal(0, 0.01, 1000)

           # Apply correction
           corrected, fitted_baseline = baseline_polynomial(
               noisy_data, x_data=x, poly_order=1
           )

           # Check results
           np.testing.assert_allclose(fitted_baseline, true_baseline, rtol=0.1)
           np.testing.assert_allclose(corrected, true_signal, atol=0.1)

       @pytest.mark.parametrize("order", [0, 1, 2, 3])
       def test_different_orders(self, order):
           """Test different polynomial orders."""
           x = np.linspace(0, 100, 500)
           y = np.ones_like(x) + 0.01 * np.random.randn(len(x))

           corrected, baseline = baseline_polynomial(y, x_data=x, poly_order=order)

           assert len(corrected) == len(y)
           assert len(baseline) == len(y)

Running Tests
~~~~~~~~~~~~~

.. code-block:: bash

   # Run all tests
   pytest

   # Run specific test file
   pytest tests/test_baseline.py

   # Run with coverage
   pytest --cov=epyr

   # Run tests for specific function
   pytest -k "test_polynomial"

Documentation Contributions
----------------------------

Types of Documentation
~~~~~~~~~~~~~~~~~~~~~~~

1. **API Documentation**: Automatically generated from docstrings
2. **User Guides**: Tutorial and how-to content
3. **Examples**: Jupyter notebooks and scripts
4. **Reference**: Technical specifications and algorithms

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Install documentation dependencies
   pip install -e .[docs]

   # Build documentation
   cd docs
   make html

   # View documentation
   open _build/html/index.html

Writing Documentation
~~~~~~~~~~~~~~~~~~~~~

* **Clear Examples**: Always include working code examples
* **Screenshots**: Add images for GUI features
* **Cross-references**: Link to related functions and concepts
* **Math**: Use LaTeX for equations when needed

**Example RST:**

.. code-block:: rst

   Advanced Baseline Correction
   ============================

   The :func:`epyr.baseline.baseline_polynomial` function supports
   advanced baseline correction with signal exclusion.

   .. math::

      y_{corrected} = y_{original} - P_n(x)

   where :math:`P_n(x)` is a polynomial of order :math:`n`.

   Example Usage
   -------------

   .. code-block:: python

      from epyr.baseline import baseline_polynomial

      # Apply quadratic correction excluding peak region
      y_corrected, baseline = baseline_polynomial(
          y_data,
          x_data=x_data,
          poly_order=2,
          exclude_regions=[(3300, 3400)]
      )

Issue Reporting
---------------

Bug Reports
~~~~~~~~~~~

When reporting bugs, please include:

1. **System Information**: OS, Python version, EPyR version
2. **Minimal Example**: Smallest code that reproduces the issue
3. **Expected vs Actual**: What should happen vs what happens
4. **Data Files**: If possible, share problem data files
5. **Error Messages**: Full traceback and error output

**Template:**

.. code-block:: markdown

   ## Bug Report

   **System:**
   - OS: macOS 12.0
   - Python: 3.9.7
   - EPyR Tools: 0.1.6

   **Issue:**
   Baseline correction fails with 2D data

   **Minimal Example:**
   ```python
   import epyr
   x, y, params, _ = epyr.eprload('2d_data.dsc')
   # Error occurs here:
   corrected, baseline = baseline_polynomial(y, x_data=x)
   ```

   **Error:**
   ```
   ValueError: y_data must be a 1D NumPy array.
   ```

   **Expected:** Should handle 2D data or give clear guidance

Feature Requests
~~~~~~~~~~~~~~~~

For feature requests:

1. **Use Case**: Describe the scientific problem
2. **Proposed Solution**: How should it work?
3. **Alternatives**: What workarounds exist?
4. **Examples**: Show expected API usage

EPR Domain Knowledge
--------------------

Contributing EPR-Specific Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding EPR-specific functionality:

* **Physical Accuracy**: Ensure equations and constants are correct
* **Units**: Be explicit about units (Gauss vs Tesla, etc.)
* **Conventions**: Follow EPR community standards
* **References**: Cite relevant papers and textbooks

**Example - g-factor calculation:**

.. code-block:: python

   def calculate_g_factor(frequency_hz: float, field_gauss: float) -> float:
       """Calculate g-factor from resonance condition.

       Uses the EPR resonance condition: hν = gμBB

       Args:
           frequency_hz: Microwave frequency in Hz
           field_gauss: Magnetic field in Gauss

       Returns:
           Dimensionless g-factor

       References:
           Weil, J. A., & Bolton, J. R. (2007). Electron paramagnetic
           resonance: elementary theory and practical applications.
           John Wiley & Sons.
       """
       from epyr.constants import PLANCK_CONSTANT, BOHR_MAGNETON

       # Convert Gauss to Tesla
       field_tesla = field_gauss * 1e-4

       # g = hν / (μB * B)
       return (PLANCK_CONSTANT * frequency_hz) / (BOHR_MAGNETON * field_tesla)

Common Contribution Areas
~~~~~~~~~~~~~~~~~~~~~~~~~

Areas where contributions are especially welcome:

1. **New File Formats**: Support for other spectrometer manufacturers
2. **Analysis Algorithms**: Advanced peak fitting, simulation
3. **Visualization**: Interactive plots, publication templates
4. **Data Processing**: Noise reduction, phase correction
5. **Integration**: Bridges to other software (Origin, MATLAB)

Code Review Process
-------------------

What to Expect
~~~~~~~~~~~~~~

1. **Automated Checks**: CI runs tests and style checks
2. **Maintainer Review**: Core team reviews code and design
3. **Community Feedback**: Other users may comment
4. **Iteration**: Expect requests for changes or improvements

Review Criteria
~~~~~~~~~~~~~~~

* **Correctness**: Does the code work as intended?
* **Testing**: Are there adequate tests?
* **Documentation**: Is it properly documented?
* **Style**: Does it follow project conventions?
* **Performance**: Is it reasonably efficient?
* **Compatibility**: Works with supported Python versions?

Getting Help
------------

If you need help contributing:

1. **GitHub Discussions**: Ask questions about implementation
2. **Issues**: Tag issues with "help wanted" or "good first issue"
3. **Email**: Contact maintainers for major changes
4. **Documentation**: Check existing docs and examples

Recognition
-----------

Contributors are recognized in:

* **CHANGELOG.md**: Major contributions noted in release notes
* **AUTHORS**: List of all contributors
* **Git History**: Detailed commit attribution
* **Documentation**: Citation in relevant sections

Thank you for contributing to EPyR Tools and supporting the EPR research community!
