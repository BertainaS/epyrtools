epyr.isotope_gui module
=======================

The ``epyr.isotope_gui`` module provides an interactive graphical user interface for exploring nuclear isotope properties relevant to EPR and NMR spectroscopy.

Overview
--------

The Isotope GUI is a Tkinter-based application that displays:

* **Periodic table interface**: Click elements to view isotope data
* **Nuclear properties**: Spin, abundance, magnetogyric ratios
* **NMR frequencies**: At specified magnetic field strengths
* **EPR relevance**: Hyperfine coupling information

Main Components
---------------

GUI Application
~~~~~~~~~~~~~~~

.. automodule:: epyr.isotope_gui.main_window
   :members: IsotopeGUI, run_gui
   :undoc-members:
   :show-inheritance:

GUI Components
~~~~~~~~~~~~~~

.. automodule:: epyr.isotope_gui.gui_components
   :members:
   :undoc-members:
   :show-inheritance:

Isotope Data Management
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: epyr.isotope_gui.isotope_data
   :members:
   :undoc-members:
   :show-inheritance:

Helper Functions
~~~~~~~~~~~~~~~~

.. automodule:: epyr.isotope_gui.gui_helpers
   :members:
   :undoc-members:
   :show-inheritance:

Usage Examples
--------------

Launch the GUI
~~~~~~~~~~~~~~

.. code-block:: python

   from epyr.isotope_gui import run_gui
   
   # Launch the interactive isotope GUI
   run_gui()

This opens a window with an interactive periodic table where you can:

1. **Click on elements** to view available isotopes
2. **Select isotopes** to see detailed nuclear properties  
3. **Set magnetic field** to calculate NMR frequencies
4. **Export data** for use in calculations

Features
--------

Interactive Periodic Table
~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Color-coded elements**: Distinguish metals, non-metals, etc.
* **Click functionality**: Access isotope data instantly
* **Hover information**: Quick property preview
* **Search function**: Find elements by name or symbol

Nuclear Property Display
~~~~~~~~~~~~~~~~~~~~~~~~

For each isotope, the GUI shows:

* **Nuclear spin** (I): Quantum number
* **Natural abundance**: Percentage in nature
* **Magnetogyric ratio** (γ): MHz/T
* **NMR frequency**: At specified field strength
* **Quadrupole moment**: For I > 1/2 nuclei

Magnetic Field Calculator
~~~~~~~~~~~~~~~~~~~~~~~~~

* **Variable field input**: Enter any magnetic field value
* **Unit selection**: Tesla, Gauss, or mT
* **Real-time calculation**: NMR frequencies update automatically
* **EPR field presets**: Common X-band, Q-band frequencies

Export Functionality
~~~~~~~~~~~~~~~~~~~~

* **Copy to clipboard**: Selected isotope data
* **Save to file**: Complete isotope database
* **Print friendly**: Formatted for documentation

EPR Applications
----------------

Hyperfine Structure Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GUI helps identify nuclear isotopes that contribute to EPR hyperfine structure:

* **I = 1/2 nuclei**: Simple doublet splitting (¹H, ¹³C, ¹⁵N, ³¹P)
* **I = 1 nuclei**: Triplet splitting (²H, ¹⁴N)  
* **I > 1 nuclei**: Complex multiplet patterns

Natural Abundance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **High abundance**: Strong hyperfine coupling (¹H: 99.98%)
* **Low abundance**: Weak but observable coupling (¹³C: 1.1%)
* **Isotope effects**: Compare different isotopes of same element

Magnetic Field Scaling
~~~~~~~~~~~~~~~~~~~~~~

Calculate how hyperfine couplings scale with magnetic field:

* **g-factor effects**: How electron resonance changes with field
* **Nuclear Zeeman**: How nuclear levels shift with field
* **Hyperfine constants**: Field-independent vs. field-dependent terms

Example Workflow
----------------

.. code-block:: python

   from epyr.isotope_gui import run_gui
   
   # 1. Launch the GUI
   run_gui()
   
   # 2. In the GUI:
   #    - Click on Hydrogen (H) element
   #    - Select ¹H isotope  
   #    - Set magnetic field to 0.335 T (X-band EPR)
   #    - Note NMR frequency: ~14.3 MHz
   
   # 3. This information helps interpret EPR spectra:
   #    - Proton hyperfine coupling ~10 G = ~28 MHz
   #    - Nuclear Zeeman splitting ~14 MHz
   #    - Expect complex hyperfine pattern

Technical Implementation
------------------------

.. automodule:: epyr.isotope_gui
   :members:
   :undoc-members:
   :show-inheritance:

The GUI is built using:

* **Tkinter**: Cross-platform GUI framework
* **Canvas widgets**: Interactive periodic table display
* **Event handling**: Mouse clicks and keyboard input
* **Data binding**: Real-time calculation updates

Complete API
------------

.. automodule:: epyr.isotope_gui
   :members:
   :undoc-members:
   :show-inheritance: