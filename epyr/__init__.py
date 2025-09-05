""" """

# Import specific, useful components from the isotopes_gui module
# The '.' means "from the current package"
from . import baseline
from .constants import *
from .eprload import *
from .fair import *
from .isotope_gui import run_gui as isotopes
from .plot import *
from .sub.utils import BrukerListFiles

__version__ = "0.1.2"


# You can add package-level initialization code here if needed in the future.
print(f"Package 'epyr' v{__version__} initialized.")  # Optional: confirmation message
