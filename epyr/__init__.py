"""

"""

# Import specific, useful components from the isotopes_gui module
# The '.' means "from the current package"
from .isotope_gui import run_gui as isotopes
from .eprload import *
from .constants import *
from .plot import *
from .fair import *
from .sub.baseline2 import *
from .sub.utils import BrukerListFiles
from . import baseline

__version__ = "0.1.0"


# You can add package-level initialization code here if needed in the future.
print(f"Package 'your_project_folder' initialized.") # Optional: confirmation message