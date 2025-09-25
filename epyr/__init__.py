"""EPyR Tools - Electron Paramagnetic Resonance Tools in Python."""

# Import configuration and logging first
from .config import config
from .logging_config import setup_logging, get_logger

# Import baseline correction functions from new modular package
from .baseline import *

# Import baseline module for compatibility and convenience
from . import baseline

# Import backend control functions for convenience
from .baseline import setup_inline_backend, setup_widget_backend, setup_notebook_backend

# Import specific, useful components from the old baseline_correction module
# Keep this for backward compatibility during transition
try:
    from .baseline_correction import *
except ImportError:
    # Old baseline_correction.py might be removed/renamed
    pass


from . import lineshapes
from . import signalprocessing
#from .constants import *
from .physics import *
from .eprload import *
from .fair import *
from .isotope_gui import run_gui as isotopes
from .lineshapes import Lineshape, gaussian, lorentzian, voigtian, pseudo_voigt
from .performance import OptimizedLoader, DataCache, get_performance_info
from .plugins import plugin_manager
#from .plot import *
from .eprplot import *
from .sub.utils import BrukerListFiles

__version__ = "0.2.0"

# Set up logging
logger = get_logger(__name__)
logger.debug(f"Package 'epyr' v{__version__} initialized.")
